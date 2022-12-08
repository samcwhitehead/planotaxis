# -*- coding: utf-8 -*-
"""
Function that extracts muscle time series signals from images. Code taken from Thad/Alysha's planotaxis/viewer.py.

Basically want to be able to run analysis in a loop, but I think most of the functionality is already there

As in my unbag.py, you should be able to give as input either a fly number (assuming Thad/Alysha dir structure) or a
filename string, corresponding to the original bagfile name (which should include a datetime, and thus keep things
relatively organized)
"""
# -----------------------------------------------------------------
import sys
import os
import h5py
import cv2
import cPickle
import numpy as np

# have to do a little work to import muscle_model
this_dir, this_filename = os.path.split(__file__)
model_path = os.path.join(this_dir, "models")
# sys.path.insert(0, os.path.join(model_path, "muscle_model"))
import muscle_model as mm

################################################################################################
################################ GENERAL INPUTS ################################################
################################################################################################
FLY_DB_PATH = '/media/sam/SamData/FlyDB'  # directory containing fly data
# FLY_DB_PATH = '/media/imager/DataExternal/FlyDB'

FLIES = None  # list of flies to analyze. If None, look for command line input
FOLDER_STR = 'Fly%04d'  # general format of folders containing fly data -- only used if FLIES contains ints
FN_STR = 'ca_cam'  # string to search for when looking for files containing image data
FN_EXT = '.hdf5'  # extension to search for when looking for files containing image data. *** Only hdf5 supported!
MIN_FILE_SIZE = 1e6  # minimum camera data file size (in bytes)
MODEL_PATH = model_path  # directory containing Thad's muscle models
DRIVER = 'GMR39E01'

IMG_KEY = 'cam_imgs'  # key name in hdf5 file for images. if None, we'll try to guess
T_KEY = 'cam_tstamps'  # key name in hdf5 file for camera time stamps. if None, we'll skip this


################################################################################################
#################### FUNCTIONS FOR DOING EXTRACTION ############################################
################################################################################################
def fit_to_model(imchunk, model, mode='pinv', fit_pix_mask=None, baseline=None):
    """
    Function to fit data to model thorax. Code taken from planotaxis/viewer.py

    INPUTS
        - imchunk : M x h x w GCaMP image array, where N is the number of images and h,w are the height and width
            NB: M here is typically a set chunk size, since these fits are done on subsets of the full set of frames
        - model :
        - mode : method for performing fits. Options are:
            - 'pinv' (pseudo-inverse with numpy)
            - 'nnls' (non-linear least squares with scipy)
        - fit_pix_mask :
        - baseline : baseline value/image, to be subtracted off from all images in imchunk

    OUTPUTS:
        - fits :
    """
    # first, subtract off baseline value if we have one
    if not (baseline is None):
        im_array = imchunk - baseline  # /baseline
    else:
        im_array = imchunk
    imshape = np.shape(im_array[0])
    im_array = im_array.reshape((-1, imshape[0] * imshape[1]))
    if mode == 'nnls':
        fits = np.empty((np.shape(model)[0], np.shape(im_array)[0]))
        for i, im2 in enumerate(im_array):
            im = im2.copy()
            im[~np.isfinite(im)] = 0
            from scipy.optimize import nnls
            if not (fit_pix_mask is None):
                fits[:, i] = nnls(model[:, fit_pix_mask].T, im[fit_pix_mask])[0]
            else:
                fits[:, i] = nnls(model.T, im)[0]
    else:
        im = im_array
        print np.shape(im_array)
        from numpy.linalg import pinv
        if not (fit_pix_mask is None):
            fits = np.dot(pinv(model[:, fit_pix_mask]).T, im[:, fit_pix_mask].T)
        else:
            fits = np.dot(pinv(model).T, im)

    print('Done with chunk')
    return fits


# ------------------------------------------------------------------------------------------------------------
def extract_gcamp_signals(imgs, fly_frame_dict, driver=DRIVER, model_type='volumetric', model_name='thorax',
                            model_dir=MODEL_PATH, mode='pinv', fit_pix_mask=None, baseline=None):
    """
    Function to run the "fit_to_model" code on a given fly. Code taken from planotaxis/viewer.py

    INPUTS:
        - imgs :  N x h x w GCaMP image array, where N is the number of images and h,w are the height and width
        - fly_frame_dict : dict containting reference frame. allows transformation from confocal data to current fly
        - driver : string giving the driver line used for imaging (e.g. 'GMR39E01' or 'GMR22H05'). Will need to expand
            to accommodate sensitivity of new GCaMP variants? Could also just directly put in muscle list...
        - model_type : either 'volumetric' or 'masks'. Not sure I understand yet, but we use 'volumetric'
        - model_name : either 'thorax' or 'muscle_model'. Tbh, not sure I understand the difference, but we use thorax
        - model_dir : where to find confocal model data (outlines.cpkl). probably looks like:
            ../planotaxis/models/thorax
        - mode : to be passed to fit_to_model (defined above)
        - fit_pix_mask : to be passed to fit_to_model (defined above)
        - baseline : to be passed to fit_to_model (defined above)

    OUTPUTS:
        -
    """
    # load the reference frame of the confocal data (will compare to 'fly_frame')
    confocal_model = mm.GeometricModel(filepath=os.path.join(model_dir, '%s'%(model_name), 'outlines.cpkl'))
    confocal_frame = confocal_model.frame

    # convert the input fly reference frame into a Frame object (just transfer over info from dict)
    fly_frame = mm.Frame()
    for key in fly_frame_dict.keys():
        fly_frame[key] = fly_frame_dict[key]

    # get the transformation matrix A (going from fly -> confocal frame) and compose with a scaling of s
    # to construct a transformation for homogeneous vectors
    s = 1  # scale
    A = fly_frame.get_transform(confocal_frame)
    Ap = np.dot([[s, 0.0, 0], [0, s, 0], [0, 0, 1]], A)

    # use the driver line name to get the list of muscles
    with open(os.path.join(model_dir, '%s'%(model_name), 'profiles', '%s.cpkl'%(driver)), 'rb') as ff:
        driver_profile = cPickle.load(ff)
    muscles = driver_profile['selected_components']

    # make sure we've converted images to a numpy array
    if not isinstance(imgs, (np.ndarray, np.generic)):
        imgs = np.array(imgs)

    # read the first image to get the output shape of the warped model
    output_shape = np.shape(imgs[0])

    # ----------------------------------------------------------------------------------------------------
    # split demixing type depeding on model_type. as far as I know, we almost exclusively use volumetric
    if model_type == 'masks':
        # get the mask of all the muscles for fit
        masks = confocal_model.get_masks(fly_frame, np.shape(imgs[0]))
        # create the model using only the muscle Childcare/Lesson Combo: that express in a given line
        model = np.vstack([masks[mask_key].T.ravel().astype(float) for mask_key in muscles])
        # construct a mask do reduce the projection to just the data within the model
        fit_pix_mask = np.sum(model, axis=0) > 0
    if model_type == 'volumetric':
        # load hdf5 file containing filter params specific to current objective
        unmixing_filters = os.path.join(model_path,
                                        'muscle_model',
                                        'unmixing_filters',
                                        'NA_0.45_200mm_Tube_FN1',
                                        'flatened_model.hdf5')
        model_data = h5py.File(unmixing_filters, 'r')
        model_muscles = [np.array(model_data[muscle]) for muscle in muscles]
        output_shapes = [output_shape for muscle in muscles]
        transforms = [Ap[:-1, :] for muscle in muscles]
        model = map(cv2.warpAffine, model_muscles, transforms, output_shapes)
        model.append(np.ones_like(model[0]))
        muscles.append('bkg')
        model = np.vstack([muscle.T.ravel() for muscle in model])
        fit_pix_mask = np.ones_like(model[0]) > 0

    # -----------------------------------------------------------------------
    # split imgs array into smaller chunks for processing
    chnk_sz = 100
    num_samps = np.shape(imgs)[0]
    print num_samps
    chunks = [slice(x, x + chnk_sz if x + chnk_sz < num_samps else num_samps) for x in range(0, num_samps, chnk_sz)]

    # create list of image chunks and also make copies of other inputs to 'fit_to_model'
    img_chunks = [np.array(imgs[chunk]) for chunk in chunks]
    models = [model for chunk in chunks]
    modes = ['pinv' for chunk in chunks]
    fit_pix_masks = [fit_pix_mask for chunk in chunks]
    baselines = [baseline for chunk in chunks]  # NB: this is a holdover from bg subtraction code that I removed for now

    # run fits on chunks
    fits = map(fit_to_model, img_chunks, models, modes, fit_pix_masks, baselines)

    # --------------------------------------------------------------------------------
    # return output
    signals_dict = dict()
    [signals_dict.update({str(mname): sig}) for sig, mname in zip(np.hstack(fits), muscles)]

    return signals_dict


# ------------------------------------------------------------------------------------------------------------
def run_gcamp_extraction(fly_id, fly_db_path=FLY_DB_PATH, fn_str=FN_STR, fn_ext=FN_EXT, save_flag=True,
                         folder_str=FOLDER_STR, min_file_size=MIN_FILE_SIZE, img_key=IMG_KEY, t_key=T_KEY):
    """
    A little wrapper function to run GCaMP muscle signal extraction on a given fly
    """
    # get path to current fly data, as well as appropriate suffix
    if isinstance(fly_id, int):
        fly_path = os.path.normpath(os.path.join(fly_db_path, folder_str % (fly_id)))
        file_suffix = ''
    else:
        fly_path = os.path.normpath(fly_db_path)
        file_suffix = fly_id
        file_suffix = os.path.splitext(file_suffix)[-1]  # make sure we don't have the '.bag' at the end

    # find cam data files containing images for
    cam_data_fns = [f for f in os.listdir(fly_path) if fn_str in f and f.endswith(file_suffix + fn_ext)]

    # remove files that are too small (this is mostly for me, since I'm still saving hdf5 files for cameras I'm not
    # using. should fix this eventually, but until then this workaround should be okay)
    cam_data_fns = [f for f in cam_data_fns if os.path.getsize(os.path.join(fly_path, f)) > min_file_size]

    # loop over camera data files (NB: we're often expecting one for left and one for right)
    for fn in cam_data_fns:
        print('Analyzing %s ...' % (fn))
        # load images for current data file
        # NB: I'm mimicking Thad's method here, but should probably explicitly look for the right h5py key
        with h5py.File(os.path.join(fly_path, fn), 'r') as h5f:
            if not img_key:
                im_key = h5f.keys()[0]
            else:
                im_key = img_key
            imgs = np.array(h5f[im_key])

            if t_key:
                tstamps = np.array(h5f[T_KEY])
            else:
                tstamps = None

        # load reference frame fits (or make them)
        rframe_fn = os.path.splitext(fn)[0] + '_rframe_fits'
        if os.path.exists(os.path.join(fly_path, rframe_fn + '.cpkl')):
            with open(os.path.join(fly_path, rframe_fn + '.cpkl'), 'rb') as pf:
                rframe = cPickle.load(pf)
        elif os.path.exists(os.path.join(fly_path, rframe_fn + '.hdf5')):
            print('Under construction!')
        else:
            print('Cannot find model fit to load -- must generate this')
            print('Under construction!')

        # run analysis on images
        signals = extract_gcamp_signals(imgs, rframe)

        # ** add time to dictionary
        signals['t'] = tstamps

        # --------------------------------------
        # save output (all)?
        if save_flag:
            # save output (cpkl)
            save_fn = os.path.splitext(fn)[0] + '_model_fits'
            with open(os.path.join(fly_path, save_fn + '.cpkl'), 'wb') as sf:
                cPickle.dump(signals, sf)

            # save output (hdf5)
            with h5py.File(os.path.join(fly_path, save_fn + '.hdf5'), 'w') as hsf:
                for key in signals.keys():
                    hsf.create_dataset(key, data=signals[key])

        # let us know we're done
        print('Completed %s' % (fn))

        # return 'signals' dict in case we want to use it
        return signals


################################################################################################
########################### RUN SCRIPT #########################################################
################################################################################################
if __name__ == '__main__':
    print(aaaaaa) 
    # read in flies from terminal input or specified list
    if not FLIES:
        flies = [int(x) for x in sys.argv[1:]]
    else:
        flies = FLIES

    # loop through flies and run analysis
    for fly in flies:
        run_gcamp_extraction(fly)
