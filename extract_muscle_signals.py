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
import glob
import h5py
import cv2
import cPickle
import numpy as np
import muscle_model as mm

# have to do a little work to import muscle_model
this_dir, this_filename = os.path.split(__file__)
model_path = os.path.join(this_dir, "models")

################################################################################################
################################ GENERAL INPUTS ################################################
################################################################################################
FLY_DB_PATH = '/media/sam/SamData/FlyDB'  # directory containing fly data
# FLY_DB_PATH = '/media/imager/DataExternal/FlyDB'

FLIES = [24]  # list of flies to analyze. If None, look for command line input
FOLDER_STR = 'Fly%04d'  # general format of folders containing fly data -- only used if FLIES contains ints
FN_STR = 'ca_cam'  # string to search for when looking for files containing image data
FN_EXT = '.hdf5'  # extension to search for when looking for files containing image data. *** Only hdf5 supported!
MIN_FILE_SIZE = 1e6  # minimum camera data file size (in bytes)
MODEL_PATH = model_path  # directory containing Thad's muscle models
DRIVER = 'GMR22H05'

IMG_KEY = 'cam_imgs'  # key name in hdf5 file for images. if None, we'll try to guess
T_KEY = 'cam_tstamps'  # key name in hdf5 file for camera time stamps. if None, we'll skip this

CHUNK_SIZE = 100  # how many images to process at a time (should also correspond to hdf5 storage)
FIT_MODE = 'pinv'  # 'pinv' or 'nnls' ; Thad's code seems to default to pinv


################################################################################################
############################ HELPER FUNCTION(S) ################################################
################################################################################################
def load_rframe_fits(data_path, side, file_suffix=''):
    """
    Quick function to load in reference frame fits saved during data collection. Ideally this
    can be altered as needed, since I'm not sure we've landed on the final way to store them

    Should output a dictionary with the rframe values
    """
    # first, get all the files in the data folder that could potentially contain rframe fits
    rframe_paths = glob.glob(os.path.join(data_path, '*%s_rframe_fits.hdf5' %(file_suffix)))

    # define list of variables we care about (later functions expect certain key/value pairings
    out_keys = ['p', 'a1', 'a2', 'A', 'Ainv', 'components']

    # intialize output
    rframe_dict = dict()

    # assume we don't need a key to access side-specific data, but we'll load one in if need be
    side_key = None

    # now we need to figure out if we saved one for each side, or if there's one file with both sides
    if (len(rframe_paths) == 1) and (side in rframe_paths):
        # in this case, we assume we have one file, but that's because we only imaged on one side
        rframe_path = rframe_paths[0]

    elif (len(rframe_paths) == 1) and (side not in rframe_paths):
        # in this case, we assume we have one hdf5 file with multiple groups corresponding to sides
        rframe_path = rframe_paths[0]
        with h5py.File(rframe_path, 'r') as h5f:
            side_keys = [key for key in h5f.keys() if (side in key)]
            if len(side_keys) != 1:
                print('Could not find unique group identifier for rframe hdf5 file, side %s' %(side))
                return rframe_dict
            else:
                side_key = side_keys[0]
    else:
        # in this case, assume we have one file per side, and read out that way
        rframe_path = [p for p in rframe_paths if (side in p)][0]

    # now that we've hopefully got the right path/key info, read data
    # NB: need to tread numeric and list data differently
    with h5py.File(rframe_path, 'r') as h5f:

        for key in out_keys:
            # read data
            if side_key:
                data = h5f[side_key][key][:]
            else:
                try:
                    data = h5f['LogRefFrame/' + key][:]
                except KeyError:
                    data = h5f[key][:]

            # write data to dict
            if isinstance(data, np.ndarray):
                rframe_dict[key] = np.unique(data)
            else:
                rframe_dict[key] = data

    # return dict
    return rframe_dict


################################################################################################
#################### FUNCTIONS FOR DOING EXTRACTION ############################################
################################################################################################
def fit_to_model(imchunk, model, fit_mode=FIT_MODE, fit_pix_mask=None, baseline=None):
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
    if fit_mode == 'nnls':
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
        print(np.shape(im_array))
        from numpy.linalg import pinv
        if not (fit_pix_mask is None):
            fits = np.dot(pinv(model[:, fit_pix_mask]).T, im[:, fit_pix_mask].T)
        else:
            fits = np.dot(pinv(model).T, im)

    print('Done with chunk')
    return fits


# ------------------------------------------------------------------------------------------------------------
def extract_gcamp_signals(imgs, fly_frame_dict, driver=DRIVER, model_type='volumetric', model_name='thorax',
                          model_dir=MODEL_PATH, fit_mode=FIT_MODE, fit_pix_mask=None, baseline=None, 
                          chunk_sz=CHUNK_SIZE):
    """
    Function to run the "fit_to_model" code on a given fly. Code taken from planotaxis/viewer.py

    INPUTS:
        - imgs :  N x h x w GCaMP image array, where N is the number of images and h,w are the height and width
        - fly_frame_dict : dict containing reference frame. allows transformation from confocal data to current fly
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
        parent_dir = os.path.abspath(os.path.join(this_dir, '..'))
        unmixing_filters = os.path.join(parent_dir, 'muscle_model', 'unmixing_filters',
                                        'NA_0.45_200mm_Tube_FN1', 'flatened_model.hdf5')
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
    num_samps = np.shape(imgs)[0]
    print(num_samps)
    chunks = [slice(x, x + chunk_sz if x + chunk_sz < num_samps else num_samps) for x in range(0, num_samps, chunk_sz)]

    # create list of image chunks and also make copies of other inputs to 'fit_to_model'
    img_chunks = [np.array(imgs[chunk]) for chunk in chunks]
    models = [model for chunk in chunks]
    modes = [fit_mode for chunk in chunks]
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
                         folder_str=FOLDER_STR, min_file_size=MIN_FILE_SIZE, img_key=IMG_KEY, t_key=T_KEY,
                         chunk_size=CHUNK_SIZE, fit_mode=FIT_MODE):
    """
    A little wrapper function to run GCaMP muscle signal extraction on a given fly
    """
    # define string we'll use for savename suffix
    save_str = '_model_fits'

    # get path to current fly data, as well as appropriate suffix
    if isinstance(fly_id, int):
        fly_path = os.path.normpath(os.path.join(fly_db_path, folder_str % (fly_id)))
        file_suffix = ''
    else:
        fly_path = os.path.normpath(fly_db_path)
        file_suffix = fly_id
        file_suffix = '_' + os.path.splitext(file_suffix)[0]  # make sure we don't have the '.bag' at the end

    # find cam data files containing images for
    cam_data_fns = [f for f in os.listdir(fly_path) if fn_str in f and f.endswith(file_suffix + fn_ext)]

    # remove files that are too small (this is mostly for me, since I'm still saving hdf5 files for cameras I'm not
    # using. should fix this eventually, but until then this workaround should be okay)
    cam_data_fns = [f for f in cam_data_fns if os.path.getsize(os.path.join(fly_path, f)) > min_file_size]

    # loop over camera data files (NB: we're often expecting one for left and one for right)
    for fn in cam_data_fns:
        # check if analysis has already been performed
        fn_no_ext = os.path.splitext(fn)[0]
        if os.path.exists(os.path.join(fly_path, fn_no_ext + save_str + '.hdf5')):
            print('Already analyzed %s -- skipping' % (fn))
            continue

        # otherwise print update for current file
        print('Analyzing %s ...' % (fn))

        # --------------------------------------------
        # load reference frame fits (or make them)
        rframe_fn = os.path.splitext(fn)[0] + '_rframe_fits'
        rframe_glob = glob.glob(os.path.join(fly_path, '*' + file_suffix[1:] + '_rframe_fits.hdf5'))

        if os.path.exists(os.path.join(fly_path, rframe_fn + '.cpkl')):
            with open(os.path.join(fly_path, rframe_fn + '.cpkl'), 'rb') as pf:
                rframe = cPickle.load(pf)
        elif os.path.exists(os.path.join(fly_path, rframe_fn + '.hdf5')):
            # for this option (which I think is ultimately better) we would need to save an extra rframe file for each
            # camera, which is not what we're doing now (based on how the rostopics are published)
            print('Under construction!')
        elif len(rframe_glob) > 0:
            # in this case, we have rframe hdf5 files stored with different filename convention from cam files
            # first figure out which side we're dealing with
            if 'left' in fn:
                side = 'left'
            elif 'right' in fn:
                side = 'right'
            else:
                print('Could not determine imaging side (needed to load rframe files)')
                return
            # then load rframe using helper function
            rframe = load_rframe_fits(fly_path, side=side, file_suffix=file_suffix)
        else:
            # in this case, we don't have anything saved. would ideally like to have this open a 'viewer' GUI window
            # (but that might be more work than is worth it)
            print('Cannot find model fit to load -- must generate this')
            print('Under construction! For now, use viewer.py GUI')
            return

        # --------------------------------------------
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

        # run analysis on images
        signals = extract_gcamp_signals(imgs, rframe, chunk_sz=chunk_size, fit_mode=fit_mode)

        # ** add time to dictionary
        signals['t'] = tstamps

        # --------------------------------------
        # save output (all)?
        if save_flag:
            # save output (cpkl)
            save_fn = os.path.splitext(fn)[0] + save_str
            with open(os.path.join(fly_path, save_fn + '.cpkl'), 'wb') as sf:
                cPickle.dump(signals, sf)

            # save output (hdf5)
            with h5py.File(os.path.join(fly_path, save_fn + '.hdf5'), 'w') as hsf:
                for key in signals.keys():
                    hsf.create_dataset(key, data=signals[key])

        # let us know we're done
        print('Completed %s' % (fn))

    # why not?
    return


################################################################################################
########################### RUN SCRIPT #########################################################
################################################################################################
if __name__ == '__main__':
    # read in flies from terminal input or specified list
    if not FLIES:
        flies = [int(x) for x in sys.argv[1:]]
    else:
        flies = FLIES

    # loop through flies and run analysis
    for fly in flies:
        run_gcamp_extraction(fly)
