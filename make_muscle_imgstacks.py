import numpy as np
import os

#create an image stack

#read images
file_path = '/home/imager/work/data/muscle_imaging_integration_bagfiles/toConvert/MI_093022/ca_camera_left_MI_2022-10-04-17-03-36.hdf5'

with h5py.File(file_path, 'r') as h5f:
    if not img_key:
        im_key = h5f.keys()[0]
    else:
        im_key = img_key
    imgs = np.array(h5f[im_key])

    if t_key:
        tstamps = np.array(h5f[T_KEY])
    else:
        tstamps = None


#self.cam_topic_names = ['/ca_camera_left/image_raw',
#                    '/ca_camera_right/image_raw',


