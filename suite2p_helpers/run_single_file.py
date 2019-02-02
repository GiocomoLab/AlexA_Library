import sys
import os
import time

from suite2p.run_s2p import run_s2p
import numpy as np

def default_ops(animal='',expid=''):
    ops = {
            'fast_disk': '', # used to store temporary binary file, defaults to save_path0 (set as a string NOT a list)
            'save_path0': '', # stores results, defaults to first item in data_path
            'delete_bin': False, # whether to delete binary file after processing
            'h5py_key': 'data',
            # main settings
            'nplanes' : 2, # each tiff has these many planes in sequence
            'nchannels' : 1, # each tiff has these many channels per plane
            'functional_chan' : 1, # this channel is used to extract functional ROIs (1-based)
            'diameter':10, # this is the main parameter for cell detection, 2-dimensional if Y and X are different (e.g. [6 12])
            'tau':  1., # this is the main parameter for deconvolution
            'fs': 30.92,  # sampling rate (total across planes)
            # output settings
            'save_mat': False, # whether to save output as matlab files
            'combined': True, # combine multiple planes into a single result /single canvas for GUI
            # parallel settings
            'num_workers': 0, # 0 to select num_cores, -1 to disable parallelism, N to enforce value
            'num_workers_roi': -1, # 0 to select number of planes, -1 to disable parallelism, N to enforce value
            # registration settings
            'do_registration': True, # whether to register data
            'nimg_init': 200, # subsampled frames for finding reference image
            'batch_size': 1000, # number of frames per batch
            'maxregshift': 0.1, # max allowed registration shift, as a fraction of frame max(width and height)
            'align_by_chan' : 1, # when multi-channel, you can align by non-functional channel (1-based)
            'reg_tif': False, # whether to save registered tiffs
            'subpixel' : 10, # precision of subpixel registration (1/subpixel steps)
            # cell detection settings
            'connected': True, # whether or not to keep ROIs fully connected (set to 0 for dendrites)
            'navg_frames_svd': 5000, # max number of binned frames for the SVD
            'nsvd_for_roi': 1000, # max number of SVD components to keep for ROI detection
            'max_iterations': 20, # maximum number of iterations to do cell detection
            'ratio_neuropil': 6., # ratio between neuropil basis size and cell radius
            'ratio_neuropil_to_cell': 3, # minimum ratio between neuropil radius and cell radius
            'tile_factor': 1., # use finer (>1) or coarser (<1) tiles for neuropil estimation during cell detection
            'threshold_scaling': 1., # adjust the automatically determined threshold by this scalar multiplier
            'max_overlap': 0.75, # cells with more overlap than this get removed during triage, before refinement
            'inner_neuropil_radius': 2, # number of pixels to keep between ROI and neuropil donut
            'outer_neuropil_radius': np.inf, # maximum neuropil radius
            'min_neuropil_pixels': 350, # minimum number of pixels in the neuropil
            # deconvolution settings
            'baseline': 'maximin', # baselining mode
            'win_baseline': 60., # window for maximin
            'sig_baseline': 10., # smoothing constant for gaussian filter
            'prctile_baseline': 8.,# optional (whether to use a percentile baseline)
            'neucoeff': .7,  # neuropil coefficient
        }
    return ops

def default_db():
    db = {
        'h5py': '', # a single h5 file path
        'h5py_key': 'data',
        'fast_disk':''
        #'fast_disk': os.path.join(os.environ['L_SCRATCH'],"s2ptmp"), # string which specifies where the binary file will be stored (should be an SSD)
    }
    return db

def get_expID(filename = 'AA_190111_029_000_001.h5'):
    h,t = filename.split('.')
    parts = h.split('_')
    return parts[-1]


if __name__ == '__main__':
    rootpath=os.path.join(os.environ['OAK'],'attialex','TEST')
    animal = 'AA_190110_023'
    animal_path = os.path.join(rootpath,animal)
    for file in os.listdir(animal_path):
        if file.endswith('.h5'):
            expID=get_expID(file)
            hf_path = os.path.join(animal_path,file)
            savepath=os.path.join(animal_path,file[:-3])
            tmp_path = os.path.join(os.environ['L_SCRATCH'],'suite2p',file[:-3])
            
            db = default_db()
            db['h5py']=hf_path
            #db['fast_disk']=tmp_path
            ops = default_ops()
            ops['save_path0'] = savepath
            ops['fast_disk'] = tmp_path
            # run one experiment
            start = time.time()

            #print('Running for: '+db['h5py']+'\n')
            #print('Saving on: '+ops['save_path0'])
            #print('Tmp saving:' + ops['fast_disk'])
            opsEnd=run_s2p(ops=ops,db=db)
            stop = time.time()
            timelog = open(os.path.join(os.environ['SCRATCH'],'suite2p_times.txt'),'a')
            timelog.write(db['h5py']+', ' +str(stop-start) +'\n')
            timelog.close()
            
            