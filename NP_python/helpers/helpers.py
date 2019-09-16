import numpy as np
import scipy as sp
def preprocess(data,nth_bin = 10):
    track_start = 0
    track_end = 400
    dx=5
    dt=0.2
    every_nth_time_bin = nth_bin
    numposbins = np.floor((track_end-track_start)/dx)
    posx_edges = np.linspace(track_start,track_end,numposbins+1)
    posx=data['posx']
    post=data['post']
    trial = data['trial']
    sp =  data['sp']
    
    # resample post, posx, and trial according to desired dt
    post_resampled = post[0::every_nth_time_bin]
    posx_resampled=posx
    posx_resampled[posx_resampled<track_start]=track_start
    posx_resampled[posx_resampled>track_end]=track_end #now happening further down
    #posx_resampled = posx[0::every_nth_time_bin]
    trial_resampled = trial[0::every_nth_time_bin]

    # get cell ids of "good" units
    good_cells = sp['cids'][sp['cgs']==2]

    # time bins for position decoding
    numtimebins = len(post_resampled)
    post_edges = np.squeeze(np.linspace(min(post)-dt/2,max(post)+dt/2,numtimebins+1))
    post_centers = post_edges[range(0,len(post_edges)-1)]+dt/2

    # posx categories for position decoding (binned)
    posx_bin = np.digitize(posx_resampled,posx_edges)
    posx_bin = posx_bin[0::every_nth_time_bin]
    posx_resampled = posx_resampled[0::every_nth_time_bin]

    # count spikes in each time bin for each cell
    spikecount = np.empty((len(good_cells),len(post_resampled),))
    spikecount[:] = np.nan
    for cell_idx in range(len(good_cells)):   
        spike_t = sp['st'][sp['clu']==good_cells[cell_idx]]
        spikecount[cell_idx,:] = np.histogram(spike_t,bins=post_edges)[0]
    data['spikecount']=np.transpose(spikecount)
    data['posx_bin']=posx_bin
    data['trial_resampled']=trial_resampled
    data['posx_resampled']=posx_resampled
    return data