import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('./helpers')
import loadmat as lm
import glob
import helpers
import tensortools as tt
import os
import scipy.ndimage as spi
from multiprocessing import Pool
from sklearn.preprocessing import normalize






def _fast_occ(occupancy,trials,bins):
    """ calculate occupancy for each position bin returns occupancy bins x trials"""
    
    for i,j in zip(trials,bins):
        if (j<0) or j>=occupancy.shape[0] or i>=occupancy.shape[1]:
            pass
        else:
            occupancy[j,i]+=1
            

#@numba.jit(nopython=True)
def _fast_bin(counts, trials, bins, neurons):
    """
    Given coordinates of spikes, compile binned spike counts. Throw away
    spikes that are outside of tmin and tmax.
    Turns into a matrix neurons x bins x trials
    """
    for i, j, k in zip(trials, bins, neurons):
        if (j < 0) or (int(j) >= counts.shape[1]) or i>=counts.shape[2]:
            pass
        else:
            counts[k, int(j), i] += 1


def prepareData(data):
    pos_edges = np.arange(0,401,5)
    posx = data['posx']
    posx[posx>=400]=399.99
    posx[posx<0]=0
    location_vec=np.digitize(posx,pos_edges)-1

    good_cells = data['sp']['cids'][data['sp']['cgs']==2]
    validSpikes = np.in1d(data['sp']['clu'],good_cells)
    spike_clu = data['sp']['clu'][validSpikes]
    (bla,spike_idx) = np.unique(spike_clu,return_inverse=True)

    spiketimes = np.digitize(data['sp']['st'][validSpikes],data['post'])
    spikelocations = location_vec[spiketimes]
    trial_idx = data['trial'][spiketimes]-1

    return (good_cells,pos_edges,trial_idx,spikelocations,spike_idx,location_vec)


def run_pipeline(filename):

    data = lm.loadmat(filename)
    session_name = os.path.basename(filename)[0:-4]
    (good_cells,pos_edges,trial_idx,spikelocations,spike_idx,location_vec) = prepareData(data)
    n_trials = 30
    n_cells = len(good_cells)
    shape = (n_cells, len(pos_edges)-1, n_trials)
    counts = np.zeros(shape, dtype=float)
    _fast_bin(counts,trial_idx,spikelocations,spike_idx)

    occupancy = np.zeros((len(pos_edges)-1,n_trials),dtype = float)
    _fast_occ(occupancy,data['trial']-1,location_vec)

    for iT in range(n_trials):
        tmp = occupancy[:,iT]
        idx_v = np.flatnonzero(tmp)
        idx_n = np.flatnonzero(tmp==0)
        tmp[idx_n]=np.interp(idx_n,idx_v,tmp[idx_v])
        occupancy[:,iT]=tmp

    spMapN = np.zeros(counts.shape)
    for iC in range(n_cells):
        spMapN[iC,:,:]=np.divide(counts[iC,:,:],occupancy)

    spMapN = spi.gaussian_filter(spMapN,(0,2,0))
    

   
    n_cells = len(good_cells)
    n_bins = len(pos_edges)-1
    spFlat = np.zeros((n_cells,n_trials*n_bins))

    for iC in range(n_cells):
        spFlat[iC,:]=spMapN[iC,:,:].ravel(order='F')
    #spFlat = spFlat-spFlat.mean(axis=1)[:,np.newaxis]
    spFlat = normalize(spFlat,axis=0,norm='l2')
    for iC in range(n_cells):
        for iT in range(n_trials):
            start = iT*n_bins
            stop = (iT+1)*n_bins
            trial_idx = np.arange(start,stop)
            tmp = spFlat[iC,trial_idx]
            spMapN[iC,:,iT]=tmp

    R=5
    # Fit CP tensor decomposition (two times).
    U = tt.ncp_bcd(spMapN, rank=R, verbose=False)
    V = tt.ncp_bcd(spMapN, rank=R, verbose=False)

    
    # Align the two fits and print a similarity score.
    sim = tt.kruskal_align(U.factors, V.factors, permute_U=True, permute_V=True)
    #print(sim)

    # Plot the results again to see alignment.
    fig, ax, po = tt.plot_factors(U.factors)
    tt.plot_factors(V.factors, fig=fig)
    fig.suptitle("aligned models")
    fig.tight_layout()
    fig.savefig('C:\\temp\\try3\\'+session_name+'_tca.png')

    ff=np.matmul(np.transpose(spFlat),spFlat)
    plt.figure()
    ax=plt.imshow(ff)
    plt.colorbar()
    plt.axvline(x=n_bins*20,color='red',ls='--',linewidth=1)
    plt.axvline(x=n_bins*21,color='green',ls='--',linewidth=1)
    plt.axhline(y=n_bins*20,color='red',ls='--',linewidth=1)
    plt.axhline(y=n_bins*21,color='green',ls='--',linewidth=1)
    plt.savefig('C:\\temp\\try3\\'+session_name+'_cov.png')
    plt.close('all')


if __name__=='__main__':
    files = glob.glob('F:/NP_DATA/np*_gain*.mat')
    #files = glob.glob('Z:/giocomo/attialex/NP_DATA/np*_gain*.mat')
    for iF in files:
        print(iF)
        run_pipeline(iF)
    
    p=Pool(processes =5)
    p.map(run_pipeline,files)

