import numpy as np
import scipy as sp
import scipy.ndimage as spi
from scipy.ndimage import gaussian_filter1d
import pandas as pd

def preprocess(data,nth_bin = 10):
    track_start = 0
    track_end = 400
    dx=5
    dt=0.2
    every_nth_time_bin = nth_bin
    numposbins = np.floor((track_end-track_start)/dx)
    posx_edges = np.linspace(track_start,track_end,numposbins+1)
    posx_centers = 0.5 * posx_edges[0:-1] + 0.5*posx_edges[1::]
    data['posx_centers']=posx_centers
    data['posx_edges']=posx_edges
    posx=data['posx']
    post=data['post']
    trial = data['trial']
    sp =  data['sp']
    
    # resample post, posx, and trial according to desired dt
    post_resampled = post[0::every_nth_time_bin]
    posx_resampled=posx
    posx_resampled[posx_resampled<track_start]=track_start
    posx_resampled[posx_resampled>=track_end]=track_end-0.001 #now happening further down
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

    #speed
    speed = calcSpeed(data['posx'])
    speed_resampled = speed[0::every_nth_time_bin]

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
    data['speed_resampled']=speed_resampled
    return data

def calcSpeed(posx):
    speed = np.diff(posx)/0.02
    speed = np.hstack((0,speed))
    speed[speed>150]=np.nan
    speed[speed<-5]=np.nan
    idx_v = np.flatnonzero(np.logical_not(np.isnan(speed)))
    idx_n = np.flatnonzero(np.isnan(speed))
    speed[idx_n]=np.interp(idx_n,idx_v,speed[~np.isnan(speed)])
    speed = spi.gaussian_filter1d(speed,10)
    return speed

def _fast_occ(occupancy,trials,bins):

    for i,j in zip(trials,bins):
        if (j<0) or j>=occupancy.shape[0] or i>=occupancy.shape[1]:
            pass
        else:
            occupancy[j,i]+=1

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
            
def calculateFiringRate(data,good_cells=None,t_edges = None):
    if good_cells is None:
        good_cells = data['sp']['cids'][data['sp']['cgs']==2]
    if t_edges is None:
        dt= 0.2;
        t_edges = np.arange(0,data['sp']['st'].max()+dt,dt)
    else:
        dt = np.mean(np.diff(t_edges))
    # count spikes in each time bin for each cell
    
    spikecount = np.full((len(good_cells),len(t_edges)-1),np.nan)
    
    for cell_idx in range(len(good_cells)):   
        spike_t = data['sp']['st'][data['sp']['clu']==good_cells[cell_idx]]
        spikecount[cell_idx,:] = np.histogram(spike_t,bins=t_edges)[0]

      
    spikecount = np.hstack((spikecount,np.zeros((spikecount.shape[0],1))))  
    spikerate = spikecount/dt
    spikes = np.transpose(spikerate)
    X = gaussian_filter1d(spikes, 2, axis=0)
    return spikes,X,t_edges

class options:
    
    def __init__(self):
        self.speed_t=0.05;
        self.extract_win = [-2,3];
        self.aux_win = [-50,50];
        self.TimeBin = 0.02;
        self.time_bins =np.arange(-2,3,0.02);
        self.extract_win = [-2,3]
        self.speedSigma = 10;
        self.smoothSigma_time = 0.2; # in sec; for smoothing fr vs time
        self.smoothSigma_dist = 2; # in cm; for smoothing fr vs distance
        self.SpatialBin = 2;
        self.TrackStart = 0
        self.TrackEnd = 400
        self.SpeedCutof = 2
        self.stab_thresh = 0.5
        self.max_lag = 30
                
    @property            
    def time_vecs(self):
        return self.time_bins[0:-1]*0.5 + self.time_bins[1:]*0.5
    @property
    def xbinedges(self):
        return np.arange(self.TrackStart,self.TrackEnd+self.SpatialBin,self.SpatialBin)
    @property
    def xbincent(self):
        return self.xbinedges[0:-1]+self.SpatialBin/2


def calculateFiringRateMap(data,trials2extract=None,good_cells = None,ops=None):
    posx=np.mod(data['posx'],ops.TrackEnd)
    post=data['post']
    trial = data['trial'] 
    sp = data['sp']
    if good_cells is None:
        good_cells = sp['cids'][sp['cgs']==2]
    if trials2extract is None:
        trials2extract = np.arange(trial.min(),trial.max()+1)
    
    
    posx_bin = np.digitize(posx,ops.xbinedges)
    validSpikes = np.in1d(data['sp']['clu'],good_cells)
    spike_clu = data['sp']['clu'][validSpikes]
    (bla,spike_idx) = np.unique(spike_clu,return_inverse=True)
    spiketimes = np.digitize(data['sp']['st'][validSpikes],data['post'])
    spikelocations = posx_bin[spiketimes]-1 # to start at 0
    spiketrials = data['trial'][spiketimes] # to start at 0
    
    valid_trialsSpike = np.in1d(spiketrials,trials2extract)
    spiketimes = spiketimes[valid_trialsSpike]
    spikelocations = spikelocations[valid_trialsSpike]
    spiketrials = spiketrials[valid_trialsSpike]
    spike_idx=spike_idx[valid_trialsSpike] 
    
    valid_trials = np.in1d(trial,trials2extract)
    occupancy = np.zeros((len(ops.xbinedges)-1,len(trials2extract)),dtype = float)
    
    _fast_occ(occupancy,trial[valid_trials]-trials2extract[0],posx_bin[valid_trials]-1)
    occupancy *=ops.TimeBin
    
    n_cells = len(good_cells)
    shape = (n_cells, len(ops.xbinedges)-1, len(trials2extract))
    counts = np.zeros(shape, dtype=float)
    _fast_bin(counts,spiketrials-spiketrials.min(),spikelocations,spike_idx)
    spMapN = np.zeros(counts.shape)
    stab =np.zeros(n_cells)
    for iC in range(n_cells):
        tmp = np.divide(counts[iC,:,:],occupancy)
        df = pd.DataFrame(tmp)
        df.interpolate(method='pchip', axis=0, limit=None, inplace=True)
        tmp = df.values
        #print((np.isnan(tmp).sum()))
        tmp_f = gaussian_filter1d(tmp,ops.smoothSigma_dist, axis=0,mode='wrap')
        spMapN[iC]=tmp_f
        cc=np.corrcoef(np.transpose(tmp_f))
        

        stab[iC]=np.nanmean(cc[np.triu(np.full(cc.shape,True),1)])
    
    return counts,spMapN,stab

def calculateSimilarityMatrix():
    pass