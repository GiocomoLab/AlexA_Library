import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.stats 
from scipy.ndimage import gaussian_filter1d
import glob
import sys
sys.path.append('./helpers')
import loadmat as lm
import helpers as helpers
import pandas as pd
from sklearn import linear_model

def _fast_occ(occupancy,trials,bins):
    """ calculate occupancy for each position bin returns occupancy bins x trials"""
    
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
            
    #return counts

from scipy.linalg import cho_factor, cho_solve
from sklearn.base import BaseEstimator

class CircularRegression(BaseEstimator):
    
    def __init__(self, alpha=0.0, tol=1e-5, max_iter=100):
        self.alpha = alpha
        self.tol = tol
        self.max_iter = max_iter
    
    def fit(self, X, y):
        """
        Parameters
        ----------
        X : array
            Independent variables, has shape (n_timepoints x n_neurons)
        y : array
            Circular dependent variable, has shape (n_timepoints x 1),
            all data should lie on the interval [-pi, +pi].
        """
        
        # Convert 1d circular variable to 2d representation
        u = np.column_stack([np.sin(y), np.cos(y)])

        # Randomly initialize weights. Ensure scaling does
        W = np.random.randn(X.shape[1], 2)
        W /= np.max(np.sum(X @ W, axis=1))

        # Cache neuron x neuron gram matrix. This is used below
        # in the M-step to solve a linear least squares problem
        # in the form inv(XtX) @ XtY. Add regularization term to
        # the diagonal.
        XtX = X.T @ X
        XtX[np.diag_indices_from(XtX)] += self.alpha
        XtX = cho_factor(XtX)

        # Compute model prediction in 2d space, and projection onto
        # each observed u.
        XW = (X @ W)
        t = np.sum(u * XW, axis=1)
        tcdf = normcdf(t)
        tpdf = normpdf(t)

        self.log_like_hist_ = [
            np.log(2 * np.pi) - 
            0.5 * np.mean(np.sum(XW * XW, axis=1), axis=0) +
            np.mean(np.log(1 + t * tcdf / tpdf))
        ]

        for itr in range(self.max_iter):

            # E-step.
            m = t + (tcdf / (tpdf + t * tcdf))
            XtY = X.T @ (m[:, None] * u)

            # M-step.
            W = cho_solve(XtX, XtY)
            
            # Recompute model prediction.
            XW = X @ W
            t = np.sum(u * XW, axis=1)
            tcdf = normcdf(t)
            tpdf = normpdf(t)

            # Store log-likelihood.
            self.log_like_hist_.append(
                np.log(2 * np.pi) - 
                0.5 * np.mean(np.sum(XW * XW, axis=1), axis=0) +
                np.mean(np.log(1 + t * tcdf / tpdf))
            )
            
            # Check convergence.
            if (self.log_like_hist_[-1] - self.log_like_hist_[-2]) < self.tol:
                break
    
        self.weights_ = W
    
    def predict(self, X):
        u_pred = X @ self.weights_
        return np.arctan2(u_pred[:, 0], u_pred[:, 1])

    def score(self, X, y):
        """
        Returns 1 minus mean angular similarity between y and model prediction.
        
        score == 1 for perfect prediction
        score == 0 in expectation for random prediction
        score == -1 if prediction is off by 180 degrees.
        """
        y_pred = self.predict(X)
        return np.mean(np.cos(y - y_pred))

normcdf = scipy.stats.norm.cdf
normpdf = scipy.stats.norm.pdf

SMOOTHNESS = 10.0
REGULARIZATION = 1e-4

track_start = 0
track_end = 400
dx=2
dt=0.02
every_nth_time_bin = 1
numposbins = np.floor((track_end-track_start)/dx)
posx_edges = np.arange(0,402,2)

def run_for_file(data,TRIALS):
    try:
        anatomy = data['anatomy']
    except:
        print('no anatomy')
        return None
    
    if 'parent_shifted' in anatomy:
        group = anatomy['parent_shifted']
    else:
        group = anatomy['cluster_parent']
    regions = ('MEC','VISp','RS')
    idx = [str(ss).startswith(regions) for ss in group]
    idx = np.array(idx)
    posx=np.mod(data['posx'],track_end)
    post=data['post']
    trial = data['trial']
    sp = data['sp']
    good_cells = sp['cids'][np.logical_and(idx,sp['cgs']==2)]


    # posx categories for position decoding (binned)
    posx_bin = np.digitize(posx,posx_edges)
    validSpikes = np.in1d(data['sp']['clu'],good_cells)
    spike_clu = data['sp']['clu'][validSpikes]
    (bla,spike_idx) = np.unique(spike_clu,return_inverse=True)
    spiketimes = np.digitize(data['sp']['st'][validSpikes],data['post'])
    spikelocations = posx_bin[spiketimes]-1
    trial_idx = data['trial'][spiketimes]-1

    occ2 = np.zeros((len(posx_edges)-1,TRIALS.max()),dtype = float)
    _fast_occ(occ2,data['trial']-1,posx_bin-1)
    n_cells = len(good_cells)
    shape = (n_cells, len(posx_edges)-1, 20)
    counts = np.zeros(shape, dtype=float)
    _fast_bin(counts,trial_idx,spikelocations,spike_idx)
    spMapN = np.zeros(counts.shape)
    stab =np.zeros(n_cells)
    for iC in range(n_cells):
        tmp = np.divide(counts[iC,:,:],occ2)
        df = pd.DataFrame(tmp)
        df.interpolate(method='pchip', axis=0, limit=None, inplace=True)
        tmp = df.values
        #print((np.isnan(tmp).sum()))
        tmp_f = gaussian_filter1d(tmp,3, axis=0,mode='wrap')
        cc=np.corrcoef(np.transpose(tmp_f[:,TRIALS-1]))
        f=cc

        stab[iC]=np.nanmean(f[np.triu(np.full((16,16),True),1)])
        #spMapN[iC,:,:]=np.divide(counts[iC,:,:],occ2)


    # count spikes in each time bin for each cell
    spikecount = np.empty((len(good_cells),len(post)-1,))
    spikecount[:] = np.nan
    for cell_idx in range(len(good_cells)):   
        spike_t = sp['st'][sp['clu']==good_cells[cell_idx]]
        spikecount[cell_idx,:] = np.histogram(spike_t,bins=post)[0]

    fr = spikecount.sum(axis=1)/post.max()
    valid_idx = np.logical_and(stab>.5,fr>=1)
    if sum(valid_idx)<5:
        return None    
    spikecount = np.hstack((spikecount,np.zeros((spikecount.shape[0],1))))  
    spikerate = spikecount/dt
    spikes = np.transpose(spikerate)
    spikes = spikes[:,valid_idx]
    X = gaussian_filter1d(spikes, SMOOTHNESS, axis=0)

    speed = helpers.calcSpeed(posx)
    position = np.mod(posx,track_end)
    n_units = spikes.shape[1]
    speed_idx = speed>5
    trial_idx = np.in1d(trial,TRIALS)
    valid_idx = np.logical_and(speed_idx,trial_idx) #only take data from trials of interest and above speed threshold
    speed_r = speed[valid_idx]
    speed_r = speed_r[0::every_nth_time_bin]
    position_r = position[valid_idx]
    position_r = position_r[0::every_nth_time_bin]
    trial_r = trial[valid_idx]
    trial_r = trial_r[0::every_nth_time_bin]
    posbin_r = posx_bin[valid_idx]
    posbin_r = posbin_r[0::every_nth_time_bin]

    X_r = X[valid_idx,:]
    X_r = scipy.stats.zscore(X_r,axis=0)
    X_r = X_r[0::every_nth_time_bin,:]
    #theta = ((position_r - position_r.min()) / position_r.max()) * 2 * np.pi - np.pi
    theta = ((position_r) / track_end) * 2 * np.pi - np.pi

    #model = CircularRegression(alpha=REGULARIZATION)
    model = linear_model.LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial', max_iter=10000, C = 0.1)
    pos_bin = posbin_r-1
    for iFold,tr in enumerate(TRIALS):
        take_idx = np.full((len(TRIALS,)), True)
        take_idx[iFold]=False
        train_TRIALS = TRIALS[take_idx]
               
        train_idx = np.in1d(trial_r,train_TRIALS)
        test_idx = np.in1d(trial_r,tr)
        
        tc = np.full((X_r.shape[1],200),np.nan);
        for ii in range(200):
            tc[:,ii] = X_r[np.logical_and(pos_bin==ii,train_idx),:].mean(axis = 0)

        df = pd.DataFrame(tc)
        df.interpolate(method='pchip', axis=1, limit=None, inplace=True)
        tc = df.values
        ff=np.dot(tc.transpose(),np.transpose(X_r[test_idx,:]))
        g=np.argmax(ff,axis=0)
        
        
        
        
        #model.fit(X_r[train_idx,:], theta[train_idx])
        #tmp = model.predict(X_r[test_idx,:])
        #pred_pos.append(tmp)
        #tmp_e = theta[test_idx]-tmp
        model.fit(X_r[train_idx,:], pos_bin[train_idx])
        tmp = model.predict(X_r[test_idx,:])
        tmp_e = theta[test_idx]-tmp
        tmp_e_mld = posx_edges[posbin_r[test_idx]] - posx_edges[g]
        if iFold == 0:
            error_list_mld = tmp_e_mld
            error_list = tmp_e
            predicted_theta = tmp
            predicted_bin = g+1
        else:
            error_list = np.hstack((error_list,tmp_e))
            error_list_mld= np.hstack((error_list_mld,tmp_e_mld))
            predicted_theta = np.hstack((predicted_theta,tmp))
            predicted_bin = np.hstack((predicted_bin,g+1))
        


    output=dict()
    output['cluID']=good_cells
    output['error_cd'] = error_list
    output['error_mld'] = error_list_mld
    output['pos_bin'] = pos_bin
    output['true_pos'] = position_r
    output['true_theta'] = theta
    output['predicted_theta'] = predicted_theta
    output['predicted_bin'] = predicted_bin
    output['region'] = group[np.logical_and(idx,sp['cgs']==2)]
    output['trial'] = trial_r
    return output

def get_gain_onsets(data,gain,contrast):
    ff = np.logical_and(data['trial_gain']==.8,data['trial_contrast']==contrast)
    onsets = np.argwhere(np.diff(np.double(ff))==1)+1
    return onsets

def run_for_file_gain(data,TRIALS):
    try:
        anatomy = data['anatomy']
    except:
        print('no anatomy')
        return None
    
    if 'parent_shifted' in anatomy:
        group = anatomy['parent_shifted']
    else:
        group = anatomy['cluster_parent']
    regions = ('MEC','VISp','RS')
    idx = [str(ss).startswith(regions) for ss in group]
    idx = np.array(idx)
    posx=np.mod(data['posx'],track_end)
    post=data['post']
    trial = data['trial']
    sp = data['sp']
    good_cells = sp['cids'][np.logical_and(idx,sp['cgs']==2)]


    # posx categories for position decoding (binned)
    posx_bin = np.digitize(posx,posx_edges)
    validSpikes = np.in1d(data['sp']['clu'],good_cells)
    spike_clu = data['sp']['clu'][validSpikes]
    (bla,spike_idx) = np.unique(spike_clu,return_inverse=True)
    spiketimes = np.digitize(data['sp']['st'][validSpikes],data['post'])
    spikelocations = posx_bin[spiketimes]-1
    trial_idx = data['trial'][spiketimes]-1

    occ2 = np.zeros((len(posx_edges)-1,TRIALS.max()),dtype = float)
    _fast_occ(occ2,data['trial']-1,posx_bin-1)
    n_cells = len(good_cells)
    shape = (n_cells, len(posx_edges)-1, TRIALS.max())
    counts = np.zeros(shape, dtype=float)
    _fast_bin(counts,trial_idx,spikelocations,spike_idx)
    spMapN = np.zeros(counts.shape)
    stab =np.zeros(n_cells)
    for iC in range(n_cells):
        tmp = np.divide(counts[iC,:,:],occ2)
        df = pd.DataFrame(tmp)
        df.interpolate(method='pchip', axis=0, limit=None, inplace=True)
        tmp = df.values
        #print((np.isnan(tmp).sum()))
        tmp_f = gaussian_filter1d(tmp,3, axis=0,mode='wrap')
        cc=np.corrcoef(np.transpose(tmp_f[:,TRIALS[0:6]-1]))
        f=cc

        stab[iC]=np.nanmean(f[np.triu(np.full((6,6),True),1)])
        #spMapN[iC,:,:]=np.divide(counts[iC,:,:],occ2)


    # count spikes in each time bin for each cell
    spikecount = np.empty((len(good_cells),len(post)-1,))
    spikecount[:] = np.nan
    for cell_idx in range(len(good_cells)):   
        spike_t = sp['st'][sp['clu']==good_cells[cell_idx]]
        spikecount[cell_idx,:] = np.histogram(spike_t,bins=post)[0]

    fr = spikecount.sum(axis=1)/post.max()
    valid_idx = np.logical_and(stab>.5,fr>=1)
    if sum(valid_idx)<5:
        return None    
    spikecount = np.hstack((spikecount,np.zeros((spikecount.shape[0],1))))  
    spikerate = spikecount/dt
    spikes = np.transpose(spikerate)
    spikes = spikes[:,valid_idx]
    X = gaussian_filter1d(spikes, SMOOTHNESS, axis=0)

    speed = helpers.calcSpeed(posx)
    position = np.mod(posx,track_end)
    n_units = spikes.shape[1]
    speed_idx = speed>5
    trial_idx = np.in1d(trial,TRIALS)
    valid_idx = np.logical_and(speed_idx,trial_idx) #only take data from trials of interest and above speed threshold
    speed_r = speed[valid_idx]
    speed_r = speed_r[0::every_nth_time_bin]
    position_r = position[valid_idx]
    position_r = position_r[0::every_nth_time_bin]
    trial_r = trial[valid_idx]
    trial_r = trial_r[0::every_nth_time_bin]
    posbin_r = posx_bin[valid_idx]
    posbin_r = posbin_r[0::every_nth_time_bin]

    X_r = X[valid_idx,:]
    X_r = scipy.stats.zscore(X_r,axis=0)
    X_r = X_r[0::every_nth_time_bin,:]
    #theta = ((position_r - position_r.min()) / position_r.max()) * 2 * np.pi - np.pi
    theta = ((position_r) / track_end) * 2 * np.pi - np.pi

    model = CircularRegression(alpha=REGULARIZATION)
    #model = linear_model.LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial', max_iter=10000, C = 0.1)
    pos_bin = posbin_r-1
 
    train_TRIALS = TRIALS[0:6]
            
    train_idx = np.in1d(trial_r,train_TRIALS)
    
    tc = np.full((X_r.shape[1],200),np.nan)
    for ii in range(200):
        tc[:,ii] = X_r[np.logical_and(pos_bin==ii,train_idx),:].mean(axis = 0)

    df = pd.DataFrame(tc)
    df.interpolate(method='pchip', axis=1, limit=None, inplace=True)
    tc = df.values
    ff=np.dot(tc.transpose(),np.transpose(X_r))
    g=np.argmax(ff,axis=0)
        
        
        
        
    model.fit(X_r[train_idx,:], theta[train_idx])
    tmp = model.predict(X_r)
    #pred_pos.append(tmp)
    #tmp_e = theta[test_idx]-tmp
    #model.fit(X_r[train_idx,:], pos_bin[train_idx])
    #tmp = model.predict(X_r[test_idx,:])
    
    predicted_theta = tmp
    predicted_bin = g+1
   

    output=dict()
    output['cluID']=good_cells
    output['pos_bin'] = pos_bin
    output['true_pos'] = position_r
    output['true_theta'] = theta
    output['predicted_theta'] = predicted_theta
    output['predicted_bin'] = predicted_bin
    output['region'] = group[np.logical_and(idx,sp['cgs']==2)]
    output['trial'] = trial_r
    return output

if __name__=='__main__':
    neuropix_folder = os.path.join('/Volumes','Samsung_T5','attialex','NP_DATA')
    files = glob.glob(os.path.join(neuropix_folder,'*.mat'))
        #files = glob.glob('/oak/stanford/groups/giocomo/attialex/NP_DATA/np*_gain*.mat'
    path = '/Volumes/Samsung_T5/attialex/python_circular_gain'
    TRIALS = np.arange(5,21)

    if not os.path.exists(path):
        os.makedirs(path)

    
    for iF in files:
        session_name = os.path.split(iF)[-1]
        print(session_name)
        if 'mismatch' in session_name or 'playback' in session_name or 'dark' in session_name:
            print('skipping {}'.format(session_name))
            continue
        data = lm.loadmat(iF)
        try:
            ons = get_gain_onsets(data,0.8,100)
        except:
            ons = []
        

        
        for nbr,iO in enumerate(ons):
            trials = iO+np.arange(-5,4)
            output = run_for_file_gain(data,trials)
            sn = session_name[0:-4]
            session_name = '{}_{}.mat'.format(sn,nbr+1)
        
            if output is not None:
                plt.subplot(211)
                fig = plt.plot(output['true_theta'])
                fig = plt.plot(output['predicted_theta'])
                plt.subplot(212)
                plt.plot(output['pos_bin'])
                plt.plot(output['predicted_bin'])
                plt.savefig(os.path.join(path,'{}.png'.format(session_name[0:-4])))
                plt.close()
                scipy.io.savemat(os.path.join(path,'{}'.format(session_name)),output)
            else:
                print('skipped {} because of too few units'.format(session_name))


'''
if __name__=='__main__':
    neuropix_folder = os.path.join('/Volumes','Samsung_T5','attialex','NP_DATA')
    files = glob.glob(os.path.join(neuropix_folder,'*.mat'))
        #files = glob.glob('/oak/stanford/groups/giocomo/attialex/NP_DATA/np*_gain*.mat'
    path = '/Volumes/Samsung_T5/attialex/python_lr'
    TRIALS = np.arange(5,21)

    if not os.path.exists(path):
        os.makedirs(path)
    
    for iF in files:
        session_name = os.path.split(iF)[-1]
        print(session_name)
        if 'mismatch' in session_name or 'playback' in session_name or 'dark' in session_name:
            print('skipping {}'.format(session_name))
            continue
        data = lm.loadmat(iF)

        output = run_for_file(data,TRIALS)

        
        if output is not None:
            plt.subplot(211)
            fig = plt.plot(output['true_theta'])
            fig = plt.plot(output['predicted_theta'])
            plt.subplot(212)
            plt.plot(output['pos_bin'])
            plt.plot(output['predicted_bin'])
            plt.savefig(os.path.join(path,'{}.png'.format(session_name[0:-4])))
            plt.close()
            scipy.io.savemat(os.path.join(path,'{}'.format(session_name)),output)
        else:
            print('skipped {} because of too few units'.format(session_name))
            '''