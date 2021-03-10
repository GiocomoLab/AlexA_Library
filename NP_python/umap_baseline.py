import umap
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.stats 
from scipy.ndimage import gaussian_filter1d
import sys
sys.path.append('./helpers')
import loadmat as lm
import helpers as helpers
from sklearn import linear_model
# %matplotlib widget
from sklearn.decomposition import PCA
import glob
from sklearn.cluster import DBSCAN,KMeans

def speedFilterFR(X,speed,speed_threshold = 0):
    X_out = X[speed>speed_threshold]
    speed_out = speed[speed>speed_threshold]
    return X_out,speed_out

def preprocess(X,speed,ds_factor = 5,gauss_win = 15,speed_threshold = 0):
    X,speed_ds = speedFilterFR(X,speed,speed_threshold = speed_threshold)
    X=gaussian_filter1d(X,gauss_win,axis=0)
    X_z = scipy.stats.zscore(X[::ds_factor,:],axis=1)
    speed_ds = speed_ds[::ds_factor]
    return X_z,speed_ds


if __name__=='__main__':

    files = glob.glob('/Volumes/T7/attialex/NP_DATA_corrected/*.mat')
    im_path ='/Volumes/T7/attialex/umap_baseline'
    if not os.path.isdir(im_path):
        os.makedirs(im_path)
    for fi in files[0:2]:
        try:
            data = lm.loadmat(fi)
            _,sn = os.path.split(fi)
        

            opt = helpers.options()
            good_cells = data['sp']['cids'][data['sp']['cgs']==2]
            counts,spMapN,stab=helpers.calculateFiringRateMap(data,good_cells=good_cells,trials2extract = None,ops=opt)
            
            spikes_gain,_,_ = helpers.calculateFiringRate(data,good_cells = good_cells,t_edges = data['post'])
            speed_gain = helpers.calcSpeed(data['posx'])

            FR = spikes_gain.mean(axis=0)
            fr_idx = FR>0.1
            spMapN=spMapN[fr_idx]
            (X_gain,speed_ds)=preprocess(spikes_gain[:,fr_idx],speed_gain,speed_threshold = 2)

            reducer = umap.UMAP(n_components=2)
            X_um=reducer.fit_transform(X_gain)

            pos,_=speedFilterFR(data['posx'],speed_gain,speed_threshold=2)
            tri,_=speedFilterFR(data['trial'],speed_gain,speed_threshold=2)
            gain,_=speedFilterFR(data['trial_gain'][data['trial']-1],speed_gain,speed_threshold = 2)
            fig,ax = plt.subplots(1,4,figsize=(15,5))
            ax[0].scatter(X_um[:,0][::2],X_um[:,1][::2],c=tri[::5][::2],s=3)
            ax[0].set_title('trial')
            ax[1].scatter(X_um[:,0][::2],X_um[:,1][::2],c=pos[::5][::2],s=3)
            ax[1].set_title('pos')
            ax[2].scatter(X_um[:,0][::2],X_um[:,1][::2],c=gain[::5][::2],s=3,cmap='Set1')
            ax[2].set_title('gain')

            sim_all = np.zeros((spMapN.shape[0],spMapN.shape[2],spMapN.shape[2]))
            for iC in range(spMapN.shape[0]):
                sim_all[iC]= np.corrcoef(spMapN[iC,:,:],rowvar=False)
            ax[3].imshow(np.nanmean(sim_all,axis=0),vmin=0,vmax=0.6)
            sn = os.path.join(im_path,sn)
            sn=sn.replace('.mat','.png')
            fig.savefig(sn)
            fig.clear()
            except:
                print('did not work for '+sn)
        

