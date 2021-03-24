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
import shutil
import traceback

def speedFilterFR(X,speed,speed_threshold = 0):
    X_out = X[speed>speed_threshold]
    speed_out = speed[speed>speed_threshold]
    return X_out,speed_out

def preprocess(X,speed,ds_factor = 5,gauss_win = 15,speed_threshold = 0):
    X,speed_ds = speedFilterFR(X,speed,speed_threshold = speed_threshold)
    X=gaussian_filter1d(X,gauss_win,axis=0)
    #X_z = scipy.stats.zscore(X[::ds_factor,:],axis=0)
    X_z = X[::ds_factor,:]
    speed_ds = speed_ds[::ds_factor]
    return X_z,speed_ds



def runUMAPForCluster(good_cells,data,trials,ds_factor=5,speed_threshold = 2):
    spikes,X,t_edges = helpers.calculateFiringRate(data,good_cells = good_cells,t_edges = data['post'])
    speed = helpers.calcSpeed(data['posx'])
    trialIDX = np.in1d(data['trial'],trials)
    FR = spikes.mean(axis=0)
    fr_idx = FR>0.1
    spikes = spikes[trialIDX]
    speed = speed[trialIDX]
    (X,speed_ds)=preprocess(spikes[:,fr_idx],speed,ds_factor=ds_factor,gauss_win=10,speed_threshold = speed_threshold)
    pca = PCA(n_components=3)
    # import pdb
    # pdb.set_trace()
    X[np.logical_not(np.isfinite(X))]=0
    #X_new=pca.fit_transform(X)
    X_new = X
    #reducer = umap.UMAP(n_components=3,metric='cosine',init='spectral',n_neighbors=1000,min_dist=0.8)
    #reducer = umap.UMAP(n_components=3,metric='cosine',init='spectral',min_dist=0.8)
    #reducer = umap.UMAP(n_components=3,metric='cosine',init='spectral',min_dist=0.8)
    #reducer = umap.UMAP(n_components = 3)
    reducer = umap.UMAP(n_components=2,metric = 'cosine',n_neighbors=15,min_dist = 0.0)
    Xu = reducer.fit_transform(X_new)
    #Xu = pca.fit_transform(X_new)
    return Xu,trialIDX

if __name__=='__main__':

    #files = glob.glob('/Volumes/T7/attialex/NP_DATA_corrected/*.mat')
    #im_path ='/Volumes/T7/attialex/umap_baseline'
    im_path = r'F:\attialex\umap_gain08_RSC_v10'
    files = glob.glob(r'F:\attialex\NP_DATA_corrected\AA*.mat')
    if not os.path.isdir(im_path):
        os.makedirs(im_path)
    speed_threshold = 2

    shutil.copy2(os.path.abspath(__file__),im_path)
    ds_factor = 5
    for fi in files:
        try:
            data = lm.loadmat(fi)

            values = (data['trial_gain']==0.8) & (data['trial_contrast']==100)
            matches =  (np.logical_not(values[:-1])) & (values[1:])
            onsets = np.where(matches)[0] +1
            if len(onsets)==0:
                continue
            trial_range = onsets[0]+np.arange(-5,11)

            try:
                anatomy = data['anatomy']
            except:
                print('no anatomy')
                continue
            
            if 'parent_shifted' in anatomy:
                group = anatomy['parent_shifted']
            else:
                group = anatomy['cluster_parent']
            #regions = ('MEC','VISp','RS')
            regions = ('RS')
            idx = [str(ss).startswith(regions) for ss in group]
            idxagl = [str(ss).startswith('RSPagl') for ss in group]

            region_idx = np.array(idx) & np.logical_not(np.array(idxagl))
            print(region_idx.sum())
            #print(region_idx.shape)
            if region_idx.sum()<5:
                continue

            _,sn = os.path.split(fi)
            # import pdb
            # pdb.set_trace()
            good_cells = data['sp']['cids'][(data['sp']['cgs']==2) & region_idx]
            opt = helpers.options()
            # _,_,stabPre=helpers.calculateFiringRateMap(data,good_cells=good_cells,trials2extract = trial_range[0:6],ops=opt)
            # good_cells = good_cells[stabPre>0.5]
            X_um,trial_idx=runUMAPForCluster(good_cells,data,trial_range,ds_factor=ds_factor,speed_threshold=speed_threshold)

            
            counts,spMapN,stab=helpers.calculateFiringRateMap(data,good_cells=good_cells,trials2extract = trial_range,ops=opt)
            
            speed = helpers.calcSpeed(data['posx'])
            twoD=True

            if twoD:
                fig,ax = plt.subplots(1,4,figsize=(15,5))
            else:
                fig = plt.figure(figsize=(15,5))
                ax = [plt.subplot(1, 4, i+1, projection='3d') for i in range(3)]
                ax.append(plt.subplot(1,4,4))
            # import pdb
            # pdb.set_trace()
            
            try:
                tri,_=speedFilterFR(data['trial'][trial_idx],speed[trial_idx],speed_threshold=speed_threshold)
                tri_tmp = tri[::ds_factor]-tri.min()
                if twoD:
                    ax[0].scatter(X_um[:,0],X_um[:,1],c=tri_tmp,s=3)
                    ax[0].set_title('trial')
                else:
                    # ax=fig.add_subplot(1,4,1,projection='3d')
                    # ax.scatter(X_um[:,0],X_um[:,1],X_um[:,2],c=tri_tmp,s=3)
                    # ax.set_title('trial')
                    ax[0].scatter(X_um[:,0],X_um[:,1],X_um[:,2],c=tri_tmp,s=3)
                    ax[0].set_title('trial')
               
            except:
                pass
            try:
                pos,_=speedFilterFR(data['posx'][trial_idx],speed[trial_idx],speed_threshold=speed_threshold)
                if twoD:
                    ax[1].scatter(X_um[:,0],X_um[:,1],c=pos[::ds_factor],s=3)
                else:
                    ax[1].scatter(X_um[:,0],X_um[:,1],X_um[:,2],c=pos[::ds_factor],s=3)

                ax[1].set_title('pos')
            except:
                pass

            try:
                
                gain,_=speedFilterFR(data['trial_gain'][data['trial']-1][trial_idx],speed[trial_idx],speed_threshold=speed_threshold)
                gain_tmp = gain[::ds_factor]
                gain_tmp[-1]=0.8
                gain_tmp[-2]=0.7
                gain_tmp[-3]=0.6
                gain_tmp[-4]=0.5
                if twoD:
                    sc=ax[2].scatter(X_um[:,0],X_um[:,1],c=gain_tmp,s=3,cmap='Set2')
                else:
                    sc=ax[2].scatter(X_um[:,0],X_um[:,1],X_um[:,2],c=gain_tmp,s=3,cmap='Set2')
                ax[2].set_title('gain')
                # import pdb 
                # pdb.set_trace()
                handles,labels = sc.legend_elements(prop='colors',alpha = 0.6)
                ax[2].legend(handles,labels,loc='best')
            except:
                pass
            sim_all = np.zeros((spMapN.shape[0],spMapN.shape[2],spMapN.shape[2]))
            for iC in range(spMapN.shape[0]):
                sim_all[iC]= np.corrcoef(spMapN[iC,:,:],rowvar=False)
            neg=ax[3].imshow(np.nanmean(sim_all,axis=0),vmin=0,vmax=0.6)
            ax[3].set_title('similarity')
            
            fig.colorbar(neg, ax=ax[3])
            sn = os.path.join(im_path,sn)
            sn=sn.replace('.mat','.png')
            fig.savefig(sn)
            fig.clear()
        except:
            print(traceback.format_exc())
            print('did not work for '+sn)
        

