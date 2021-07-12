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
from distance_tuning_clustering import cluster_plotXCorrs,speedFilterFR,preprocess

def runForFile(good_cells,sn_this,labels,umap_save_path,good_cells_orig,xcorr=None):
    data = lm.loadmat(sn_this)
    summary = []
    ds_factor = 5
    for iClu,cluID in enumerate(np.unique(labels)):
        n=np.sum(labels==cluID)
        pwd_this = mean_pwd[iClu]
        if n>=10 and pwd_this>0.88:
            if xcorr is not None:
                xcorr_this = xcorr[labels==cluID]
            else:
                xcorr_this = None
            good_cells_this = good_cells[labels==cluID]
            (Xu,X_pca)=runUMAPForCluster(good_cells_this,data,ds_factor=ds_factor)
            #fig=plotResults(Xu,data['trial'],data['posx'],speed)
            _,sn=os.path.split(fi)
            sn_new = sn.replace('.mat','_clu{}.png'.format(cluID))
            savepath = os.path.join('/Volumes/T7/attialex/umap_dark',sn_new)
            #fig.savefig(savepath)
            #plt.close(fig)
            summary.append((Xu,cluID,X_pca,xcorr_this))
        else:
            print('skipping clu {}, n: {},pwd: {:.2f}'.format(cluID,n,pwd_this))
    # import pdb
    # pdb.set_trace()
    (Xu,X_pca)=runUMAPForCluster(good_cells_orig,data,ds_factor=ds_factor) # run once for all cells as sanity check
    summary.append((Xu,cluID,X_pca,None))
    # import pdb
    # pdb.set_trace()
    if len(summary)>0:        
        fig = plotSummary(summary,data,ds_factor)
        _,sn = os.path.split(sn_this)
        sn = sn.replace('.mat','_UMAPSummary.png')
#         import pdb
#         pdb.set_trace()
        fig.savefig(os.path.join(umap_save_path,sn))
    return    

def plotSummary(summary_this,data,ds_factor):
    fig,axs = plt.subplots(len(summary_this),4,figsize=(12,12))
    speed = helpers.calcSpeed(data['posx'])
    speed_ds,_ = speedFilterFR(speed,speed,speed_threshold = 2)
    for ii,(X_um,cluID,Xpc,xcorr) in enumerate(summary_this):
        if len(summary_this)>1:
            ax = axs[ii]
        else:
            ax = axs #for when there is only 1 cluster
        try:
            tri,_=speedFilterFR(data['trial'],speed,speed_threshold=2)
            # import pdb
            # pdb.set_trace()
            ax[0].scatter(X_um[:,0][::2],X_um[:,1][::2],c=tri[::ds_factor][::2],s=3,marker='.')
            ax[0].set_title('trial')
        except:
            pass
        try:
            pos,_=speedFilterFR(data['posx'],speed,speed_threshold=2)

            ax[1].scatter(X_um[:,0][::2],X_um[:,1][::2],c=pos[::ds_factor][::2],s=3,marker='.')
            ax[1].set_title('pos')
        except:
            pass

        try:
            gain,_=speedFilterFR(data['trial_gain'][data['trial']-1],speed,speed_threshold = 2)

            sc=ax[2].scatter(X_um[:,0][::2],X_um[:,1][::2],c=gain[::ds_factor][::2],s=3,cmap='Set1',marker='.')
            ax[2].set_title('gain')
            # import pdb 
            # pdb.set_trace()
            handles,labels = sc.legend_elements(prop='colors',alpha = 0.6)
            ax[2].legend(handles,labels,loc='best')
        except:
            pass
        try:
            ax[3].imshow(xcorr,aspect='auto')
        except:
            pass
        for ax_this in ax:
            ax_this.set_xticks([])
            ax_this.set_yticks([])
    return fig

def runUMAPForCluster(good_cells,data,ds_factor=5):
    spikes,X,t_edges = helpers.calculateFiringRate(data,good_cells = good_cells,t_edges = data['post'])
    speed = helpers.calcSpeed(data['posx'])
    FR = spikes.mean(axis=0)
    fr_idx = FR>0.1
    (X,speed_ds)=preprocess(spikes[:,fr_idx],speed,ds_factor=ds_factor,gauss_win=10,speed_threshold = 2)
    pca = PCA(n_components=6)
    # import pdb
    # pdb.set_trace()
    X[np.logical_not(np.isfinite(X))]=0
    X_new=pca.fit_transform(X)
    #X_new = X
    #reducer = umap.UMAP(n_components=3,metric='cosine',init='spectral',n_neighbors=1000,min_dist=0.8)
    #reducer = umap.UMAP(n_components=3,metric='cosine',init='spectral',min_dist=0.8)
    #reducer = umap.UMAP(n_components=3,metric='cosine',init='spectral',min_dist=0.8)
    #reducer = umap.UMAP(n_components = 3)
    reducer = umap.UMAP(n_components=2)
    Xu = reducer.fit_transform(X_new)
    return Xu,X_new[:,0]




if __name__=='__main__':
    root = '/Users/attialex/distance_tuning'
    root =r'C:\Users\attialex\Documents\distance_tuning'
    data_dir = r'F:\attialex\NP_DATA_corrected'
    umap_version = 'Vanilla_otherData_pcaumap'
    files = glob.glob(os.path.join(root,'*.mat'))
    umap_save_path = os.path.join(root,umap_version)
    

    if not os.path.isdir(umap_save_path):
        os.makedirs(umap_save_path)
    # import pdb
    # pdb.set_trace()
    shutil.copy2(os.path.abspath(__file__),umap_save_path)

    for fi in files:
        print(fi)
        data_out = lm.loadmat(fi)
        data_out = data_out['data_out']
        idx = data_out['pvals']<0.05
        if sum(idx)<30:
            continue

        _,sn_darkData=os.path.split(fi)
    #     data = lm.loadmat(os.path.join(data_path,sn_darkData))

        xcorrs = data_out['xcorrs'][idx]

        reducer = umap.UMAP(n_components=2)
        #reducer = PCA(n_components=2)
        X_new = reducer.fit_transform(xcorrs)

        labels,fig,mean_pwd = cluster_plotXCorrs(X_new,data_out['peak_loc_all'][idx]/5,xcorrs)
        sn = fi.replace('.mat','_DBSCAN.png')
        fig.savefig(sn)
        good_cells = data_out['good_cells'][idx]
        
        elements = sn_darkData.split('_')
        stem = elements[0]+'_'+elements[1]+'*.mat'
        pot_files = glob.glob(os.path.join(data_dir,stem))
        
        for sn_this in pot_files:
            runForFile(good_cells,sn_this,labels,umap_save_path,data_out['good_cells'],xcorrs)

        
