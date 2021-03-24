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

def cluster_plotXCorrs(X_umap,peak_locs,xcorrs):
    db = DBSCAN(min_samples=5,leaf_size=10).fit(X_umap)
    sidx = np.argsort(db.labels_)
    (fig,axs)=plt.subplots(1,3,figsize=(12,5))
    axs=axs.flatten()
    axs[0].scatter(X_umap[:,0],X_umap[:,1],marker='.',s=4,c=db.labels_,cmap='Set2')
    axs[1].imshow(xcorrs[sidx],aspect='auto')
    axs[1].scatter(peak_locs[sidx],np.arange(len(sidx)),s=4,c='r')
    labs = np.unique(db.labels_)
    sep=[]
    m_xcorrs = np.zeros((len(labs),xcorrs.shape[1]))
    mean_pwd = np.zeros((len(labs,)))
    for ii,iL in enumerate(labs):
        iidx = db.labels_==iL
        mu = X_umap[iidx].mean(axis=0)
        axs[0].plot(mu[0],mu[1],marker='x',c='k')
        mean_pwd[ii]=1-scipy.spatial.distance.pdist(xcorrs[db.labels_==iL],metric='correlation').mean()
        lab = '{}: {:.2f}'.format(iL,mean_pwd[ii])
        axs[0].annotate(lab,(mu[0],mu[1]),(mu[0]+.5,mu[1]+.5),size=10)
        sep.append(iidx.sum())
        m_xcorrs[ii]=xcorrs[iidx].mean(axis=0)
        #axs[2].plot(xcorrs[iidx].mean(axis=0))
    axs[2].imshow(m_xcorrs,aspect='auto',vmax=.7)
    locs = np.cumsum(sep)
    axs[1].hlines(locs[0:-1],xmin = 0,xmax = xcorrs.shape[1]-1,colors='r')
    
    return db.labels_,fig,mean_pwd

def cluster_plotXCorrsKM(X_umap,peak_locs,xcorrs):
    clu_list= np.arange(1,23)
    wcss=calculate_wcss(xcorrs,clu_list)
    nclu,d=optimal_number_of_clusters(wcss,clu_list)
    db = KMeans(n_clusters=nclu).fit(X_umap)
    sidx = np.argsort(db.labels_)
    (fig,axs)=plt.subplots(1,3,figsize=(12,5))
    axs=axs.flatten()
    axs[0].scatter(X_umap[:,0],X_umap[:,1],marker='.',s=4,c=db.labels_,cmap='Set2')
    axs[1].imshow(xcorrs[sidx],aspect='auto')
    axs[1].scatter(peak_locs[sidx],np.arange(np.sum(idx)),s=4,c='r')
    labs = np.unique(db.labels_)
    sep=[]
    m_xcorrs = np.zeros((len(labs),xcorrs.shape[1]))
    for ii,iL in enumerate(labs):
        iidx = db.labels_==iL
        mu = X_umap[iidx].mean(axis=0)
        axs[0].plot(mu[0],mu[1],marker='x',c='k')
        axs[0].annotate(iL,(mu[0],mu[1]),(mu[0]+.5,mu[1]+.5),size=12)
        sep.append(iidx.sum())
        m_xcorrs[ii]=xcorrs[iidx].mean(axis=0)
        #axs[2].plot(xcorrs[iidx].mean(axis=0))
    axs[2].imshow(m_xcorrs,aspect='auto',vmax=.7)
    locs = np.cumsum(sep)
    axs[1].hlines(locs[0:-1],xmin = 0,xmax = xcorrs.shape[1]-1,colors='r')
    return db.labels_,fig

def calculate_wcss(xcorrs,clu_list):
    wcss = []
    for n in clu_list:
        kmeans = KMeans(n_clusters=n)
        kmeans.fit(X=xcorrs)
        wcss.append(kmeans.inertia_)

    return wcss


def optimal_number_of_clusters(wcss,clu_list):
    x1, y1 = clu_list[0], wcss[0]
    x2, y2 = clu_list[-1], wcss[-1]

    distances = []
    for i in range(len(wcss)):
        x0 = clu_list[i]
        y0 = wcss[i]
        numerator = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
        denominator = np.sqrt((y2 - y1)**2 + (x2 - x1)**2)
        distances.append(numerator/denominator)
    
    return clu_list[distances.index(max(distances))],distances

def speedFilterFR(X,speed,speed_threshold = 0):
    X_out = X[speed>speed_threshold]
    speed_out = speed[speed>speed_threshold]
    return X_out,speed_out

def preprocess(X,speed,ds_factor = 5,gauss_win = 15,speed_threshold = 0):
    X,speed_ds = speedFilterFR(X,speed,speed_threshold = speed_threshold)
    X=gaussian_filter1d(X,gauss_win,axis=0)
    fr = X.mean(axis=0)
    X=X[:,fr>0.1]
    X_z = scipy.stats.zscore(X,axis=0)
    X_z = X_z[::ds_factor,:]
    speed_ds = speed_ds[::ds_factor]
    return X_z,speed_ds

def runUMAPForCluster(good_cells,data,ds_factor=5):
    spikes,X,t_edges = helpers.calculateFiringRate(data,good_cells = good_cells,t_edges = data['post'])
    speed = helpers.calcSpeed(data['posx'])
    (X,speed_ds)=preprocess(spikes,speed,ds_factor=ds_factor,gauss_win=10,speed_threshold = 2)
    pca = PCA(n_components=6)
    # import pdb
    # pdb.set_trace()
    X[np.logical_not(np.isfinite(X))]=0
    X_new=pca.fit_transform(X)
    #X_new = X
    #reducer = umap.UMAP(n_components=3,metric='cosine',init='spectral',n_neighbors=1000,min_dist=0.8)
    #reducer = umap.UMAP(n_components=3,metric='cosine',init='spectral',min_dist=0.8)
    #reducer = umap.UMAP(n_components=3,metric='cosine',init='spectral',min_dist=0.8)
    reducer = umap.UMAP(n_components = 3)
    #reducer = umap.UMAP(n_components=3,metric='cosine',min_dist = 1)
    Xu = reducer.fit_transform(X_new)
    return Xu,X_new[:,0]
    
def plotResults(Xu,trial,posx,speed,speed_threshold=2,ds_factor=5):
    fig,ax = plt.subplots(1,3,figsize=(12,5))

    tri,_=speedFilterFR(trial,speed,speed_threshold=2)

    ax[0].scatter(Xu[:,0][::2],Xu[:,1][::2],c=tri[::ds_factor][::2],s=3)
    ax[0].set_title('trial')

    pos,_=speedFilterFR(posx,speed,speed_threshold=2)

    ax[1].scatter(Xu[:,0][::2],Xu[:,1][::2],c=pos[::ds_factor][::2],s=3)
    ax[1].set_title('pos')
    ax[2].scatter(Xu[:,0],Xu[:,1],marker='.',s=4,c='gray')
    plot_idx = np.arange(800,1000)
    ax[2].plot(Xu[plot_idx,0],Xu[plot_idx,1],'r') #
    ax[2].scatter(Xu[plot_idx,0],Xu[plot_idx,1],c=speed_ds[2800:3000],cmap='plasma') #
    return fig



if __name__=='__main__':
    root = '/Users/attialex/distance_tuning'
    umap_version = 'Vanilla'
    files = glob.glob(os.path.join('/Users/attialex/distance_tuning','*.mat'))
    umap_save_path = os.path.join(root,umap_version)
    shutil.copy2('/Users/attialex/code/AlexA_Library/NP_python/distance_tuning_clustering.py',umap_save_path)

    if not os.path.isdir(umap_save_path):
        os.makedirs(umap_save_path)

    for fi in files:
        print(fi)
        data_out = lm.loadmat(fi)
        data_out = data_out['data_out']
        idx = data_out['pvals']<0.05
        if sum(idx)<30:
            continue

        _,sn_darkData=os.path.split(fi)
        data_path = '/Volumes/T7/attialex/NP_DATA_corrected'
        data = lm.loadmat(os.path.join(data_path,sn_darkData))

        xcorrs = data_out['xcorrs'][idx]
        
        reducer = umap.UMAP(n_components=2)
        #reducer = PCA(n_components=2)
        X_new = reducer.fit_transform(xcorrs)

        labels,fig,mean_pwd = cluster_plotXCorrs(X_new,data_out['peak_loc_all'][idx]/5,xcorrs)
        sn = fi.replace('.mat','_DBSCAN.png')
        fig.savefig(sn)

        # dbScan,fig = cluster_plotXCorrsKM(X_new,data_out['peak_loc_all'][idx]/5,xcorrs)
        # sn = fi.replace('.mat','_KMeans.png')
        # fig.savefig(sn)
        good_cells = data_out['good_cells'][idx]
        summary = []
        ds_factor = 10
        for iClu,cluID in enumerate(np.unique(labels)):
            n=np.sum(labels==cluID)
            pwd_this = mean_pwd[iClu]
            if n>=10 and pwd_this>0.88:
                good_cells_this = good_cells[labels==cluID]
                (Xu,X_pca)=runUMAPForCluster(good_cells_this,data,ds_factor=ds_factor)
                #fig=plotResults(Xu,data['trial'],data['posx'],speed)
                _,sn=os.path.split(fi)
                sn_new = sn.replace('.mat','_clu{}.png'.format(cluID))
                savepath = os.path.join('/Volumes/T7/attialex/umap_dark',sn_new)
                #fig.savefig(savepath)
                #plt.close(fig)
                summary.append((Xu,cluID,X_pca))
            else:
                print('skipping clu {}, n: {},pwd: {:.2f}'.format(cluID,n,pwd_this))
        try:
            if len(summary)>0:
                fig,ax=plt.subplots(len(summary),2,figsize=(5,12))
                for iS,(Xu,CluID,X_pca) in enumerate(summary):
                    ax[iS,0].scatter(Xu[:,0],Xu[:,1],marker='.',s=4,c='gray')
                    plot_idx = np.arange(np.round(10)/0.02/ds_factor,np.round(15)/0.02/ds_factor)
                    # import pdb
                    # pdb.set_trace() 
                    plot_idx = plot_idx.astype('int')
                    ax[iS,0].plot(Xu[plot_idx,0],Xu[plot_idx,1],'r') #
                    plot_idx = np.arange(Xu.shape[0]-plot_idx[-1],Xu.shape[0])
                    plot_idx = plot_idx.astype('int')
                    ax[iS,0].plot(Xu[plot_idx,0],Xu[plot_idx,1],'c') #
                    ax[iS,1].imshow(xcorrs[labels==CluID],aspect='auto')

                

                _,sn = os.path.split(fi)
                sn = sn.replace('.mat','_UMAPSummary.png')
                fig.savefig(os.path.join(umap_save_path,sn))

                fig = plt.figure(figsize=(12,12))
                for iS,(Xu,CluID,X_pca) in enumerate(summary):
                    
                    ax = fig.add_subplot(int(np.ceil(len(summary)/3)), 3, iS+1, projection='3d')
                    #import pdb
                    #pdb.set_trace()
                    ax.scatter(Xu[:,0],Xu[:,1],Xu[:,2],marker='.',s=1,c=X_pca,alpha=.3)
                    ax.grid(False)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.set_zticks([])
                
                _,sn = os.path.split(fi)
                sn = sn.replace('.mat','_UMAPSummary_3D.png')
                fig.savefig(os.path.join(umap_save_path,sn))

        except:
            print('error!!')
            pass

        plt.close('all')
        
