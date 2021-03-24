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



if __name__=='__main__':

    #files = glob.glob('/Volumes/T7/attialex/NP_DATA_corrected/*.mat')
    #im_path ='/Volumes/T7/attialex/umap_baseline'
    im_path = r'F:\attialex\umap_gain05AverageSpatialMap_V1_v1'
    files = glob.glob(r'F:\attialex\NP_DATA_corrected\AA*.mat')
    if not os.path.isdir(im_path):
        os.makedirs(im_path)
    

    shutil.copy2(os.path.abspath(__file__),im_path)
    ds_factor = 5
    for fi in files:
        try:
            data = lm.loadmat(fi)
            gain_val = 0.5
            values = (data['trial_gain']==gain_val) & (data['trial_contrast']==100)
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
            regions = ('VISp')
            idx = [str(ss).startswith(regions) for ss in group]
            region_idx = np.array(idx)

            _,sn = os.path.split(fi)
            # import pdb
            # pdb.set_trace()
            good_cells = data['sp']['cids'][(data['sp']['cgs']==2) & region_idx]
            if len(good_cells)<5:
                continue

            opt = helpers.options()
            counts,spMapN,stab=helpers.calculateFiringRateMap(data,good_cells=good_cells,trials2extract = trial_range,ops=opt)
            
            
            tc = spMapN.mean(axis=2).T
            tcPre = spMapN[:,:,0:6].mean(axis = 2).T
            tcGain = spMapN[:,:,6:10].mean(axis=2).T
            tcPost = spMapN[:,:,10::].mean(axis=2).T
            fr = tc.mean(axis=0)
            tcComb = np.vstack((tcPre,tcGain,tcPost))
            # import pdb
            # pdb.set_trace()
            tcZ = scipy.stats.zscore(tcComb[:,fr>0.1],axis=0)
            tcZ = tcComb[:,fr>0.1]

            reducer = PCA(n_components=3)

      
            pos =  np.vstack((np.arange(200),np.arange(200),np.arange(200)))
            gain = np.vstack((np.ones((200,1)),np.ones((200,1))*gain_val,np.ones((200,1))))
            tri = np.vstack((np.ones((200,1)),np.ones((200,1))*2,np.ones((200,1))*3))

           
            #reducer = umap.UMAP(n_components=2,metric='cosine')
            X_um = reducer.fit_transform(tcZ)
            #reducer  = PCA(n_components=3)
            #X_um = reducer.fit_transform(Xz)

            twoD=False
            # import pdb
            # pdb.set_trace()
            if twoD:
                fig,ax = plt.subplots(1,4,figsize=(15,5))
            else:
                fig = plt.figure(figsize=(15,5))
                ax = [plt.subplot(1, 4, i+1, projection='3d') for i in range(3)]
                ax.append(plt.subplot(1,4,4))
            # import pdb
            # pdb.set_trace()
            
            try:
                
                if twoD:
                    ax[0].scatter(X_um[:,0],X_um[:,1],c=tri,s=3)
                    ax[0].set_title('trial')
                else:
                    # ax=fig.add_subplot(1,4,1,projection='3d')
                    # ax.scatter(X_um[:,0],X_um[:,1],X_um[:,2],c=tri_tmp,s=3)
                    # ax.set_title('trial')
                    ax[0].scatter(X_um[:,0],X_um[:,1],X_um[:,2],c=tri,s=3)
                    ax[0].set_title('trial')
               
            except:
                pass
            try:
                if twoD:
                    ax[1].scatter(X_um[:,0],X_um[:,1],c=pos,s=3)
                else:
                    ax[1].scatter(X_um[:,0],X_um[:,1],X_um[:,2],c=pos,s=3)

                ax[1].set_title('pos')
            except:
                pass

            try:
                
                gain_tmp = gain
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
        except ZeroDivisionError:
            print(traceback.format_exc())
            print('did not work for '+sn)
            

