import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import os
# Helper function for generating shifted data
from affinewarp import SpikeData, ShiftWarping
from affinewarp import visualization as vis
import glob


def load_data(fn):
    data = sio.loadmat(matfile)
    count_vec=data['count_vec']
    trial_vec = data['trial_vec']
    trial_vec_random = data['trial_vec_random']
    good_cells = np.squeeze(data['good_cells'])
    theta_response = np.squeeze(data['theta'])
    return count_vec,trial_vec,trial_vec_random,good_cells,theta_response

def getDataset(trial_vec, tmin = -2, tmax=3):
    return SpikeData(trials = np.uint16(trial_vec[:,2]),spiketimes = trial_vec[:,0],neurons = np.uint16(trial_vec[:,1]),tmin = tmin, tmax = tmax)

def getSubSet(dataset,good_cells,sortidx,nCells=12,tmin=0., tmax = 1.):
    subset = dataset.select_neurons(np.sort(good_cells[sortidx[range(nCells)]])).squeeze_neurons()
    subset_time = subset.crop_spiketimes(tmin,tmax)
    nbins = subset_time.tmax-subset_time.tmin
    nbins=int(nbins*50)
    binned = subset_time.bin_spikes(nbins)
    return subset,subset_time,binned

def fitModel(binned,cells2fit = [0]):
    model = ShiftWarping(maxlag=.15, smoothness_reg_scale=10.)
    data = binned[:,:,cells2fit]
    print(data.shape)
    model.fit(data, iterations=20)
    return model

def plotFitted(model,binned):
    fig,axes = plt.subplots(binned.shape[2],2,figsize=(15, 5))
    transformed = model.transform(binned)
    print(axes.shape)
    if binned.shape[2]>1:
        for ii in range(binned.shape[2]):
            axes[ii,0].imshow(binned[:,:,ii].squeeze(),aspect='auto')
            axes[ii,1].imshow(transformed[:,:,ii].squeeze(),aspect='auto')

    else:
        axes[0].imshow(binned[:,:,0].squeeze(),aspect='auto')
        axes[1].imshow(transformed[:,:,0].squeeze(),aspect='auto')

    return fig

    #fig.savefig('/Users/attialex/testfig.png')
    
    #plt.show()
def plotTuningCurves(tc,tc_random,tc_baseline,x_vec=np.linspace(-2,3,250)):
    n_fig = tc.shape[1]
    n_cols = 3
    n_rows = int(np.ceil(n_fig/n_cols))
    fig,axes = plt.subplots(n_rows,n_cols,figsize=(15, 5),gridspec_kw={'hspace': 0},sharex=True, sharey=False)
    axes = axes.flatten()
    for iF in range(n_fig):
        l=axes[iF].axvline(0)
        l.set_color('k')
        l.set_linewidth(1)
        l.set_linestyle(':')
        axes[iF].plot(x_vec,tc[:,iF],label='R fit')
        axes[iF].plot(x_vec,tc_random[:,iF],label='Random')
        axes[iF].plot(x_vec,tc_baseline[:,iF],label = 'R orig',color='k',linewidth=1)
        axes[iF].set_xlim(x_vec.min(),x_vec.max())
    axes[-1].legend()
    


    return fig

def runPipeline(dataset,good_cells,sortidx,imprefix='MM_MMasc'):
    nCells = 12
    subset,subset_time,binned = getSubSet(dataset,good_cells,sortidx,nCells=nCells,tmin=0.1, tmax = 1.1)
    figs = dict()
    cells2fit = np.arange(2)
    order = np.arange(subset.n_neurons)
    order[np.argsort(good_cells[sortidx[0:nCells]])]=order
    cells2fit=order[0:2]
    mod = fitModel(binned,cells2fit=cells2fit)
    fig_fitted = plotFitted(mod,binned[:,:,cells2fit])
    figs[imprefix +'fitted_neurons']=fig_fitted
    ordered = subset.shift_each_trial_by_constant(mod.fractional_shifts) #this works as long as fit period is 1s
    (fig_rast,ax) = vis.rasters(ordered,order=order,subplots=(6,2))
    figs[imprefix+'rasterplot']=fig_rast
    nBins = (subset.tmax-subset.tmin)*50
    binned = ordered.bin_spikes(nBins)
    tc = binned.mean(axis=0)
    tc = tc[:,order]
    tc_bl = subset.bin_spikes(nBins).mean(axis=0)
    tc_bl = tc_bl[:,order]
    return figs,tc,tc_bl

def save_figs(fig_dict,im_save_dir):
    for k in fig_dict.keys():
        fig = fig_dict[k]
        name = k +'.png'
        path = os.path.join(im_save_dir,name)
        fig.savefig(path)




if __name__== '__main__':
    #matfile = r'/Users/attialex/temp/AA_200122_2_200206_mismatch_2.mat'
    matfiles=glob.glob(os.path.join(r'/Users/attialex/temp/','*.mat'))
    for matfile in matfiles:
        sn = os.path.splitext(os.path.split(matfile)[1])[0]
        imdir = os.path.join('/Users/attialex/images/warp_fit/MM',sn)
        #print(sn)
        os.makedirs(imdir,exist_ok=True)
        count_vec,trial_vec,trial_vec_random,good_cells,theta_response= load_data(matfile)

        resp = count_vec[:,101:150].mean(axis=1)-count_vec[:,75:100].mean(axis=1)
        method = 'desc'
        if method == 'asc':
            sortidx=np.argsort(resp)
        elif method == 'desc':
            sortidx=np.argsort(resp)[::-1]
            sortidx = np.argsort(-resp)
        elif method == 'theta':
            sortidx = np.argsort(theta_response)[::-1]


        
        dataset = getDataset(trial_vec)
        (figs,tc,tc_baseline)=runPipeline(dataset,good_cells,sortidx,'MM_MM'+method+'_')
        dataset_random = getDataset(trial_vec_random)
        (figsR,tcR,tc_rbaseline)=runPipeline(dataset_random,good_cells,sortidx,'MMR_MM'+method+'_')
        tc_fig = dict()
        tc_fig['MM'+method+'_tuning']=plotTuningCurves(tc,tcR,tc_baseline)
        #save_figs({**figs,**figsR,**tc_fig},imdir)
        save_figs(tc_fig,imdir)
        plt.close('all')
        #plt.show()
        #fig_rast.savefig(os.path.join(im_dir,sn,'raster.png'))

