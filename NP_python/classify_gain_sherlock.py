import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import sys
sys.path.append('./helpers')
import loadmat as lm
from sklearn import linear_model
from scipy import signal
from sklearn.model_selection import cross_val_score, cross_validate
import os
from sklearn.metrics import make_scorer, confusion_matrix, mean_absolute_error
import glob
from helpers import preprocess, calcSpeed
from multiprocessing import Pool
#import multiprocessing
#from joblib import Parallel, delayed
# ml python/3.6.1
# ml py-scipystack/1.0_py36

def eval_and_train(dataset,trialrange=[]):
    if len(trialrange)>0:
        bl_trials = trialrange
    else:
        bl_trials = np.nonzero(np.all([dataset['trial_contrast']==100, dataset['trial_gain']==1],axis = 0))
    trialidx = np.in1d(dataset['trial_resampled'],bl_trials)
    speed = dataset['speed_resampled']
    speedidx = speed>2
    trialidx = np.logical_and(trialidx,speedidx)
    #scoring = {'prec_macro': 'precision_macro','mae_macro': make_scorer(mean_absolute_error)}
    #model = linear_model.LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial', max_iter=10000, C = 0.1)
    #scores = cross_val_score(model, , cv=5)
    #scores = cross_validate(model,dataset['spikecount'][trialidx,:],dataset['posx_bin'][trialidx],scoring=scoring, cv=5,return_estimator=True)
    scores = np.nan
    model = linear_model.LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial', max_iter=10000, C = 0.1)
    model.fit(dataset['spikecount'][trialidx,:],dataset['posx_bin'][trialidx])
    return model,scores

def score_baseline_model(model,test_data):
    ma_errors = []
    m_errors = []
    precision = []
    conf_matrix = []
    for gain in [1, 0.8, 0.7, 0.6, 0.5]:
        

        bl_trials = np.nonzero(np.all([test_data['trial_contrast']==100, test_data['trial_gain']==gain],axis = 0))
        trialidx = np.in1d(test_data['trial_resampled'],bl_trials)
        if trialidx.sum()>0:
            #print(model.score(test_data['spikecount'][trialidx,:],test_data['posx_bin'][trialidx]))
            Yhat = model.predict(test_data['spikecount'][trialidx,:])
            ma_errors.append(mean_absolute_error(Yhat,test_data['posx_bin'][trialidx]))
            precision.append(model.score(test_data['spikecount'][trialidx,:],test_data['posx_bin'][trialidx]))
            err = Yhat-test_data['posx_bin'][trialidx]
            m_errors.append(err.mean())
            c_matrix=confusion_matrix(test_data['posx_bin'][trialidx],Yhat)
            conf_matrix.append(c_matrix)
        else: 
            ma_errors.append(np.nan)
            m_errors.append(np.nan)
            precision.append(np.nan)
            conf_matrix.append(np.nan)
        
        
    return (ma_errors,m_errors,precision,conf_matrix,Yhat)

def score_gain_model(model,test_data):
    ma_errors = []
    m_errors = []
    precision = []
    conf_matrix = []  
    posx_centers = test_data['posx_centers']
    trialidx = np.in1d(test_data['trial_resampled'],np.arange(21,26))
    if trialidx.sum()>0:
        #print(model.score(test_data['spikecount'][trialidx,:],test_data['posx_bin'][trialidx]))
        Yhat = model.predict(test_data['spikecount'][trialidx,:])
        #ma_errors.append(mean_absolute_error(Yhat,test_data['posx_bin'][trialidx]))
        #precision.append(model.score(test_data['spikecount'][trialidx,:],test_data['posx_bin'][trialidx]))
        #err = Yhat-test_data['posx_bin'][trialidx]
        #m_errors.append(err.mean())
        c_matrix=confusion_matrix(test_data['posx_bin'][trialidx],Yhat)
        Ypred = Yhat
        Ytrue = test_data['posx_resampled'][trialidx]
        speed = test_data['speed_resampled'][trialidx]
        trial = test_data['trial_resampled'][trialidx]
        #conf_matrix.append(c_matrix)
    else: 
        ma_errors.append(np.nan)
        m_errors.append(np.nan)
        precision.append(np.nan)
        conf_matrix.append(np.nan)
        
        
    return (Ypred,Ytrue,speed,trial,c_matrix)

def run_pipeline(iF):
    try:
        
        print('Now working on '+ iF)
        dataset = lm.loadmat(iF)
        dataset = preprocess(dataset)
        if 'anatomy' not in dataset.keys():
            return
        else:
            anatomy = dataset['anatomy']
            if 'parent_shifted' in anatomy:
                group = anatomy['parent_shifted']
            else:
                group = anatomy['cluster_parent']
        region = 'MEC'
        idx = [region in ss for ss in group]
        idx = np.array(idx)
        idx = idx[dataset['sp']['cgs']==2]

        if idx.sum()==0:
            return
        
        dataset['spikecount']=dataset['spikecount'][:,idx]

        (model, bl_scores) = eval_and_train(dataset)
        (Ypred,Ytrue,speed,trial,c_matrix) = score_gain_model(model,dataset)
        plt.plot(Ytrue)

        plt.plot(dataset['posx_centers'][Ypred-1])
        name = os.path.basename(iF)[0:-4]
        plt.savefig('F:\\temp\\classifier_out\\'+region +'_'+ name + '.png')
        plt.close()
        tmp_array = np.array([Ypred,Ytrue,speed,trial,dataset['posx_edges']])
        np.save('F:\\temp\\classifier_out\\'+region +'_'+ name + '_scores.npy',tmp_array)
        #np.save('/oak/stanford/groups/giocomo/attialex/processed_data/classifier_output1/'+region +'_'+ name + '_scores.npy',tmp_array)
        #np.save('/oak/stanford/groups/giocomo/attialex/processed_data/classifier_output1/'+region +'_'+ name + '_confMatrix.npy',conf_matrix)
    except Exception as e:
        print(str(e))
        print('not working')
        pass

if __name__=='__main__':

    files = glob.glob('F:/NP_DATA/np*_gain*.mat')
    #files = glob.glob('/oak/stanford/groups/giocomo/attialex/NP_DATA/np*_gain*.mat')
    gain_scores = []
    baseline_scores = []
    for iF in files:
        run_pipeline(iF)    
    #Parallel(n_jobs=1)(delayed(run_pipeline)(i) for i in files)
    #p=Pool(processes =1)
    #p.map(run_pipeline,files)
    #p.close()
    #p.join()
