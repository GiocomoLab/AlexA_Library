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
from helpers import preprocess

def eval_and_train(dataset):
    bl_trials = np.nonzero(np.all([dataset['trial_contrast']==100, dataset['trial_gain']==1],axis = 0))
    trialidx = np.in1d(dataset['trial_resampled'],bl_trials)
    scoring = {'prec_macro': 'precision_macro','mae_macro': make_scorer(mean_absolute_error)}
    model = linear_model.LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial', max_iter=10000, C = 0.1)
    #scores = cross_val_score(model, , cv=5)
    scores = cross_validate(model,dataset['spikecount'][trialidx,:],dataset['posx_bin'][trialidx],scoring=scoring, cv=5,return_estimator=True)
    
    model = linear_model.LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial', max_iter=10000, C = 0.1)
    model.fit(dataset['spikecount'][trialidx,:],dataset['posx_bin'][trialidx])
    return model,scores

def score_baseline_model(model,test_data):
    ma_errors = []
    m_errors = []
    precision = []
    conf_matrix = []
    for gain in [1, 0.8, 0.7, 0.6, 0.5]:
        

        bl_trials = np.nonzero(np.all([test_data['trial_contrast']==gain, test_data['trial_gain']==1],axis = 0))
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
        
        
    return (ma_errors,m_errors,precision,conf_matrix)


files = glob.glob('/oak/stanford/groups/giocomo/attialex/NP_DATA/np*_gain_*.mat')
gain_scores = []
baseline_scores = []
for iF in files[0:3]:
    if not 'baseline' in iF:
        print('Now working on '+ iF)
        dataset = lm.loadmat(iF)
        dataset = preprocess(dataset)
        (model, bl_scores) = eval_and_train(dataset)
        tmp_scores = score_baseline_model(model,dataset)
        gain_scores.append(tmp_scores)
        baseline_scores.append(bl_scores)