# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 21:26:27 2019

@author: hill103
"""

'''
Using logistic regression
'''

from sklearn.metrics import confusion_matrix, accuracy_score, matthews_corrcoef, roc_auc_score, roc_curve
from scikitplot.metrics import plot_roc, plot_calibration_curve
from sklearn.externals import joblib
from sklearn import metrics
from pylab import mpl
from matplotlib.pyplot import plot,savefig
import numpy as np
import pandas as pd
import xgboost as xgb
import os
import argparse
import math
#parser the argument
parser = argparse.ArgumentParser(
  description='extract the complex region' )

parser.add_argument('--train_TP_file', type=str, default="training_tp_tensor.record", 
            help="path to the tensor file contains the variant information, default:training_tp_tensor.record")

parser.add_argument('--train_FP_file', type=str, default="training_fp_tensor.record", 
            help="path to the tensor file contains the variant information, default:training_fp_tensor.record")

parser.add_argument('--test_TP_file', type=str, default="testing_tp_file.record", 
            help="path to the tensor file contains the variant information, default:testing_tp_file.record")

parser.add_argument('--test_FP_file', type=str, default="testing_fp_file.record", 
            help="path to the tensor file contains the variant information, default:testing_fp_file.record")

parser.add_argument('--out_model', type=str, default="RandomForestClassifier.pkl", 
            help="the output training model name with path, default:./RandomForestClassifier.pkl")

parser.add_argument('--out_path', type=str, default="./", 
            help="path to the tensor file contains the variant information, default:./")
args = parser.parse_args()

train_TP_file=args.train_TP_file
train_FP_file=args.train_FP_file
test_TP_file=args.test_TP_file
test_FP_file=args.test_FP_file
out_model=args.out_model
out_path=args.out_path



def readData(file_path1, file_path2):
    # Read file and combine
    data1 = pd.read_csv(file_path1, header=None, sep=' ', low_memory=False)
    data2 = pd.read_csv(file_path2, header=None, sep=' ', low_memory=False)
    
    combined = pd.concat([data1, data2])
    # Return X, y
    return combined.iloc[:, 6:].values, combined.iloc[:, 0].values

# read training data
X_train, y_train = readData(train_TP_file,train_FP_file)
X_test, y_test = readData(test_TP_file,test_FP_file)
dtrain=xgb.DMatrix(X_train,label=y_train)
dtest=xgb.DMatrix(X_test)

#booster:
params={'booster':'gbtree',
        'objective': 'binary:logistic',
        'eval_metric': 'auc',
        'max_depth':8,
        'n_estimators':500,
        'lambda':15,
        'subsample':0.75,
        'colsample_bytree':0.75,
        'min_child_weight':1,
        'eta': 0.025,
        'seed':0,
        'silent':0,
        'gamma':0.15,
        'nthread':12,
        'learning_rate' : 0.01}

watchlist = [(dtrain,'train')]

bst=xgb.train(params,dtrain,num_boost_round=500, evals=watchlist)
bst.save_model(out_model)
ypred=bst.predict(dtest)
y_pred = (ypred >= 0.5)*1
mcc_test = matthews_corrcoef(y_test, y_pred)
confusion_data=metrics.confusion_matrix(y_test, y_pred)
tn=int(confusion_data[0][0])
fp=int(confusion_data[0][1])
fn=int(confusion_data[1][0])
tp=int(confusion_data[1][1])
fnr=fn/(fn+tp)
fpr=fp/(tn+fp)
fdr=fp/(fp+tp)
tnr=tn/(tn+fp)

summary_file='xgboost_merge_gatk_our_summary.txt'
f=open(os.path.join(out_path,summary_file),'w')
f.write(str('AUC: %.4f' % roc_auc_score(y_test,ypred))+'\n')
f.write(str('MCC: %.4f' % mcc_test)+'\n')
f.write(str('F1-score: %.4f' %metrics.f1_score(y_test,y_pred))+'\n')
f.write(str('ACC: %.4f' % metrics.accuracy_score(y_test,y_pred))+'\n')
f.write(str('FDR: %.4f' % fdr)+'\n')
f.write(str('FNR: %.4f' % fnr)+'\n')
f.write(str('FPR: %.4f' % fpr)+'\n')
f.write(str('TNR: %.4f' % tnr)+'\n')
f.write(str('Recall: %.4f' % metrics.recall_score(y_test,y_pred))+'\n')

f.write(str('Precesion: %.4f' %metrics.precision_score(y_test,y_pred))+'\n')
f.write(str('compound matrix:')+'\n'+str(metrics.confusion_matrix(y_test,y_pred))+'\n')
f.close()
 

mpl.rcParams['font.sans-serif'] = ['SimHei']  
xgb.plot_importance(bst)
savefig(os.path.join(out_path, 'feature_importance.png'))

