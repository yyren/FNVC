# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 21:26:27 2019

@author: yyren
"""


from sklearn.metrics import confusion_matrix, accuracy_score, matthews_corrcoef, roc_auc_score, roc_curve, precision_recall_curve, auc
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
import pylab as pl
import seaborn as sns
import matplotlib.pyplot as plt

#parser the argument
parser = argparse.ArgumentParser(
  description='extract the complex region' )

parser.add_argument('--test_TP_file', type=str, default="testing_tp_file.record", 
            help="path to the tensor file contains the true positive variants, default:testing_tp_file.record")

parser.add_argument('--test_FP_file', type=str, default="testing_fp_file.record", 
            help="path to the tensor file contains the false positive variants, default:testing_fp_file.record")

parser.add_argument('--in_model', type=str, default="FNVC_gatk_model.bin", 
            help="the output training model name with path, default:FNVC_gatk_model.bin")

parser.add_argument('--out_file', type=str, default="FNVC_gatk_summary.txt", 
            help="statistic result file, default:FNVC_gatk_summary.txt")
            
parser.add_argument('--acc_cutoff', type=float, default="0.5", 
            help="regard as true if larger than the value, default:0.5")
            
parser.add_argument('--mcc_cutoff', type=str, default="0.5", 
            help="regard as true if larger than the value, default:0.5")
args = parser.parse_args()

test_TP_file=args.test_TP_file
test_FP_file=args.test_FP_file
in_model=args.in_model
out_file=args.out_file
acc_cutoff=float(args.acc_cutoff)
f1_cutoff=float(args.mcc_cutoff)

def readData(file_path1, file_path2):
    # Read file and combine
    data1 = pd.read_csv(file_path1, header=None, sep=' ', low_memory=False)
    data2 = pd.read_csv(file_path2, header=None, sep=' ', low_memory=False)
    
    combined = pd.concat([data1, data2])
    # Return X, y
    return combined.iloc[:, 6:].values, combined.iloc[:, 0].values

# read training data
X_test, y_test = readData(test_TP_file,test_FP_file)
dtest=xgb.DMatrix(X_test)

bst = xgb.Booster(model_file=in_model)
ypred=bst.predict(dtest)
y_pred = (ypred >= acc_cutoff)*1
mcc_test = matthews_corrcoef(y_test, y_pred)
confusion_data=metrics.confusion_matrix(y_test, y_pred)
tn=int(confusion_data[0][0])
fp=int(confusion_data[0][1])
fn=int(confusion_data[1][0])
tp=int(confusion_data[1][1])
fnr=fn/(fn+tp)
fpr=fp/(tn+fp)
fdr=fp/(fp+tp)
info=tn/fn

precision, recall, thresholds = precision_recall_curve(y_test, ypred)
# fig = plt.figure(figsize=(15,5))
# ax = fig.add_subplot(1,2,1)
# sns.heatmap(confusion_data,cmap='coolwarm_r',linewidths=0.5,annot=True,ax=ax,fmt ='d')
# plt.title('Confusion Matrix')
# plt.ylabel('Real Classes')
# plt.xlabel('Predicted Classes')
# savefig('/lustre/home/acct-clslh/clslh/renyongyong/project/paper_data/script/get_tensor/script/gatk_feature/confusion.png')
#summary_file='summary.txt'
f=open(out_file,'w')
f.write('max_ACC max_MCC'+'\n')
f.write(str(acc_cutoff)+' '+str(f1_cutoff)+'\n')
f.write(str('AUC: %.6f' % roc_auc_score(y_test,ypred))+'\n')
f.write(str('AUPRC: %.6f' % auc(recall,precision))+'\n\n')

f.write(str('max_acc_MCC: %.5f' % mcc_test)+'\n')
f.write(str('max_acc_F1-score: %.5f' %metrics.f1_score(y_test,y_pred))+'\n')
f.write(str('max_acc_ACC: %.5f' % metrics.accuracy_score(y_test,y_pred))+'\n')
f.write(str('max_acc_FDR: %.5f' % fdr)+'\n')
f.write(str('max_acc_FNR: %.5f' % fnr)+'\n')
f.write(str('max_acc_info: %.5f' % info)+'\n')
f.write(str('max_acc_Recall: %.5f' % metrics.recall_score(y_test,y_pred))+'\n')
f.write(str('max_acc_Precesion: %.5f' %metrics.precision_score(y_test,y_pred))+'\n')
f.write(str('max_acc_FPR: %.5f' % fpr)+'\n')

f.write(str('max_acc_compound matrix:')+'\n'+str(metrics.confusion_matrix(y_test,y_pred))+'\n\n')

#max MCC
y_pred = (ypred >= f1_cutoff)*1
mcc_test = matthews_corrcoef(y_test, y_pred)
confusion_data=metrics.confusion_matrix(y_test, y_pred)
tn=int(confusion_data[0][0])
fp=int(confusion_data[0][1])
fn=int(confusion_data[1][0])
tp=int(confusion_data[1][1])
fnr=fn/(fn+tp)
fpr=fp/(tn+fp)
fdr=fp/(fp+tp)
info=tn/fn

#summary_file='summary.txt'
f.write(str('max_mcc_MCC: %.5f' % mcc_test)+'\n')
f.write(str('max_mcc_F1-score: %.5f' %metrics.f1_score(y_test,y_pred))+'\n')
f.write(str('max_mcc_ACC: %.5f' % metrics.accuracy_score(y_test,y_pred))+'\n')
f.write(str('max_mcc_FDR: %.5f' % fdr)+'\n')
f.write(str('max_mcc_FNR: %.5f' % fnr)+'\n')
f.write(str('max_mcc_info: %.5f' % info)+'\n')
f.write(str('max_mcc_Recall: %.5f' % metrics.recall_score(y_test,y_pred))+'\n')
f.write(str('max_mcc_Precesion: %.5f' %metrics.precision_score(y_test,y_pred))+'\n')
f.write(str('max_mcc_FPR: %.5f' % fpr)+'\n')

f.write(str('max_mcc_compound matrix:')+'\n'+str(metrics.confusion_matrix(y_test,y_pred))+'\n')
f.close()