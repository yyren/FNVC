# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 21:26:27 2019

@author: yyren
"""


from sklearn.metrics import confusion_matrix, accuracy_score, matthews_corrcoef, roc_auc_score, roc_curve,precision_recall_curve, auc
from scikitplot.metrics import plot_roc, plot_calibration_curve
from sklearn.externals import joblib
from sklearn import metrics
from pylab import mpl
from matplotlib.pyplot import plot,savefig
import numpy as np
import pandas as pd
import xgboost as xgb
import argparse
import math
import pylab as pl
import re
import os

#parser the argument
parser = argparse.ArgumentParser(
  description='extract the complex region' )

parser.add_argument('--test_TP_file', type=str, default="testing_tp_file.record", 
            help="path to the tensor file contains the true positive variants, default:testing_tp_file.record")

parser.add_argument('--test_FP_file', type=str, default="testing_fp_file.record", 
            help="path to the tensor file contains the false positive variants, default:testing_fp_file.record")
            
parser.add_argument('--out_file', type=str, default="garfield_results_statistic.txt", 
            help="the statistic results file, default: garfield_results_statistic.txt")
            
parser.add_argument('--sample_tag', type=str, default="HG001", 
            help="sample name in the 10th col of vcf, default: HG001")
parser.add_argument('--model_type', type=str, default="snv", 
            help="snv/indel, default:snv")
args = parser.parse_args()

tp_garanno_out=args.test_TP_file
fp_garanno_out=args.test_FP_file
out_file=args.out_file
sample_tag=args.sample_tag
model_type=args.model_type

cutoff=0.8369
if model_type == 'indel':
    cutoff=0.6301

# read training data
def readData(file,sample,tag):
    gar_tp_pred=[]
    gar_tp_y=[]
    with open(file) as f:
        for row in f.readlines():
            if not (re.findall(r'^#',row)):
                recs = row.strip().split('\t')
                infors = recs[7].split(';')
                gar_value=0
                for rec in infors:
                    if sample in rec:
                        gar_value=rec.split('=')[1]
                gar_tp_pred.append(gar_value)
                gar_tp_y.append(tag)
    f.close()
    gar_tp_pred=np.array(gar_tp_pred).astype(np.float64)
    gar_tp_y=np.array(gar_tp_y).astype(np.int32)
    return gar_tp_y, gar_tp_pred


y_tp_tag, ypred_tp = readData(tp_garanno_out,sample_tag,1)
y_fp_tag, ypred_fp = readData(fp_garanno_out,sample_tag,0)
ypred=np.concatenate((ypred_tp,ypred_fp),axis=0)
y_test=np.concatenate((y_tp_tag,y_fp_tag),axis=0)

#max ACC and F1 in illumina indel(0.6301) and snp(0.8369)
y_pred = (ypred >= cutoff)*1
mcc_test = matthews_corrcoef(y_test, y_pred)
confusion_data=metrics.confusion_matrix(y_test, y_pred)
tn=int(confusion_data[0][0])
fp=int(confusion_data[0][1])
fn=int(confusion_data[1][0])
tp=int(confusion_data[1][1])
fnr=fn/(fn+tp)
fpr=fp/(tn+fp)
fdr=fp/(tp+fp)
acc_info=fn/tn

precision, recall, thresholds = precision_recall_curve(y_test, ypred)
f=open(out_file,'w')
f.write(str('AUC: %.4f' % roc_auc_score(y_test,ypred))+'\n\n')
f.write(str('AUPRC: %.6f' % auc(recall,precision))+'\n\n')
f.write(str('max_acc_MCC: %.4f' % mcc_test)+'\n')
f.write(str('max_acc_ACC: %.4f' % metrics.accuracy_score(y_test,y_pred))+'\n')
f.write(str('max_acc_Recall: %.4f' % metrics.recall_score(y_test,y_pred))+'\n')
f.write(str('max_acc_F1-score: %.4f' %metrics.f1_score(y_test,y_pred))+'\n')
f.write(str('max_acc_Precesion: %.4f' %metrics.precision_score(y_test,y_pred))+'\n')
f.write(str('max_acc_FNR: %.4f' % fnr)+'\n')
f.write(str('max_acc_FPR: %.4f' % fpr)+'\n')
f.write(str('max_acc_FDR: %.4f' % fdr)+'\n')
f.write(str('acc_info: %.4f' % acc_info)+'\n')
f.write(str('max_acc_compound matrix:')+'\n'+str(metrics.confusion_matrix(y_test,y_pred))+'\n\n')


#suggest value in github
y_pred = (ypred >= cutoff)*1
mcc_test = matthews_corrcoef(y_test, y_pred)
confusion_data=metrics.confusion_matrix(y_test, y_pred)
tn=int(confusion_data[0][0])
fp=int(confusion_data[0][1])
fn=int(confusion_data[1][0])
tp=int(confusion_data[1][1])
fnr=fn/(fn+tp)
fpr=fp/(tn+fp)
fdr=fp/(tp+fp)
acc_info=fn/tn

#summary_file='summary.txt'
f.write(str('suggest_MCC: %.4f' % mcc_test)+'\n')
f.write(str('suggest_ACC: %.4f' % metrics.accuracy_score(y_test,y_pred))+'\n')
f.write(str('suggest_Recall: %.4f' % metrics.recall_score(y_test,y_pred))+'\n')
f.write(str('suggest_F1-score: %.4f' %metrics.f1_score(y_test,y_pred))+'\n')
f.write(str('suggest_Precesion: %.4f' %metrics.precision_score(y_test,y_pred))+'\n')
f.write(str('suggest_FNR: %.6f' % fnr)+'\n')
f.write(str('suggest_FPR: %.6f' % fpr)+'\n')
f.write(str('suggest_FDR: %.6f' % fdr)+'\n')
f.write(str('acc_info: %.4f' % acc_info)+'\n')
f.write(str('suggest_compound matrix:')+'\n'+str(metrics.confusion_matrix(y_test,y_pred))+'\n')
f.close()

