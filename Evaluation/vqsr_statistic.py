# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 21:26:27 2019

@author: yyren
"""


from sklearn.metrics import confusion_matrix, accuracy_score, matthews_corrcoef, roc_auc_score, roc_curve,precision_recall_curve, auc
#from scikitplot.metrics import plot_roc, plot_calibration_curve
from sklearn.externals import joblib
from sklearn import metrics
from pylab import mpl
from matplotlib.pyplot import plot,savefig
import numpy as np
import pandas as pd
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

parser.add_argument('--out_file', type=str, default="vqsr_statistic_results.txt", 
            help="statistic result files, default:vqsr_statistic_results.txt")
            
args = parser.parse_args()

test_TP_file=args.test_TP_file
test_FP_file=args.test_FP_file
out_file=args.out_file

vqsr_score_key='VQSLOD'
# read training data
def readData(file,tag,count):
    y_score=[]
    y_tag=[]
    y_pred=[]
    right_numb=0
    wrong_numb=0
    idx=1
    with open(file) as f:
        for row in f.readlines():
            if (not (re.findall(r'^#',row))):
                recs = row.strip().split('\t')
                if recs[6] == 'PASS':
                    y_pred.append(1)
                else:
                    y_pred.append(0)
                infors = recs[7].split(';')
                target_value=0
                for rec in infors:
                    if vqsr_score_key in rec:
                        target_value=rec.split('=')[1]
                y_score.append(target_value)
                y_tag.append(tag)
                idx +=1
    f.close()
    y_tag=np.array(y_tag).astype(np.int32)
    y_score=np.array(y_score).astype(np.float64)
    y_pred=np.array(y_pred).astype(np.int32)
    return y_tag, y_score, y_pred

def get_count(file):
    idx=0
    with open(file) as f:
        for row in f.readlines():
            if not (re.findall(r'^#',row)):
                idx +=1
    f.close()
    return int(idx)

sample_number=get_count(test_FP_file)

y_tp_tag, ypred_tp, y_pred_tp = readData(test_TP_file,1,sample_number)
y_fp_tag, ypred_fp, y_pred_fp = readData(test_FP_file,0,sample_number)
ypred=np.concatenate((ypred_tp,ypred_fp),axis=0)
y_test=np.concatenate((y_tp_tag,y_fp_tag),axis=0)
y_pred=np.concatenate((y_pred_tp,y_pred_fp),axis=0)

#max ACC and F1 in illumina snp
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
f.write(str('AUC: %.4f' % roc_auc_score(y_test,ypred))+'\n')
f.write(str('AUPRC: %.6f' % auc(recall,precision))+'\n\n')
f.write(str('acc_info: %.4f' % acc_info)+'\n')
f.write(str('max_f1_MCC: %.4f' % mcc_test)+'\n')
f.write(str('max_f1_F1-score: %.4f' %metrics.f1_score(y_test,y_pred))+'\n')
f.write(str('max_f1_ACC: %.4f' % metrics.accuracy_score(y_test,y_pred))+'\n')
f.write(str('max_f1_FDR: %.4f' % fdr)+'\n')
f.write(str('max_f1_FNR: %.4f' % fnr)+'\n')
f.write(str('max_f1_Recall: %.4f' % metrics.recall_score(y_test,y_pred))+'\n')
f.write(str('max_f1_Precesion: %.4f' %metrics.precision_score(y_test,y_pred))+'\n')
f.write(str('max_f1_FPR: %.4f' % fpr)+'\n')
f.write(str('max_f1_compound matrix:')+'\n'+str(metrics.confusion_matrix(y_test,y_pred))+'\n\n')


f.close()

