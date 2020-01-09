# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 2019

@author: ryy
"""

'''
predict variant confidence
'''


import numpy as np
import pandas as pd
import xgboost as xgb
import os
import argparse
import math
import re
import change_format

#parser the argument
parser = argparse.ArgumentParser(
  description='filtering the false variants in NGS variant calling results' )

parser.add_argument('--in_file', type=str, default="input.vcf", 
            help="the vcf file need to be quality control, default:input.vcf")
            
parser.add_argument('--out_file', type=str, default="./", 
            help="the quality controled vcf file, default:./")
            
parser.add_argument('--bam_file', type=str, default="./file.bam", 
            help="required bam file for FVC-bam model, default:./file.bam")
            
parser.add_argument('--reference', type=str, default="./hg19.fa", 
            help="the reference fasta file used for alignement, default:./hg19.fa")
            
parser.add_argument('--snv_model', type=str, default="FVC_gatk_snv.bin", 
            help="the trained snv model:./FVC_gatk_snv.bin")
            
parser.add_argument('--indel_model', type=str, default="FVC_gatk_indel.bin", 
            help="the trained indel model, default:./FVC_gatk_indel.bin")
            
parser.add_argument('--snv_cutoff', type=float, default="0.5126", 
            help="regard as true snv variant if above the setted value, suggest:FVC-gatk:0.5126; FVC-bam:0.5062; FVC-mf:0.4453")
            
parser.add_argument('--indel_cutoff', type=float, default="0.5824", 
            help="regard as true indel variant if above the setted value, suggest:FVC-gatk:0.5824; FVC-bam:0.4784; FVC-mf:0.4673")
            
parser.add_argument('--model', type=str, default="gatk", 
            help="Model to be selected: gatk,bam,mf, default: gatk")
            
args = parser.parse_args()

vcf_file=args.in_file
out_file=args.out_file
bam_file=args.bam_file
reference=args.reference
snv_model=args.snv_model
indel_model=args.indel_model
snv_cutoff=float(args.snv_cutoff)
indel_cutoff=float(args.indel_cutoff)
model=args.model


### define functions ####
class MagicDict(dict):
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def readData(file_path, model):
    data = pd.read_csv(file_path, header=None, sep=' ', low_memory=False)
    if model == 'gatk' or model == 'mf':
        return data.iloc[:, 6:].values
    else:
        return data.iloc[:, 11:-2].values

def separate_snv_indel(data):
    data_path=os.path.abspath(os.path.dirname(data))
    snv_file=os.path.join(data_path,'snv_temp.vcf')
    indel_file=os.path.join(data_path,'indel_temp.vcf')
    snv_out=open(snv_file,'w')
    indel_out=open(indel_file,'w')
    with open(data) as infile:
        for row in infile.readlines():
            row=row.strip()
            if not re.findall(r'^#',row):
                recs = row.strip().split('\t')
                ref=recs[3]
                alt=recs[4]
                if re.findall(r'\*',recs[4]) or re.findall(r'\,',recs[4]):
                    alts=alt.split(",")
                    if (len(ref) != len(alts[0]) and re.match("[ATCGN]",alts[0])) or (len(ref) != len(alts[1]) and re.match("[ATCGN]",alts[1])):
                        indel_out.write(str(row)+'\n')
                    else:
                        snv_out.write(str(row)+'\n')
                else:
                    if len(ref) == len(alt):
                        snv_out.write(str(row)+'\n')
                    else:
                        indel_out.write(str(row)+'\n')
    infile.close()
    snv_out.close()
    indel_out.close()
    return(snv_file,indel_file)



#### seperate the data to SNV and INDEL ######################
snv_file,indel_file=separate_snv_indel(vcf_file)

#### convert the vcf format to FVC record format #############
snv_record=change_format.make_tensor(snv_file,model,'snv',bam_file,reference)
indel_record=change_format.make_tensor(indel_file,model,'indel',bam_file,reference)

#### load the model ##########
SNV_MODEL = xgb.Booster(model_file=snv_model)
INDEL_MODEL = xgb.Booster(model_file=indel_model)

#### load data ###############
predict_conf={}
predict_pass={}
pass_number=0
variant_types=['snv','indel']
for variant_type in variant_types:
    record_file=snv_record
    FVC_model=SNV_MODEL
    conf_cutoff=snv_cutoff
    filter_tag='FVC'+str(conf_cutoff)
    if variant_type == 'indel':
        record_file=indel_record
        FVC_model=INDEL_MODEL
        conf_cutoff=indel_cutoff
        filter_tag='FVC'+str(conf_cutoff)
    if os.path.getsize(record_file)>0:
        FVC_record= xgb.DMatrix(readData(record_file,model))
        ypred_prob=FVC_model.predict(FVC_record)
        idx=0
        with open(record_file) as infile:
            for row in infile.readlines():
                if not(re.findall(r'^#',row)):
                    recs = row.strip().split(' ')
                    #chrom=str(re.sub('chr','',recs[1]))#chr1->1
                    chrom=recs[1]
                    key='_'.join([chrom,recs[3],recs[4],recs[5]])
                    predict_conf[key]=ypred_prob[idx]
                    if float(predict_conf[key])>= conf_cutoff:
                        predict_pass[key]='PASS'
                        pass_number +=1
                    else:
                        predict_pass[key]=str(filter_tag)
                    idx +=1
        infile.close()
    os.remove(record_file)

#### output the new vcf file with predict confidence value #######
predict_out=open(out_file,'w')
with open(vcf_file) as infile:
    for row in infile.readlines():
        row=row.strip()
        if re.findall(r'^#',row):
            predict_out.write(str(row)+'\n')
        else:
            recs = row.split('\t')
            #chrom=str(re.sub('chr','',recs[0]))#chr1->1
            chrom=recs[0]
            key='_'.join([chrom,recs[1],recs[3],recs[4]])
            predict_value=('%.4f' % float(predict_conf[key]))
            recs[6]=str(recs[6])+';'+str(predict_pass[key])
            recs[7]=str(recs[7])+';FVC_score='+str(predict_value)
            sep_str="\t"
            new_line=sep_str.join(recs)
            predict_out.write(str(new_line)+'\n')
infile.close()
predict_out.close()
filtered_numb=len(predict_conf)-pass_number

#### remove the temp files #########
os.remove(snv_file)
os.remove(indel_file)

print("Total: %d\nFiltered: %d\nPass: %d\n" %(len(predict_conf),filtered_numb,pass_number))