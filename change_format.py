# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 2019

@author: ryy
"""

'''
convert the vcf record to FVC tensor format
'''

import get_gatk_features
import get_bam_features
import re
import os

def merge_gatk_bam_features(gatk_tensor,bam_tensor,out_merged):
    script_path=os.path.abspath(os.path.dirname(__file__))
    script=os.path.join(script_path,'combine_features.pl')
    merge_cmd=' '.join(['perl',str(script),gatk_tensor,bam_tensor,out_merged])
    os.system(merge_cmd)

def make_tensor(vcf_file,model,variant_type,bam_file="input.bam",ref_fasta="hg19.fa",tag='0'):
    data_path=os.path.abspath(os.path.dirname(vcf_file))
    fvc_tensor=os.path.join(data_path,str(variant_type)+'_feature.record')
    tensor_seqcontext_abs_path=os.path.join(data_path,str(variant_type)+'_seqcontext.record')
    if model == 'gatk':
        get_gatk_features.out_record(vcf_file,fvc_tensor)
    elif model == 'bam':
        get_bam_features.ourt_record(vcf_file,bam_file,ref_fasta,fvc_tensor,tensor_seqcontext_abs_path)
        os.remove(tensor_seqcontext_abs_path)
    elif model == 'mf':
        gatk_feature_tensor=os.path.join(data_path,str(variant_type)+'_gatk_feature.record')
        bam_feature_tensor=os.path.join(data_path,str(variant_type)+'_bam_feature.record')
        bam_seqContext_tensor=os.path.join(data_path,str(variant_type)+'_bam_seqcontext.record')
        merged_tensor=os.path.join(data_path,str(variant_type)+'_feature.record')
        
        get_gatk_features.out_record(vcf_file,gatk_feature_tensor)
        get_bam_features.ourt_record(vcf_file,bam_file,ref_fasta,bam_feature_tensor,bam_seqContext_tensor)
        merge_gatk_bam_features(gatk_feature_tensor,bam_feature_tensor,fvc_tensor)
        os.remove(bam_seqContext_tensor)
        os.remove(gatk_feature_tensor)
        os.remove(bam_feature_tensor)
    return(str(fvc_tensor))
