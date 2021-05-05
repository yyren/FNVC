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
            
parser.add_argument('--caller', type=str, default="gatk", 
            help="variant caller to be selected: gatk,bam,mf, default: gatk")
            
parser.add_argument('--tag', type=str, default="0", 
            help="0 for fp and 1 for tp, default: 0")
            
args = parser.parse_args()

vcf_file=args.in_file
out_file=args.out_file
tag=args.tag
caller=args.caller



### define functions ####
class MagicDict(dict):
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

# def readData(file_path, caller):
    # data = pd.read_csv(file_path, header=None, sep=' ', low_memory=False)
    # if caller == 'gatk' or caller == 'mf':
        # return data.iloc[:, 6:].values
    # else:
        # return data.iloc[:, 11:-2].values



#### convert the vcf format to FVC record format #############
print(vcf_file+':\n')
FNVC_feature_file=change_format.make_tensor(vcf_file,caller,out_file,'NA','NA',tag)
