import argparse
import numpy as np
import collections
import math
import re
import sys


class MagicDict(dict):
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

result_record=MagicDict()
chrom2num = {"X":'23',"Y":'24',"M":'25',"MT":'25'}
def sort_vcf(args):
    in_file=args.in_file
    out_file=args.out_file
    out_records=open(out_file,'w+')
    with open(in_file) as f:
        for row in f.readlines():
            if re.findall(r'^#',row):
                out_records.write(row)
            else:
                records = row.strip().split()
                chrom=re.sub('chr','',records[0])#chr1->1
                if chrom in ('X','Y','M','MT'):
                    chrom=chrom2num[chrom]
                if chrom in ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25'):
                    chrom=int(chrom)
                    start=int(records[1])
                    result_record[chrom][start][records[3]][records[4]]=row
    f.close()
    for chrom in sorted(result_record.keys()):
        #print(str(chrom)+'\n')
        for start in sorted(result_record[chrom].keys()):
            for ref in sorted(result_record[chrom][start].keys()):
                for alt in sorted(result_record[chrom][start][ref].keys()):
                    out_records.write(result_record[chrom][start][ref][alt])
    out_records.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
      description='sort the vcf file' )
      
    parser.add_argument('--in_file', type=str, default="input.vcf", 
            help="path to the input file, default:input.vcf")
            
    parser.add_argument('--out_file', type=str, default="out_file.vcf", 
            help="path to the output file, default:out_file.vcf")
    args = parser.parse_args()
    sort_vcf(args)
    