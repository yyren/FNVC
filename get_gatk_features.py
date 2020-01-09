from __future__ import (absolute_import, division, print_function,
                        unicode_literals, with_statement)
import argparse
import numpy as np
import collections
import math
import re
import itertools
import datetime
import sys
def group(n, iterable):
    """group items to iterables of size n.
 
    the last part can have less than n elements.
 
    Args:
        n: group by this number
        iterable: any iterable
 
    """
    if n < 1:
        raise ValueError("group by N, N should be at least 1")
    one_element = []
    #for index, e in itertools.izip(itertools.cycle(range(n)), iterable):
    for index, e in zip(itertools.cycle(range(n)), iterable):
        one_element.append(e)
        if index == n - 1:
            yield one_element[:]
            one_element = []
    if one_element:
        yield one_element
 
 
def find_minimum_repeat_unit(text):
    """Find minimum repeat unit and how many times it repeats.
 
    Args:
        text: the string to test.
 
    Return:
        (unit, repeat_times)
 
    """
    l = len(text)
    for i in range(l):
        unit_length = i + 1
        if l % unit_length != 0:
            continue
        sequences = list(group(unit_length, text))
        for e in sequences[1:]:
            # print("comparing %s with %s" % (e, sequences[0]))
            if e != sequences[0]:
                break
        else:
            return "".join(sequences[0]), l // unit_length
    return text,1
    #assert False    # never reach

class MagicDict(dict):
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

base2num = dict(zip("ACGTNM", (0,1,2,3,4,4)))
# Base calling ambiguities, See https://www.bioinformatics.org/sms/iupac.html
AMBIGUITY_CODES = {
    'K':[0, 0, 0.5, 0.5], 'M':[0.5, 0.5, 0, 0], 'R':[0.5, 0, 0, 0.5], 'Y':[0, 0.5, 0.5, 0], 'S':[0, 0.5, 0, 0.5],
    'W':[0.5, 0, 0.5, 0], 'B':[0,0.333,0.333,0.334], 'V':[0.333,0.333,0,0.334], 'H':[0.333,0.333,0.334,0],
    'D':[0.333,0,0.333,0.334], 'X':[0.25,0.25,0.25,0.25], 'N':[0.25,0.25,0.25,0.25]
}
DNA_SYMBOLS = {'A':0, 'C':1, 'G':2, 'T':3}


def get_seq(chr,pos,refSeq):
    seqBases=(refSeq.fetch(reference=chr, start=pos-21, end=pos+20)).upper()
    seq_tensor=np.zeros( (41, 4) )
    for i in range(len(seqBases)):
        if seqBases[i] in DNA_SYMBOLS:
            seq_tensor[i][DNA_SYMBOLS[seqBases[i]]]=1
        elif seqBases[i] in AMBIGUITY_CODES:
            seq_tensor[i]=AMBIGUITY_CODES[seqBases[i]]
    return(seq_tensor)


def get_homolen_and_distance(refSeq,variant_chr,pos,center):
    ref=(refSeq.fetch(reference=variant_chr, start=center-21, end=center+20)).upper()
    pos_idx=20+int(pos)-int(center)
    index=0
    max_len=1
    homo_start_idx=0
    while index <=len(ref)-1:
        #print(index)
        currLen=1
        currBase=ref[index]
        next_index=index
        if index+1 < len(ref):
            for j in range(index+1,len(ref)):
                #print(j)
                next_index +=1
                if currBase==ref[j]:
                    currLen +=1
                else:
                    break
        else:
            next_index=index+1
        if currLen>=max_len:
            max_len=currLen
            homo_start_idx=index
        index =next_index
        #print(index,len(ref))
    in_homo=0
    homo_len=0
    dist_to_homo=20
    if max_len>3:
        dist_to_homo=min(abs(pos_idx-homo_start_idx),abs(pos_idx-(homo_start_idx+max_len-1)))
        if homo_start_idx <=pos_idx and (int(homo_start_idx) +int(max_len)-1>=pos_idx):
            in_homo=1
            dist_to_homo=0
        homo_len=max_len
    #print(pos_idx,homo_start_idx,int(homo_start_idx) +int(max_len)-1,homo_len,dist_to_homo)
    return(in_homo,homo_len,dist_to_homo)

def get_GC_freq(refSeq,variant_chr,pos):
    ref=(refSeq.fetch(reference=variant_chr, start=pos-21, end=pos+21)).upper()
    gc_count=0
    gc_count_around_5bp=0
    gc_count_around_10bp=0
    at_count=0
    at_count_around_5bp=0
    at_count_around_10bp=0
    for i in range(0,len(ref)):
        if str(ref[i]) == 'G' or str(ref[i]) == 'C':
            gc_count +=1
            if i >=15 and i <=25:
                gc_count_around_5bp +=1
            if i >=10 and i <=30:
                gc_count_around_10bp +=1
        elif str(ref[i]) == 'A' or str(ref[i]) == 'T':
            at_count +=1
            if i >=15 and i <=25:
                at_count_around_5bp +=1
            if i >=10 and i <=30:
                at_count_around_10bp +=1
    gc_freq=('%.2f' %(float(gc_count)/40.0))
    gc_around_5bp_freq=('%.2f' %(float(gc_count_around_5bp)/40.0))
    gc_around_10bp_freq=('%.2f' %(float(gc_count_around_10bp)/40.0))
    at_freq=('%.2f' %(float(at_count)/40.0))
    at_around_5bp_freq=('%.2f' %(float(at_count_around_5bp)/40.0))
    at_around_10bp_freq=('%.2f' %(float(at_count_around_10bp)/40.0))
    return(gc_around_5bp_freq,gc_around_10bp_freq,gc_freq,at_around_5bp_freq,at_around_10bp_freq,at_freq)


def get_repeat_unit(seq):
    unit_seq, times = find_minimum_repeat_unit(seq)
    return(unit_seq)


def get_next_cigar(cigar_index,cigar_array,query_cigar_index,cigarLength):
    nextCigar = 'other'
    if (cigar_index < len(cigar_array)-1):#not the last cigar
        if ((cigar_array[cigar_index+1] == 1) or (cigar_array[cigar_index+1] == 4)) and (query_cigar_index == cigarLength-1):#the last cigar is ins and the pos is final_base-1 of the cigar
            nextCigar = 'ins'
        else:
            nextCigar = 'other'
    return(nextCigar)

def get_center_index_base(chr,center,read_start,cigarLine):
    pos = read_start
    center_idx_in_fragment = 0
    for ( cigarType, cigarLength ) in cigarLine:
        if cigarType in (0, 7, 8):
            for i in range(cigarLength):
                if pos == center:
                    center_type = 'snv'
                    return(center_type,center_idx_in_fragment)
                    break
                center_idx_in_fragment +=1
                pos +=1
        elif cigarType in (1, 4):# soft clip or Insertion
            for i in range(cigarLength):
                if pos == center:
                    center_type = 'ins'
                    return(center_type,center_idx_in_fragment)
                    break
                center_idx_in_fragment +=1
        elif cigarType == 2:# deletion
            for i in range(cigarLength):
                if pos == center:
                    center_type = 'del'
                    return(center_type,center_idx_in_fragment)
                    break
                pos +=1
        else:
            print('Uknown')
            raise ValueError('Not include center base')

def on_unit_region(ref_seq,alt_unit_seq):
    ref_unit_seq = get_repeat_unit(ref_seq)
    if ref_unit_seq == alt_unit_seq:
        in_unit_region=1
    else:
        in_unit_region=0
    return(in_unit_region)

Infor_feature = {'MQ':0, 'DP':1, 'SOR':2, 'FS':3, 'QD':4, 'MQRankSum':5, 'QUAL':6, 'ReadPosRankSum':7, 'BaseQRankSum':8}
def get_annotations(QUAL,INFO,FORMAT,Value):
    annotation_feature=np.zeros((1,10))
    if QUAL != '.':
        annotation_feature[0][Infor_feature['QUAL']]=QUAL
    infors=INFO.strip().split(';')
    for i in range(len(infors)):
        if re.findall(r'=',infors[i]):
            infor_recs=infors[i].strip().split('=')
            if infor_recs[0] in Infor_feature:
                annotation_feature[0][Infor_feature[infor_recs[0]]]=infor_recs[1]
    formats=FORMAT.strip().split(':')
    values=Value.strip().split(':')
    for i in range(len(formats)):
        if formats[i] == 'GQ':
            annotation_feature[0][9]=values[i]
    return annotation_feature

def out_record(vcf_file,out_tensor,tag='0',max_sample_numb=9999999):
    begin2end = MagicDict()
    sample_numb=1
    annotation_data = np.zeros((1,10))
    feature = open(out_tensor,'w')
    time_stamp = datetime.datetime.now()
    print('Starting get features:',time_stamp.strftime('%Y.%m.%d-%H:%M:%S.%f')[:-3],'\n')
    with open(vcf_file) as f:
        for row in f.readlines():
            if not (re.findall(r'^#',row)):
                if sample_numb <= max_sample_numb:
                    row = row.strip().split()
                    chrom=re.sub('chr','',row[0])#chr1->1
                    #chrom = row[0]
                    pos = int(row[1])
                    alt=row[4]
                    GT=(row[9].strip().split(':'))[0]
                    #print(GT)
                    if GT == '0/0':
                        gt_tag=10
                    elif GT == '0/1' or GT == '1/0':
                        gt_tag=11
                    elif GT == '1/1':
                        gt_tag=12
                    elif GT == '2/1' or GT == '1/2':
                        gt_tag=13
                    if (re.findall(r'\*',row[4])) or (re.findall(r',',row[4])):
                        alts=alt.split(',')
                        if ((len(row[3])>len(alts[0])) and alts[0] != '*') or ((len(row[3])>len(alts[1])) and alts[1] != '*'):
                            pos = int(row[1])+1
                    else:
                        if(len(row[3])>len(row[4])):
                            pos = int(row[1])+1
                    output_line = []
                    #print(int(tag),chrom,pos,row[1],str(row[3]),str(row[4]))
                    output_line.append( "%d %s %d %d %s %s" % (int(tag),chrom,int(pos),int(row[1]),str(row[3]),str(row[4])) )
                    annotation_feature = get_annotations(row[5],row[7],row[8],row[9])
                    for record in np.reshape(annotation_feature, 1*10):
                        #output_line.append("%0.2f" % record)
                        output_line.append("%.4f" % record)
                    out_line=" ".join(output_line)
                    feature.write(out_line+'\n')
                else:
                    break
                sample_numb +=1
    f.close()
    feature.close()
    time_stamp = datetime.datetime.now()
    print('Finish get features:',time_stamp.strftime('%Y.%m.%d-%H:%M:%S.%f')[:-3],'\n')
