from __future__ import (absolute_import, division, print_function,
                        unicode_literals, with_statement)
import pysam
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
    ref=(refSeq.fetch(reference=variant_chr, start=center-11, end=center+10)).upper()
    pos_idx=10+int(pos)-int(center)
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
    dist_to_homo=1.00
    if max_len>=3:
        dist_to_homo=min(abs(pos_idx-homo_start_idx),abs(pos_idx-(homo_start_idx+max_len-1)))
        if homo_start_idx <=pos_idx and (int(homo_start_idx) +int(max_len)-1>=pos_idx):
            in_homo=1
            dist_to_homo=0
        homo_len=max_len
    #print(pos_idx,homo_start_idx,int(homo_start_idx) +int(max_len)-1,homo_len,dist_to_homo)
    if homo_len<3:
        dist_to_homo=('%.3f' %(float(1.0)))
    else:
        dist_to_homo=('%.3f' %(float(dist_to_homo)/float(10.0)))
    if homo_len <3:
        homo_len=('%.3f' %(float(0.0)))
    elif 3<=homo_len<13:
        homo_len=('%.3f' %(float(homo_len-2)/float(10.0)))
    elif homo_len>=13:
        homo_len=('%.3f' %(float(1.0)))
    
    return(in_homo,homo_len,dist_to_homo)

def get_GC_freq(refSeq,variant_chr,pos):
    ref=(refSeq.fetch(reference=variant_chr, start=pos-21, end=pos+20)).upper()
    gc_count=0
    gc_count_around_5bp=0
    gc_count_around_10bp=0
    at_count=0
    at_count_around_5bp=0
    at_count_around_10bp=0
    for i in range(0,len(ref)):
        if str(ref[i]) == 'G' or str(ref[i]) == 'C':
            gc_count +=1
            if i >=16 and i <=26:
                gc_count_around_5bp +=1
            if i >=11 and i <=31:
                gc_count_around_10bp +=1
        elif str(ref[i]) == 'A' or str(ref[i]) == 'T':
            at_count +=1
            if i >=16 and i <=26:
                at_count_around_5bp +=1
            if i >=11 and i <=31:
                at_count_around_10bp +=1
    #print(len(ref),gc_count_around_5bp)
    gc_freq=('%.3f' %(float(gc_count)/41.0))
    gc_around_5bp_freq=('%.3f' %(float(gc_count_around_5bp)/11.0))
    gc_around_10bp_freq=('%.3f' %(float(gc_count_around_10bp)/21.0))
    at_freq=('%.3f' %(float(at_count)/41.0))
    at_around_5bp_freq=('%.3f' %(float(at_count_around_5bp)/11.0))
    at_around_10bp_freq=('%.3f' %(float(at_count_around_10bp)/21.0))
    return(gc_around_5bp_freq,gc_around_10bp_freq,gc_freq,at_around_5bp_freq,at_around_10bp_freq,at_freq)


def get_repeat_unit(seq):
    unit_seq, times = find_minimum_repeat_unit(seq)
    return(unit_seq,times)


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

# def on_unit_region(ref_seq,alt_seq):
    # base=['A','T','C','G']
    # next_base=ref_seq[11]
    # last_base=ref_seq[9]
    # break_step='no'
    # units_freq=0
    # change_to_repeat=0
    # for i in range(0,4):
        # if break_step == 'yes':
            # break
        # for j in range(0,4):
            # ref_double_unit=(base[i]+base[j]) *2
            # ref_trible_unit=(base[i]+base[j]) *3
            # #print(double_unit,trible_unit)
            # exist_double=re.findall(ref_double_unit,ref_seq)
            # exist_trible=re.findall(ref_trible_unit,ref_seq)
            # if len(exist_double)>=2 or exist_trible:
                # units=re.findall(base[i]+base[j],ref_seq)
                # units_freq=float(len(units))*2.0/21.0
                # break_step='yes'
                # break
    # alt_double_unit=(alt_seq+next_base) *2
    # alt_trible_unit=(alt_seq+next_base) *3
    # alt_exist_double=re.findall(alt_double_unit,ref_seq)
    # alt_exist_trible=re.findall(alt_trible_unit,ref_seq)
    # last_alt_double_unit=(last_base+alt_seq) *2
    # last_alt_trible_unit=(last_base+alt_seq) *3
    # last_alt_exist_double=re.findall(last_alt_double_unit,ref_seq)
    # last_alt_exist_trible=re.findall(last_alt_trible_unit,ref_seq)
    # if (len(alt_exist_double)>=2 or alt_exist_trible) or (len(last_alt_exist_double)>=2 or last_alt_exist_trible):
        # change_to_repeat=1
    # return(units_freq,change_to_repeat)

def on_unit_region(ref_seq,alt_seq):
    base=['A','T','C','G']
    next_base=ref_seq[2]
    last_base=ref_seq[0]#snv or ins
    break_step='no'
    units_freq=0
    change_to_repeat=0
    unit_length=0
    #whether exist repeat with more than two bases
    unit_seq, times=find_minimum_repeat_unit(alt_seq)
    exist_unit=re.findall(unit_seq,ref_seq)
    if len(exist_unit) >=3 and len(unit_seq)>=2:
        break_step='yes'
        unit_length=len(unit_seq)
        units_freq=round(float(len(exist_unit)*unit_length)/21.0,3)
        #print(units_freq)
    #one or two bases repeat
    if break_step == 'no':
        for i in range(0,4):
            if break_step == 'yes':
                break
            for j in range(0,4):
                ref_double_unit=(base[i]+base[j]) *2
                ref_trible_unit=(base[i]+base[j]) *3
                #print(double_unit,trible_unit)
                exist_double=re.findall(ref_double_unit,ref_seq)
                exist_trible=re.findall(ref_trible_unit,ref_seq)
                if len(exist_double)>=2 or exist_trible:
                    units=re.findall(base[i]+base[j],ref_seq)
                    units_freq=round(float(len(units))*2.0/21.0,3)
                    break_step='yes'
                    if i == j:
                        unit_length=1
                    else:
                        unit_length=2
                    break
    alt_with_nextbase_double_unit=(alt_seq+next_base) *2
    alt_with_nextbase_trible_unit=(alt_seq+next_base) *3
    alt_with_nextbase_exist_double=re.findall(alt_with_nextbase_double_unit,ref_seq)
    alt_with_nextbase_exist_trible=re.findall(alt_with_nextbase_trible_unit,ref_seq)
    last_alt_double_unit=(last_base+alt_seq) *2
    last_alt_trible_unit=(last_base+alt_seq) *3
    last_alt_exist_double=re.findall(last_alt_double_unit,ref_seq)
    last_alt_exist_trible=re.findall(last_alt_trible_unit,ref_seq)
    alt_with_double=(alt_seq)*2
    alt_exist_double=re.findall(alt_with_double,ref_seq)
    alt_with_trible=(alt_seq)*3
    alt_exist_trible=re.findall(alt_with_trible,ref_seq)
    if (len(alt_with_nextbase_exist_double)>=2 or alt_with_nextbase_exist_trible) or (len(alt_exist_double)>=2 or alt_exist_trible) or (len(last_alt_exist_double)>=2 or last_alt_exist_trible):
        change_to_repeat=1
    #print('sub:',units_freq,change_to_repeat,unit_length)
    return(units_freq,change_to_repeat,unit_length)
    
    
def dist_to_end(read_length,pos):
    distance=0;
    read_len=int(read_length)
    if read_len % 2 == 1:
        mid_idx=(read_len+1)/2
        if pos<=mid_idx:
            distance=(pos-1)/(mid_idx-1)
        else:
            distance=(read_len-pos)/(mid_idx-1)
    else:
        mid_idx=(read_len)/2
        if pos<=mid_idx:
            distance=(pos-1)/(mid_idx-1)
        else:
            distance=(read_len-pos)/(mid_idx-1)
    return(distance)
        
        
def store_vcf_record(tag,vcf_file,max_sample_numb):
    begin2end = MagicDict()
    sample_numb=1
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
                    begin2end[ chrom ][ pos ][ row[3] ][ row[4] ] = (pos - 20, pos + 20, tag, row[1], row[3],row[4],gt_tag)
                else:
                    break
                sample_numb +=1
    f.close()
    return(begin2end)

def ourt_record(vcf_file,bam_file,ref_fasta,out_aln_tensor,out_refBases,max_sample_numb=100000000,tag='0'):
    max_sample_numb=int(max_sample_numb)
    time_stamp = datetime.datetime.now()
    print('Starting get variant features:',time_stamp.strftime('%Y.%m.%d-%H:%M:%S.%f')[:-3],'\n')
    begin2end = store_vcf_record(tag,vcf_file,max_sample_numb)
    time_stamp = datetime.datetime.now()
    #print('finish store the vcf information:',time_stamp.strftime('%Y.%m.%d-%H:%M:%S.%f')[:-3],'\n')
    f = open(out_aln_tensor,'w')
    refBases_out = open(out_refBases,'w')
    refSeq=pysam.Fastafile(ref_fasta)    
    bamFP = pysam.AlignmentFile(bam_file, "rb")
    read_numb=0
    time_stamp = datetime.datetime.now()
    #print('finish open bam:',time_stamp.strftime('%Y.%m.%d-%H:%M:%S.%f')[:-3],'\n')
    #tensor_record = MagicDict()
    for variant_chr in sorted(begin2end.keys()):
        tensor_record = MagicDict()
        tensor_record_snp_count = MagicDict()
        read_count = MagicDict()
        exit_flag = 'false'
        for center in sorted(begin2end[ variant_chr ].keys()):
            for vcf_ref in begin2end[ variant_chr ][ center ].keys():
                for vcf_alt in begin2end[ variant_chr ][ center ][ vcf_ref ].keys():
                    exit_flag = 'false'
                    time_stamp = datetime.datetime.now()
                    #print('start center:',time_stamp.strftime('%Y.%m.%d-%H:%M:%S.%f')[:-3],'\n')
                    repeat_range_ref_seq=(refSeq.fetch(reference=variant_chr, start=center-2, end=center+19)).upper()
                    #print(repeat_range_ref_seq)
                    for read in bamFP.fetch(variant_chr, center-1, center +1):
                        calculate_snp_numb = 'no'
                        MAQ=read.mapping_quality
                        cigarLine = read.cigar
                        ReadStart = read.pos+1
                        ReadSeq = read.seq
                        ReadSeq = ReadSeq.upper()#some base in lowercase format
                        time_stamp = datetime.datetime.now()
                        #print(MAQ)
                        #print('finish read infor:',time_stamp.strftime('%Y.%m.%d-%H:%M:%S.%f')[:-3],'\n')
                        fragment_index=0
                        ReadEnd=0
                        read_length=len(ReadSeq)
                        pos=ReadStart
                        cigar_array = []
                        center_base = ""
                        snp_numb_in_region=0
                        snp_numb_in_abnormal_read=0
                        center_base_temp=''
                        abnormal_cigar=0
                        if read.is_proper_pair:
                            Proper_Paire=1
                        else:
                            Proper_Paire=0
                            #print("no",MAQ)
                        if ReadStart <= center:
                            for ( cigarType, cigarLength ) in cigarLine:
                                cigar_array.append(cigarType)
                                if cigarType in (0, 2, 7, 8):# M,D,=,X
                                    ReadEnd += cigarLength
                                elif cigarType not in (0, 1, 2, 4, 7, 8):
                                    abnormal_cigar=1
                            ReadEnd = ReadEnd+pos-1
                            #print(abnormal_cigar)
                            #if( (not( read.is_unmapped ))& (MAQ>=25) & (abnormal_cigar ==0) ):
                            if( (not( read.is_unmapped )) & (abnormal_cigar ==0) & (MAQ>=20)):
                                time_stamp = datetime.datetime.now()
                                #print(MAQ)
                                #print('query read infor:',time_stamp.strftime('%Y.%m.%d-%H:%M:%S.%f')[:-3],'\n')
                                cigar_index = 0
                                pos_temp=pos
                                base_idx=0
                                #print('query:',variant_chr,center,ReadStart,cigarLine,ReadSeq)
                                center_type,center_idx_in_fragment = get_center_index_base(variant_chr,center,ReadStart,cigarLine)
                                
                                for ( cigarType, cigarLength ) in cigarLine:
                                    if(pos>center+20):# exit if the pos out of the range [center-20,center+20]
                                        break
                                    #if pos==28470449:
                                    #print('pos, count:',pos, read_count[ variant_chr ][ center ][ pos ])
                                    if cigarType in (0, 7, 8):# M,=,X
                                        time_stamp = datetime.datetime.now()
                                        #print('snv start:',cigarLength,' ',time_stamp.strftime('%Y.%m.%d-%H:%M:%S.%f')[:-3],'\n')
                                        for j in range(cigarLength):
                                            if(pos>center+2):# exit if the pos out of the range [center-2,center+2]
                                                break 
                                            if (pos == center):# only query pos located on the region [center-20,center+20]
                                                ref=(refSeq.fetch(reference=variant_chr, start=pos-1, end=pos)).upper()
                                                if ref == 'N':
                                                    if tensor_record[ variant_chr ][ center ] != {}:#all the value will be zero in this pos
                                                        del tensor_record[ variant_chr ][ center ]
                                                        exit_flag = 'true'
                                                        break
                                                
                                                if read_count[ variant_chr ][ center ][ pos ]=={}:
                                                    read_count[ variant_chr ][ center ][ pos ]=1
                                                else:
                                                    read_count[ variant_chr ][ center ][ pos ] +=1
                                                
                                                read_base = base2num[(ReadSeq[ fragment_index ]).upper()]
                                                ref = base2num[ref]
                                                nextCigar = get_next_cigar(cigar_index,cigar_array,j,cigarLength)
                                                if nextCigar != 'ins':
                                                    if pos == center:
                                                        # if pos == center and read.is_read1 and read_base == 0:
                                                            # read_numb +=1
                                                            # print(read.qname,':',read_numb)
                                                            # if read.qname == 'D00360:94:H2YT5BCXX:1:1101:13256:71812' and pos ==center:
                                                                # print('base,cigarLine,ReadStart,center:',read_base,cigarLine,ReadStart,center)
                                                                # print('seq: ',ReadSeq)
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'type' ] = 'snv'
                                                        if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'RC' ] == {}:
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'RC' ] = 0
                                                        if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'FC' ] == {}:
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'FC' ] = 0
                                                        if read.is_read1:
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'FC' ] +=1
                                                        else:
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'RC' ] +=1
                                                        if pos == center:# store the homopolymer infor, average of the MAQ of the read support alt, number of forward and reverse reads support the snv
                                                            center_base_temp=read_base
                                                        if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'MAQ' ] == {}:
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'MAQ' ] = float(MAQ)
                                                        else:
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'MAQ' ] += float(MAQ)
                                                            #print('MAQ:',pos,read_base,MAQ,tensor_record[ variant_chr ][ center ][ pos ][ read_base ][ 'MAQ' ],'\n')
                                                        if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'BAQ' ] == {}:
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'BAQ' ] = float(read.query_qualities[fragment_index])
                                                        else:
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'BAQ' ] += float(read.query_qualities[fragment_index])
                                                            #print('BAQ:',pos,read_base,read.query_qualities[fragment_index],tensor_record[ variant_chr ][ center ][ pos ][ read_base ][ 'BAQ' ],'\n')
                                                        if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'Proper_Paire' ] == {}:
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'Proper_Paire' ] = Proper_Paire
                                                        else:
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'Proper_Paire' ] += Proper_Paire
                                                        base_in_homo,homo_len,dist_to_homo=get_homolen_and_distance(refSeq,variant_chr,pos,center)
                                                        gc_freq_around_5bp,gc_freq_around_10bp,gc_freq_around_20bp,at_freq_around_5bp,at_freq_around_10bp,at_freq_around_20bp=get_GC_freq(refSeq,variant_chr,pos)
                                                        if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'ambigurous_alt_count' ] == {}:
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'ambigurous_alt_count' ] = 0
                                                        #print(repeat_range_ref_seq,(ReadSeq[ fragment_index ]).upper())
                                                        in_unit,change_to_repeat,unit_length=on_unit_region(repeat_range_ref_seq,(ReadSeq[ fragment_index ]).upper())
                                                        if read_base == ref:
                                                            change_to_repeat=0
                                                        #if(read_base != 3):
                                                        #print(read.qname,cigarLine,read_base,MAQ,Proper_Paire,read.query_qualities[fragment_index],ReadStart,fragment_index,ReadSeq[fragment_index-5:fragment_index+5])
                                                            #print(read)
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ]['base_in_homo']=base_in_homo
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ]['homo_len']=homo_len
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ]['dist_to_homo']=dist_to_homo
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ]['gc_freq_around_5bp']=gc_freq_around_5bp
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ]['gc_freq_around_10bp']=gc_freq_around_10bp
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ]['gc_freq_around_20bp']=gc_freq_around_20bp
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ]['at_freq_around_5bp']=at_freq_around_5bp
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ]['at_freq_around_10bp']=at_freq_around_10bp
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ]['at_freq_around_20bp']=at_freq_around_20bp
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ]['in_unit']=in_unit
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ]['change_to_repeat']=change_to_repeat
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ]['unit_length']=unit_length
                                                        distance=dist_to_end(read_length,fragment_index)
                                                        if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'distance' ] == {}:
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'distance' ] = distance
                                                        else:
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'distance' ] += distance
                                                        if (pos-ReadStart<=3) or (ReadEnd-pos<=3):
                                                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ read_base ][ 'ambigurous_alt_count' ] += 1
                                                    #if read_base ==0:
                                                    #    print(center,ref,read_base,calculate_snp_numb,snp_numb_in_region,cigarLine,cigarType,cigar_index,fragment_index,j,ReadStart)
                                                else:
                                                    pos += 1
                                                    fragment_index += 1
                                                    break
                                            #end if
                                            j += 1        
                                            pos += 1
                                            fragment_index += 1
                                        #end for j
                                        time_stamp = datetime.datetime.now()
                                        #print('snv end:',cigarLength,' ',time_stamp.strftime('%Y.%m.%d-%H:%M:%S.%f')[:-3],'\n')
                                    elif cigarType in (1, 4):# soft clip or Insertion
                                        if cigar_index == 0 and cigarType == 4:
                                            for i in range(cigarLength):                                     
                                                fragment_index += 1
                                            cigar_index += 1
                                            continue
                                        
                                        pos = pos -1
                                        cigar_final_idx=len(cigar_array)-1
                                        if(cigar_index == cigar_final_idx) and((cigar_array[cigar_final_idx-1] in(1,4)) and (cigar_array[cigar_final_idx] in(1,4))):
                                            break
                                        if (pos>=center-20) and (pos<=center+20):
                                            alt_without_ancor_base=''
                                            alt_without_ancor = 0
                                            in_unit=0
                                            readbase_in_homo=0
                                            ref=(refSeq.fetch(reference=variant_chr, start=pos-1, end=pos)).upper()
                                            if ref == 'N':
                                                if tensor_record[ variant_chr ][ center ] != {}:
                                                    del tensor_record[ variant_chr ][ center ]
                                                    exit_flag = 'true'
                                                    break
                                            if cigarLength == 1:
                                                alt_without_ancor_base = (ReadSeq[fragment_index]).upper()
                                            else:
                                                alt_without_ancor_base = (ReadSeq[fragment_index:(fragment_index+cigarLength)]).upper()
                                            for i in range(cigarLength):
                                                read_base = base2num[(alt_without_ancor_base[i]).upper()]
                                                alt_without_ancor += read_base+10
                                                #ins_index += 1
                                            ref=base2num[ref]
                                            if pos == center:
                                                if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'MAQ' ] == {}:
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'MAQ' ] = float(MAQ)
                                                else:
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'MAQ' ] += float(MAQ)
                                                if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'BAQ' ] == {}:
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'BAQ' ] = float(read.query_qualities[fragment_index])#use the quality of the first inser base as the base quality of the insert bases
                                                else:
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'BAQ' ] += float(read.query_qualities[fragment_index])
                                                if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'Proper_Paire' ] == {}:
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'Proper_Paire' ] = Proper_Paire
                                                else:
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'Proper_Paire' ] += Proper_Paire
                                                    
                                                tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'type' ] = 'ins'
                                                if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'FC' ] == {}:
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'FC' ] = 0
                                                if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'RC' ] == {}:
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'RC' ] = 0
                                                if read.is_read1:
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'FC' ] += 1
                                                else:
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'RC' ] += 1
                                                if pos == center:
                                                    center_base_temp=ref+alt_without_ancor+10
                                                    calculate_snp_numb = 'yes'
                                                in_unit,change_to_repeat,unit_length=on_unit_region(repeat_range_ref_seq,alt_without_ancor_base)
                                                base_in_homo,homo_len,dist_to_homo=get_homolen_and_distance(refSeq,variant_chr,pos,center)
                                                gc_freq_around_5bp,gc_freq_around_10bp,gc_freq_around_20bp,at_freq_around_5bp,at_freq_around_10bp,at_freq_around_20bp=get_GC_freq(refSeq,variant_chr,pos)
                                                tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ]['base_in_homo']=base_in_homo
                                                tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ]['homo_len']=homo_len
                                                tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ]['dist_to_homo']=dist_to_homo
                                                tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ]['in_unit']=in_unit
                                                tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ]['change_to_repeat']=change_to_repeat
                                                tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ]['unit_length']=unit_length
                                                tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ]['gc_freq_around_5bp']=gc_freq_around_5bp
                                                tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ]['gc_freq_around_10bp']=gc_freq_around_10bp
                                                tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ]['gc_freq_around_20bp']=gc_freq_around_20bp
                                                tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ]['at_freq_around_5bp']=at_freq_around_5bp
                                                tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ]['at_freq_around_10bp']=at_freq_around_10bp
                                                tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ]['at_freq_around_20bp']=at_freq_around_20bp
                                                distance=dist_to_end(read_length,fragment_index)
                                                if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'distance' ] == {}:
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'distance' ] = distance
                                                else:
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'distance' ] += distance
                                                if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'ambigurous_alt_count' ] == {}:
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'ambigurous_alt_count' ] = 0
                                                if (pos-ReadStart<=3) or (ReadEnd-pos<=3):
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref+alt_without_ancor+10 ][ 'ambigurous_alt_count' ] += 1
                                            #end if
                                            
                                            for i in range(cigarLength):                                     
                                                fragment_index += 1
                                        pos += 1
                                        #end if
                                    elif cigarType == 2:# deletion
                                        first_base_pos = pos
                                        for i in range(cigarLength):
                                            if(pos>center+20):#query pos located on the region [center+10,ReadEnd]
                                                break
                                            if (pos>=center-20) and (pos<=center+20):
                                                if pos == center:
                                                    ref=refSeq.fetch(reference=variant_chr, start=pos-1, end=pos)
                                                    ref = ref.upper()
                                                    if ref == 'N':
                                                        if tensor_record[ variant_chr ][ center ] != {}:
                                                            del tensor_record[ variant_chr ][ center ]
                                                            exit_flag = 'true'
                                                            break
                                                    if read_count[ variant_chr ][ center ][ pos ]=={}:
                                                        read_count[ variant_chr ][ center ][ pos ]=1
                                                    else:
                                                        read_count[ variant_chr ][ center ][ pos ] +=1
                                                    ref = base2num[ref]
                                                    #in_unit=0
                                                    readbase_in_homo_from_first=0
                                                    readbase_in_homo_from_second=0
                                                    homo_len_from_first=0
                                                    homo_len_from_second=0
                                                    del_bases_2_number=0
                                                    #whether ref is same with the one base around ref(both side)
                                                    #whether ref is same with the three bases around ref(both side)
                                                    del_bases = refSeq.fetch(reference=variant_chr, start=first_base_pos-1, end=first_base_pos+cigarLength-1)
                                                    for i in range(cigarLength):
                                                        read_base = base2num[(del_bases[i]).upper()]
                                                        del_bases_2_number += read_base+10
                                                    #print(center,pos,del_bases,unit_seq)
                                                    if pos == center:
                                                        center_base_temp=del_bases_2_number
                                                        #if (ReadStart<=center-10) and (ReadEnd>=center+10):
                                                        calculate_snp_numb = 'yes'
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'type' ] = 'del'
                                                    if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'BAQ' ] == {}:
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'BAQ' ] = float(read.query_qualities[fragment_index])#use the quality of the base before the deltion as the base quality of the delete bases
                                                    else:
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'BAQ' ] += float(read.query_qualities[fragment_index])
                                                    if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'MAQ' ] == {}:
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'MAQ' ] = float(MAQ)
                                                    else:
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'MAQ' ] += float(MAQ)
                                                    if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'Proper_Paire' ] == {}:
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'Proper_Paire' ] = Proper_Paire
                                                    else:
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'Proper_Paire' ] += Proper_Paire
                                                    if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'FC' ] == {}:
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'FC' ] = 0
                                                    if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'RC' ] == {}:
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'RC' ] = 0
                                                    if read.is_read1:
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'FC' ] += 1
                                                    else:
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'RC' ] += 1
                                                    in_unit,change_to_repeat,unit_length=on_unit_region(repeat_range_ref_seq,del_bases)
                                                    base_in_homo,homo_len,dist_to_homo=get_homolen_and_distance(refSeq,variant_chr,pos,center)
                                                    gc_freq_around_5bp,gc_freq_around_10bp,gc_freq_around_20bp,at_freq_around_5bp,at_freq_around_10bp,at_freq_around_20bp=get_GC_freq(refSeq,variant_chr,pos)
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ]['in_unit']=in_unit
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ]['change_to_repeat']=change_to_repeat
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ]['unit_length']=unit_length
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ]['base_in_homo']=base_in_homo
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ]['homo_len']=homo_len
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ]['dist_to_homo']=dist_to_homo
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ]['gc_freq_around_5bp']=gc_freq_around_5bp
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ]['gc_freq_around_10bp']=gc_freq_around_10bp
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ]['gc_freq_around_20bp']=gc_freq_around_20bp
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ]['at_freq_around_5bp']=at_freq_around_5bp
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ]['at_freq_around_10bp']=at_freq_around_10bp
                                                    tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ]['at_freq_around_20bp']=at_freq_around_20bp
                                                    distance=dist_to_end(read_length,fragment_index)
                                                    if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'distance' ] == {}:
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'distance' ] = distance
                                                    else:
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'distance' ] += distance
                                                    
                                                    #print('xxx:',center,pos,del_bases_2_number,readbase_in_homo_from_first,homo_len_from_first,readbase_in_homo_from_second,homo_len_from_second)
                                                    
                                                    #print('xxx2:',center,pos,tensor_record[ variant_chr ][ center ][ pos ][ del_bases_2_number ]['in_homo_from_first'])
                                                    if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'ambigurous_alt_count' ] == {}:
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'ambigurous_alt_count' ] = 0
                                                    if (pos-ReadStart<=3) or (ReadEnd-pos<=3):
                                                        tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ del_bases_2_number ][ 'ambigurous_alt_count' ] += 1
                                                #end if pos==center
                                            #end if pos in [-20,20]
                                            pos += 1
                                        #end for i
                                        if exit_flag == 'true':
                                            break#exit to the for read line loop
                                    if exit_flag == 'true':
                                        break#exit to the for read line loop
                                    #end if
                                    cigar_index += 1
                            #end for cigar
                            if exit_flag == 'true':
                                break
                    #end for read
                #end for ref
            #end for alt
        #end for center
        for center in sorted((tensor_record[ variant_chr ].keys())):
            for vcf_ref in list(tensor_record[ variant_chr ][ center ].keys()):
                for vcf_alt in list(tensor_record[ variant_chr ][ center ][ vcf_ref ].keys()):
                    tensor_record_out = np.zeros( (1, 19) )
                    tensor_record_out[0][6]=1.00
                    ReadPos=center
                    ref=base2num[(refSeq.fetch(reference=variant_chr, start=ReadPos-1, end=ReadPos)).upper()]
                    for ReadBase in list(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ].keys()):
                        #print(ReadBase)
                        if ReadBase != ref and tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ReadBase]['type'] !={}:
                            temp_AltBC = tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'FC' ]+tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'RC' ]
                            for ReadBase2 in list(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ].keys()):
                                ReadBase_type = tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ReadBase]['type']
                                ReadBase2_type = tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ReadBase2]['type']
                                if (ReadBase != ReadBase2) and (ReadBase2 != ref) and (ReadBase_type == ReadBase2_type): 
                                    
                                    temp_AltBC2 = tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase2 ][ 'FC' ]+tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase2 ][ 'RC' ]
                                    if temp_AltBC >= temp_AltBC2:
                                        del tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase2 ]
                                    else:
                                        del tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ]
                                        ReadBase = ReadBase2
                                        temp_AltBC = tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'FC' ]+tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'RC' ]
                            alt_count = tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'FC' ] + tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'RC' ]
                            if alt_count <= 2:# exclude the variant site with the reads not more than 2
                                del tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ]
                                read_count[ variant_chr ][ center ][ ReadPos]=read_count[ variant_chr ][ center ][ ReadPos]-alt_count
                            #break
                        #end if
                    #end for ReadBase
                    
                    k = 0
                    exist_ref_allele=0
                    for ReadBase in list(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ].keys()):
                        k = k+1#allele number
                        if ReadBase == ref:
                            exist_ref_allele=1
                        if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ReadBase]['type']=={}:
                            del tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ReadBase]
                            k = k-1
                            if ReadBase == ref:
                                exist_ref_allele=0
                    for ReadBase in list(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ].keys()):
                        variant_type = 0
                        if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ReadBase]['type'] == 'del':#deletion
                            variant_type=3
                        elif tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ReadBase]['type'] == 'ins':#insertion
                            variant_type=2
                        elif ReadBase != ref:#snv
                            variant_type=1
                        if ref == 4:#ambigurous base(M,N...)
                            variant_type=1
                        tensor_record_out[0][ref]=1#ref:0->4
                        if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'FC' ]=={}:
                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'FC' ]=0
                        if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'RC' ]=={}:
                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'RC' ]=0
                        if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'FC' ]=={}:
                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'FC' ]=0
                        if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'RC' ]=={}:
                            tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'RC' ]=0
                        base_count=float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'FC' ])+float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'RC' ])
                        ref_base_count=float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'FC' ])+float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'RC' ])
                        RefBaq=41.0
                        RefMaq=60.0
                        if variant_type !=0:
                            AF=float(base_count)/float(read_count[ variant_chr ][ center ][ ReadPos])
                            AF=1.0/(1.0+math.exp(-AF))*(float(base_count-2)/float(base_count))
                            
                            AD_FC=float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'FC' ])
                            AD_RC=float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'RC' ])
                            strand_bias_AD=('%.3f' %(float(abs(AD_FC-AD_RC))/(float(AD_FC+AD_RC)/2.0)))
                            AltBaq=float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'BAQ' ])/float(base_count)/float(41.0)
                            AltMaq=float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'MAQ' ])/float(base_count)/float(60.0)
                            Alt_proper_paire=('%.3f' %(float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'Proper_Paire' ])/float(base_count)))
                            if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'BAQ' ] != {}:
                                RefBaq=float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'BAQ' ])/float(ref_base_count)/float(41.0)
                            if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'MAQ' ] != {}:
                                RefMaq=float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'MAQ' ])/float(ref_base_count)/float(60.0)
                            #print(RefBaq,AltBaq)
                            Baq_bias=float(RefBaq-AltBaq)/(float(RefBaq+AltBaq)/2.0)
                            Maq_bias=float(RefMaq-AltMaq)/(float(RefMaq+AltMaq)/2.0)
                            change_to_repeat=tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'change_to_repeat' ]
                            unit_length=tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'unit_length' ]
                            Dist_bias=0
                            if tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'distance' ] != {}:
                                Dist_sum=float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'distance' ])+float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'distance' ])
                                Dist_bias=(float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'distance' ])-float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'distance' ]))/(float(Dist_sum)/2.0)
                            
                            #print(center,AD_FC,AD_RC,tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'MAQ' ],base_count,read_count[ variant_chr ][ center ][ ReadPos])
                            #print(center,ReadBase,AD_FC,AD_RC,strand_bias_AD)
                            #print(Baq_bias,Maq_bias,RefMaq,AltMaq,RefBaq,AltBaq)
                            #print(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ref ][ 'dist_to_homo' ],tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'dist_to_homo' ])
                            tensor_record_out[0][5]=AF
                            tensor_record_out[0][6]=strand_bias_AD
                            tensor_record_out[0][7]=('%.3f' % AltBaq)
                            tensor_record_out[0][8]=('%.3f' % AltMaq)
                            tensor_record_out[0][9]=('%.3f' % Baq_bias)
                            tensor_record_out[0][10]=('%.3f' % Maq_bias)
                            tensor_record_out[0][11]=Alt_proper_paire
                            tensor_record_out[0][13]=Dist_bias
                            tensor_record_out[0][15]=change_to_repeat
                            tensor_record_out[0][18]=unit_length
                        #end if alt base
                        GC_5bp_around2side_freq=tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'gc_freq_around_5bp' ]
                        GC_10bp_around2side_freq=tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'gc_freq_around_10bp' ]
                        homo_len=tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'homo_len' ]
                        dist_to_homo=tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'dist_to_homo' ]
                        
                        tandem_repeat_region=tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'in_unit' ]
                        distance_2_end=('%.3f' %(float(tensor_record[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][ ReadBase ][ 'distance' ])/float(base_count)))
                        tensor_record_out[0][12]=GC_10bp_around2side_freq
                        tensor_record_out[0][14]=distance_2_end
                        tensor_record_out[0][16]=tandem_repeat_region
                        tensor_record_out[0][17]=read_count[ variant_chr ][ center ][ ReadPos]
                        
                    #end for
                    if k ==1 and exist_ref_allele == 1:
                        tensor_record_out[0][6]=('%.3f' %(float(1.0)))
                    #print(tensor_record_out[0][6])
                    #i +=1
                    #end for ReadPos
                    output_line = []
                    output_refBases= []
                    output_line.append( "%d %s %d %d %s %s" % (int(tag),variant_chr,center,int(begin2end[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][3]),str(vcf_ref),str(vcf_alt)) )
                    output_refBases.append( "%d %s %d %d %s %s" % (int(tag),variant_chr,center,int(begin2end[ variant_chr ][ center ][ vcf_ref ][ vcf_alt ][3]),str(vcf_ref),str(vcf_alt)) )
                    #print(tensor_record_out[21][0][8])
                    for record in np.reshape(tensor_record_out, 1*19):
                        #output_line.append("%0.2f" % record)
                        output_line.append("%.3f" % record)
                    out_line=" ".join(output_line)
                    refBases=get_seq(variant_chr,center,refSeq)
                    for record in np.reshape(refBases, 41*4):
                        output_refBases.append("%.3f" % record)
                    out_refBases=" ".join(output_refBases)
                    f.write(out_line)
                    f.write('\n')
                    refBases_out.write(out_refBases+'\n')
                #end for vcf_alt
            #end for vcf_ref
        #end for center
    #end for chr
    f.close()
    refBases_out.close()
    time_stamp = datetime.datetime.now()
    print('Finish get features:',time_stamp.strftime('%Y.%m.%d-%H:%M:%S.%f')[:-3],'\n')
                       
