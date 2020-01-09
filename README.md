# FNVC
FNVC: Filtering for Variant Calls in NGS<br>

This software is developed by Yongyong Ren and Yan Kong.<br>

It is freely available for academic use, but the copyrights of all the data in this paper belong to the original copyright owners or organizations.<br>

Introduction
------------
As costs continue to decline, whole genome sequencing is increasingly being used in the large-scale population genetic research and clinical diagnosis of rare congenital disorders. 
However, a significant proportion of identified genetic variants are false positives which may mislead the diag-nosis or the genetic study. 
Nearly half of them can be eliminated by using the post-processing filtering methods such as VQSR, Hard-Filter (HF) or GARFIELD-NGS, but comparing with the number of the eliminated false variants, ~5-23 times (up to ~110,000 SNVs and ~19,000 INDELs) detected true variants were filtered out at the same time when using these widely used methods. 
Here, we present an effective tool named FNVC (Filtering for NGS Variant Calls), which relies on the gradient boosting decision tree method to eliminate ~80% false SNVs and ~47% false INDELs, and re-called >97% true genetic variants that had been filtered out before. 
FNVC is the only method that eliminates more false variants than the loss of true, and showed a significant improvement of OFO, MCC and F1-score (P<0.005). 
Furthermore, different from the existing post-processing methods developed only for GATK variant calls, FNVC extends the application to all variant calls from the Illumina sequencing data by using the features in the read aligned BAM file.<br>

![AUC](https://github.com/yyren/FNVC/raw/master/Picture/AUC_performance.bmp)<br>

![](https://github.com/yyren/FNVC/raw/master/Picture/performance_for_patho_variants.bmp)<br>

![](https://github.com/yyren/FNVC/raw/master/Picture/performance_accross_genome.bmp)<br>

Prerequisites
------------
* python v3.4 or later
* perl v5.0 or later

Installation
------------
Uncompress the installation zip:<br>
		$ cd /my/install/dir/<br>
		$ unzip /path/to/FNVC-tools/dist/FNVC-tools-VERSION.zip<br>
		$ ./INSTALL.sh<br>

Usage
------------
FNVC-GATK model:<br>
		python FNVC_predict.py \
			--in_file input.vcf --out_file out.vcf \
			--snv_model ./models/FNVC-gatk_snv_model.bin \
			--indel_model ./models/FNVC-gatk_indel_model.bin \
			--snv_cutoff 0.5126 --indel_cutoff 0.5824 \
			--model gatk<br>

FNVC-MF model:<br>
		python FNVC_predict.py \
			--in_file input.vcf --out_file out.vcf \
			--snv_model ./models/FNVC-mf_snv_model.bin \
			--indel_model ./models/FNVC-mf_indel_model.bin \
			--snv_cutoff 0.4453 --indel_cutoff 0.4673 \
			--bam_file input.bam --reference human_g1k_v37_decoy.fasta \
			--model mf

FNVC-BAM model:<br>
		python FNVC_predict.py \
			--in_file input.vcf --out_file out.vcf \
			--snv_model ./models/FNVC-bam_snv_model.bin \
			--indel_model ./models/FNVC-bam_indel_model.bin \
			--snv_cutoff 0.5062 --indel_cutoff 0.4784 \
			--bam_file input.bam --reference human_g1k_v37_decoy.fasta \
			--model bam

Training
------------
FNVC-gatk model:<br>
		python FNVC_gatk_train.py
			--train_TP_file train_tp_features.txt --train_FP_file train_fp_features.txt \
			--test_TP_file test_tp_features.txt --test_FP_file test_fp_features.txt \
			--out_model FNVC_gatk_model.bin --out_path ~/FNVC

Data
------------
### Leave one individual out cross validation study data
The training and testing data used for the leave one individual out cross validation study are available in the: <br>

https://drive.google.com/open?id=17AOz1SChj22PTlHGPTLNlBrgjk5ExYHV <br>

The high confidence variants used in the paper for training and testing can be download from: <br>

ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release <br>

The raw sequencing data used in the paper are available in the: <br>

https://github.com/genome-in-a-bottle/giab_data_indexes/tree/master/ <br>
### Additional Test Data
The independent test data (NA12877) is described in this [paper](https://genome.cshlp.org/content/27/1/157.full). <br>

The high confidence [variants](https://github.com/Illumina/PlatinumGenomes/blob/master/files/2017-1.0.files) and raw Illumina sequencing [data](https://www.ebi.ac.uk/ena/browser/view/PRJEB3381) are also available.<br>

The features used for the VQSR, GARFIELD-NGS, Hard-Filter and FNVC can be download from: <br>

https://drive.google.com/open?id=17AOz1SChj22PTlHGPTLNlBrgjk5ExYHV <br>