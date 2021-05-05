#!/bin/bash

#get arguments
docker_image_file=$1 # full path
bam_file=$2 #
vcf_file=$3 # please use the same sampleId both in bam and vcf file
reference=$4 # eg. hg19.fa
sampleID=$5 #the SM tag in the bam file, eg. HG001
out_file=$6 # the reannotated output file


# SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)

#######################################################
#################### Checking environment ################

if ! type singularity >/dev/null 2>&1; then
    echo 'please install: singularity'
    exit
else
    echo 'already install: singularity'
fi

#######################################################
#################### Checking arguments ##################

file_list=($docker_image_file $bam_file $vcf_file $reference)
for file in ${file_list[*]}
do
    if [[ ! -f $docker_image_file ]]; then
        echo "Not exist: $file"
        exit
    else
        echo "Checking: $file exist"
    fi
done

########################################################
################## Define new argument #################

GATK="singularity exec $docker_image_file /home/bmap/software/GATK/gatk-4.1.9.0/gatk"

########################################################
################## Reannotation #################
$GATK --java-options "-Xmx30G -Djava.io.tmpdir=./" VariantAnnotator \
-R $reference -I $bam_file -V $vcf_file -O $out_file \
--enable-all-annotations true


