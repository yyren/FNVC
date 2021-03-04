#!/usr/bin/perl
BEGIN {push @INC, "please type the path of this script in here"};

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use vcf_modules;
use base_modules;
use merge_prediction;

my $Time_Start= base_modules::get_datetime();
my $version='1.0';

##------------------------------------------------------------------------------------
##Getoptions
my %opts;
GetOptions(\%opts, "i=s", "t=s", "snpm=s", "im=s", "c=s", "tp=s", "fp=s", "o=s", "h");
if(!defined($opts{"t"}) || !defined($opts{"o"}) || !defined($opts{"c"}) || defined($opts{"h"}))
{
	print <<"	Usage End.";
	Start Time :[$Time_Start]
	Description v$version: extract qc information from html file.
	
	Usage
	Forced parameter:
		-i          vcf file for filtering                                   <infile> must be given
		-t          [annotation] annotation with the bam file or [filter] annotation with the exist features in vcf file
		            or [retrain] retaining on existent model with additional data      <character> must be given
		-snpm        pre-trained snp model, eg. gatk_snp_hg001.model
		-im          pre-trained indel model, eg. gatk_indel_hg001.model
		-c          variant caller: [gatk], [mutect], [others]
		-tp         vcf file with true-positive variants for retraining
		-fp         vcf file with false-positive variants for retraining
		-o          output annotated vcf file or the prefix name of re-trained model (with full path)                 <infile> must be given
		
	Other parameter
		-h      help document
	Example:
	filtering:
	perl $Script -i input.vcf -t annotation -snpm gatk_snp_hg001.model -im gatk_indel_hg001.model -c gatk -o /home/ubuntu/output.vcf
	
	retrain:
	perl $Script -tp hg002_tp.vcf -fp hg002_fp.vcf -t retrain -snpm gatk_snp_hg001.model -im gatk_indel_hg001.model
	     -c gatk -o /home/ubuntu/hg001_gatk_retrained
	
	Usage End.
	exit;

}

##------------------------------------------------------------------------------------
##get parameters
print "Start: $Time_Start\n";
my $infile=$opts{"i"};
my $tool=$opts{"t"};
my $snp_model=$opts{"snpm"};
my $indel_model=$opts{"im"};
my $tp_file=$opts{"tp"};
my $fp_file=$opts{"fp"};
my $variant_caller=$opts{"c"};
my $outfile=$opts{"o"};

my $infile_basename= base_modules::get_basename($infile);
my $out_dir= base_modules::get_parent_path($outfile);
my $script_path= base_modules::get_parent_path($Script);
my $snp_file=$out_dir.'/'.$infile_basename.'.snp';
my $indel_file=$out_dir.'/'.$infile_basename.'.indel';
my $snp_features_file=$out_dir.'/'.$infile_basename.'.snp.features';
my $indel_features_file=$out_dir.'/'.$infile_basename.'.indel.features';
my $prediction_snp=$out_dir.'/'.$infile_basename.'.snp.prediction';
my $prediction_indel=$out_dir.'/'.$infile_basename.'.indel.prediction';

##------------------------------------------------------------------------------------
##
if($tool eq 'annotation' or $tool eq 'filter'){
	my $out_temp_file=$out_file.'.temp';
	### divide the vcf into snp file and indel file
	vcf_modules::separate_snv_indel($infile,$snp_file,$indel_file);

	### get features
	`python '$script_path'/get_features_main.py --in_file '$snp_file' --out_file '$snp_features_file' --caller '$variant_caller' --tag 0`;
	`python '$script_path'/get_features_main.py --in_file '$indel_file' --out_file '$indel_features_file' --caller '$variant_caller' --tag 0`;
	
	### FNVC prediction
	python $script_path/FNVC_Prediction.py --in_file $snp_features_file --model $snp_model --out_file $prediction_snp
	python $script_path/FNVC_Prediction.py --in_file $indel_features_file --model $indel_model --out_file $prediction_indel
	
	### merge and sort the vcf file
	merge_predict_results::merge_prediction($prediction_snp,$prediction_indel,$snp_file,$indel_file,$out_temp_file);
	`python '$script_path'/sort_vcf.py --in_file '$out_temp_file' --out_file '$out_file'`;
	
	### delete temp files
	my @temp_files=($snp_features_file,$indel_features_file,$prediction_snp,$prediction_indel,$snp_file,$indel_file,$out_temp_file);
	unlink @temp_files;
}

if($tool eq 'retrain' ){
	my $out_snp_model=$out_file.'_retain_snp.model';
	my $out_indel_model=$out_file.'_retrain_indel.model';
	
	my $tp_snp_file=$out_dir.'/'.$tp_file.'.snp';
	my $tp_indel_file=$out_dir.'/'.$tp_file.'.indel';
	my $fp_snp_file=$out_dir.'/'.$fp_file.'.snp';
	my $fp_indel_file=$out_dir.'/'.$fp_file.'.indel';
	
	my $tp_snp_feature_file=$out_dir.'/'.$tp_file.'.snp.features';
	my $tp_indel_feature_file=$out_dir.'/'.$tp_file.'.indel.features';
	my $fp_snp_feature_file=$out_dir.'/'.$fp_file.'.snp.features';
	my $fp_indel_feature_file=$out_dir.'/'.$fp_file.'.indel.features';
	
	### divide the vcf into snp file and indel file
	vcf_modules::separate_snv_indel($tp_file,$tp_snp_file,$tp_indel_file);
	vcf_modules::separate_snv_indel($fp_file,$fp_snp_file,$fp_indel_file);
	
	### get features
	`python '$script_path'/get_features_main.py --in_file '$fp_snp_file' --out_file '$fp_snp_feature_file' --caller '$variant_caller' --tag 0 &`;
	`python '$script_path'/get_features_main.py --in_file '$tp_indel_file' --out_file '$tp_indel_feature_file' --caller '$variant_caller' --tag 1 &`;
	`python '$script_path'/get_features_main.py --in_file '$fp_indel_file' --out_file '$fp_indel_feature_file' --caller '$variant_caller' --tag 0 &`;
	`python '$script_path'/get_features_main.py --in_file '$tp_snp_file' --out_file '$tp_snp_feature_file' --caller '$variant_caller' --tag 1`;
	`wait`;
	
	### retrain new models
	`python '$script_path'/FNVC_Retrain.py --in_tp '$tp_snp_feature_file' --in_fp '$fp_snp_feature_file' --model '$snp_model' --out_file '$out_snp_model'`;
	`python '$script_path'/FNVC_Retrain.py --in_tp '$tp_indel_feature_file' --in_fp '$fp_indel_feature_file' --model '$indel_model' --out_file '$out_indel_model'`;
	
	my @temp_files=($tp_snp_file,$tp_indel_file,$fp_snp_file,$fp_indel_file,$tp_snp_feature_file,$tp_indel_feature_file,$fp_snp_feature_file,$fp_indel_feature_file);
	unlink @temp_files;
}

my $Time_End= base_modules::get_datetime();
print "Finish: $Time_End\n";

