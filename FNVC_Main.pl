#!/usr/bin/perl
BEGIN {push @INC, "/home/ubuntu/project/FNVC"};

use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use FindBin qw($Bin $Script);
use vcf_modules;
use base_modules;
use merge_predict_results;

my $Time_Start= base_modules::get_datetime();
my $version='1.0';

##------------------------------------------------------------------------------------
##Getoptions
my %opts;
GetOptions(\%opts, "i=s", "t=s", "snpm=s", "im=s", "c=s", "tp=s", "fp=s", "o=s", "h");
if(!defined($opts{"t"}) || !defined($opts{"c"}) || defined($opts{"h"}))
{
	print <<"	Usage End.";
	Start Time :[$Time_Start]
	Description v$version: extract qc information from html file.
	
	Usage
	Forced parameter:
		-i          vcf file for filtering                                   <file with path>
		-t          [filter] annotation with the exist features in vcf file
		            or [retrain] retaining with variants from different variant caller
		            or [incremental] incremental learning on existent model with additional data <character> must be given
		-snpm        pre-trained snp model, eg. gatk_snp_hg001.model <character>
		-im          pre-trained indel model, eg. gatk_indel_hg001.model <character>
		-c          variant caller: [gatk], [mutect], [others] <character>
		-tp         vcf file with true-positive variants for retraining  <character>
		-fp         vcf file with false-positive variants for retraining  <character>
		-o          output annotated vcf file or the prefix name of re-trained model (with full path)
		
	Other parameter
		-h      help document
	Example:
	filtering:
	perl $Script -i input.vcf -t annotation -snpm gatk_snp_hg001.model -im gatk_indel_hg001.model -c gatk -o /home/ubuntu/output.vcf
	
	incremental:
	perl $Script -tp hg002_tp.vcf -fp hg002_fp.vcf -t incremental -snpm gatk_snp_hg001.model -im gatk_indel_hg001.model
	     -c gatk
	
	retrain:
	perl $Script -tp hg002_tp.vcf -fp hg002_fp.vcf -t retrain -c gatk -o /home/ubuntu/new_model
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
my $out_file=$opts{"o"};


#my $infile_basename= base_modules::get_basename($infile);
#my $out_dir= base_modules::get_parent_path($out_file);
#my $script_path= base_modules::get_parent_path($Script);
my $script_path=$Bin;
my ($vol, $infile_dir, $infile_basename) = File::Spec->splitpath(File::Spec->rel2abs($infile));
my ($vol2, $out_dir, $outfile_basename) = File::Spec->splitpath(File::Spec->rel2abs($out_file));
my $snp_file=$out_dir.$infile_basename.'.snp';
my $indel_file=$out_dir.$infile_basename.'.indel';
my $snp_features_file_temp=$out_dir.$infile_basename.'.snp.features.temp.record';
my $indel_features_file_temp=$out_dir.$infile_basename.'.indel.features.temp.record';

my $snp_features_file=$out_dir.$infile_basename.'.snp.features.record';
my $indel_features_file=$out_dir.$infile_basename.'.indel.features.record';


my $prediction_snp=$out_dir.$infile_basename.'.snp.prediction';
my $prediction_indel=$out_dir.$infile_basename.'.indel.prediction';

#print"$snp_model\n";
##------------------------------------------------------------------------------------
##
if($tool eq 'filter'){
	my $out_temp_file=$out_file.'.temp';
	### divide the vcf into snp file and indel file
	vcf_modules::separate_snv_indel($infile,$snp_file,$indel_file);
	
	print"Getting features...\n";
	### get features
	`python '$script_path'/get_features_main.py --in_file '$snp_file' --out_file '$snp_features_file_temp' --caller '$variant_caller' --tag 0`;
	`python '$script_path'/get_features_main.py --in_file '$indel_file' --out_file '$indel_features_file_temp' --caller '$variant_caller' --tag 0`;
	
	### add feature 'CT' 
	`perl '$script_path'/add_region_feature.pl '$snp_features_file_temp' '$indel_features_file_temp' '$snp_features_file' '$indel_features_file'`;
	
	print"Predicting...\n";
	### FNVC prediction
	`python $script_path/FNVC_Prediction.py --in_file $snp_features_file --model $snp_model --out_file $prediction_snp --workpath $out_dir`;
	`python $script_path/FNVC_Prediction.py --in_file $indel_features_file --model $indel_model --out_file $prediction_indel --workpath $out_dir`;
	
	### merge and sort the vcf file
	merge_predict_results::merge_prediction($prediction_snp,$prediction_indel,$snp_file,$indel_file,$out_temp_file);
	`python '$script_path'/sort_vcf.py --in_file '$out_temp_file' --out_file '$out_file'`;
	
	### delete temp files
	my @temp_files=($snp_features_file,$indel_features_file,$snp_features_file_temp,$indel_features_file_temp,$prediction_snp,$prediction_indel,$snp_file,$indel_file,$out_temp_file);
	unlink @temp_files;

}elsif(($tool eq 'incremental') or ($tool eq 'retrain') ){
	
	my ($temp1, $tpfile_dir, $tpfile_basename) = File::Spec->splitpath(File::Spec->rel2abs($tp_file));
	my ($temp2, $fpfile_dir, $fpfile_basename) = File::Spec->splitpath(File::Spec->rel2abs($fp_file));
	
	### creat temp folder
	my $temp_folder=$tpfile_dir.'fnvc_temp/';
	if(-e $temp_folder){
		print "rm folder: $temp_folder\n";
		`rm -r '$temp_folder'`;
		print "creat folder: $temp_folder\n";
		`mkdir '$temp_folder'`;
	}else{
		print "creat folder: $temp_folder\n";
		`mkdir '$temp_folder'`;
	}
	
	my $tp_snp_file=$temp_folder.$tpfile_basename.'.snp';
	my $tp_indel_file=$temp_folder.$tpfile_basename.'.indel';
	my $fp_snp_file=$temp_folder.$fpfile_basename.'.snp';
	my $fp_indel_file=$temp_folder.$fpfile_basename.'.indel';
	
	
	my $tp_snp_feature_file=$temp_folder.$tpfile_basename.'.snp.features';
	my $tp_indel_feature_file=$temp_folder.$tpfile_basename.'.indel.features';
	my $fp_snp_feature_file=$temp_folder.$fpfile_basename.'.snp.features';
	my $fp_indel_feature_file=$temp_folder.$fpfile_basename.'.indel.features';
	
	my $tp_snp_feature_file_temp=$temp_folder.$tpfile_basename.'.snp.features.temp';
	my $tp_indel_feature_file_temp=$temp_folder.$tpfile_basename.'.indel.features.temp';
	my $fp_snp_feature_file_temp=$temp_folder.$fpfile_basename.'.snp.features.temp';
	my $fp_indel_feature_file_temp=$temp_folder.$fpfile_basename.'.indel.features.temp';
	
	### divide the vcf into snp file and indel file
	vcf_modules::separate_snv_indel($tp_file,$tp_snp_file,$tp_indel_file);
	vcf_modules::separate_snv_indel($fp_file,$fp_snp_file,$fp_indel_file);
	
	print"Getting features...\n";
	### get features
	`python '$script_path'/get_features_main.py --in_file '$fp_snp_file' --out_file '$fp_snp_feature_file_temp' --caller '$variant_caller' --tag 0`;
	`python '$script_path'/get_features_main.py --in_file '$tp_indel_file' --out_file '$tp_indel_feature_file_temp' --caller '$variant_caller' --tag 1`;
	`python '$script_path'/get_features_main.py --in_file '$fp_indel_file' --out_file '$fp_indel_feature_file_temp' --caller '$variant_caller' --tag 0`;
	`python '$script_path'/get_features_main.py --in_file '$tp_snp_file' --out_file '$tp_snp_feature_file_temp' --caller '$variant_caller' --tag 1`;
	
	### add feature 'CT' 
	`perl '$script_path'/retrain_add_region_feature.pl '$tp_snp_feature_file_temp' '$fp_snp_feature_file_temp' '$tp_indel_feature_file_temp' '$fp_indel_feature_file_temp' '$tp_snp_feature_file' '$fp_snp_feature_file' '$tp_indel_feature_file' '$fp_indel_feature_file'`;
	
	### incremental learning with additional WGS data
	if($tool eq 'retrain'){
		my $out_snp_model=$out_file.'.snp.model';
		my $out_indel_model=$out_file.'.indel.model';
		
		print"Start training\n";
		`python '$script_path'/FNVC_Retrain.py --in_tp '$tp_snp_feature_file' --in_fp '$fp_snp_feature_file' --out_file '$out_snp_model'`;
		`python '$script_path'/FNVC_Retrain.py --in_tp '$tp_indel_feature_file' --in_fp '$fp_indel_feature_file' --out_file '$out_indel_model'`;
	
	}elsif($tool eq 'incremental'){
		my ($temp3, $snpmodel_dir, $snpmodel_basename) = File::Spec->splitpath(File::Spec->rel2abs($snp_model));
		my ($temp4, $indelmodel_dir, $indelmodel_basename) = File::Spec->splitpath(File::Spec->rel2abs($indel_model));
		
		my $out_snp_model=$snp_model.'.incremental';
		my $out_indel_model=$indel_model.'.incremental';
		print"Start incremental learning\n";
		`python '$script_path'/FNVC_Retrain.py --in_tp '$tp_snp_feature_file' --in_fp '$fp_snp_feature_file' --model '$snp_model' --out_file '$out_snp_model'`;
		`python '$script_path'/FNVC_Retrain.py --in_tp '$tp_indel_feature_file' --in_fp '$fp_indel_feature_file' --model '$indel_model' --out_file '$out_indel_model'`;
	}
	`rm -r '$temp_folder'`;
	# #my @temp_files=($tp_snp_file,$tp_indel_file,$fp_snp_file,$fp_indel_file,$tp_snp_feature_file,$tp_indel_feature_file,$fp_snp_feature_file,$fp_indel_feature_file,$tp_snp_feature_file_temp,$tp_indel_feature_file_temp,$fp_snp_feature_file_temp,$fp_indel_feature_file_temp);
	# #unlink @temp_files;
}

my $Time_End= base_modules::get_datetime();
print "Finish: $Time_End\n";

