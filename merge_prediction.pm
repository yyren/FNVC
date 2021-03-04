#!/usr/bin/perl
use strict;
use warnings;


#merge_prediction.pm

package merge_prediction;
require Exporter;

@ISA=qw(Exporter);
@EXPORT=qw(merge_prediction);

sub merge_prediction {
	my @args=@_;
	my $snp_predict_results=$args[0];
	my $indel_predict_results=$args[1];
	my $snp_vcf_file=$args[2];
	my $indel_vcf_file=$args[3];
	my $out_file=$args[4];

	my $numb=0;
	my (%fnvc_predict_snp,%fnvc_predict_indel)=();

	open(IN,"<$snp_predict_results") or die "can not open $snp_predict_results\n";
	while(my $line=<IN>){
		chomp $line;
		my @recs=split(/\t/,$line);
		if($numb!=0){
			my $id=join("\_",@recs[1..4]);
			$fnvc_predict_snp{$id}=$recs[6];
		}
		$numb++;
	}
	close IN;

	$numb=0;
	open(IN,"<$indel_predict_results") or die "can not open $snp_predict_results\n";
	while(my $line=<IN>){
		chomp $line;
		my @recs=split(/\t/,$line);
		if($numb!=0){
			my $id=join("\_",@recs[1..4]);
			$fnvc_predict_indel{$id}=$recs[6];
		}
		$numb++;
	}
	close IN;

	my $infor_endline='F';
	my $fnvc_infor='##INFO=<ID=FNVC,Number=1,Type=Float,Description="FNVC prediction for variant"';
	open(IN, "<$indel_vcf_file") or die "can not open $indel_vcf_file";
	open(OUT,">$out_file") or die "can not open $out_file";
	while(my $line=<IN>){
		chomp $line;
		if($line=~m/^#/){
			if($line=~m/^##contig/ and $infor_endline eq 'F'){
				$infor_endline='T';
				print OUT "$fnvc_infor\n$line\n";
			}else{
				print OUT "$line\n";
			}
		}else{
			my @recs=split(/\t/,$line);
			my $id=join("\_",@recs[0..1,3..4]);
			my $predict_value=1.0;
			if(exists $indel_predict_results{$id}){
				$predict_value=$indel_predict_results{$id};
			}
			my $recs[7]=$recs[7].';FNVC='.$predict_value;
			my $new_line=join("\t",@recs);
			print OUT "$new_line\n";
		}
	}
	close IN;

	open(IN, "<$snp_vcf_file") or die "can not open $snp_vcf_file";
	while(my $line=<IN>){
		chomp $line;
		if($line!~m/^#/){
			my @recs=split(/\t/,$line);
			my $id=join("\_",@recs[0..1,3..4]);
			my $predict_value=1.0;
			if(exists $snp_predict_results{$id}){
				$predict_value=$snp_predict_results{$id};
			}
			my $recs[7]=$recs[7].';FNVC='.$predict_value;
			my $new_line=join("\t",@recs);
			print OUT "$new_line\n";
		}
	}
	close IN;
	close OUT;
}

1;


