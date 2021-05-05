#!/usr/bin/perl
use strict;
use warnings;

my $tp_snp = $ARGV[0];
my $fp_snp = $ARGV[1];
my $tp_indel = $ARGV[2];
my $fp_indel = $ARGV[3];
my $out_tp_snp=$ARGV[4];
my $out_fp_snp=$ARGV[5];
my $out_tp_indel=$ARGV[6];
my $out_fp_indel=$ARGV[7];


my $seq_len = 10;#1000000000
my $hash_seed = 7;

my (%all_rec,%snp_rec,%indel_rec) = ();
print"Start store tpsnp infor...\n";
open(TSNP,"<$tp_snp") or die "can not open $tp_snp";
while(my $line = <TSNP>)
{
	chomp $line;
	my @recs = split/ /,$line;
	if($recs[1]=~/^chr/){
		$recs[1]=~s/chr//;
	}
	if($recs[1]=~/^MT/){
		$recs[1]=~s/MT/M/;
	}
	$snp_rec{$recs[1]}{$recs[2]}='snp';
	$all_rec{$recs[1]}{$recs[2]}=0;
}
close TSNP;

print"Start store fpsnp infor...\n";
open(Fsnp,"<$fp_snp") or die "can not open $fp_snp";
while(my $line = <Fsnp>)
{
	chomp $line;
	my @recs = split/ /,$line;
	if($recs[1]=~/^chr/){
		$recs[1]=~s/chr//;
	}
	if($recs[1]=~/^MT/){
		$recs[1]=~s/MT/M/;
	}
	$snp_rec{$recs[1]}{$recs[2]}='snp';
	$all_rec{$recs[1]}{$recs[2]}=0;
}
close Fsnp;

print"Start store tpindel infor...\n";
open(Tindel,"<$tp_indel") or die "can not open $tp_indel";
while(my $line = <Tindel>)
{
	chomp $line;
	my @recs = split/ /,$line;
	if($recs[1]=~/^chr/){
		$recs[1]=~s/chr//;
	}
	if($recs[1]=~/^MT/){
		$recs[1]=~s/MT/M/;
	}
	$indel_rec{$recs[1]}{$recs[2]}='indel';
	$all_rec{$recs[1]}{$recs[2]}=0;
}
close Tindel;

print"Start store fpindel infor...\n";
my %group_rec=();
open(Findel,"<$fp_indel") or die "can not open $fp_indel";
while(my $line = <Findel>)
{
	chomp $line;
	my @recs = split/ /,$line;
	if($recs[1]=~/^chr/){
		$recs[1]=~s/chr//;
	}
	if($recs[1]=~/^MT/){
		$recs[1]=~s/MT/M/;
	}
	$indel_rec{$recs[1]}{$recs[2]}='indel';
	$all_rec{$recs[1]}{$recs[2]}=0;
}
close Findel;

print"Start merge snp and indel infor...\n";
foreach my $chr(keys %all_rec){
	foreach my $pos(keys %{$all_rec{$chr}}){
		my $group = &getIndex($seq_len,$hash_seed,$pos);
		my @pos_arr=();
		$pos_arr[0]=$pos;
		if(exists $group_rec{$chr}{$group}){
			push @{$group_rec{$chr}{$group}},@pos_arr;
		}else{
			@{$group_rec{$chr}{$group}}=@pos_arr;
		}
	}
}

print"Start out tpsnp infor...\n";
open(Tsnp,"<$tp_snp") or die "can not open $tp_snp";
open(Otsnp,">$out_tp_snp") or die "can not open $out_tp_snp";
while(my $line = <Tsnp>)
{
	chomp $line;
	my @recs = split/ /,$line;
	if($recs[1]=~/^chr/){
		$recs[1]=~s/chr//;
	}
	if($recs[1]=~/^MT/){
		$recs[1]=~s/MT/M/;
	}
	my $group=&getIndex($seq_len,$hash_seed,$recs[2]);
	
	my @positions=@{$group_rec{$recs[1]}{$group}};
	my $tag_label=0;
	if($#positions>=1){
		my $type_temp='snp';
		for(my $i=0;$i<=$#positions;$i++){
			if(($positions[$i]>=$recs[2]-10) and ($positions[$i]<=$recs[2]+10)){
				if(exists $indel_rec{$recs[1]}{$positions[$i]}){
					$tag_label=3;#snp+indel
					last;
				}elsif($positions[$i] !=$recs[2]){
					$tag_label=1;#snp+snp
				}
			}
		}
	}
	
	print Otsnp "$line $tag_label\n";
}
close Tsnp;
close Otsnp;

print"Start out fpsnp infor...\n";
open(Fsnp,"<$fp_snp") or die "can not open $fp_snp";
open(Ofsnp,">$out_fp_snp") or die "can not open $out_fp_snp";
while(my $line = <Fsnp>)
{
	chomp $line;
	my @recs = split/ /,$line;
	if($recs[1]=~/^chr/){
		$recs[1]=~s/chr//;
	}
	if($recs[1]=~/^MT/){
		$recs[1]=~s/MT/M/;
	}
	my $group=&getIndex($seq_len,$hash_seed,$recs[2]);
	
	my @positions=@{$group_rec{$recs[1]}{$group}};
	my $tag_label=0;
	if($#positions>=1){
		my $type_temp='snp';
		for(my $i=0;$i<=$#positions;$i++){
			if(($positions[$i]>=$recs[2]-10) and ($positions[$i]<=$recs[2]+10)){
				if(exists $indel_rec{$recs[1]}{$positions[$i]}){
					$tag_label=3;#snp+indel
					last;
				}elsif($positions[$i] !=$recs[2]){
					$tag_label=1;#snp+snp
				}
			}
		}
	}
	
	print Ofsnp "$line $tag_label\n";
}
close Fsnp;
close Ofsnp;

print"Start out tpindel infor...\n";
open(Tindel,"<$tp_indel") or die "can not open $tp_indel";
open(Otindel,">$out_tp_indel") or die "can not open $out_tp_indel";
while(my $line = <Tindel>)
{
	chomp $line;
	my @recs = split/ /,$line;
	if($recs[1]=~/^chr/){
		$recs[1]=~s/chr//;
	}
	if($recs[1]=~/^MT/){
		$recs[1]=~s/MT/M/;
	}
	my $group=&getIndex($seq_len,$hash_seed,$recs[2]);
	
	my @positions=@{$group_rec{$recs[1]}{$group}};
	my $tag_label=0;
	if($#positions>=1){
		my $type_temp='indel';
		for(my $i=0;$i<=$#positions;$i++){
			if(($positions[$i]>=$recs[2]-10) and ($positions[$i]<=$recs[2]+10)){
				if($positions[$i] !=$recs[2]){
					if(exists $indel_rec{$recs[1]}{$positions[$i]}){
						$tag_label=2;#indel+indel
						last;
					}else{
						$tag_label=3;#indel+snp
					}
				}elsif(exists $snp_rec{$recs[1]}{$positions[$i]}){
					$tag_label=3;#indel+snp
				}
			}
		}
	}
	
	print Otindel "$line $tag_label\n";
}
close Tindel;
close Otindel;

print"Start out fpindel infor...\n";
open(Findel,"<$fp_indel") or die "can not open $fp_indel";
open(Ofindel,">$out_fp_indel") or die "can not open $out_fp_indel";
while(my $line = <Findel>)
{
	chomp $line;
	my @recs = split/ /,$line;
	if($recs[1]=~/^chr/){
		$recs[1]=~s/chr//;
	}
	if($recs[1]=~/^MT/){
		$recs[1]=~s/MT/M/;
	}
	my $group=&getIndex($seq_len,$hash_seed,$recs[2]);
	
	my @positions=@{$group_rec{$recs[1]}{$group}};
	my $tag_label=0;
	if($#positions>=1){
		my $type_temp='indel';
		for(my $i=0;$i<=$#positions;$i++){
			if(($positions[$i]>=$recs[2]-10) and ($positions[$i]<=$recs[2]+10)){
				if($positions[$i] !=$recs[2]){
					if(exists $indel_rec{$recs[1]}{$positions[$i]}){
						$tag_label=2;#indel+indel
						last;
					}else{
						$tag_label=3;#indel+snp
					}
				}elsif(exists $snp_rec{$recs[1]}{$positions[$i]}){
					$tag_label=3;#indel+snp
				}
			}
		}
	}
	
	print Ofindel "$line $tag_label\n";
}
close Findel;
close Ofindel;


sub getIndex()
{
	my $len = shift;
	my $seed = shift;
	my $num = shift;
	my $complete = substr('0'x($len - length($num)).$num,0,$seed);
	return int($complete);
}