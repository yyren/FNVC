#!/usr/bin/perl

use strict;
use warnings;

my $in_gatk=$ARGV[0];
my $in_bam=$ARGV[1];
my $out_file=$ARGV[2];

my %hash_rec=();
open(IN,"<$in_gatk") or die "can not open $in_gatk:$!";
while(my $line=<IN>){
    chomp $line;
    my @recs=split(' ',$line);
    my $id=join("_",@recs[1..5]);
    my $new_line=join(" ",@recs[0..6,8..$#recs]);
    $hash_rec{$id}=$new_line;
}
close IN;

# combine the features from tensor file
open(IN,"<$in_bam") or die "can not open $in_bam:$!";
open(OUT,">$out_file") or die "can not open $out_file:$!";
while(my $line=<IN>){
    chomp $line;
    my @recs=split(' ',$line);
    my $id=join("_",@recs[1..5]);
    if(exists $hash_rec{$id}){
        my $gatk_line=$hash_rec{$id};
        my $new_line=join(" ",($gatk_line,@recs[11..$#recs]));
        print OUT "$new_line\n";
    }
}
close IN;
close OUT;