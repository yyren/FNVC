#!/usr/bin/perl

#vcf_modules.pm

package vcf_modules;
require Exporter;

@ISA=qw(Exporter);
@EXPORT=qw(get_datetime,get_parent_path,get_basename);

##------------------------------------------------------------------------------------
## get time
sub get_datetime{
	my($sec,$min,$hour,$day,$mon,$year,$wday, $yday, $isdst)=localtime(time());
	$wday = $yday = $isdst = 0;
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

## get parent path
sub get_parent_path {
	my $file=$_[0];
	my @files=split(/\//,$file);
	my $parent_path=join("\/",@files[0..$#files-1]);
	return $parent_path;

}

## get file name without path
sub get_basename {
	my $file=$_[0];
	my $file_basename=$file;
	if($file=~m/\//){
		my @files=split(/\//,$file);
		$file_basename=$files[$#files];
	}
	return $file_basename;

}

1;
