#!/usr/bin/perl
use strict;
use warnings;

#---Usage---
#perl $0 meta tn5_bc

# save good barcodes
my $tn5 = $ARGV[1];
# iterate over metadata and save barcode
open (M1, '<', $ARGV[0]) or die;
my $line = 0;
while(<M1>){
	chomp;
	$line ++;
	next if $line == 1; #skip first line
	my @tmp = split;
	my $cell_bc = substr $tmp[0],5;
	my $tn5_bc = substr $tmp[0], 26, 5;
	
	if ($tn5_bc eq $tn5)
	{
	 print "$cell_bc\n";
	}
}
close M1;
