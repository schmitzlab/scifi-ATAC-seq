#!/usr/bin/perl
use strict;
use warnings;

my $usage = <<_EOUSAGE_;
#------Usage----------
#perl $0 meta test.vcf ref_geno_code_in_vcf min_depth
#

#------Description---------
#make comparision vcf for ref geno type.
# ref geno mean the ture geno, and the alt geno means the other 7 genotypes other than ref geno.
#output vcf with normal vcf with 2 genotype, REF, ALT

_EOUSAGE_
    ;
die $usage if (scalar @ARGV == 1);


# filter BAM file
open (VCF, '<', $ARGV[0]) or die;

my $min_dp = 0;
if($ARGV[2])
{
 $min_dp = $ARGV[2]; 
}
my $ref = $ARGV[1];
my $ref_index = 0;
while(<VCF>){
    if (/^##/)
	{
	 print;
	 next;
	}elsif(/^#/){
	 chomp;
         my @tmp = split /\t/;
	 for (my $i = 9; $i <= $#tmp; $i++)
	  {
	   $ref_index = $i if $tmp[$i] eq $ref;
	  }
	 die "Error: Can not find genotype $ref in $ARGV[0]\n" if $ref_index == 0;
	 print join("\t", @tmp[0..8]);
	 print "\t$ref\tAlt\n";
	}else{
	 chomp;
	 my @tmp1 = split /\t/;
         #ignore low depth ref;
         my @ref_info = split /:/, $tmp1[$ref_index];
         next if $ref_info[1] < $min_dp;
	 my $alt;
	 my $dp = 0;
	 if ($tmp1[$ref_index] =~ /^0\/0(:\S+)/){
	   for (my $i = 9; $i <= $#tmp1; $i++)
	   {
	    if($tmp1[$i] =~ /^1\/1:\S+/){
		 my @tmp2 = split /:/,$tmp1[$i];
		 #save the max dp as alt;
		 if($tmp2[1] > $dp){
		  $dp = $tmp2[1];
		  $alt = $tmp1[$i];
		 }
		}
	   }
	   next if $dp < $min_dp;
	   print join("\t", @tmp1[0..8]);
	   print "\t$tmp1[$ref_index]\t$alt\n";
	 }elsif($tmp1[$ref_index] =~ /^1\/1(:\S+)/){
	   for (my $i = 9; $i <= $#tmp1; $i++)
	   {
	    if($tmp1[$i] =~ /^0\/0:\S+/){
		 my @tmp2 = split /:/,$tmp1[$i];
		 #save the max dp as alt;
		 if($tmp2[1] > $dp){
		  $dp = $tmp2[1];
		  $alt = $tmp1[$i];
		 }
		}
	   }
	   next if $dp < $min_dp;
	   print join("\t", @tmp1[0..8]);
	   print "\t$tmp1[$ref_index]\t$alt\n";
	 }
	}
}
close VCF;
