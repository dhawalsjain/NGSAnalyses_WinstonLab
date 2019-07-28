#!/usr/bin/perl
use strict;
use warnings;

#################### This script removes PCR duplicates using following approach
#### 5'end of the sequenced fragment bears random hexameric molecular barcodes to suppress PCR duplicates.
#### Removal of duplicates can be done by selecting uniq 5'end_chr_barcode locations and extracting the reads corresponding to them
#### The script needs initialization of samtools (or simply load it in Orchestra before running the code)


my ($bam, $out) = @ARGV;
my (%hashp, %hashn, %ids);
my $count_total;
my $count_noduplicates;
open(REPORT,">>PCR_duplicate removal report.txt") or die $!;
 
if (not defined $bam) {
  die "Provide a bamfile\n";
}
if (not defined $out) {
  die "Specify output file\n";
}


if (defined $bam) {

### Working on + strand
open(BAM,"samtools view -F 16 -h $bam |") or die $!;
while(<BAM>){
	 if(/^(\@)/) {next;}
	 $count_total++;
 	 my @sam = split("\t",$_);
	 $sam[9] =~ s/\s*//g;
	 $sam[2] =~ s/\s*//g;
	 $sam[3] =~ s/\s*//g;
	 my ($barcode) = ($sam[0] =~ m/_MolecularBarcode:(\S*)/g);
 	 $hashp{$barcode."_Pos_".$sam[2]."_".$sam[3]} = $sam[0];
}
close BAM;

### Working on - strand
open(BAM1,"samtools view -f 16 $bam |") or die $!;
while(<BAM1>){
	 if(/^(\@)/) {next;}
	 $count_total++;
 	 my @sam = split("\t",$_);
	 my ($barcode) = ($sam[0] =~ m/_MolecularBarcode:(\S*)/g);
	 $barcode =~ s/\s*//g;	
     $sam[9] =~ s/\s*//g;
	 $sam[2] =~ s/\s*//g;
	 $sam[3] =~ s/\s*//g;
  	 $sam[3] = $sam[3] + length($sam[9]);                  ## This step should give a 5' end of the antisense strand mapped reads
 	 $hashn{$barcode."_Neg_".$sam[2]."_".$sam[3]} = $sam[0];
}
close BAM1;

### Swapping the index values to identify names of unique reads
my ($pos,$neg);
while (my ($key, $value) = each %hashp) { $ids{$value}++; } 
undef%harshp;
while (my ($key, $value) = each %hashn) { $ids{$value}++; }
undef%harshn;

### Extractinig  unique reads from a bam file into a new output bam formatted file
open(BAM2,"samtools view -h $bam |") or die $!;
open(BOUT,"| samtools view -Sb - > $out") or die $!;
while(<BAM2>){
	 if(/^(\@)/) {
		 print BOUT "$_"; 
	     next;
	 }
	 my @sam = split(/\t+/);
    $sam[0] =~ s/\s*//g;
     if($ids{$sam[0]}) {
		  print BOUT "$_"; 
		  $count_noduplicates++;
	 }
}
 
 print REPORT "File= $bam\t Total-counted reads=$count_total \t PCR duplicate removed reads= $count_noduplicates\n";
 
 close BAM2;
 close BOUT;
 close REPORT;
}

exit;
