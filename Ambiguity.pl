#!/bin/env perl
use strict;
use warnings;
use Data::Dumper;
$Data::Dumper::Purity = 1;


#########################################################################################################
# This program checks peptide ambiguity, by analyzing the alternative peptides explaining the spectra.  #
# In case an alternative peptide exists with a similar score (less than 20), it is checked in uniprot   #
#########################################################################################################

my $wdir=$ARGV[0];
my $up="oneline.UP000005640.nr.fa"; #uniprot one line tab-delimited file 
#Input format:
#>uniprot header	sequence

my $i=$wdir."/best.msms.for.mimcry.txt"; #from MaxQuant msms.txt file ("All sequences" column without decoys)
#This input can be generater also by "Arrange.msms.inputs.py" script

#input file should be arranged in following format:
#["RID", "Sequence","Score", "All scores", "All sequences"]
#freetext    Sequence        Score   All scores      All sequences

#input example:
#QEP2_13046_YSA_921_1_301220@1   AAAAAAAAATN     67.032  67.03204;32.32016;29.14901      AAAAAAAAATN;TGPAAVNATN;PSAAAPAAYP
#QEP2_13046_YSA_917_1_301220@2   AAAAAAAAAW      75.109  75.10919;41.50215;36.39517      AAAAAAAAAW;AAGQIQAW;SGQKPAAW
#QEP2_13046_YSA_917_1_301220@3   AAAAAAAAAW      63.565  63.56507;37.55836;37.55836      AAAAAAAAAW;SGQKPAAW;FQAATGHL
#QEP2_13046_YSA_917_1_301220@4   AAAAAAAAAW      101.23  101.2281;49.58986;46.01568      AAAAAAAAAW;SGQKPAAW;ASGKPGAAW

my %o=();

#read all uniprot sequences
my $line;
open (UP, "<$up") or die "Couldn't open $up, $!";
while ($line=<UP>) {
	chomp($line);
	my @sq=split(/\t/,$line);
	$o{$sq[0]}=$sq[1];
	}
close(UP);

#print column header
print"Sequence	Alternative-peptide	Delta-score	Altentive Uniprot description	MS-file	RID	Score\n";


open (MF, "<$i") or die "Couldn't open $i, $!";
while ($line=<MF>) {
	chomp($line);

	my @mq=split(/\t/,$line);
	my $group=$mq[0]; #aureus
	my $pep=$mq[1]; #peptide sequence
	my $pep_score=$mq[2];
	my $all_scores=$mq[3];
	my $alter=$mq[4]; #all peptide alternatves. For example: "ANGGTYNNAHSIQKVVT;ANNTMINNAAAGGTLDII;ANNGQKPQNIASFVATN"
	my $pep_PEP=$mq[5];

	my @scores=split(/\;/,$all_scores);
	my @alts=split(/\;/,$alter);

	#go over all alternatives peptides of this peptide
	for (my $i=0; $i<scalar(@alts);$i++) {
		my $a=$alts[$i]; #alternative pep
		my $s=$scores[$i]; #alternative score

		my $delta_score=$pep_score-$s;
		if($a ne $pep) {
			#check further only if the score of the alternative peptide is close (less than 20) to that of the original peptide
			if($delta_score<20) { 
			foreach my $gdesc (keys %o) { #check if alternative peptide exists in uniprot
					my $seq=$o{$gdesc};
					if($seq=~/$a/){
						$group=~/^(.*?)\@/;
						printf "%s\t%s\t%3.2f\t%s\t%s\t%s\n",$pep,$a,$delta_score,$gdesc,$1,$group;
						last; #no need to find all matching WT sequences
						}
					}
				}
			}
		}
	}
