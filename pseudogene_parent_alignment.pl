############# pseudogene_parent_alignment.pl #############################
# Usage: perl pseudogene_parent_alignment.pl input                       #
# Args: input = pefix of the output from pseudogene_expression_parser.py #
# Author: Chelsea Ju                                                     #
# Date: 2013-05-31                                                       #
# Function: parse the sequence from pseudogene reads and parents         #
#           read, compare the two sequences to get an alignment          #
#           score. 2 output files are generated: .score and .align       #
#           .align = sequence alignment .score = preserve the info       #
#           from input, plus the alignment score                         #
##########################################################################

#!/usr/bin/perl

use strict;
use warnings;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlast;

## check user input
if($#ARGV != 0){
	print "Usage: perl pseudogene_parent_alignment.pl input\n";
	exit;
}

my $arg = $ARGV[0];
my $in_txt = $arg.".txt";
my $out_score = $arg.".score";
my $out_algn = $arg.".algn";

my @data;

#read in data
open(FILE, $in_txt) || die "Can't open $in_txt";
my $line = <FILE>;
@data = <FILE>;
close FILE;

## prepare blast

my @params = (program  => 'blastn');
my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

## retrieve data
open(SCORE, ">$out_score") || die "Can't write to $out_score";
open(ALIGN, ">$out_algn") || die "Can't write to $out_algn";

foreach $line (@data){
	my @mini_data = split("\t", $line);

	my $pseudo_id = $mini_data[0];
	my $pseudo_read_start = $mini_data[5];
	my $pseudo_read_end = $mini_data[6];
	my $pseudo_start_diff = $mini_data[9];
	my $pseudo_end_diff = $mini_data[11];
	my $pseudo_read = $mini_data[12];
	my $parent_id = $mini_data[13];
	my $parent_read_start = $mini_data[17];
	my $parent_read_end = $mini_data[18];
	my $parent_start_diff = $mini_data[21];
	my $parent_end_diff = $mini_data[23];
	my $parent_read = $mini_data[24];

	my $title = "==".$pseudo_id."::".$pseudo_read_start."-".$pseudo_read_end." vs ".$parent_id."::".$parent_read_start."-".$parent_read_end."==";

	my $pseudo_diff = abs($pseudo_start_diff) + abs($pseudo_end_diff);
	my $parent_diff = abs($parent_start_diff) + abs($parent_end_diff);

	my $pseudo_input = Bio::Seq->new(-id => "pseudo", -seq => $pseudo_read);
	my $parent_input = Bio::Seq->new(-id => "parent", -seq => $parent_read);
	my $blast_report = $factory -> bl2seq ($pseudo_input, $parent_input);

	my $result_obj = $blast_report->next_result;
	my $hit_obj = $result_obj -> next_hit;
	my $hit_score = 0;

	if($hit_obj){
		my $hsp_obj = $hit_obj -> next_hsp;

		my $alignment = $title."\n".($hsp_obj->query_string)."\n".($hsp_obj->homology_string)."\n".($hsp_obj->hit_string)."\n\n";
		print ALIGN $alignment;
		$hit_score = $hsp_obj->score;
	}
	$mini_data[12] = $pseudo_diff;
	$mini_data[24] = $parent_diff;

	my $score = join("\t",@mini_data);
	$score .= "\t".$hit_score."\n";

	print SCORE $score;
 }


close SCORE;
close ALIGN;

