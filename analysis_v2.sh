#!/bin/bash
## sh analysis_v2.sh reference human72_filter 100 100
DIR=$1
PREFIX=$2
GENOME="/home/chelseaju/Database/Homo_sapiens_UCSC/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
GENOMEINDEX="/home/chelseaju/Database/Homo_sapiens_UCSC/UCSC/hg19/Sequence/Bowtie2Index/genome"
COVERAGE=$3
READLEN=$4
ABUNDANCE=$5
INSERT=200
INSERTSTD=20
CUFFLINKS="cufflinks_out"
TOPOUT="tophat_out"
SUBDIR=$COVERAGE"X_"$READLEN"L_"$ABUNDANCE"A"

echo "Data Analysis Version 2:"

# step 1: count number of fragment to each transcript
CMD1a="python fragments_counter.py -p $DIR/$PREFIX""_parents.bed -s $DIR/$PREFIX""_pseudo.bed -e $DIR/$SUBDIR/$TOPOUT/accepted_hits.bam -o $DIR/$SUBDIR/$PREFIX""_transcript_none -d 0 -c 0"

CMD1b="python fragments_counter.py -p $DIR/$PREFIX""_parents.bed -s $DIR/$PREFIX""_pseudo.bed -e $DIR/$SUBDIR/$TOPOUT/accepted_hits.bam -o $DIR/$SUBDIR/$PREFIX""_transcript_random -d 0 -c 1"

CMD1c="python fragments_counter.py -p $DIR/$PREFIX""_parents.bed -s $DIR/$PREFIX""_pseudo.bed -e $DIR/$SUBDIR/$TOPOUT/accepted_hits.bam -o $DIR/$SUBDIR/$PREFIX""_transcript_all -d 0 -c 2"

echo "\n\nStep 1: Counting Mapped Fragments For Transcript"
#echo "$CMD1a"
#$CMD1a
#echo "$CMD1b"
#$CMD1b
echo "$CMD1c"
$CMD1c


# step 2: count number of fragment to each transcription locus (region with overlapped genes)
CMD2a="python fragments_counter.py -p $DIR/$PREFIX""_parents.bed -s $DIR/$PREFIX""_pseudo.bed -e $DIR/$SUBDIR/$TOPOUT/accepted_hits.bam -o $DIR/$SUBDIR/$PREFIX""_locus_none -d 1 -c 0"
CMD2b="python fragments_counter.py -p $DIR/$PREFIX""_parents.bed -s $DIR/$PREFIX""_pseudo.bed -e $DIR/$SUBDIR/$TOPOUT/accepted_hits.bam -o $DIR/$SUBDIR/$PREFIX""_locus_random -d 1 -c 1"
CMD2c="python fragments_counter.py -p $DIR/$PREFIX""_parents.bed -s $DIR/$PREFIX""_pseudo.bed -e $DIR/$SUBDIR/$TOPOUT/accepted_hits.bam -o $DIR/$SUBDIR/$PREFIX""_locus_all -d 1 -c 2"


echo "\n\nStep 2: Counting Mapped Fragments For Locus"
#$CMD2a
#$CMD2b
$CMD2c


# step 3: extract expression distribution to formula for transcripts
CMD3a="python extract_distribution_matrix.py -i $DIR/$SUBDIR/$PREFIX""_transcript_all.fragment -o $DIR/$SUBDIR/$PREFIX""_transcript_all_fragment.eqn"

CMD3b="python extract_distribution_matrix.py -i $DIR/$SUBDIR/$PREFIX""_transcript_none.fragment -o $DIR/$SUBDIR/$PREFIX""_transcript_none_fragment.eqn"

CMD3c="python extract_distribution_matrix.py -i $DIR/$SUBDIR/$PREFIX""_transcript_random.fragment -o $DIR/$SUBDIR/$PREFIX""_transcript_random_fragment.eqn"

CMD3d="python extract_distribution_matrix.py -i $DIR/$SUBDIR/$PREFIX""_transcript_all.fpkm -o $DIR/$SUBDIR/$PREFIX""_transcript_all_fpkm.eqn"

CMD3e="python extract_distribution_matrix.py -i $DIR/$SUBDIR/$PREFIX""_transcript_none.fpkm -o $DIR/$SUBDIR/$PREFIX""_transcript_none_fpkm.eqn"

CMD3f="python extract_distribution_matrix.py -i $DIR/$SUBDIR/$PREFIX""_transcript_random.fpkm -o $DIR/$SUBDIR/$PREFIX""_transcript_random_fpkm.eqn"

echo "\n\nStep 3: Extract Expression Distribution For Transcript"

$CMD3a
#$CMD3b
#$CMD3c
$CMD3d
#$CMD3e
#$CMD3f


# Step 4: extract expression distribution to formula for locus
CMD4a="python extract_distribution_matrix.py -i $DIR/$SUBDIR/$PREFIX""_locus_all.fragment -o $DIR/$SUBDIR/$PREFIX""_locus_all_fragment.eqn"

CMD4b="python extract_distribution_matrix.py -i $DIR/$SUBDIR/$PREFIX""_locus_none.fragment -o $DIR/$SUBDIR/$PREFIX""_locus_none_fragment.eqn"

CMD4c="python extract_distribution_matrix.py -i $DIR/$SUBDIR/$PREFIX""_locus_random.fragment -o $DIR/$SUBDIR/$PREFIX""_locus_random_fragment.eqn"

CMD4d="python extract_distribution_matrix.py -i $DIR/$SUBDIR/$PREFIX""_locus_all.fpkm -o $DIR/$SUBDIR/$PREFIX""_locus_all_fpkm.eqn"

CMD4e="python extract_distribution_matrix.py -i $DIR/$SUBDIR/$PREFIX""_locus_none.fpkm -o $DIR/$SUBDIR/$PREFIX""_locus_none_fpkm.eqn"

CMD4f="python extract_distribution_matrix.py -i $DIR/$SUBDIR/$PREFIX""_locus_random.fpkm -o $DIR/$SUBDIR/$PREFIX""_locus_random_fpkm.eqn"

echo "\n\nStep 4: Extract Expression Distribution For Locus"

#$CMD4a
#$CMD4b
#$CMD4c
#$CMD4d
#$CMD4e
#$CMD4f


## step 5: Extract the pseudogene equation to tab-delimited file for transcript
CMD5a="python extract_pseudogene_equation.py -i $DIR/$SUBDIR/$PREFIX""_transcript_all_fragment.eqn -o $DIR/$SUBDIR/$PREFIX""_transcript_all_fragment_pseudo_distribution.xls"
CMD5b="python extract_pseudogene_equation.py -i $DIR/$SUBDIR/$PREFIX""_transcript_none_fragment.eqn -o $DIR/$SUBDIR/$PREFIX""_transcript_none_fragment_pseudo_distribution.xls"
CMD5c="python extract_pseudogene_equation.py -i $DIR/$SUBDIR/$PREFIX""_transcript_random_fragment.eqn -o $DIR/$SUBDIR/$PREFIX""_transcript_random_fragment_pseudo_distribution.xls"

CMD5d="python extract_pseudogene_equation.py -i $DIR/$SUBDIR/$PREFIX""_transcript_all_fpkm.eqn -o $DIR/$SUBDIR/$PREFIX""_transcript_all_fpkm_pseudo_distribution.xls"
CMD5e="python extract_pseudogene_equation.py -i $DIR/$SUBDIR/$PREFIX""_transcript_none_fpkm.eqn -o $DIR/$SUBDIR/$PREFIX""_transcript_none_fpkm_pseudo_distribution.xls"
CMD5f="python extract_pseudogene_equation.py -i $DIR/$SUBDIR/$PREFIX""_transcript_random_fpkm.eqn -o $DIR/$SUBDIR/$PREFIX""_transcript_random_fpkm_pseudo_distribution.xls"

echo "Step 5:Extract pseudogene equation for transcript\n"
$CMD5a
#$CMD5b
#$CMD5c
$CMD5d
#$CMD5e
#$CMD5f


## step 6: Extract the pseudogene equation to tab-delimited file for locus
CMD6a="python extract_pseudogene_equation.py -i $DIR/$SUBDIR/$PREFIX""_locus_all_fragment.eqn -o $DIR/$SUBDIR/$PREFIX""_locus_all_fragment_pseudo_distribution.xls"
CMD6b="python extract_pseudogene_equation.py -i $DIR/$SUBDIR/$PREFIX""_locus_none_fragment.eqn -o $DIR/$SUBDIR/$PREFIX""_locus_none_fragment_pseudo_distribution.xls"
CMD6c="python extract_pseudogene_equation.py -i $DIR/$SUBDIR/$PREFIX""_locus_random_fragment.eqn -o $DIR/$SUBDIR/$PREFIX""_locus_random_fragment_pseudo_distribution.xls"

CMD6d="python extract_pseudogene_equation.py -i $DIR/$SUBDIR/$PREFIX""_locus_all_fpkm.eqn -o $DIR/$SUBDIR/$PREFIX""_locus_all_fpkm_pseudo_distribution.xls"
CMD6e="python extract_pseudogene_equation.py -i $DIR/$SUBDIR/$PREFIX""_locus_none_fpkm.eqn -o $DIR/$SUBDIR/$PREFIX""_locus_none_fpkm_pseudo_distribution.xls"
CMD6f="python extract_pseudogene_equation.py -i $DIR/$SUBDIR/$PREFIX""_locus_random_fpkm.eqn -o $DIR/$SUBDIR/$PREFIX""_locus_random_fpkm_pseudo_distribution.xls"

echo "Step 6: Extract pseudogene equation for locus\n"
#$CMD6a
#$CMD6b
#$CMD6c
#$CMD6d
#$CMD6e
#$CMD6f







