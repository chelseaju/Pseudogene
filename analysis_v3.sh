#!/bin/bash

## bash analysis_v3.sh DIR COVERAGE READLEN ABUNDANCE
## ex bash analysis_v3.sh select_one_pseudogene_110 1 100 10

DIR=$1
COVERAGE=$2
READLEN=$3
ABUNDANCE=$4
TOPOUT="tophat_out"
SUBDIR=$COVERAGE"X_"$READLEN"L_"$ABUNDANCE"A"
BAM="accepted_hits.bam"
CHROMOSOME=('chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8'
	'chr9' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16'
	'chr17' 'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chrX' 'chrY' 'chrM')

#CHROMOSOME=('chr1' 'chr2')

echo "Data Analysis Version 3:"

echo "Step 0: Folders and Files Preparation"
CMD0="mkdir -p $DIR/$SUBDIR/$TOPOUT/mapping"
CMD1="python bamfile_sorter.py -i $BAM -d $DIR/$SUBDIR/$TOPOUT"
echo ">$CMD0"
echo ">$CMD1"
$CMD0
$CMD1

# STEP 1: Read Counts - a series of scripts
## Compute the expected number of read for each gene
## Map each read to known gene or identify the mapped region
## Compute the observed number of read for each region
## Compute the distribution of read generated from parent genes

echo ""
echo "Step 1: Read Counts For Transcripts"

# expect read count
python expected_counter.py -d $DIR/$SUBDIR/$TOPOUT 

for chr in "${CHROMOSOME[@]}"
do
  echo ""
  echo "CHROMOSOME ${chr} =========================================="
  python exon_identifier.py -d $DIR/$SUBDIR/$TOPOUT -c ${chr}
  python transcript_identifier.py -d $DIR/$SUBDIR/$TOPOUT -c ${chr}
  python observed_counter.py  -d $DIR/$SUBDIR/$TOPOUT -c ${chr} -t transcripts
  python observed_separator.py -d $DIR/$SUBDIR/$TOPOUT -c ${chr} -t transcripts

done

# observed read distribution
python distribution_equation.py -d $DIR/$SUBDIR/$TOPOUT -t transcripts
python distribution_matrix.py -d $DIR/$SUBDIR/$TOPOUT -t transcripts
#python distribution_coefficient_merger.py -d $DIR/$SUBDIR/$TOPOUT/transcripts_distribution.matrix -c $DIR/transcripts_lasso_coefficient.xls

echo ""
echo "Step 2: Read Counts for Genes"

# expect read count
python expected_counter.py -d $DIR/$SUBDIR/$TOPOUT
python ENSP2ENSG.py -d $DIR/$SUBDIR/$TOPOUT -i transcripts_expected_read_count.txt -o genes_expected_read_count.txt

for chr in "${CHROMOSOME[@]}"
do
  echo ""
  echo "CHROMOSOME ${chr} =========================================="
  
  python exon_identifier.py -d $DIR/$SUBDIR/$TOPOUT -c ${chr}
  python gene_identifier.py -d $DIR/$SUBDIR/$TOPOUT -c ${chr}
  python observed_counter.py -d $DIR/$SUBDIR/$TOPOUT -c ${chr} -t genes
  python observed_separator.py -d $DIR/$SUBDIR/$TOPOUT -c ${chr} -t genes
  python ENSP2ENSG.py -d $DIR/$SUBDIR/$TOPOUT/mapping -i ${chr}_genes_distribution.txt -o ${chr}_genes_distribution.txt.new
  mv $DIR/$SUBDIR/$TOPOUT/mapping/${chr}_genes_distribution.txt.new $DIR/$SUBDIR/$TOPOUT/mapping/${chr}_genes_distribution.txt

done

# observed read distribution
python distribution_equation.py -d $DIR/$SUBDIR/$TOPOUT -t genes
mv $DIR/$SUBDIR/$TOPOUT/genes_distribution.eqn $DIR/$SUBDIR/$TOPOUT/genes_distribution.eqn.backup
python ENSP2ENSG.py -d $DIR/$SUBDIR/$TOPOUT -i genes_distribution.eqn.backup -o genes_distribution.eqn
python distribution_matrix.py -d $DIR/$SUBDIR/$TOPOUT -t genes



echo ""
echo "COMPLETE !! "
