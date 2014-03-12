#!/bin/bash
## sh preparation.sh DIR COVERAGE READLEN ABUNDANCE  

DIR=$1
COVERAGE=$2
READLEN=$3
ABUNDANCE=$4

GENOME="/u/home/c/chelseaj/project/database/Ensembl/Genome/Homo_sapiens_GRCh37_74_genome.fa"
#GENOME="/home/chelseaju/Database/Homo_sapiens_UCSC/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
GENOMEINDEX="/u/home/c/chelseaj/project/database/Bowtie2Index/"
#ERROR_1="/u/scratch/c/chelseaj/database/SequenceQuality/hiseq_1.fastq"
#ERROR_2="/u/scratch/c/chelseaj/database/SequenceQuality/hiseq_2.fastq"
CUFFLINKS="cufflinks_out"
TOPOUT="tophat_out"
SUBDIR=$COVERAGE"X_"$READLEN"L_"$ABUNDANCE"A"
GTF="SELECT_GENES.gtf"
SIMCONFIG="$DIR/$SUBDIR/config.txt"

## Step 0: Create directory
echo ""
echo "Step 0 - Creating Directories"

mkdir -p $DIR/$SUBDIR/$TOPOUT
mkdir -p $DIR/$SUBDIR/$CUFFLINKS

echo "$DIR/$SUBDIR/$TOPOUT created"
echo "$DIR/$SUBDIR/$CUFFLINKS created"

echo ""
echo "Step 1 - Generating Reads"
CMD1="RNAseqSim-CJ.jar $SIMCONFIG -GTF_File=$DIR/$GTF -FASTA_File=$GENOME -Chromosome=chr[0-9XYM]* -Chromosome_Matching=Fuzzy -Abundance_File=$PREFIX""_parents_abundance_$ABUNDANCE.txt -Abundance_Overwritten=No -Expressed_Transcript_Percentage=1 -SV_Allowed=No -Read_ID_Prefix=ENSP- -Read_Length=$READLEN -Coverage_Factor=$COVERAGE -Quality_Generator=Prefect -Flip_And_Reverse=Yes -Unknown_Factor=0.0001 -Output_Fastq_1=$DIR/$SUBDIR/1.fq -Output_Fastq_2=$DIR/$SUBDIR/2.fq"

