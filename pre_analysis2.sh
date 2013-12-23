#!/bin/bash
## sh pre_analysis2.sh reference human72 0.02 100 100 
DIR=$1
PREFIX=$2
GENOME="/u/scratch/c/chelseaj/database/genome.fa"
GENOMEINDEX="/u/scratch/c/chelseaj/database/Bowtie2Index/genome"
ERROR_1="/u/scratch/c/chelseaj/database/SequenceQuality/hiseq_1.fastq"
ERROR_2="/u/scratch/c/chelseaj/database/SequenceQuality/hiseq_2.fastq"
#ERROR=$3
COVERAGE=$3
READLEN=$4
ABUNDANCE=$5
#INSERT=200
#INSERTSTD=20
CUFFLINKS="cufflinks_out"
TOPOUT="tophat_out"
SUBDIR=$COVERAGE"X_"$READLEN"L_"$ABUNDANCE"A"
SIMCONFIG="$DIR/$SUBDIR/config.txt"

## Step 0: Create directory
CMD0="mkdir -p $DIR/$SUBDIR/$TOPOUT"
CMD0a="mkdir -p $DIR/$SUBDIR/$CUFFLINKS"

$CMD0
$CMD0a

## Step 1: Generate Reads
#CMD1="dwgsim $DIR/$PREFIX""_parents.fa $DIR/$SUBDIR/$PREFIX -e $ERROR -E $ERROR -d $INSERT -s $INSERTSTD -C $COVERAGE -1 $READLEN -2 $READLEN -r 0 -F 0 -c 0 -R 0 -X 0 -y 0"


CMD1="java -Xmx1G -Xms1000m -jar /u/home/c/chelseaj/project/software/RNAseqSim-CJ.jar $SIMCONFIG -GTF_File=$DIR/$PREFIX""_parents.gtf -FASTA_File=$GENOME -Chromosome=chr[0-9XY]* -Chromosome_Matching=Fuzzy -Abundance_File=$DIR/$PREFIX""_parents_abundance_$ABUNDANCE.txt -Abundance_Overwritten=No -Expressed_Transcript_Percentage=1 -SV_Allowed=No -Read_ID_Prefix=ENSP- -Read_Length=$READLEN -Coverage_Factor=$COVERAGE -Quality_Generator=Prefect -Flip_And_Reverse=Yes -Unknown_Factor=0.0001 -Output_Fastq_1=$DIR/$SUBDIR/1.fq -Output_Fastq_2=$DIR/$SUBDIR/2.fq"

#CMD1="java -Xmx1G -Xms1000m -jar /u/home/c/chelseaj/project/software/RNAseqSim-CJ.jar $SIMCONFIG -GTF_File=$DIR/$PREFIX""_parents.gtf -FASTA_File=$GENOME -Chromosome=chr[0-9XY]* -Chromosome_Matching=Fuzzy -Abundance_File=$DIR/$PREFIX""_parents_abundance_$ABUNDANCE.txt -Abundance_Overwritten=No -Expressed_Transcript_Percentage=1 -SV_Allowed=No -Read_ID_Prefix=ENSP- -Read_Length=$READLEN -Coverage_Factor=$COVERAGE -Quality_Generator=Real -Real_Quality_Score_Fastq_1=$ERROR_1 -Real_Quality_Score_Fastq_2=$ERROR_2 -Flip_And_Reverse=Yes -Unknown_Factor=0.0001 -Output_Fastq_1=$DIR/$SUBDIR/1.fq -Output_Fastq_2=$DIR/$SUBDIR/2.fq"



echo "\n\nStep 1: Generate Reads: read_length= $READLEN, coverage= $COVERAGE"
echo ">$CMD1"
$CMD1

# step 2: mapping the reads
CMD2="tophat2 -o $DIR/$SUBDIR/$TOPOUT $GENOMEINDEX $DIR/$SUBDIR/1.fq $DIR/$SUBDIR/2.fq"

echo "\n\nStep 2: Mapping Reads with Tophat2"
echo ">$CMD2"
$CMD2

# step 3: run cufflinks to confirm FPKM
CMD3="cufflinks -o $DIR/$SUBDIR/$CUFFLINKS $DIR/$SUBDIR/$TOPOUT/accepted_hits.bam"
echo "Step 3: Assemble the Transcripts"
echo ">$CMD3"
#$CMD3


