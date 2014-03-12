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
GTF="$DIR/SELECT_GENES.gtf"
SIMCONFIG="$DIR/$SUBDIR/config.txt"
GENELIST="$DIR/SELECT_GENES.txt"
EXPRESSION="$DIR/$SUBDIR/abundance.txt"

## Step 0: Create directory
echo ""
echo "Step 0 - Creating Directories"

mkdir -p $DIR/$SUBDIR/$TOPOUT
mkdir -p $DIR/$SUBDIR/$CUFFLINKS

echo ""
echo "$DIR/$SUBDIR/$TOPOUT created"
echo "$DIR/$SUBDIR/$CUFFLINKS created"

## Step 1: Generate gene abundance file, and simulate reads
echo ""
echo "Step 1 - Generating Reads"
python abundance_assigner.py -i $GENELIST -a $ABUNDANCE -o $EXPRESSION
java -Xmx1G -Xms1000m -jar /u/home/c/chelseaj/project/software/RNAseqSim-CJ.jar $SIMCONFIG -GTF_File=$GTF -FASTA_File=$GENOME -Chromosome=chr[0-9XYM]* -Chromosome_Matching=Fuzzy -Abundance_File=$EXPRESSION -Abundance_Overwritten=No -Expressed_Transcript_Percentage=1 -SV_Allowed=No -Read_ID_Prefix=ENST- -Read_Length=$READLEN -Coverage_Factor=$COVERAGE -Quality_Generator=Prefect -Flip_And_Reverse=Yes -Unknown_Factor=0.0001 -Output_Fastq_1=$DIR/$SUBDIR/1.fq -Output_Fastq_2=$DIR/$SUBDIR/2.fq


## Step 2: Mapping the reads using tophat
echo ""
echo "Step 2 - Mapping Reads"
tophat2 -o $DIR/$SUBDIR/$TOPOUT $GENOMEINDEX $DIR/$SUBDIR/1.fq $DIR/$SUBDIR/2.fq


# step 3: run cufflinks to compute FPKM
echo ""
echo "Step 3 - Assembling Transcripts"
cufflinks -o $DIR/$SUBDIR/$CUFFLINKS $DIR/$SUBDIR/$TOPOUT/accepted_hits.bam










