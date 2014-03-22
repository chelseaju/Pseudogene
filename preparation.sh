#!/bin/bash
## sh preparation.sh DIR COVERAGE READLEN ABUNDANCE  

DIR=$1
COVERAGE=$2
READLEN=$3
ABUNDANCE=$4

SIMULATOR="/u/home/c/chelseaj/project/software/RNAseqSim-CJ.jar"
GENOME="/u/home/c/chelseaj/project/database/Ensembl/Genome/Homo_sapiens_GRCh37_74_genome.fa"
GENOMEINDEX="/u/home/c/chelseaj/project/database/Ensembl/Bowtie2Index/genome"
ERROR1="/u/home/c/chelseaj/project/database/SequenceQuality/paired_101/hiseq_1.fastq"
ERROR2="/u/home/c/chelseaj/project/database/SequenceQuality/paired_101/hiseq_2.fastq"

#SIMULATOR="RNAseqSim-CJ.jar"
#GENOME="/home/chelseaju/Database/Ensembl/Genome/Homo_sapiens_GRCh37_74_genome.fa"
#GENOMEINDEX="/home/chelseaju/Database/Ensembl/Bowtie2Index/genome"
#VCF="/home/chelseaju/Database/Ensembl/VCF/ASW-12156-01.vcf"
#ERROR1="/home/chelseaju/Database/SequenceQuality/paired_101/hiseq_1.fastq"
#ERROR2="/home/chelseaju/Database/SequenceQuality/paired_101/hiseq_2.fastq"

CUFFLINKS="cufflinks_before"
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
java -Xmx1G -Xms1000m -jar $SIMULATOR $SIMCONFIG -GTF_File=$GTF -FASTA_File=$GENOME -Chromosome=[0-9XYM]* -Chromosome_Matching=Fuzzy -Abundance_File=$EXPRESSION -Abundance_Overwritten=No -Expressed_Transcript_Percentage=1 -Read_ID_Prefix=ENST- -Read_Length=$READLEN -Coverage_Factor=$COVERAGE -Quality_Generator=Real -Real_Quality_Score_Fastq_1=$ERROR1 -Real_Quality_Score_Fastq_2=$ERROR2 -Flip_And_Reverse=Yes -Unknown_Factor=0.0001 -Output_Fastq_1=$DIR/$SUBDIR/1.fq -Output_Fastq_2=$DIR/$SUBDIR/2.fq


#$SIMULATOR $SIMCONFIG -GTF_File=$GTF -FASTA_File=$GENOME -Chromosome=[0-9XYMT]* -Chromosome_Matching=Fuzzy -Abundance_File=$EXPRESSION -Abundance_Overwritten=No -Expressed_Transcript_Percentage=1 -SV_Allowed=No -Read_ID_Prefix=ENST- -Read_Length=$READLEN -Coverage_Factor=$COVERAGE -Quality_Generator=Real -Real_Quality_Score_Fastq_1=$ERROR1 -Real_Quality_Score_Fastq_2=$ERROR2 -Flip_And_Reverse=Yes -Unknown_Factor=0.0001 -Output_Fastq_1=$DIR/$SUBDIR/1.fq -Output_Fastq_2=$DIR/$SUBDIR/2.fq




## Step 2: Mapping the reads using tophat
echo ""
echo "Step 2 - Mapping Reads"
tophat2 -o $DIR/$SUBDIR/$TOPOUT $GENOMEINDEX $DIR/$SUBDIR/1.fq $DIR/$SUBDIR/2.fq


# step 3: run cufflinks to compute FPKM
echo ""
echo "Step 3 - Assembling Transcripts"
cufflinks -o $DIR/$SUBDIR/$CUFFLINKS $DIR/$SUBDIR/$TOPOUT/accepted_hits.bam










