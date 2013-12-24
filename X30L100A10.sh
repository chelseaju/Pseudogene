#! /bin/bash
#$ -S /bin/bash

. /u/local/Modules/default/init/modules.sh
module load samtools

export PATH=$PATH:/u/home/c/chelseaj/project/software/tophat-2.0.7.Linux_x86_64
export PATH=$PATH:/u/home/c/chelseaj/project/software/bowtie2-2.0.6


cd /u/scratch/c/chelseaj/Pseudogene
echo $HOME
sh pre_analysis2.sh reference_top500 human72_top500_filter 30 100 10
sh analysis_v2.sh reference_top500 human72_top500_filter 30 100 10

