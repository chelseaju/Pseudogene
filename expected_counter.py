"""
Usage: python expected_counter.py -i accepted_hit.bam -o outputPrefix
Input: -i bam file with mapped reads, -o the directory and prefix of output files
Output: A list of gene with the number of reads. {gene_name \t number_of_read}
Function: Assuming read with the same prefix (gene name) comes from the same gene.
        The script iterates through the bam file and count expected number of read for each gene.
Date: 2013-12-29
Author: Chelsea Ju
Note: this script is adopted from fragments_counter
"""

