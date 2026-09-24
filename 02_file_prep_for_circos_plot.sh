#!/bin/bash
# Laura Dean
# 23/9/26

# script to prep files for and draw circos plot

# setup env
srun --partition defq --cpus-per-task 4 --mem 20g --time 08:00:00 --pty bash

# setup software
source $HOME/.bash_profile
conda activate circos

# set variables
PROJECT=Klebsiella_pneumoniae
ASSEMBLY=/gpfs01/home/mbzlld/data/circos_tradis_plot/CP008827.1.fa
ANNOTATION=/gpfs01/home/mbzlld/data/circos_tradis_plot/CP008827.1.gff
TRADIS=/gpfs01/home/mbzlld/data/circos_tradis_plot/trimmed.fq.ENA_CP008827_CP008827.1.insert_site_plot_combined.gz

PROJECT=Escherichia_coli
ASSEMBLY=/gpfs01/home/mbzlld/data/circos_tradis_plot/CP009273.1.fa
ANNOTATION=/gpfs01/home/mbzlld/data/circos_tradis_plot/CP009273.1.gff
TRADIS=/gpfs01/home/mbzlld/data/circos_tradis_plot/trimmed.fq.ENA_CP009273_CP009273.1.insert_site_plot.gz

# setup wkdir
cd /gpfs01/home/mbzlld/data/circos_tradis_plot
mkdir -p /gpfs01/home/mbzlld/data/circos_tradis_plot/$PROJECT
cd /gpfs01/home/mbzlld/data/circos_tradis_plot/$PROJECT



#####################
### PREP ASSEMBLY ###
#####################

# convert assembly to the right format
# first index it
conda activate samtools1.24
samtools faidx $ASSEMBLY
conda deactivate
# without naming the chrs with their chr names
awk '{print "chr - " $1 " " $1 " 0 " $2 " chr1"}' $ASSEMBLY.fai > karyotype.txt

##############################
### PREP GENOME ANNOTATION ###
##############################

# convert annotation to circos format
awk '$3=="CDS"' $ANNOTATION | awk '{print $1, $4, $5}' OFS="\t" > genes.txt
sort -k1,1 -k2,2n genes.txt > genes.txt.tmp && mv genes.txt.tmp genes.txt

# make separate annotation files for genes on fwd and rev strands (strand info is 7th column)
# fwd strand
awk '$3=="CDS" && $7=="+"' $ANNOTATION | awk '{print $1, $4, $5}' OFS="\t" > genes_fwd_strand.txt
sort -k1,1 -k2,2n genes_fwd_strand.txt > genes_fwd_strand.txt.tmp && mv genes_fwd_strand.txt.tmp genes_fwd_strand.txt
# rev strand
awk '$3=="CDS" && $7=="-"' $ANNOTATION | awk '{print $1, $4, $5}' OFS="\t" > genes_rev_strand.txt
sort -k1,1 -k2,2n genes_rev_strand.txt > genes_rev_strand.txt.tmp && mv genes_rev_strand.txt.tmp genes_rev_strand.txt

################################
### PREP INSERTION SITE DATA ###
################################

# convert the tradis insertion site output to bed format
rm insertions_fwd_strand.bed insertions_rev_strand.bed
zcat "$TRADIS" | awk '
BEGIN {
    while ((getline < "'$ASSEMBLY.fai'") > 0) {
        chr[++n] = $1
        len[n] = $2
    }
    c = 1
    pos = 0
}
{
    print chr[c], pos, pos+1, $1 >> "insertions_fwd_strand.bed"
    print chr[c], pos, pos+1, $2 >> "insertions_rev_strand.bed"
    pos++
    if (pos >= len[c]) {
        c++
        pos = 0
    }
}
' OFS='\t'

# make genome windows to count insertion sites in
bedtools makewindows -g $ASSEMBLY.fai -w 20000 > windows_20kb.bed
bedtools makewindows -g $ASSEMBLY.fai -w 5000 > windows_5kb.bed
bedtools makewindows -g $ASSEMBLY.fai -w 1000 > windows_1kb.bed

# count the insertions per window
bedtools map -a windows_20kb.bed -b insertions_fwd_strand.bed -c 4 -o sum -null 0 > insertions_fwd_strand_20kb.bed
bedtools map -a windows_20kb.bed -b insertions_rev_strand.bed -c 4 -o sum -null 0 > insertions_rev_strand_20kb.bed

bedtools map -a windows_5kb.bed -b insertions_fwd_strand.bed -c 4 -o sum -null 0 > insertions_fwd_strand_5kb.bed
bedtools map -a windows_5kb.bed -b insertions_rev_strand.bed -c 4 -o sum -null 0 > insertions_rev_strand_5kb.bed

bedtools map -a windows_1kb.bed -b insertions_fwd_strand.bed -c 4 -o sum -null 0 > insertions_fwd_strand_1kb.bed
bedtools map -a windows_1kb.bed -b insertions_rev_strand.bed -c 4 -o sum -null 0 > insertions_rev_strand_1kb.bed

#####################
### TO RUN CIRCOS ###
#####################
# make the circos.conf file stipulating how you want the plot to be then
cp /gpfs01/home/mbzlld/github/circos_tradis_plot/circos.conf ./
cp /gpfs01/home/mbzlld/github/circos_tradis_plot/ticks.conf ./
# in the dir with the circos.conf file run:
circos

# cleanup env
conda deactivate

