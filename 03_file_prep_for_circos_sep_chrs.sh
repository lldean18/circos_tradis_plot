#!/bin/bash
# Laura Dean
# 23/9/26

# script to prep files for and draw circos plot for each contig

# setup env
#srun --partition defq --cpus-per-task 4 --mem 20g --time 08:00:00 --pty bash

# setup software
source $HOME/.bash_profile
conda activate circos

# set variables
PROJECT=Photorhabdus_khanii
#ASSEMBLY=/gpfs01/home/mbzlld/data/bryant/11d3a3246d_20251024_Bryant1L/assembly/323630L_Photorhabduskhanii_reordered.fasta
#ANNOTATION=/gpfs01/home/mbzlld/data/bryant/11d3a3246d_20251024_Bryant1L/my_assembly_edits/323630L_Photorhabduskhanii_for_fasta_rev_order.gff
ASSEMBLY=/share/bryant_lab/reference_genomes/323630L_Photorhabduskhanii.fna
ANNOTATION=/share/bryant_lab/reference_genomes/323630L_Photorhabduskhanii.gff
# set the window size in kb for plotting insertion sites
window_size=2

# setup wkdir
cd /gpfs01/home/mbzlld/data/circos_tradis_plot
mkdir -p /gpfs01/home/mbzlld/data/circos_tradis_plot/$PROJECT
cd /gpfs01/home/mbzlld/data/circos_tradis_plot/$PROJECT

# copy the ticks file to the current working dir
cp /gpfs01/home/mbzlld/github/circos_tradis_plot/ticks.conf ./

CONTIG_LIST=( contig_1 contig_2 )

#TRADIS_LIST=(/gpfs01/home/mbzlld/photorhabdus_test_PK123_WITH_newcutadapt/biotradis/PK123_combined.out.contig_1.insert_site_plot.gz /gpfs01/home/mbzlld/photorhabdus_test_PK123_WITH_newcutadapt/biotradis/PK123_combined.out.contig_2.insert_site_plot.gz)

#####################
### PREP ASSEMBLY ###
#####################

# convert assembly to the right format
# first index it
conda activate samtools1.24
samtools faidx $ASSEMBLY
conda deactivate

for CONTIG in ${CONTIG_LIST[@]}
do
# without naming the chrs with their chr names
grep "$CONTIG" $ASSEMBLY.fai |  awk '{print "chr - " $1 " " $1 " 0 " $2 " chr1"}' > karyotype_$CONTIG.txt
done

##############################
### PREP GENOME ANNOTATION ###
##############################

for CONTIG in ${CONTIG_LIST[@]}
do
# convert annotation to circos format
awk 'match($1, "'$CONTIG'") && $3=="CDS"' $ANNOTATION | awk '{print $1, $4, $5}' OFS="\t" > genes_$CONTIG.txt
sort -k1,1 -k2,2n genes_$CONTIG.txt > genes_$CONTIG.txt.tmp && mv genes_$CONTIG.txt.tmp genes_$CONTIG.txt

# make separate annotation files for genes on fwd and rev strands (strand info is 7th column)
# fwd strand
awk 'match($1, "'$CONTIG'") && $3=="CDS" && $7=="+"' $ANNOTATION | awk '{print $1, $4, $5}' OFS="\t" > genes_fwd_strand_$CONTIG.txt
sort -k1,1 -k2,2n genes_fwd_strand_$CONTIG.txt > genes_fwd_strand_$CONTIG.txt.tmp && mv genes_fwd_strand_$CONTIG.txt.tmp genes_fwd_strand_$CONTIG.txt
# rev strand
awk 'match($1, "'$CONTIG'") && $3=="CDS" && $7=="-"' $ANNOTATION | awk '{print $1, $4, $5}' OFS="\t" > genes_rev_strand_$CONTIG.txt
sort -k1,1 -k2,2n genes_rev_strand_$CONTIG.txt > genes_rev_strand_$CONTIG.txt.tmp && mv genes_rev_strand_$CONTIG.txt.tmp genes_rev_strand_$CONTIG.txt
done

################################
### PREP INSERTION SITE DATA ###
################################

for CONTIG in ${CONTIG_LIST[@]}
do
# reformat the tradis insertion site data to have the contig name and location at the start
grep "$CONTIG" $ASSEMBLY.fai > $CONTIG.info.txt
rm insertions_fwd_strand_$CONTIG.bed insertions_rev_strand_$CONTIG.bed
zcat /gpfs01/home/mbzlld/photorhabdus_test_PK123_WITH_newcutadapt/biotradis/PK123_combined.out.$CONTIG.insert_site_plot.gz |
awk -v OFS='\t' '
NR == FNR {
    chr = $1
    pos = 0
    next
}
{
    print chr, pos, pos+1, $1 >> "insertions_fwd_strand_'$CONTIG'.bed"
    print chr, pos, pos+1, $2 >> "insertions_rev_strand_'$CONTIG'.bed"
    pos++
}
' "$CONTIG.info.txt" -
# make genome windows to count insertion sites in
bedtools makewindows -g $CONTIG.info.txt -w ${window_size}000 > windows_${window_size}kb_$CONTIG.bed
# count the insertions per window
bedtools map -a windows_${window_size}kb_$CONTIG.bed -b insertions_fwd_strand_$CONTIG.bed -c 4 -o sum -null 0 > insertions_fwd_strand_${window_size}kb_$CONTIG.bed
bedtools map -a windows_${window_size}kb_$CONTIG.bed -b insertions_rev_strand_$CONTIG.bed -c 4 -o sum -null 0 > insertions_rev_strand_${window_size}kb_$CONTIG.bed
# cleanup
rm $CONTIG.info.txt

##  ######################################
##  # prev version
##  # convert the tradis insertion site output to bed format
##  rm insertions_fwd_strand.bed insertions_rev_strand.bed
##  touch insertions_fwd_strand.bed insertions_rev_strand.bed
##  zcat /gpfs01/home/mbzlld/photorhabdus_test_PK123_WITH_newcutadapt/biotradis/PK123_combined.out.$CONTIG.insert_site_plot.gz | awk '
##  BEGIN {
##      while ((getline < "'$ASSEMBLY.fai'") > 0) {
##          chr[++n] = $1
##          len[n] = $2
##      }
##      c = 1
##      pos = 0
##  }
##  {
##      print chr[c], pos, pos+1, $1 >> "'insertions_fwd_strand_$CONTIG.bed'"
##      print chr[c], pos, pos+1, $2 >> "'insertions_rev_strand_$CONTIG.bed'"
##      pos++
##      if (pos >= len[c]) {
##          c++
##          pos = 0
##      }
##  }
##  ' OFS='\t'
##  
##  # make genome windows to count insertion sites in
##  bedtools makewindows -g $ASSEMBLY.fai -w 5000 > windows_5kb.bed
##  
##  # count the insertions per window
##  bedtools map -a windows_5kb.bed -b insertions_fwd_strand_$CONTIG.bed -c 4 -o sum -null 0 > insertions_fwd_strand_5kb_$CONTIG.bed
##  bedtools map -a windows_5kb.bed -b insertions_rev_strand_$CONTIG.bed -c 4 -o sum -null 0 > insertions_rev_strand_5kb_$CONTIG.bed

##########################
### PREP THE ORIC FILE ###
##########################

awk 'match($1, "'$CONTIG'") && $3=="oriC" {
    midpoint=int(($4+$5)/2)
    print $1, midpoint, midpoint, 1
}' OFS="\t" "$ANNOTATION" > oriC_$CONTIG.txt

awk 'match($1, "'$CONTIG'") && $3=="oriC" {
    midpoint=int(($4+$5)/2)
    print $1, midpoint, midpoint, "oriC"
}' OFS="\t" "$ANNOTATION" > oriC_label_$CONTIG.txt
done

#############################
# PREP THE CIRCOS CONF FILE #
#############################

for CONTIG in ${CONTIG_LIST[@]}
do
echo "
<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
<<include ticks.conf>>

<image>
background = white
dir   = .
#dir  = conf(configdir)
file  = circos_$CONTIG.png
png   = yes
svg   = yes

# radius of inscribed circle in image
radius         = 1500p
# by default angle=0 is at 3 o'clock position
angle_offset      = -90
#angle_orientation = counterclockwise
auto_alpha_colors = yes
auto_alpha_steps  = 5
</image>

karyotype = karyotype_$CONTIG.txt

########################################

<ideogram>

show_label       = yes
label_font       = default
#label_radius     = 1r + 110p
label_radius     = 1.17r
label_size       = 40
label_parallel   = yes

<spacing>
default = 0.005r
</spacing>

radius*    = 0.85r
thickness = 20p
fill      = yes
color = black
</ideogram>

########################################

<plots>

########################
# gene annotation ring
########################

<plot>
type = tile
file = genes_fwd_strand_$CONTIG.txt
r1   = 0.97r
r0   = 0.92r
color = dgrey
layers = 1
margin      = 0.05u
orientation = center
stroke_thickness = 1
stroke_color     = dgrey
thickness = 50
padding = 8
</plot>

<plot>
type = tile
file = genes_rev_strand_$CONTIG.txt
r1   = 0.90r
r0   = 0.85r
color = dgrey
layers = 1
margin = 0.05u
orientation = center
stroke_thickness = 1
stroke_color = dgrey
thickness = 50
padding = 8
</plot>

############################
### INSERTION SITE PLOTS ###
############################

<plot>
type = histogram
file = insertions_fwd_strand_${window_size}kb_$CONTIG.bed
r1   = 0.80r
r0   = 0.55r
color = vdred
fill_color = vdred
thickness = 0.5
</plot>

<plot>
type = histogram
file = insertions_rev_strand_${window_size}kb_$CONTIG.bed
r1   = 0.55r
r0   = 0.30r
orientation = in
color = orange
fill_color = orange
thickness = 0.5
</plot>

##################
### oriC marker ###
##################

<plot>
type = scatter
file = oriC_$CONTIG.txt
r1 = 0.83r
r0 = 0.81r
color = dblue
stroke_color = blue
stroke_thickness = 2
glyph = circle
glyph_size = 22
</plot>

##################
### oriC label ###
##################

<plot>
type = text
file = oriC_label_$CONTIG.txt
r0 = 1.05r
r1 = 1.15r
label_size = 35p
label_font = bold
color = dblue
orientation = out
</plot>

</plots>
" > circos_$CONTIG.conf

##################
### RUN CIRCOS ###
##################

circos -conf circos_$CONTIG.conf

done




# cleanup env
conda deactivate

