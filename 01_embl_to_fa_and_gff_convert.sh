#!/bin/bash
# Laura Dean
# 24/9/26

# Info on how I converted embl flatfile format to fasta and gff


## I used this website to convert the .embl file to fasta: https://www.bioinformatics.org/sms2/embl_fasta.html
## but it didn't retain the contig names so I used this awk code instead:
## then I used this to convert embl file to fasta
EMBL=/gpfs01/home/mbzlld/data/circos_tradis_plot/CP008827.1.embl
EMBL=/gpfs01/home/mbzlld/data/circos_tradis_plot/CP009273.1.embl

# awk conversion command
awk '
/^ID   / {
    id = $2
    sub(/;/, "", id)
    sv = $4
    sub(/;/, "", sv)
    print ">" id "." sv
    inseq = 0
}

/^SQ   / {
    inseq = 1
    next
}

/^\/\// {
    inseq = 0
    next
}

inseq {
    gsub(/[0-9 ]/, "")
    print
}
' $EMBL > ${EMBL%.*}.fa

# and this website to convert it to gff3: https://www.ebi.ac.uk/ena/gff3/converter/
# but its out of order managed to get round that by sorting the downstream files


