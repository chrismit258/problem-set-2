 #! /usr/bin/env bash

 datasets='/vol2/home/smithc/data-sets'

 # Use BEDtools intersect to identify the size of the largest overlap
 # between CTCF and H3K4me3 locations. 
 CTCF="$datasets/bed/encode.tfbs.chr22.bed.gz"
 H3K4me3="$datasets/bed/encode.h3k4me3.hela.chr22.bed.gz"

answer_1=$(zcat $CTCF \
    | awk '$4=="CTCF"' \
    | bedtools intersect -wo -a - -b $H3K4me3 \
    | sort -k15nr \
    | head -n1 \
    | awk '{print $15}')

echo "answer-1: $answer_1"

#Use BEDtools to calculate the GC content of nucleotides 19,000,000 to
#19,000,500 on chr22 of hg19 genome build. Report the GC content as a
#fraction (e.g., 0.50). 
fasta="$datasets/fasta/hg19.chr22.fa"
interval="$datasets/bed/interval.bed"

answer_2=$(bedtools nuc -fi $fasta -bed $interval \
    | grep -v '^#' \
    | awk '{print $5}')

echo "answer-2: $answer_2"

#Use BEDtools to identify the length of the CTCF ChIP-seq peak (i.e.,
#interval) that has the largest mean signal in ctcf.hela.chr22.bg.gz.
HeLa="$datasets/bedtools/ctcf.hela.chr22.bg.gz"

answer_3=$(zcat $CTCF \
    | awk '$4=="CTCF"' \
    | bedtools map -c 4 -o mean -a - -b $HeLa \
    | sort -k5nr \
    | awk 'BEGIN {OFS="\t"} {print $0, $3 - $2}' \
    | head -n1 \
    | awk '{print $6}')

echo "answer-3: $answer_3"

#Use BEDtools to identify the gene promoter (defined as 1000 bp upstream
#of a TSS) with the highest median signal in ctcf.hela.chr22.bg.gz. Report
#the gene name (e.g., 'ABC123').
TSS="$datasets/bed/tss.hg19.chr22.bed.gz"

answer_4=$(bedtools map -c 4 -o median -a $TSS -b $HeLa \
    | sort -k7nr \
    | head -n1 \
    | awk '{print $4}')

echo "answer-4: $answer_4"

#Use BEDtools to identify the longest interval on chr22 that is not
#covered by genes.hg19.bed.gz. Report the interval like chr1:100-500.
genes=$datasets/bed/genes.hg19.bed.gz
genome=$datasets/genome/hg19.genome

answer_5=$(bedtools complement -i $genes -g $genome \
    | awk '$1 == "chr22"' \
    | awk 'BEGIN {OFS="\t"} {print $0, $3-$2}' \
    | sort -k4nr \
    | head -n1 \
    | awk '{print$1 ":" $2 "-" $3}')

echo "answer-5: $answer_5"





