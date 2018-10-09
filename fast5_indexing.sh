#!/bin/bash
awk '{print $1 "\t" $2}'  $1 | sort -k1 > summary_awk_sort.txt & 

P1=$!

find $2 -name '*.fast5' > fast5_find.txt;

cat fast5_find.txt | sed 's/.*\///' | paste fast5_find.txt - | sort -k2 - > fast5_find_sort.txt &

P2=$!

wait $P1 $P2

join -1 1 -2 2 -t $'\t' summary_awk_sort.txt fast5_find_sort.txt | awk '{print $2 "\t" $3}' | sort -k1 - | gzip - > ${3}_fast5_index.gz
rm summary_awk_sort.txt
rm fast5_find_sort.txt
rm fast5_find.txt