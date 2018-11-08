#!/bin/bash
if [[ -z $1 || -z $2 || -z $3 ]]; then
  echo 'Error: one or more mandatory variables are undefined'
  exit 1
fi

if [ ! -f $1 ]; then
    echo "Error: summary file not found"
    exit 1
fi

if [ ! -d $2 ]; then
    echo "Error: fast5 directory not found"
    exit 1
fi

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
