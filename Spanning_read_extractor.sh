#!/bin/bash

if [[ -z $1 || -z $2 || -z $3 ]]; then
  echo 'Error: one or more mandatory variables are undefined'
  exit 1
fi

if [ ! -f $1 ]; then
    echo "Error: Bam file not found"
    exit 1
fi

if [ ! -f $1.bai ]; then
    echo "Error: Bam index not found"
    exit 1
fi

if [ ! -f $3 ]; then
    echo "Error: fast5 index not found"
    exit 1
fi

bam=$1
region=$2
index=$3
filebase=`echo ${1##*/} | cut -d '.' -f1`
extension=1000

chromosome=`echo $region | cut -d ':' -f1`
TR_left=`echo $region | cut -d ':' -f2 | cut -d '-' -f1`
TR_right=`echo $region | cut -d ':' -f2 | cut -d '-' -f2`

####### Positive reads #######

positive_left=(`samtools view -F 0x10 $1 $chromosome:$((TR_left - $extension))-$TR_left | cut -f1`)

positive_right=(`samtools view -F 0x10 $1 $chromosome:$TR_right-$((TR_right + $extension)) | cut -f1`)

l2=" ${positive_right[*]} "
for item in ${positive_left[@]}; do
  if [[ $l2 =~ " $item " ]] ; then
    positive_result+=($item)
  fi
done

####### Negative reads #######

negative_left=(`samtools view -f 0x10 $1 $chromosome:$((TR_left - $extension))-$TR_left | cut -f1`)

negative_right=(`samtools view -f 0x10 $1 $chromosome:$TR_right-$((TR_right + $extension)) | cut -f1`)

l2=" ${negative_right[*]} "
for item in ${negative_left[@]}; do
  if [[ $l2 =~ " $item " ]] ; then
    negative_result+=($item)
  fi
done

####### Write temporary output file #######

for read in ${positive_result[@]};
do echo -e "$read\tpositive" >> ${filebase}_tmp_spanning_reads.tsv;
done

for read in ${negative_result[@]};
do echo -e "$read\tnegative" >> ${filebase}_tmp_spanning_reads.tsv;
done

sort -k1 ${filebase}_tmp_spanning_reads.tsv > ${filebase}_tmp_spanning_reads_sort.tsv

####### Join with index #######

echo -e "name\tstrand\tpath" > ${filebase}_spanning_reads.tsv

zcat $3 | join -t $'\t' -1 1 -2 1 ${filebase}_tmp_spanning_reads_sort.tsv - >> ${filebase}_spanning_reads.tsv
rm ${filebase}_tmp_spanning_reads.tsv
rm ${filebase}_tmp_spanning_reads_sort.tsv
