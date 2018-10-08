#!/bin/bash

######## Get parameters from input ########
reference_file=''
region=''
flank_size=250
TR_unit=''
TR_number=10
prefix="output"
report="yes"

while getopts 'f:r:l:u:n:p:e:' flag; do
  case "${flag}" in
    f) reference_file="${OPTARG}" ;;
    r) region="${OPTARG}" ;;
    l) flank_size="${OPTARG}" ;;
    u) TR_unit="${OPTARG}" ;;
	n) TR_number="${OPTARG}" ;;
	p) prefix="${OPTARG}" ;;
	e) report="${OPTARG}" ;;
    *) error "Unexpected option ${flag}" ;;
  esac
done

######## Obtain coordinates from input ########
chromosome=`echo $region | cut -d ':' -f1`

TR_left=`echo $region | cut -d ':' -f2 | cut -d '-' -f1`
left_flank_start=$((TR_left - flank_size))

TR_right=`echo $region | cut -d ':' -f2 | cut -d '-' -f2`
right_flank_stop=$((TR_right + flank_size))

######## Change coordinates to 1-based (samtools faidx) #######

TR_left=$(($TR_left -1))
TR_right=$(($TR_right+1))

######## Obtain sequences #######

left_flank_fasta=`samtools faidx $reference_file $chromosome:$left_flank_start-$TR_left | tr '[:lower:]' '[:upper:]'`
right_flank_fasta=`samtools faidx $reference_file $chromosome:$TR_right-$right_flank_stop | tr '[:lower:]' '[:upper:]'`

TR_rep_fasta=`yes $TR_unit | tr '[:upper:]' '[:lower:]' | head -n $TR_number`
TR_rep_fasta_plus=`yes $TR_unit | tr '[:upper:]' '[:lower:]' | head -n $((TR_number + 2))`

########Generate clean positive and negative sequences########

positive_left_flank=$(echo $left_flank_fasta | cut -d " " -f2- | tr -d ' ')
negative_right_flank=$(echo "$positive_left_flank" | rev | tr "ATGCatgc" "TACGtacg")

positive_right_flank=$(echo $right_flank_fasta | cut -d " " -f2- | tr -d ' ')
negative_left_flank=$(echo "$positive_right_flank" | rev | tr "ATGCatgc" "TACGtacg")

positive_left_flank_plus=$(echo `echo $left_flank_fasta | cut -d " " -f2-`$TR_rep_fasta_plus| tr -d ' ')
negative_right_flank_plus=$(echo "$positive_left_flank_plus" | rev | tr "ATGCatgc" "TACGtacg")

positive_right_flank_plus=$(echo $TR_rep_fasta_plus`echo $right_flank_fasta | cut -d " " -f2-`| tr -d ' ')
negative_left_flank_plus=$(echo "$positive_right_flank_plus" | rev | tr "ATGCatgc" "TACGtacg")

positive_TR_sequence=$(echo $TR_rep_fasta_plus| tr -d ' ')
negative_TR_sequence=$(echo "$positive_TR_sequence" | rev | tr "ATGCatgc" "TACGtacg")

########Write fasta files########
#Left flank
echo ">left_flank_positive;coordinates=$chromosome:$left_flank_start-$TR_left;flank_size=$flank_size" > ${prefix}.fasta
echo $positive_left_flank >> ${prefix}.fasta

echo ">left_flank_negative;coordinates=$chromosome:$left_flank_start-$TR_left;flank_size=$flank_size" >> ${prefix}.fasta
echo $negative_left_flank >> ${prefix}.fasta

#Right flank
echo ">right_flank_positive;coordinates=$chromosome:$TR_right-$right_flank_stop;flank_size=$flank_size" >> ${prefix}.fasta
echo $positive_right_flank >> ${prefix}.fasta

echo ">right_flank_negative;coordinates=$chromosome:$TR_right-$right_flank_stop;flank_size=$flank_size" >> ${prefix}.fasta
echo $negative_right_flank >> ${prefix}.fasta

#Left flank plus TR
echo ">left_flank_positive_plus;coordinates=$chromosome:$left_flank_start-$TR_left;flank_size=$flank_size;TR_number=$TR_number" >> ${prefix}.fasta
echo $positive_left_flank_plus >> ${prefix}.fasta

echo ">left_flank_negative_plus;coordinates=$chromosome:$left_flank_start-$TR_left;flank_size=$flank_size;TR_number=$TR_number" >> ${prefix}.fasta
echo $negative_left_flank_plus >> ${prefix}.fasta

#Right flank plus TR
echo ">right_flank_positive_plus;coordinates=$chromosome:$TR_right-$right_flank_stop;flank_size=$flank_size;TR_number=$TR_number" >> ${prefix}.fasta
echo $positive_right_flank_plus >> ${prefix}.fasta

echo ">right_flank_negative_plus;coordinates=$chromosome:$TR_right-$right_flank_stop;flank_size=$flank_size;TR_number=$TR_number" >> ${prefix}.fasta
echo $negative_right_flank_plus >> ${prefix}.fasta

#TR sequence
echo ">TR_sequence_positive;TR_number=$TR_number;TR_number_plus=$((TR_number + 2))" >> ${prefix}.fasta
echo $positive_TR_sequence >> ${prefix}.fasta

echo ">TR_sequence_negative;TR_number=$TR_number;TR_number_plus=$((TR_number + 2))" >> ${prefix}.fasta
echo $negative_TR_sequence >> ${prefix}.fasta

########Generate squiggles########
scrappie squiggle ${prefix}.fasta > ${prefix}.squiggle

########Homemade index file for easier handling in R########
grep -n "#" ${prefix}.squiggle > ${prefix}.squiggle.index

########Create report if requested########
if [[ $report = "yes" ]]
then
	path=${BASH_SOURCE[0]%/*}
	Rscript $path/scripts/Generated_squiggles_exec.R ${prefix}.squiggle ${prefix} $path 2>&1 >/dev/null
fi
