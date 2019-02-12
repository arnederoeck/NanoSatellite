set -e

wget -O- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz | gunzip -c > chr19_hg19.fa

bash Squiggle_generator.sh -f chr19_hg19.fa -r chr19:1049437-1050028 -u GTGAGCCCCCCACCACTCCCTCCCC -p ABCA7_VNTR
bash no_summary_indices.sh test-NA19240-fast5 example

# bash Spanning_read_extractor.sh \
#  /storage/bams/example.bam \
#  chr19:1049437-1050028 \
#  example_fast5_index.gz

Rscript Signal2chunk.R  test-NA19240-fast5/test-spanning-reads.txt ABCA7_VNTR.squiggle
Rscript scripts/example.R
