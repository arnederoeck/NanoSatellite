set -e

wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR286/ERR2864496/NA19240.tar.gz && tar -xzf NA19240.tar.gz && rm NA19240.tar.gz
wget -O- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz | gunzip -c > chr19_hg19.fa

bash Squiggle_generator.sh -f chr19_hg19.fa -r chr19:1049437-1050028 -u GTGAGCCCCCCACCACTCCCTCCCC -p ABCA7_VNTR
bash no_summary_indices.sh NA19240 example

# bash Spanning_read_extractor.sh \
#  /storage/bams/example.bam \
#  chr19:1049437-1050028 \
#  example_fast5_index.gz

sed "s|your_fast5_directory/||" ABCA7_VNTR_example/NA19240_spanning_reads.txt | head -n 5 > test_spanning_reads.txt
Rscript Signal2chunk.R  test_spanning_reads.txt ABCA7_VNTR.squiggle
Rscript scripts/example.R
