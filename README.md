# NanoSatellite
Dynamic time warping of Oxford Nanopore squiggle data to characterize tandem repeats.

## Rationale

Several tools exist to analyze tandem repeats (e.g. [tandem-genotypes](https://github.com/mcfrith/tandem-genotypes) and [RepeatHMM](https://github.com/WGLab/RepeatHMM)). While they do a great job for many tandem repeats in a relatively fast fashion, their quality depends on base calling and alignment. For some tandem repeats (in particular expanded and/or GC-rich tandem repeats), base calling and alignment perform poorly with suboptimal (and DNA strand biased) tandem repeat length and sequence estimations. To overcome these issues, we developed NanoSatellite, a dynamic time warping based algorithm to analyze tandem repeats on raw Oxford Nanopore squiggle data. More information can be found in our [preprint](https://www.biorxiv.org/content/early/2018/10/09/439026).


## UNDER CONSTRUCTION!
* NanoSatellite is tested on the ABCA7 VNTR (a tandem repeat with a 25bp motif) and PromethION data.
* I am updating NanoSatellite to be more user friendly and efficient.


### Dependencies
* [Scrappie](https://github.com/nanoporetech/scrappie "Scrappie") needs to be in path
* Samtools in path
* R and following R packages: 
  * rhdf5
  * dtw
  * ggplot2
  * dplyr
  * tidyr
  * rmarkdown
  * dtwclust
  * knitr
```
#To install these packages, open R and run the following command:
install.packages(c("rhdf5","dtw","ggplot2","dplyr","tidyr","rmarkdown","dtwclust","knitr"))
```

### Usage

#### Generate reference squiggles for the ABCA7 VNTR

```
sh Squiggle_generator.sh -f genome_hg19.fa -r chr19:1049437-1050028 -u GTGAGCCCCCCACCACTCCCTCCCC -p ABCA7_VNTR
```

#### Delineate and segment tandem repeat spanning reads

This script expects the reference squiggles generated above and tab a separated input file (e.g. *spanning_reads.txt*) containing a *name*, *strand* and *path* to fast5. 

Example:

```
name strand	path
02eca7e3-ec99-4ade-bef6-d1a2d0a0ab0f	positive	/storage/fast5/read_7237_ch_2597_strand.fast5
162df1f1-bd5d-4cd1-a558-e6dd1525c6cc	negative	/storage/fast5/read_8421_ch_1457_strand.fast5
40b2e343-cc45-4c7c-9714-bede62e39748	negative	/storage/fast5/read_4311_ch_1669_strand.fast5
```
One possible way to obtain this input file is by first indexing the fast5 containing folder and summary files produced by Albacore or Guppy and then using the temporary *Spanning_read_extractor.sh* script

```
sh fast5_indexing.sh /storage/albacore/summary.tsv /storage/fast5 example

sh Spanning_read_extractor.sh \
/storage/bams/example.bam \
chr19:1049437-1050028 \
example_index.gz

```

Another option is [tandem-genotypes](https://github.com/mcfrith/tandem-genotypes) with the "-v" (verbose) option and then processing the output file.

Run delineation and segmentation:

```
Rscript Signal2chunk.R spanning_reads.txt ABCA7_VNTR.squiggle
```

The output consists of:
* several *.chunk* files, each corresponding to a tandem repeat squiggle unit.
* A *read_metadata.table* containing the number of repeat units per read
* An html report: *chunk_report.html*
