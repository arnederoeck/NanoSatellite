# NanoSatellite
Dynamic time warping of Oxford Nanopore squiggle data to characterize tandem repeats.

## UNDER CONSTRUCTION!
* NanoSatellite is tested on the ABCA7 VNTR (a tandem repeat with a 25bp motif) and PromethION data.
* The scripts are being updated to become more user friendly and efficient.

### Dependencies
* [Scrappie](https://github.com/nanoporetech/scrappie "Scrappie") needs to be in path
* R and following R packages:
  * rhdf5
  * dtw
  * ggplot2
  * dplyr
  * tidyr
  * rmarkdown
  * dtwclust
  * knitr

### Usage

#### Generate reference squiggles for the ABCA7 VNTR

```
sh Squiggle_generator.sh -f genome_hg19.fa -r chr19:1049437-1050028 -u GTGAGCCCCCCACCACTCCCTCCCC -p ABCA7_VNTR
```

#### Delineate and segment tandem repeat spanning reads

This script expects the reference squiggles generated above and tab a separated input file (e.g. *spanning_reads.txt*) containing a *name*, *strand* and *path* to fast5. 
One possible way to obtain this input file is by using [tandem-genotypes](https://github.com/mcfrith/tandem-genotypes) with the "-v" (verbose) option and then processing the output file.
Another is to 

Example:

```
name strand	path
02eca7e3-ec99-4ade-bef6-d1a2d0a0ab0f	positive	/storage/fast5/read_7237_ch_2597_strand.fast5
162df1f1-bd5d-4cd1-a558-e6dd1525c6cc	negative	/storage/fast5/read_8421_ch_1457_strand.fast5
40b2e343-cc45-4c7c-9714-bede62e39748	negative	/storage/fast5/read_4311_ch_1669_strand.fast5
```

Run Singal2chunk.R:

```
Rscript Signal2chunk.R spanning_reads.txt ABCA7_VNTR_squiggle_analysis.squiggle
```
