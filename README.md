# NanoSatellite
Dynamic time warping of Oxford Nanopore squiggle data to characterize tandem repeats.

# UNDER CONSTRUCTION!

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

This script expects the reference squiggles generated above and tab separated input file containing a *name*, *strand* and *path to fast5* header.
