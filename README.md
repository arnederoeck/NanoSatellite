# NanoSatellite
Dynamic time warping of Oxford Nanopore squiggle data to characterize tandem repeats.

## Rationale

Several tools exist to analyze tandem repeats (e.g. [tandem-genotypes](https://github.com/mcfrith/tandem-genotypes) and [RepeatHMM](https://github.com/WGLab/RepeatHMM)). While they do a great job for many tandem repeats in a relatively fast fashion, their quality depends on base calling and alignment. For some tandem repeats (in particular expanded and/or GC-rich tandem repeats), base calling and alignment perform poorly with suboptimal (and DNA strand biased) tandem repeat length and sequence estimations. To overcome these issues, we developed NanoSatellite, a dynamic time warping based algorithm to analyze tandem repeats on raw Oxford Nanopore squiggle data. The figure below illustrates how NanoSatellite can delineate a tandem repeat in squiggle data and subsequently segment that tandem repeat in seperate units. More information can be found in our [preprint](https://www.biorxiv.org/content/early/2018/10/09/439026).

![NanoSatellite delineation and segmentation](https://github.com/arnederoeck/NanoSatellite/blob/master/raw_positive_squiggle_plot_chunk_colorized20180813.png)

## UNDER CONSTRUCTION!
* NanoSatellite is tested on the ABCA7 VNTR (a tandem repeat with a 25bp motif) and PromethION data.
* I am updating NanoSatellite to be more user friendly and efficient.


## Dependencies
* [Scrappie](https://github.com/nanoporetech/scrappie "Scrappie") (>= 1.3.1) needs to be in path. (installation through [conda](https://anaconda.org/bioconda/scrappie) is also possible)
* Samtools (>= 1.3) in path
* R (>= 3.4.2) and following R packages: 
  * rhdf5 (>= 2.22.0)
  * dtw (>= 1.18-1)
  * ggplot2 (>= 2.2.1)
  * dplyr (>= 0.7.7)
  * tidyr (>= 0.7.2)
  * rmarkdown (>= 1.8)
  * dtwclust (>= 5.5.0)
  * knitr (>= 1.18)
  * devtools (>=2.0.0)
```
#To install these packages, open R and run the following commands:
install.packages(c("dtw","ggplot2","dplyr","tidyr","rmarkdown","dtwclust","knitr","devtools"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rhdf5", version = "3.8")
```
* [NanoSatelliteR](https://github.com/arnederoeck/NanoSatelliteR) (>= 0.1.0) for downstream quality control, plotting, and clustering of NanoSatellite results 
```
#Install in R
devtools::install_github("arnederoeck/NanoSatelliteR")

#If you get an "Installation failed: error in running command" error (sometimes happens when R is installed via conda) use:
options(unzip = "internal")
devtools::install_github("arnederoeck/NanoSatelliteR")
```


## Usage

### Clone this repository

```
git clone https://github.com/arnederoeck/NanoSatellite

# Optionally, make shortcuts to the scripts:
ln -s $PWD/NanoSatellite/{*.sh,*.R} ~/bin
```

### Generate reference squiggles (Squiggle_generator.sh)

```
sh Squiggle_generator.sh [arguments]

Mandatory arguments
    -f Fasta file with your genome of interest (necessary in the algorithm to determine the flanking sequences)
    -r Coordinates of your tandem repeat of interest (e.g. you can obtain these from the Tandem Repeats Finder / Simple Repeats track in the UCSC browser)
    -u Consensus motif of your tandem repeat of interest (e.g. you can obtain these from the Tandem Repeats Finder / Simple Repeats track in the UCSC browser)
    
Additional arguments
    -p Prefix added to the output files
    -e If you add "-e no" then generation of an optional markdown report is suppressed
    -l Change the size of the flanking sequence in the reference squiggles
    -n Change the number of tandem repeat units in the reference squiggles
```

#### ABCA7 VNTR example:
```
sh Squiggle_generator.sh -f genome_hg19.fa -r chr19:1049437-1050028 -u GTGAGCCCCCCACCACTCCCTCCCC -p ABCA7_VNTR
```

### Index fast5 directory (fast5_indexing.sh)
```
sh fast5_indexing.sh <summary_file> <fast5_directory> <prefix>

Example:
sh fast5_indexing.sh /storage/albacore/summary.tsv /storage/fast5 example
```

### Extract tandem repeat spanning reads (Spanning_read_extractor.sh)

Reads spanning the tandem repeat of interest are extracted from an aligned BAM file and coupled to the corresponding fast5 file.
In general, commonly used alignment tools (e.g. [minimap2](https://github.com/lh3/minimap2)) should suffice. However, some tandem repeats are more difficult to align. In such a case, an approach as described [here](https://github.com/mcfrith/last-rna/blob/master/last-long-reads.md) could potentially yield more spanning reads.

```
sh Spanning_read_extractor.sh <in.bam> <region> <index>

Example:
sh Spanning_read_extractor.sh \
/storage/bams/example.bam \
chr19:1049437-1050028 \
example_index.gz
```

### Delineate and segment tandem repeat spanning reads on the squiggle level (Signal2chunk.R)

```
Rscript Signal2chunk.R <spanning_reads file> <reference_squiggles.squiggle>

Example:
Rscript Signal2chunk.R spanning_reads.txt ABCA7_VNTR.squiggle
```

The output consists of:
* several *.chunk* files, each corresponding to a tandem repeat squiggle unit.
* A *read_metadata.table* containing the number of repeat units per read
* An html report: *chunk_report.html*

multiple samples and/or tandem repeats can be run in parallel using external tools.
[GNU Parallel](https://www.gnu.org/software/parallel/) example:

```
#Multiple samples, each with a separate "spanning_reads file" in the "spanning_reads/" directory:
ls spanning_reads/* | parallel 'Rscript Signal2chunk.R {} ABCA7_VNTR.squiggle'
```

### Downstream processing in R with NanoSatelliteR

#### Quality control

```
library(NanoSatelliteR)
chunk_dir="/storage/NanoSatellite_chunks/"
df <- load_summary(chunk_dir)
qc <- summary_qc(df)
```
The output consists of plots displaying normalized "flank" and "center" dynamic time warping distance, respectively corresponding to delineation of tandem repeat squiggles from flanking squiggles, and the segmentation of the tandem repeat squiggle. In addition, cutoff values (red lines in the plots) corresponding to 1.5 times the interquartile range from the 75th percentile are returned.

![NanoSatelliteR QC](https://github.com/arnederoeck/NanoSatellite/blob/master/NanoSatelliteR_qc.png)

#### Tandem repeat length plotting

```
df2 <- qual_reads(df,qc$center_cutoff)
plot_lengths(df2)
```
Each sample is depicted in a separate panel, the number of tandem repeat units is shown on the y-axis, and colored dots correspond to individual sequencing reads originating from positive (red) and negative (blue) DNA strands.

![NanoSatelliteR_plot_lengths](https://github.com/arnederoeck/NanoSatellite/blob/master/NanoSatelliteR_plot_lengths.png)

#### Clustering to identify alternative motifs

```
#Load the tandem repeat unit squiggles
squiggles <- load_squiggles(chunk_dir,df2)

# Clustering with the dtwclust package
#      When clustering a large amount of squiggles, parallelization with doParallel is advised

library(doParallel)
library(dtwclust)

#Number of clusters
k_clusters=2

#Number of cores your computing system has available
registerDoParallel(cores=8)

#The following command clusters the tandem repeat squiggle units originating from the positive DNA strand. The same can be done for the negative strand
positive_clustering <- tsclust(squiggles$positive,type="h",k=k_clusters,trace=TRUE,distance = "dtw_basic", control=hierarchical_control(method="ward.D",symmetric = T))

#To visualy inspect clustering a heatmap can be generated:
ns_heatmap(positive_clustering@distmat,"example.png",max_dist=200,rm0=T)

```

##### Centroid and cluster extraction

```
cent <- extract_centroids(positive_clustering)

library(ggplot2)
ggplot(cent,aes(x=pos,y=signal,colour=factor(cluster)))+geom_point()+geom_line()+theme_minimal()+facet_grid(. ~ cluster)+guides(colour=guide_legend(title="cluster"))
```

![NanoSatelliteR_centroid_plot](https://github.com/arnederoeck/NanoSatellite/blob/master/NanoSatelliteR_centroid_plot.png)
