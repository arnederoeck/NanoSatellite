# NanoSatellite
Dynamic time warping of Oxford Nanopore squiggle data to characterize tandem repeats.

## Rationale

Several tools exist to analyze tandem repeats (e.g. [tandem-genotypes](https://github.com/mcfrith/tandem-genotypes) and [RepeatHMM](https://github.com/WGLab/RepeatHMM)). While they do a great job for many tandem repeats in a relatively fast fashion, their quality depends on base calling and alignment. For some tandem repeats (in particular expanded and/or GC-rich tandem repeats), base calling and alignment perform poorly with suboptimal (and DNA strand biased) tandem repeat length and sequence estimations. To overcome these issues, we developed NanoSatellite, a dynamic time warping based algorithm to analyze tandem repeats on raw Oxford Nanopore squiggle data. The figure below illustrates how NanoSatellite can delineate a tandem repeat in squiggle data and subsequently segment that tandem repeat in seperate units. More information can be found in our [preprint](https://www.biorxiv.org/content/early/2018/10/09/439026).

![NanoSatellite delineation and segmentation](https://github.com/arnederoeck/NanoSatellite/blob/master/figures/raw_positive_squiggle_plot_chunk_colorized20180813.png)

Dynamic time warping (embedded in NanoSatellite) requires a minimum "pattern size" for good pattern recognition. Hence for tandem repeats with a very short tandem repeat motif, or low numbers of repeating motifs NanoSatellite will produce suboptimal results.

## Citation
If this tool is useful for your work, please consider citing our [preprint](https://www.biorxiv.org/content/early/2018/10/09/439026).

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

### *ABCA7* VNTR PromethION example data

The NanoSatellite algorithm was originally tested on fast5 reads spanning the *ABCA7* VNTR, originating from whole genome sequencing from 11 individuals on the Oxford Nanopore PromethION platform as described [here](https://www.biorxiv.org/content/early/2018/10/09/439026). The *ABCA7* VNTR fast5's are publicly accessible from [ENA](https://www.ebi.ac.uk/ena/data/view/PRJEB29458). These can be downloaded and unpacked from command line with the following code:

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA164/ERA1642169/oxfordnanopore_native/*
for f in `ls *.tar.gz`; do tar -xzvf $f; done
```

The corresponding annotation files (&ast;_spanning_reads.txt) can be found in NanoSatellite/ABCA7_VNTR_example. These files replace steps 2. and 3. below. If you'd like to use this example data, it's important to set the path to the directory where the fast5 files were unpacked. The following code can be used:

```
cd NanoSatellite/ABCA7_VNTR_example/
sed -i "s|your_fast5_directory|/your/fast5/directory|" *_spanning_reads.txt
```



### 1. Generate reference squiggles (Squiggle_generator.sh)

```
bash Squiggle_generator.sh [arguments]

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

#### *ABCA7* VNTR example:
```
bash Squiggle_generator.sh -f genome_hg19.fa -r chr19:1049437-1050028 -u GTGAGCCCCCCACCACTCCCTCCCC -p ABCA7_VNTR
```

### 2. Index fast5 directory (fast5_indexing.sh)
```
bash fast5_indexing.sh <summary_file> <fast5_directory> <prefix>

Example:
bash fast5_indexing.sh /storage/albacore/summary.tsv /storage/fast5 example
```

### 3. Extract tandem repeat spanning reads (Spanning_read_extractor.sh)

Reads spanning the tandem repeat of interest are extracted from an aligned BAM file and coupled to the corresponding fast5 file.
In general, commonly used alignment tools (e.g. [minimap2](https://github.com/lh3/minimap2)) should suffice. However, some tandem repeats are more difficult to align. In such a case, an approach as described [here](https://github.com/mcfrith/last-rna/blob/master/last-long-reads.md) could potentially yield more spanning reads.

```
bash Spanning_read_extractor.sh <in.bam> <region> <index>

Example:
bash Spanning_read_extractor.sh \
/storage/bams/example.bam \
chr19:1049437-1050028 \
example_index.gz
```

### 4. Delineate and segment tandem repeat spanning reads on the squiggle level (Signal2chunk.R)

#### 4.1 The simplest way (needs a functioning rhdf5 installation in R)

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
#### 4.2 Alternative approach in case of rhdf5 error using approach 4.1

Some versions of R and rhdf5 don't work well together and can result in an HDF5 error when running `Signal2chunk.R`. If this occurs we made a workaround using `fast5_extract.py` (with thanks to [Wouter](https://github.com/wdecoster)), which depends on [ont_fast5_api](https://github.com/nanoporetech/ont_fast5_api). `fast5_extract.py` extracts the squiggle data from fast5, which can be written away to a file that can subsequently be supplied to `Signal2chunk.R` as the third argument.

usage:

```
fast5_extract.py [-h] [-d DIR] [-f [FILE [FILE ...]]] [-r]

Extract signal level from fast5 files

optional arguments:
  -h, --help            show this help message and exit
  -d DIR, --dir DIR     directory with fast5 file(s) to extract signal from
  -f [FILE [FILE ...]], --file [FILE [FILE ...]]
                        fast5 file(s) to extract signal from
  -r, --recursive       recursively go through directories
```

Either supply a directory of fast5 files using -d/--dir, searched optionally recursive with -r/--recursive OR supply on or more fast5 files with -f/--file

Example:

```
python fast5_extract.py -r -d /storage/fast5/ > fast5_pre-extracted.txt

Rscript Signal2chunk.R spanning_reads.txt ABCA7_VNTR.squiggle fast5_pre-extracted.txt
```



### 5. Downstream processing in R with NanoSatelliteR

#### 5.1. Quality control

```
library(NanoSatelliteR)
chunk_dir="/storage/NanoSatellite_chunks/"
df <- load_summary(chunk_dir)
qc <- summary_qc(df)
```
The output consists of plots displaying normalized "flank" and "center" dynamic time warping distance, respectively corresponding to delineation of tandem repeat squiggles from flanking squiggles, and the segmentation of the tandem repeat squiggle. In addition, cutoff values (red lines in the plots) corresponding to 1.5 times the interquartile range from the 75th percentile are returned.

![NanoSatelliteR QC](https://github.com/arnederoeck/NanoSatellite/blob/master/figures/NanoSatelliteR_qc.png)

#### 5.2. Tandem repeat length plotting

```
df2 <- qual_reads(df,qc$center_cutoff)
plot_lengths(df2)
```
Each sample is depicted in a separate panel, the number of tandem repeat units is shown on the y-axis, and colored dots correspond to individual sequencing reads originating from positive (red) and negative (blue) DNA strands.

![NanoSatelliteR_plot_lengths](https://github.com/arnederoeck/NanoSatellite/blob/master/figures/NanoSatelliteR_plot_lengths.png)

#### 5.3. Clustering to identify alternative motifs

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

Clustering in more than 2 groups is also possible. However, by increasing the number of clusters, differences between clusters become smaller which can impair accuracy. Hence supervision of clustering by the user is warranted. Re-clustering can be done by re-running `tsclust()` and supplying the previously generated distance matrix to the `distmat` paramater in `hierarchical_control()`.


##### 5.4. Centroid extraction

```
cent <- extract_centroids(positive_clustering)

library(ggplot2)
ggplot(cent,aes(x=pos,y=signal,colour=factor(cluster)))+geom_point()+geom_line()+theme_minimal()+facet_grid(. ~ cluster)+guides(colour=guide_legend(title="cluster"))
```

![NanoSatelliteR_centroid_plot](https://github.com/arnederoeck/NanoSatellite/blob/master/figures/NanoSatelliteR_centroid_plot.png)

For each cluster, a centroid is assigned. In this particular case of positive *ABCA7* VNTR tandem repeat unit squiggles, the differences observed in cluster 1, are caused by a guanine insertion, or cytosine to adenine substitution at nucleotide ten of the *ABCA7* VNTR consensus motif (inferred from comparison to Scrappie reference squiggles).

##### 5.5. Cluster sequence per sequencing read

```
cpr <- clusters_per_read(positive_clustering)
```

This provides a data.frame with a comma-separated string of clusters per sequencing read, ordered from start to end according to the original sequencing read squiggle.
