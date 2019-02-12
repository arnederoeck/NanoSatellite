set -e

conda install -y -c bioconda scrappie samtools bioconductor-rhdf5 ont-fast5-api
conda install -y -c conda-forge r-dtw r-devtools r-knitr r-ggplot2 r-dplyr r-tidyr r-rmarkdown r-dtwclust
R -e 'Sys.setenv(TAR = "/bin/tar") ; options(unzip = "internal") ; devtools::install_github("arnederoeck/NanoSatelliteR")'
