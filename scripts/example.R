library(NanoSatelliteR)
library(dtwclust)


chunk_dir="chunks"
df <- load_summary(chunk_dir)
qc <- summary_qc(df)

df2 <- qual_reads(df,qc$center_cutoff)
plot_lengths(df2)

#Load the tandem repeat unit squiggles
squiggles <- load_squiggles(chunk_dir,df2)

#The following command clusters the tandem repeat squiggle units originating from the positive DNA strand. The same can be done for the negative strand
positive_clustering <- tsclust(squiggles$positive,
                               type="h",
                               k=2,
                               trace=TRUE,
                               distance = "dtw_basic",
                               control=hierarchical_control(method="ward.D",symmetric = T))

#To visualy inspect clustering a heatmap can be generated:
ns_heatmap(positive_clustering@distmat,
           "example.png",
           max_dist=200,
           rm0=T)

cent <- extract_centroids(positive_clustering)

library(ggplot2)
ggplot(cent, aes(x=pos,y=signal,colour=factor(cluster))) +
    geom_point() +
    geom_line() +
    theme_minimal() +
    facet_grid(. ~ cluster) +
    guides(colour=guide_legend(title="cluster"))

cpr <- clusters_per_read(positive_clustering)
