
source("lib.R")

###

input_file <- paste0(DATA_DIR, 'Statistics_clusters_with_TTA.txt')
df <- read.table(input_file, sep="\t", header=TRUE, as.is = TRUE)
head(df)
tail(df)


df %>%
  mutate(
    log_pvalue = -log10(p_value_min),
    TTA_peak_perc = 100*TTA_peak.num_TTA/Number_all_TTA_genes) %>%
  ggplot() +
  aes(x = Cluster_size, y = log_pvalue, size = TTA_peak_perc) %>%
  geom_point(alpha=0.45) +
  theme(legend.position="bottom") +
  # Rename legend: http://www.cookbook-r.com/Graphs/Legends_(ggplot2)/
  scale_size_continuous(name  ="% aligned TTA codons: ") +
  xlab("Cluster size") +
  ylab("-log10(min SynPlot2 p-value)") +
  ggtitle(paste0("Total number of clusters = ", nrow(df)))
ggsave('cofs_buble_plot.pdf', path = OUT_DIR, height = 10, width = 10)
