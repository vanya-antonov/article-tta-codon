#!/usr/bin/env Rscript

source("lib.R")

###

input_file <- paste0(DATA_DIR, 'Clusters_X_with_names.csv')
df <- read.table(input_file, sep="\t", header=TRUE, comment.char="#", as.is = TRUE)
head(df)
# str(df)

# cofID	ORGS	NUM_FS_GENES	NUM_COF_ORGS	NUM_FS_and_TTA_GENES	NUM_TTA_GENES	NUM_TTA_GENE_ORGS	FULL_COF_GENES	FULL_COF_ORGS	FULL_COF_GENES_GBK	COF_NAME
# 1004907	Bacteria	73	35	16	58	27	131	49	66	IS5 family transposase
#...
# 1009299	Bacteria+Phages	2	2	0	1	1	3	3	1	3'-5' exonuclease
# 1009320	Phages	2	2	0	2	2	4	4	2	chitosanase

# NUM_FS_GENES --- кол-во СРС-генов в кластере (это множество СРС-генов называем "ядром", seed)
# NUM_COF_ORGS --- кол-во уникальных организмов, к которым относятся эти СРС-гены из кластера
# NUM_FS_and_TTA_GENES --- кол-во СРС-генов кластера, содержащих и ТТА кодон
# NUM_TTA_GENES --- кол-во ТТА-генов (обычных, не СРС-генов), добавленых к кластеру by BLASTp


#------- 3rd subfigure: COF_genes(X) vs. TTA_genes(Y) vs. ORGs (Z) -------

# cofdata <- subset(df, NUM_TTA_GENES > 0)	# отфильтровываем кластеры
cofdata <- df
cofdata$TTAtoFS <- cofdata$NUM_TTA_GENES / cofdata$NUM_FS_GENES

#svg(file=paste0(OUT_DIR, 'cofs_tta_genes_buble_plot.svg'))

cofdata %>%
  ggplot() +
    aes(x=FULL_COF_GENES, y=NUM_TTA_GENES, size=FULL_COF_ORGS, colour=ORGS) +
#    coord_trans(x="log10") + # , y="log10"  scale_y_continuous(trans='log10') + # scale_y_log10() +
    geom_point(alpha=0.45) +
    # geom_text(
    #   # aes(label=ifelse((TTAtoFS > 14 | NUM_FS_and_TTA_GENES >2 | ORGS != 'Bacteria') & COF_NAME != 'hypothetical protein', as.character(COF_NAME),'')),
    #   position=position_dodge(width=0.9),
    #   vjust=0,
    #   size=5) +	#  hjust=0.5, подписать точки
    scale_size_binned(range = c(.1, 10), name="Number of organisms") +
    xlab("Cluster size") +
    ylab("Number of TTA-containing genes in cluster") +
    ggtitle(paste0("All clusters\nNumber of COFs = ", nrow(cofdata))) + 
    scale_fill_brewer(palette="Set2")
ggsave('cofs_buble_plot.pdf', path = OUT_DIR, height = 13, width = 12)

#dev.off()
