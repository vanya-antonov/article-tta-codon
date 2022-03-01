
source("lib.R")

###

TYPE_COLORS <- c(
  'num_tta_genes' = 'red',
  'num_tta_fs_genes' = 'blue',
  'num_fshifts' = 'green')
TYPE_LABELS <- c(
  'num_tta_genes' = 'Genes with TTA codon(s)',
  'num_tta_fs_genes' = 'Frameshifted genes with TTA codon(s)',
  'num_fshifts' = 'Frameshifted genes')

###

orgs_df <- read.delim(paste0(DATA_DIR, "all_orgs.tsv"), as.is=TRUE, na.strings = "NULL") %>%
  filter(num_annotated_genes > 0) %>%
  # Adjust two categories
  mutate(
    num_fshifts = num_fshifts - num_tta_fs_genes,
    num_tta_genes = num_tta_genes - num_tta_fs_genes)
head(orgs_df)

get_org_hist <- function(org_df, cur_kingdom)
{
  #cur_kingdom <- 'Bacteria'
  #cur_kingdom <- 'Viruses'
  gg_df <- orgs_df %>% filter(kingdom == cur_kingdom)
  sorted_names_v <- gg_df %>% arrange(-num_tta_genes) %>% pull(name)
  #head(gg_df)
  
  title <- sprintf('Number of genomes (%s) = %d', cur_kingdom, nrow(gg_df))
  subT <- sprintf(
    'TTA-genes = %d, frameshifted genes = %d, frameshifted TTA-genes = %d',
    sum(gg_df$num_tta_genes), sum(gg_df$num_fshifts), sum(gg_df$num_tta_fs_genes))
  gg_df %>%
    select(name, num_fshifts, num_tta_genes, num_tta_fs_genes) %>%
    gather(key = 'type', value = 'num', -name) %>%
    # Define the order of bars
    mutate(name = factor(name, levels=sorted_names_v)) %>%
    ggplot() +
    aes(x = name, y = num, fill = type) +
    geom_bar(stat = "identity", position = 'stack') +
    ylab('Number of genes') +
    xlab('Genomes') +
    scale_fill_manual(
      name = NULL,
      values = TYPE_COLORS,
      breaks = names(TYPE_COLORS),
      labels = TYPE_LABELS) +
    ggtitle(title, subtitle = subT) +
    # Legend items in one column
    guides(fill=guide_legend(ncol=1)) +
    # Remove org names: https://stackoverflow.com/a/35090981/310453
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position='bottom')
}

bact_plot <- get_org_hist(org_df, 'Bacteria')
virus_plot <- get_org_hist(org_df, 'Viruses')

plot_grid(bact_plot, virus_plot, ncol=1, labels = "AUTO", label_size = 18)
ggsave('all_orgs_hist.pdf', path = OUT_DIR, height = 13, width = 12)

