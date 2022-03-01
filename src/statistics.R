
source("lib.R")

###

# "First, we analyzed all the genes annotated in 240 phage and 98 Streptomyces genomes"
orgs_df <- read.delim(paste0(DATA_DIR, "all_orgs.tsv"), as.is=TRUE, na.strings = "NULL") %>%
  filter(num_annotated_genes > 0 )
orgs_df %>% filter(kingdom == 'Bacteria') %>% nrow()
orgs_df %>% filter(kingdom == 'Viruses') %>% nrow()

# "In total, among all the 759,420 annotated genes there were 21,487 genes (2.8\%) with TTA codon(s)."
sum(orgs_df$num_annotated_genes)
sum(orgs_df$num_tta_genes)
100 * sum(orgs_df$num_tta_genes) / sum(orgs_df$num_annotated_genes)

# "Interestingly, the fraction of the TTA-containing genes was higher in the bacterial phage rather than in the viral genomes (2.86\% vs 1.91\%, respectively)."
bacteria_df <- orgs_df %>% filter(kingdom == 'Bacteria')
virus_df <- orgs_df %>% filter(kingdom == 'Viruses')
100 * sum(bacteria_df$num_tta_genes) / sum(bacteria_df$num_annotated_genes)
100 * sum(virus_df$num_tta_genes) / sum(virus_df$num_annotated_genes)
