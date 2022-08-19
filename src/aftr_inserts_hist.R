
source("lib.R")

###

input_file <- paste0(DATA_DIR, 'AFTR_INTERGENIC_eA-sB.tsv')
df <- read.table(input_file, sep="\t", header=TRUE, as.is = TRUE)
colnames(df) <- c('strand', 'len', 'num', 'seq')
head(df)

df <- df %>% mutate(seq = paste0(seq, ' (', len, ')'))

seqs_sorted <- df %>% arrange(-num, len) %>% pull(seq)

df %>%
  # Define the order
  mutate(seq = factor(seq, levels = rev(seqs_sorted))) %>%
  ggplot() +
  aes(x = seq, y = num) +
  geom_bar(stat = 'identity') +
  ylab('Number of AFTR genes') +
  xlab('Insert sequence') +
  coord_flip()
ggsave('aftr_inserts_hist.pdf', path = OUT_DIR, height = 10, width = 20)
