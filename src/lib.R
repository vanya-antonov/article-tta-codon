
# install.packages('cowplot')

library(ggplot2)
library(cowplot)   # To have ggplots side-by-side: plot_grid()

library(dplyr)
library(tidyr)     # separate(), gather() and spread() functions
library(tibble)    # for rownames_to_column() and column_to_rownames()

# library(ComplexHeatmap)
# library(circlize)   # colorRamp2()


###

DATA_DIR <- "../data/"
OUT_DIR <- "../images/"

# FS_COLORS <- c('Yes' = 'black', 'No' = 'white')
# EVALUE_COLORS <- colorRamp2(c(0, 6, 100), c("white", "yellow", "red"))
# EVALUE_COLORS_BW <- colorRamp2(c(0, 6, 100), c("white", "gray", "black"))

theme_set(theme_bw(base_size = 19))  # increase the font size: https://stackoverflow.com/a/11955412/310453

