
# source("lib.R")

###

dat <- data.frame(
  "internal_TTA_yes" = c(177, 7),
  "internal_TTA_no" = c(19, 15),
  row.names = c("uORF_TTA_yes", "uORF_TTA_no"),
  stringsAsFactors = FALSE
)
# colnames(dat) <- c("Non-smoker", "Smoker")
dat

test <- fisher.test(dat)
test