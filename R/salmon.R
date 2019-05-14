parseSalmonCount <- function(dir, sample, column="TPM") {
  file <- paste0(dir, sample, "/quant.sf")
  stopifnot(file.exists(file))
  f <- read.table(file, sep="\t", header=TRUE)
  res <- data.frame(sample=f[, column])
  rownames(res) <- f$Name
  names(res) <- sample
  res
}
