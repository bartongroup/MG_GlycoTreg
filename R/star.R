parseStarLog <- function(dir, sample) {
  file <- paste0(dir, sample, "_Log.final.out")
  fs <- readLines(file)
  res <- list()
  for(f in fs) {
    if(grepl("|", f, fixed=TRUE)) {
      s <- unlist(strsplit(f, "\t",  perl=TRUE))
      key <- s[1]
      val <- s[2]
      key <- gsub(" \\|", "", key)
      key <- gsub("^\\s+", "", key, perl=TRUE)
      res[key] <- val
    }
  }
  res <- data.frame(sample=unlist(res))
  names(res) <- sample
  res
}

parseStarCount <- function(dir, sample, column=2, suffix=".txt") {
  file <- paste0(dir, sample, suffix)
  stopifnot(file.exists(file))
  f <- read.table(file, sep="\t", header=FALSE, skip=4)
  # second column - unstranded data
  res <- data.frame(sample=f[, column])
  rownames(res) <- f$V1
  names(res) <- sample
  res
}

readBedgraphs <- function(dir, metadata, windowSize=100000) {
  P <- lapply(1:nrow(metadata), function(i) {
    file <- paste0(dir, metadata$name[i], ".", sprintf("%d", windowSize), ".bedgraph")
    b <- setNames(read.table(file, header=FALSE, sep="\t"), c("chr", "start", "end", "score"))
    b <- cbind(sample=metadata$sample[i], b, cumpos=(seq_along(b$start) - 0.5) * windowSize / 1e6)
  })
  bg <- do.call(rbind, P)
}

plotStarSense <- function(file) {
  
  oneplot <- function(d, col1, col2) {
    d[[col1]] <- log10(d[[col1]])
    d[[col2]] <- log10(d[[col2]])
    ggplot(d, aes_string(col1, col2)) +
      theme_classic() +
      geom_point(size=0.1) +
      geom_abline(slope=1, intercept=0, colour="red") +
      labs(x=col1, y=col2)
  }
  
  stopifnot(file.exists(file))
  f <- read.table(file, sep="\t", header=FALSE, skip=4)
  names(f) <- c("gene", "Unstranded", "First", "Second")
  g1 <- oneplot(f, "Unstranded", "First")
  g2 <- oneplot(f, "Unstranded", "Second")
  g3 <- oneplot(f, "First", "Second")
  
  plot_grid(g1, g2, g3, ncol=3)
}

