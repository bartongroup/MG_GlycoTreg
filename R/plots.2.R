funCondition <- function(dat, metadata, FUN) {
  P <- lapply(unique(metadata$condition), function(cond) {
    sel <- metadata[metadata$condition == cond, "sample"]
    setNames(data.frame(apply(dat[, sel], 1, function(x) FUN(x, na.rm=TRUE))), cond)
  })
  do.call(cbind, P)
}


plotMV <- function(dat, metadata, title="", mid.gradient=0.3, bins=80) {
  P <- lapply(unique(metadata$condition), function(cond) {
    sel <- metadata[metadata$condition == cond, "sample"]
    ds <- dat[, sel]
    d <- data.frame(
      condition = cond,
      M = apply(ds, 1, function(x) mean(x, na.rm=TRUE)),
      V = apply(ds, 1, function(x) sd(x, na.rm=TRUE)^2)
    )
  })
  d <- do.call(rbind, P)
  
  ggplot(d, aes(log10(M), log10(V))) + 
    labs(x="log10 Mean", y="log10 Variance", title=title) +
    facet_grid(. ~ condition) +
    stat_binhex(bins=bins) +
    scale_fill_gradientn(colours=c("seagreen","yellow", "red"), values=c(0, mid.gradient, 1), name="count", na.value=NA) +
    geom_abline(slope=1, intercept=0, colour="red")
    
}
  

plotMA <- function(M, conditions, title="", gene=NULL) {
  stopifnot(length(conditions) == 2)
  stopifnot(all(conditions %in% names(M)))
  M <- M[, conditions]
  M$x <- log10(M[[1]] + M[[2]])
  M$y <- log2(M[[2]] / M[[1]])
  g <- ggplot(M) +
    theme_classic() +
    geom_point(aes(x, y), size=0.3) +
    labs(
      x = paste("log10", paste(conditions, collapse=" + ")),
      y = paste("log2", paste(conditions, collapse=" / ")),
      title=title
    ) +
    geom_hline(yintercept=c(-1,1), colour="red", alpha=0.5)
  if(!is.null(gene)) {
    g <- g + geom_point(data=M[gene,], aes(x, y), colour="red", size=2)
  }
  g
}



plotManhattan <- function(bg) {
  #bg$chr <- as.factor(bg$chr)
  
  # colour map of alternating chromosomes
  r <- rle(as.character(bg$chr))$values
  cmap <- setNames(seq_along(r) %% 2, r)
  bg$fill <- as.factor(cmap[as.character(bg$chr)])
  
  width <- bg$cumpos[2] - bg$cumpos[1]
  ggplot(bg) +
    theme_classic() +
    geom_col(aes(x=cumpos, y=score, fill=fill), width=width) +
    facet_grid(sample ~ .) +
    labs(x="Position (Mb)", y="Count") +
    theme(legend.position = "none") +
    scale_fill_manual(values=c("#E69F00", "#56B4E9"))
}