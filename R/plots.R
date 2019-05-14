plotMA <- function(dat, title="", gene=NULL) {
  dat$x <- log10(dat$CTRL + dat$NC)
  dat$y <- log2(dat$NC / dat$CTRL)
  g <- ggplot(dat) +
    theme_classic() +
    geom_point(aes(x, y), size=0.3) +
    labs(x="log10 CTRL + NC", y="log2 NC / CTRL", title=title) +
    geom_hline(yintercept=c(-1,1), colour="red", alpha=0.5)
  if(!is.null(gene)) {
    g <- g + geom_point(data=dat[gene,], aes(x, y), colour="red", size=2)
  }
  g
}

plotXY <- function(dat, title="", gene=NULL) {
  dat$x <- log10(dat$CTRL)
  dat$y <- log10(dat$NC)
  g <- ggplot(dat) +
    theme_classic() +
    geom_point(aes(x, y), size=0.3) +
    labs(x="log10 CTRL", y="log10 NC", title=title)
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