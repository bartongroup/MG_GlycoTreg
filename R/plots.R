cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

funCondition <- function(dat, metadata, FUN) {
  P <- lapply(unique(metadata$condition), function(cond) {
    sel <- metadata[metadata$condition == cond, "name"]
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



plotManhattan <- function(bg, ymax=NA, yscale=1) {
  #bg$chr <- as.factor(bg$chr)
  
  # colour map of alternating chromosomes
  r <- rle(as.character(bg$chr))$values
  cmap <- setNames(seq_along(r) %% 2, r)
  bg$colour <- as.factor(cmap[as.character(bg$chr)])
  
  ylab <- "Count"
  if(yscale != 1) {
    bg$score <- bg$score / yscale
    ylab <- paste("Count x", yscale)
  }
  
  ggplot(bg) +
    theme_classic() +
    geom_segment(aes(x=cumpos, xend=cumpos, y=0, yend=score, colour=colour)) +
    facet_grid(sample ~ .) +
    labs(x="Position (Mb)", y=ylab) +
    theme(legend.position = "none") +
    scale_colour_manual(values=c("#E69F00", "#56B4E9")) +
    scale_y_continuous(limits=c(0, ymax), expand=c(0,0)) +
    theme(strip.background = element_blank())
}


plotDistanceMatrix <- function(cnt, metadata, distance=c("correlation"), text.size=10) {
  distance <- match.arg(distance)
  
  corr.mat <- cor(cnt, use="complete.obs")
  m <- reshape2::melt(corr.mat, varnames=c("Sample1", "Sample2"))
  m$Sample1 <- factor(m$Sample1, levels=metadata$name)
  m$Sample2 <- factor(m$Sample2, levels=metadata$name)
  ggplot(m, aes_(x=~Sample1, y=~Sample2)) +
    geom_tile(aes_(fill=~value)) +
    viridis::scale_fill_viridis(option="cividis") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=text.size),
      axis.text.y = element_text(size=text.size)
    ) +
    labs(x='Sample', y='Sample', fill="Correlation")
}


plotClustering <- function(cnt, x.text.size=10) {
  corr.mat <- cor(cnt, use="complete.obs")
  dis <- as.dist(1 - corr.mat)  # dissimilarity matrix
  hc <- hclust(dis)
  dendr <- ggdendro::dendro_data(hc)
  dat <- ggdendro::segment(dendr)
  theme.d <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust=0.5, size=x.text.size),
    axis.line.x = ggplot2::element_blank(),
    axis.line.y = ggplot2::element_line(size=0.5),
    axis.ticks.x = ggplot2::element_blank()
  )
  ggplot() +
    theme.d +
    geom_segment(data=dat, aes_(x=~x, y=~y, xend=~xend, yend=~yend)) +
    scale_x_continuous(breaks = seq_along(dendr$labels$label), labels = dendr$labels$label) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(dat$y) * 1.03)) +
    labs(x="Sample", y="Distance")
}


plotPCA <- function(rc, metadata, with.legend=TRUE, text.size=12, legend.size=10) {
  # remove rows with zero variance
  rc <- rc[apply(rc, 1, sd) > 0, ]
  pca <- prcomp(t(rc), scale.=TRUE, center=TRUE)
  
  p <- data.frame(x = pca$x[, 1], y = pca$x[, 2], label=metadata$replicate, condition=metadata$condition)
  var.perc <- 100 * (pca$sdev)^2 / sum((pca$sdev)^2)
  pca1 <- sprintf("PCA1 (%5.1f%%)", var.perc[1])
  pca2 <- sprintf("PCA2 (%5.1f%%)", var.perc[2])
  
  g <- ggplot(p, aes(x, y)) + 
    theme_classic() +
    #coord_cartesian(xlim=exlim(p$x), ylim=exlim(p$y)) +
    geom_point(aes(colour=condition, shape=condition), size=2) +
    geom_text(aes(label=label), nudge_x=5, nudge_y=4, size=2.5, colour="grey40") +
    #scale_color_viridis(discrete=TRUE, option="cividis") +
    scale_colour_manual(values=cbPalette) +
    theme(text = element_text(size=text.size)) +
    labs(x=pca1, y=pca2)
  #theme(legend.key.size = unit(0.2, "cm"), legend.text = element_text(size=legend.size))
  if(!with.legend) g <- g + theme(legend.position="none")
  g
}