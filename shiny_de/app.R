libDir <- "/cluster/gjb_lab/mgierlinski/R_shiny/library/3.6/"
if(dir.exists(libDir)) .libPaths(libDir)

dirs <- c("/home/mgierlinski/projects/glycotreg", "/cluster/gjb_lab/mgierlinski/projects/glycotreg")
for(d in dirs){
  if(dir.exists(d)) topDir <- d
}

library(shiny)
library(ggplot2)
library(dplyr)
library(DT)
library(gplots)

css <- "table{font-size: 11px; background-color: #EAF5FF}"

### Read data ###

sfile <- file.path(topDir, "RData/shiny_data.RData")
load(sfile)
pairs <- names(DE)
pair <- pairs[1]
expr <- expr[DE[[pair]]$gene_id, ] # make sure the order of genes is the same
annotation <- DE[[pair]][, c("gene_id", "gene_name", "description")]
ens2name <- setNames(annotation$gene_name, annotation$gene_id)

DE <- lapply(DE, function(tab) {
  tab$logp <- -log10(tab$PValue)
  tab$logCPM <- tab$logCPM / log2(10)
  tab$sel <- tab$gene_id %in% glyco.ens
  tab
})

names(gene.go)[2] <- "term_id"
names(go.terms)[1] <- "term_id"
go.terms$go_domain <- NULL
names(gene.reactome)[2] <- "term_id"
names(reactome.terms)[1] <- "term_id"

max_points <- 500

#######################################################################


plotVolcano <- function(tab, alpha=0.05, title="") {
  lim <- -log10(max(tab[tab$FDR < alpha, "PValue"]))
  ggplot(tab, aes(x, y)) +
    theme_classic() +
    theme(text = element_text(size=18)) +
    geom_point(size=0.3, colour="grey60") +
    geom_point(data=tab[tab$sel, ], size=0.6, colour="black") +
    geom_vline(xintercept = 0, colour="grey70") +
    geom_hline(yintercept = lim, colour="red", linetype="dashed") +
    labs(x="log FC", y="log10 P", title=title) +
    scale_y_continuous(expand=c(0,0), limits=c(0, max(tab$y)*1.03))
}

plotMA <- function(tab, alpha=0.05, title="") {
  ggplot(tab, aes(x, y)) +
    theme_classic() +
    theme(text = element_text(size=18)) +
    geom_point(size=0.3, colour="grey60") +
    geom_point(data=tab[tab$sel, ], size=0.6, colour="black") +
    geom_hline(yintercept = 0, colour="grey70") +
    labs(x="log CPM", y="log FC", title=title) 
    #scale_y_continuous(expand=c(0,0), limits=c(0, max(tab$y)*1.03))
}

allGeneTable <- function(tab) {
  cols <- c("gene_id", "gene_name", "logFC", "FDR", "description")
  d <- tab[, cols]
  d[, 3:4] <- sapply(d[, 3:4], function(x) signif(x, 3))
  d <- DT::datatable(d, class = 'cell-border strip hover')
  DT::formatStyle(d, 0, cursor = 'pointer')
}

jitterPlot <- function(dat, scale, max_points=100) {
  n <- nrow(dat)
  if(n == 0) return(NULL)
  
  if(scale == "log"){
    dat <- log10(dat)
  }
  m <- colMeans(dat, na.rm = TRUE)
  s <- apply(dat, 2, function(x) sd(x, na.rm = TRUE) / sqrt(na.omit(length(x))))
  p <- data.frame(
    intensity = m,
    lo = m - s,
    up = m + s,
    condition = metadata$condition,
    replicate = as.factor(metadata$replicate)
  )
  p$shape <- rep(21, length(p$intensity))
  p$shape[which(p$intensity==0)] <- 24
  pd <- ggplot2::position_dodge(width = 0.3)
  
  th <- theme_classic() + theme(
    text = element_text(size=18),
    axis.text.x = element_text(angle=90, hjust=0.5),
    panel.grid.major.x = element_line(colour="grey70"),
    legend.position = "none"
  )
  
  if(n == 1) {
    mx <- max(dat) * 1.03
  } else {
    mx <- max(p$up) * 1.03
  }
  
  g <- ggplot(p, aes(x=condition, y=intensity, ymin=lo, ymax=up, colour=replicate, shape=shape, fill=replicate)) +
    th +
    geom_point(position=pd, size=4) +
    scale_shape_identity() +  # necessary for shape mapping
    viridis::scale_fill_viridis(discrete=TRUE, option="cividis") +
    xlab("")
  if(scale == "log") {
    g <- g + labs(y = "Log Intensity")
  } else {
    g <- g +
      labs(y = "Normalized count") +
      scale_y_continuous(limits=c(0, mx), expand=c(0, 0))
  }
  if(n > 1) {
    g <- g + geom_errorbar(position=pd, width = 0.1)
  }
  g
}

# Enrichment function
#
# Performs hypergeometric test.
# Input: all.genes and sel.genes are vectors with gene IDs.
# gene2term: data frame with columns gene_id and term_id
# term.info: data frame with columns term_id and other columns with name, description and so on. These additional columns will be inluded in the output.

functionalEnrichment <- function(all.genes, sel.genes, gene2term, term.info, min.count=3, sig.limit=0.05) {
  
  # select only terms represented in our gene set
  gene2term <- gene2term[gene2term$gene_id %in% all.genes,]
  
  # all terms present in the selection
  terms <- unique(gene2term[gene2term$gene_id %in% sel.genes,]$term_id)
  terms <- terms[terms != ""]  # remove empty
  
  # number of selected genes
  Nsel <- length(sel.genes)
  # size of the universe
  Nuni <- length(all.genes)
  
  # empty line for missing terms
  na.term <- term.info[1, ]
  na.term[, !(names(na.term) == "term_id")] <- NA
  
  P <- lapply(terms, function(term) {
    info <- term.info[term.info$term_id == term, ]
    if(nrow(info) == 0) {
      na.term$term_id <- term
      info <- na.term
    }
    
    term.uni <- gene2term[gene2term$term_id == term, ]$gene_id
    term.sel <- term.uni[term.uni %in% sel.genes]
    # number of selected genes with GO-term
    nsel <- length(term.sel)
    # all genes with GO-term
    nuni <- length(term.uni)
    expected <- nuni * Nsel / Nuni
    fish <- matrix(c(nsel, nuni - nsel, Nsel - nsel, Nuni + nsel - Nsel - nuni), nrow=2)
    ft <- fisher.test(fish, alternative = 'greater')
    p <- as.numeric(ft$p.value)
    
    row <- cbind(
      info,
      data.frame(
        tot = nuni,
        sel = nsel,
        expect = expected,
        enrich = nsel / expected,
        genes = paste(ens2name[term.sel], collapse=","),
        P = p
      )
    )
  })
  tab <- do.call(rbind, P)
  tab$P <- p.adjust(tab$P, method="BH")
  tab <- tab[tab$sel >= min.count & tab$P <= sig.limit, ]
  tab <- tab[order(-tab$enrich), ]

  tab$enrich <- round(tab$enrich, 1)
  tab$expect <- round(tab$expect, 2)
  tab$P <- NULL
  tab
}



#######################################################################

ui <- shinyUI(fluidPage(
  
  tags$style(css),
  
  titlePanel("Differential expression rest vs active"),

  fluidRow(
    column(7,
      fluidRow(
        column(8, 
          #radioButtons("condPair", "Condition pair:", choices = pairs, inline = TRUE),
          radioButtons("selGenes", "Selection:", choices = c("All", "Glyco"), inline = TRUE),
          radioButtons("plotType", "Plot type:", choices = c("Volcano", "MA"), inline=TRUE),
          plotOutput("mainPlot", height="480px", width="100%", brush="plot_brush", hover="plot_hover")
        ),
        column(4,
          radioButtons("intensityScale", "Intesity scale:", choices = c("lin" = "", "log"="log"), inline = TRUE),
          plotOutput("jitterPlot", height = "400px",width = "100%")
        )
      ),
      fluidRow(
        DT::dataTableOutput("allGeneTable")
      )
    ),
#    column(2, 
#      radioButtons("intensityScale", "Intesity scale:", choices = c("lin" = "", "log"="log"), inline = TRUE),
#      plotOutput("jitterPlot", height = "400px",width = "100%")
#    ),
    column(5,
      p("Gene list"),
      div(style = 'height: 250px; overflow-y: scroll', tableOutput("geneInfo")),
      br(),
      p("GO enrichment"),
      div(style = 'height: 250px; overflow-y: scroll', tableOutput("GOEnrichment")),
      br(),
      p("Pathway enrichment"),
      div(style = 'height: 250px; overflow-y: scroll', tableOutput("ReactomeEnrichment"))
    )
  )

  # Show main gene table
  #fluidRow(
  #  column(width = 12,
  #         DT::dataTableOutput("allGeneTable"))
  #)
)
)


########################################################################################


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  getData <- function() {
    #pair <- input$condPair
    #tab <- DE[[pair]]
    tab <- DE[[1]]
    if(input$selGenes == "Glyco") {
      sel <- tab$sel
      tab <- tab[sel, ]
    }
    if(input$plotType == "Volcano") {
      tab$x <- tab$logFC
      tab$y <- tab$logp
    } else {
      tab$x <- tab$logCPM
      tab$y <- tab$logFC
    }
    tab
  }
  
  selectGene <- function(max_hover=1) {
    tab <- getData()
    sel <- NULL
    tab_idx <- as.numeric(input$allGeneTable_rows_selected)
    if(!is.null(input$plot_brush)){
      brushed <- na.omit(brushedPoints(tab, input$plot_brush))
      sel <- as.numeric(rownames(brushed))
    } else if(!is.null(input$plot_hover)) {
      near <- nearPoints(tab, input$plot_hover, threshold = 20, maxpoints = max_hover)
      sel <- as.numeric(rownames(near))
    } else if(length(tab_idx) > 0) {
      sel <- tab_idx
    }
    return(sel)  # for rowname selection
  }
  
  output$geneInfo <- renderTable({
    sel <- selectGene()
    if (!is.null(sel) && length(sel) >= 1 && length(sel) <= max_points) {
      df <- data.frame(annotation[sel,])
      df[order(df$gene_name), ]
    } else if (length(sel) > max_points) {
      data.frame(Error=paste0('only ',max_points,' points can be selected.'))
    }
  })
  
  output$GOEnrichment <- renderTable({
    tab <- getData()
    sel <- NULL
    fe <- NULL
    if(!is.null(input$plot_brush)){
      brushed <- na.omit(brushedPoints(tab, input$plot_brush))
      sel <- as.numeric(rownames(brushed))
      n <- length(sel)
      if(n > 0 && n <= max_points) {
        all.genes <- tab$gene_id
        sel.genes <- tab[as.character(sel), "gene_id"]  # selection by row name
        fe <- functionalEnrichment(all.genes, sel.genes, gene.go, go.terms)
      } else if(n > 0) {
        fe <- data.frame(Error=paste0('only ',max_points,' points can be selected.'))
      }
    }
    fe
  })
  
  output$ReactomeEnrichment <- renderTable({
    tab <- getData()
    sel <- NULL
    fe <- NULL
    if(!is.null(input$plot_brush)){
      brushed <- na.omit(brushedPoints(tab, input$plot_brush))
      sel <- as.numeric(rownames(brushed))
      n <- length(sel)
      if(n > 0 && n <= max_points) {
        all.genes <- tab$gene_id
        sel.genes <- tab[as.character(sel), "gene_id"]  # selection by row name
        fe <- functionalEnrichment(all.genes , sel.genes, gene.reactome, reactome.terms)
      } else if(n > 0) {
        fe <- data.frame(Error=paste0('only ',max_points,' points can be selected.'))
      }
    }
    fe
  })

  output$jitterPlot <- renderPlot({
    tab <- getData()
    sel <- selectGene()
    if(!is.null(sel) && length(sel) <= max_points) {
      d <- expr[sel,, drop=FALSE]
      jitterPlot(d, input$intensityScale)
    }
  })
  
  output$mainPlot <- renderPlot({
    tab <- getData()
    tab_idx <- as.numeric(input$allGeneTable_rows_selected)
    
    if(input$plotType == "Volcano") {
      g <- plotVolcano(tab)
    } else {
      g <- plotMA(tab)
    }
    if(length(tab_idx) >= 1) {
      g <- g + geom_point(data=tab[tab_idx, ], colour="red", size=2)
    }
    g
  })

  output$allGeneTable <- DT::renderDataTable({
    #allGeneTable(DE[[input$condPair]])
    tab <- getData()
    allGeneTable(tab)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
