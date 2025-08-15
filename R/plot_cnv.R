#' Plot Heatmap
#'
#' Expects an annotated PennCNV output. The `gene` and `sample` columns are used
#' to create a matrix with the number/type of CNVs. The matrix is plotted as a
#' heatmap with the number of CNVs in each gene on the right as a bar plot and
#' the number of affected genes in each sample on the top.
#'
#' @param cnv A GRanges object such as the one created using `read_cnv`
#' @param ... Other arguments passed to `Heatmap`
#'
#' @importFrom tidyr unnest
#' @importFrom dplyr filter mutate
#' @importFrom reshape2 acast
#' @importFrom grid gpar
#' @importFrom ComplexHeatmap columnAnnotation rowAnnotation anno_barplot Heatmap
#'
#' @examples
#' # locate and read cnv file
#' fl <- system.file('extdata/family.gene.cnv', package = 'cnvr')
#' col_names <- c('region', 'numsnp', 'length', 'cn', 'sample', 'startsnp', 'endsnp', 'conf', 'gene', 'exon')
#' cnv <- read_cnv(fl, col_names)
#'
#' # plot heatmap
#' plot_heatmap(cnv, column_labels = c('father', 'mother', 'offspring'))
#'
#' @export
plot_heatmap <- function(cnv, ...) {
  # tidy cnv
  cnv <- as.data.frame(cnv)
  cnv <- transform(cnv, gene = strsplit(gene, ','))
  cnv <- unnest(cnv, gene)
  cnv <- filter(cnv, gene != 'NOT_FOUND')
  cnv <- mutate(cnv, cn = ifelse(cn > 2, 'loss', 'gain'))

  # heatmap
  # columns: sample_count
  sample_count <- table(cnv$sample)
  ca <- columnAnnotation(
    '# Samples' = anno_barplot(as.integer(sample_count))
  )

  # rows: gene_count
  gene_count <- table(cnv$gene)
  ra <- rowAnnotation(
    '# Variants' = anno_barplot(as.integer(gene_count))
  )

  # heatmap: cnvs
  heat_matrix <- acast(
    cnv,
    gene ~ sample,
    value.var = 'cn',
    fill = 'none',
    fun.aggregate = function(x) paste(sort(unique(x)), collapse = ',')
  )

  # colors
  colors <- c("darkblue", "white", "darkred", 'black')
  names(colors) <- c('loss', 'none', 'gain', 'gain,loss')


  # plot
  Heatmap(
    heat_matrix,
    col = colors,
    left_annotation = ra,
    top_annotation = ca,
    name = 'CNV',
    rect_gp = grid::gpar(col = "white", lwd = 2),
    ...
  )
}

#' Plot Signal
#'
#' Plots the position in the `gr` object at the x-axis and a metadata column at
#' the y-axis: lrr or baf. When `plot_gene = TRUE` a GRanges object with the
#' gene model is expected and plotted using `plot_gene` underneath the signal
#' panel.
#'
#' @param gr A `GRanges` object.
#' @param type Default `'LRR'`. Other accepted signals are `'BAF'`
#' @param gene_model A `GRanges` object.
#' @param plot_gene Default TRUE.
#' @param ... Other arguemnts passed to `plot`
#'
#' @examples
#' # get gene model
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' org <- org.Hs.eg.db::org.Hs.eg.db
#' gene_model <- get_genemodel(txdb, org, 'OR4C11')
#'
#' # signal
#' fl <- system.file('extdata/offspring.txt', package = 'cnvr')
#' col_names <- c('name', 'chr', 'pos', 'gt', 'lrr', 'baf')
#' signal <- read_signal(fl, col_names)
#'
#' # cnv
#' fl <- system.file('extdata/family.gene.cnv', package = 'cnvr')
#' col_names <- c('region', 'numsnp', 'length', 'cn', 'sample', 'startsnp', 'endsnp', 'conf', 'gene', 'exon')
#' cnv <- read_cnv(fl, col_names)
#' cnv <- cnv[grepl('OR4C11', cnv$gene) & cnv$cn == 0]
#'
#' # get the overlap
#' ol <- get_overlap(cnv, signal, flank = 5000)
#'
#' plot_signal(ol, type = 'LRR', gene_model)
#' plot_signal(ol, type = 'BAF', gene_model)
#'
#' @export
plot_signal <- function(gr, type = 'LRR', gene_model = NULL, plot_gene = TRUE, ...) {
  # extract position, x
  d <- as.data.frame(gr)
  x <- as.integer(d$start)

  # extract signal, y
  if (type == 'LRR') {
    y <- d$lrr
  } else if (type == 'BAF') {
    y <- d$baf
  } else {
    message("type should be one of LRR and BAF.")
  }

  if ( plot_gene ) {
    # stop if gene model is not provided
    if (is.null(gene_model)) stop("gene_model is required when plot_gene = TRUE.")

    # create a 1 x 2 layout
    layout(matrix(c(1,2)), heights = c(2, 1))

    # plot signal
    par(mar = c(1,5,1,1))
    plot(x, y, xlab = '', xaxt = 'n', ...)

    # plot gene
    par(mar = c(5,5,1,1))
    xlim <- c(min(d$start), max(d$end))
    plot_gene(gene_model, xlim)
  } else {
    # plot signal only when plot_gene = FALSE
    plot(x, y, ...)
  }
}

#' Get Overlap
#'
#' Expands the `cnv` GRanges object with amount `flank` on both sides. Then the
#' overlap between the two objects is extracted.
#'
#' @param cnv A GRanges object
#' @param signal A GRanges object
#'
#' @examples
#' # signal
#' fl <- system.file('extdata/offspring.txt', package = 'cnvr')
#' col_names <- c('name', 'chr', 'pos', 'gt', 'lrr', 'baf')
#' signal <- read_signal(fl, col_names)
#'
#' # cnv
#' fl <- system.file('extdata/family.gene.cnv', package = 'cnvr')
#' col_names <- c('region', 'numsnp', 'length', 'cn', 'sample', 'startsnp', 'endsnp', 'conf', 'gene', 'exon')
#' cnv <- read_cnv(fl, col_names)
#' cnv <- cnv[grepl('OR4C11', cnv$gene) & cnv$cn == 0]
#' cnv$gene <- 'OR4C11'
#'
#' get_overlap(cnv, signal)
#'
#' @importFrom IRanges subsetByOverlaps
#'
#' @export
get_overlap <- function(cnv, signal, flank = 1000) {
  # expand the cnv by flank on both sides
  cnv <- cnv + flank

  # get overlap betwee the two GRanges objects
  ol <- subsetByOverlaps(signal, cnv)

  return(ol)
}

#' Get Gene Model
#'
#' Extracts the coding regions (exons), and the gaps between them (introns) from
#' a TxDb object.
#'
#' @param txdb A TxDb object
#' @param org An OrgDb object
#' @param symbol A character gene symbol.
#'
#' @examples
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' org <- org.Hs.eg.db::org.Hs.eg.db
#' gr <- get_genemodel(txdb, org, c('OR4C11', 'SLX4IP'))
#'
#' @importFrom AnnotationDbi select
#' @importFrom GenomicFeatures exons
#' @importFrom GenomeInfoDb seqlevels `seqlevels<-` seqlevelsInUse
#' @importFrom GenomicRanges gaps GRangesList
#'
#' @export
get_genemodel <- function(txdb, org, symbol) {
  res <- GRangesList()
  for ( i in seq_along(symbol)) {
    # map the gene symbol to a gene_id
    gene_id <- select(org, symbol[i], 'ENTREZID', 'SYMBOL')$ENTREZID

    # extract coding regions, exons
    exons <- exons(txdb, filter = list(gene_id = gene_id))
    exons$feature <- 'coding_region'
    seqlevels(exons) <- seqlevelsInUse(exons)

    # extract gaps, introns
    introns <- gaps(exons, ignore.strand = TRUE)
    introns <- introns[2:(length(introns)-1)]
    seqlevels(introns) <- seqlevelsInUse(introns)
    introns$feature <- 'intron'

    # merge gr objects
    combined <- c(exons, introns)
    combined$gene <- symbol[i]

    res[[i]] <- combined
  }

  res <- unlist(res)
  return(res)
}

#' Plot Gene Model
#'
#' Plot the exons and introns as rectangles and lines.
#'
#' @param gr A GRanges object
#' @param xlim A numeric vector of length 2.
#'
#' @examples
#' d <- data.frame(
#'   seqnames = rep('chr1', 7),
#'   start = c(10, 40, 60, 30, 50, 05, 80),
#'   end   = c(30, 50, 80, 40, 60, 10, 90),
#'   feature  = rep(c('coding_region', 'intron', "5' utr", "3' utr"), c(3,2,1,1)),
#'   gene = rep('gene1', 7)
#' )
#'
#' gr1 <- GenomicRanges::makeGRangesFromDataFrame(d, keep.extra.columns = TRUE)
#' gr2 <- GenomicRanges::shift(gr1, 10)
#' gr2$gene <- 'gene2'
#'
#' gr <- c(gr1, gr2)
#'
#' plot_gene(gr, c(0, 100))
#'
#' @export
plot_gene <- function(gr, xlim, ...) {
  # convert to a data.frame
  d <- as.data.frame(gr)
  d$gene <- as.factor(d$gene)
  d$y <- as.numeric(as.factor(d$gene))

  # create empty plot frame
  plot(1,
       xlim = xlim,
       ylim = c(0, length(levels(d$gene)) + 1),
       type = "l",
       yaxt = 'n',
       xlab = unique(d$seqnames),
       ylab = '',
       frame.plot = FALSE,
       ...)

  axis(2, unique(d$y) - 0.5, labels = levels(d$gene), las = 2)

  # loop over features and plot
  apply(d, 1, function(x) {
    feature = x['feature']
    start = as.integer(x['start'])
    end = as.integer(x['end'])
    gene <- as.integer(x['gene'])
    y <- as.integer(x['y']) - 1

    if ( feature == 'coding_region' ) {
      rect(start, y + .1, end, y + .9, col = 'gray')
    } else if ( feature %in% c("5' utr", "3' utr") ) {
      rect(start, y + .1, end, y + .9, col = 'darkgray')
    } else if  ( feature == 'intron' ) {
      mid <- c(start + end) / 2
      segments(x0 = start, y0 = y + .3, x1 = mid, y1 = y + .7)
      segments(x0 = mid, y0 = y + .7, x1 = end, y1 = y + .3)
    }
  })
}
