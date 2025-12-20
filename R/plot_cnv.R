#' Plot Heatmap
#'
#' Expects an annotated PennCNV output. The `gene` and `sample` columns are used
#' to create a matrix with the number/type of CNVs. The matrix is plotted as a
#' heatmap with the number of CNVs in each gene on the right as a bar plot and
#' the number of affected genes in each sample on the top.
#'
#' @param cnv A GRanges object such as the one created using `read_cnv`
#' @param top_samples An integer. The number of top samples to plot
#' @param top_genes An integer. The number of top genes to plot
#' @param ... Other arguments passed to `Heatmap`
#'
#' @importFrom tidyr unnest separate
#' @importFrom dplyr filter mutate
#' @importFrom reshape2 acast
#' @importFrom grid gpar unit
#' @importFrom ComplexHeatmap columnAnnotation rowAnnotation anno_barplot Heatmap
#' @importFrom circlize colorRamp2
#'
#' @examples
#' # locate and read cnv file
#' fl <- system.file('extdata/family.gene.cnv', package = 'cnvr')
#' col_names <- c('region', 'numsnp', 'length', 'cn', 'sample', 'startsnp', 'endsnp', 'conf', 'gene', 'distance')
#' cnv <- read_cnv(fl, col_names)
#' cnv <- cnv[cnv$gene != 'NOT_FOUND']
#'
#' # plot heatmap
#' plot_heatmap(cnv)
#'
#' @export
plot_heatmap <- function(cnv, top_samples = NULL, top_genes = NULL, ...) {
  # tidy cnv
  cnv <- as.data.frame(cnv)
  cnv <- transform(cnv, gene = strsplit(gene, ','))
  cnv <- unnest(cnv, gene)

  # subset to top_genes
  if (is.null(top_genes)) top_genes <- length(unique(cnv$gene))
  top_genes <- head(names(sort(table(cnv$gene), decreasing = TRUE)), top_genes)
  cnv <- filter(cnv, gene %in% top_genes)

  # subset to top_samples
  if (is.null(top_samples)) top_samples <- length(unique(cnv$sample))
  top_samples <- head(names(sort(table(cnv$sample), decreasing = TRUE)), top_samples)
  cnv <- filter(cnv, sample %in% top_samples)

  # columns: sample_count
  sample_count <- table(cnv$sample)
  ca <- columnAnnotation(
    '# Samples' = anno_barplot(as.integer(sample_count)),
    annotation_height = unit(1, 'in')
  )

  # rows: gene_count
  gene_count <- table(cnv$gene)

  ra <- rowAnnotation(
    '# Variants' = anno_barplot(as.integer(gene_count)),
    annotation_width = unit(1, 'in')
  )

  # heatmap: cnvs
  heat_matrix <- acast(
    cnv,
    gene ~ sample,
    value.var = 'cn',
    fill = 2,
  )

  # colors
  colors <- colorRamp2(
    breaks = c(1,2,3),
    colors = c("darkblue", "white", "darkred")
  )

  # plot
  Heatmap(
    heat_matrix,
    col = colors,
    right_annotation = ra,
    top_annotation = ca,
    rect_gp = grid::gpar(col = "white", lwd = 2),
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    heatmap_legend_param = list(
      title = "CNV",
      at = c(1,2,3),
      labels = c(1,2,3)
    ),
    ...
  )
}

#' Plot Signal
#'
#' Plots the position in the `signal` `GRanges` object at the x-axis and a
#' metadata column at the y-axis: lrr or baf. When `plot_gene = TRUE` a `GRanges`
#' object with the gene model when provided.
#'
#' @param signal A `GRanges` object.
#' @param cnv A `GRanges` object.
#' @param type Default `'LRR'`. Other accepted signals are `'BAF'`
#' @param gene_model A `GRanges` object.
#' @param bands A `GRanges` object.
#' @param ... Other arguemnts passed to `plot`
#'
#' @examples
#' # get gene model
#' gene <- 'SLX4IP'
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' org <- org.Hs.eg.db::org.Hs.eg.db
#' gene_model <- get_genemodel(txdb, org, gene)
#'
#' # signal
#' fl <- system.file('extdata/offspring.txt', package = 'cnvr')
#' col_names <- c('name', 'chr', 'pos', 'gt', 'lrr', 'baf')
#' signal <- read_signal(fl, col_names)
#'
#' # cnv
#' fl <- system.file('extdata/family.gene.cnv', package = 'cnvr')
#' col_names <- c('region', 'numsnp', 'length', 'cn', 'sample', 'startsnp', 'endsnp', 'conf', 'gene', 'distance')
#' cnv <- read_cnv(fl, col_names)
#'
#' sample <- "data-raw/PennCNV/example/offspring.txt"
#' cn <- 1
#'
#' cnv <- cnv[cnv$sample == sample & cnv$cn == cn & grepl(gene, cnv$gene)]
#' cnv$gene <- gene
#'
#' flank <- 200000
#' # read cytobands file
#' fl <- system.file('extdata/hg38.cytoBand.txt', package = 'cnvr')
#' cytobands <- read_cytobands(fl)
#' cytobands <- IRanges::subsetByOverlaps(cytobands, cnv)
#'
#' # plot
#' plot_signal(signal, cnv, flank, type = 'LRR', gene_model, bands = cytobands)
#' plot_signal(signal, cnv, flank, type = 'BAF', gene_model, bands = cytobands)
#'
#' @export
plot_signal <- function(signal, cnv, flank = 200000, type = 'LRR',
                        gene_model = NULL, bands = NULL, ...) {
  # get the overlap
  gr <- get_overlap(cnv, signal, flank)

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

  if ( !is.null(gene_model) ) {
    # create a 1 x 2 layout
    layout(matrix(c(1,2)), heights = c(2, 1))
    xlim <- c(min(d$start), max(d$end))

    # plot signal
    par(mar = c(.5,5,5,1))
    plot(x, y, xlim = xlim, xlab = '', xaxt = 'n', type = 'n', ...)
    if ( !is.null(bands) ) {
      plot_bands(bands)
    }

    points(x, y, ...)

    region <- as.data.frame(cnv)
    abline(h = 0, lty = 2, col = 'red')
    abline(h = 0, lty = 2, col = 'red')
    abline(v = c(region$start, region$end), lty = 2, col = 'red')

    # plot gene
    par(mar = c(5,5,.5,1))
    plot_gene(gene_model, xlim)
  } else {
    # plot signal only when plot_gene = FALSE
    plot(x, y, type = 'n', ...)
    if ( !is.null(bands) ) {
      plot_bands(bands)
    }

    points(x, y, ...)
    abline(h = 0, lty = 2, col = 'red')
    abline(v = c(GenomicRanges::start(cnv), GenomicRanges::end(cnv)), lty = 2, col = 'red')

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
#' col_names <- c('region', 'numsnp', 'length', 'cn', 'sample', 'startsnp', 'endsnp', 'conf', 'gene', 'distance')
#' cnv <- read_cnv(fl, col_names)
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

  d2 <- subset(d, select = c('start', 'y'))
  d2 <- aggregate(d2, list(gene = d$gene), min)
  text(d2$start, d2$y - 0.5, labels = d2$gene, xpd = NA, adj = 1.2)

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

#' Plot cytobands
#'
#' @param gr A GRanges object
#' @param shade A Logical, default TRUE. Shade bands
#' @param names A Logical, default TRUE. Plot band names
#'
#' @return A plot
#'
#' @examples
#' x <- 1:10
#' set.seed(1234)
#' y <- rnorm(10)
#' b <- tibble::tibble(
#'   chrom = rep('chr1', 3),
#'   start = c(1, 4, 8),
#'   end = c(4, 8, 10),
#'   band = c('p1', 'p2', 'p3'),
#'   stain = c('gneg', 'gpos', 'gneg')
#' )
#' plot(x, y, type = 'n')
#' plot_bands(b)
#' points(x, y)
#'
#' @export
plot_bands <- function(gr, shade = TRUE, names = TRUE) {
  b <- as.data.frame(gr)
  if ( shade ) {
    shade <- c('white', 'gray')[as.factor(b$stain)]
  } else {
    shade <- NULL
  }

  if ( names ) {
    center <- (b$start+b$end) / 2
    nms <- b$band
  } else {
    center <- (b$start+b$end) / 2
    nms <- ''
  }

  rect(b$start, -10, b$end, 10, col = shade)
  mtext(nms, at = center)
}
