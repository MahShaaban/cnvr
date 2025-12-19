#' Read Signal
#'
#' Expects a tab separate file with at 5 columns: name, chr, pos, lrr, and baf.
#' Column names can be changed using the `col_names` argument.
#'
#' @param file A tab separated file
#' @param col_names A character vector with length equal to the `file` columns.
#' @param as.granges Default TRUE.
#' @param style Default 'UCSC'
#' @param ... Other arguments passed to `read_tsv`
#'
#' @examples
#' # locate signal file
#' fl <- system.file('extdata/offspring.txt', package = 'cnvr')
#'
#' # read signal file
#' col_names <- c('name', 'chr', 'pos', 'gt', 'lrr', 'baf')
#' read_signal(fl, col_names)
#'
#' # no pos
#' fl <- system.file('extdata/example.pfb', package = 'cnvr')
#' col_names <- c('name', 'chr', 'pos', 'pfb')
#' pfb <- read_pfb(fl, col_names)
#'
#' # read signal file, no pos
#' fl <- system.file('extdata/offspring.nopos.txt', package = 'cnvr')
#' col_names <- c('name', 'gt', 'lrr', 'baf')
#' read_signal(fl, col_names, pfb)
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr inner_join
#' @importFrom GenomicRanges makeGRangesFromDataFrame sort
#' @importFrom GenomeInfoDb `seqlevelsStyle<-` seqlevelsStyle
#'
#' @export
read_signal <- function(file, col_names = NULL, pfb = NULL, as.granges = TRUE,
                        style = 'UCSC', ...) {
  # read tsv file
  d <- read_tsv(file, ...)

  # change name when given
  if (!is.null(col_names)) {
    names(d) <- col_names
  }

  # merge pfb when provided
  if (!is.null(pfb)) {
    d <- inner_join(d, pfb)
  }

  if ( as.granges ) {
    # create a GRanges object from signal and sort
    d <- makeGRangesFromDataFrame(
      d,
      start.field = 'pos',
      end.field = 'pos',
      keep.extra.columns = TRUE
    )
    d <- sort(d)
    seqlevelsStyle(d) <- style
  }

  return(d)
}

#' Read PFB file
#'
#' Expects a tab separate file with at 4 columns: name, chr, pos, and pfb.
#' Column names can be changed using the `col_names` argument.
#'
#' @inheritParams read_signal
#' @param as.granges Default FALSE.
#'
#' @examples
#' # locate pfb file
#' fl <- system.file('extdata/example.pfb', package = 'cnvr')
#'
#' # read pfb file
#' col_names <- c('name', 'chr', 'pos', 'pfb')
#' read_pfb(fl, col_names)
#' read_pfb(fl, col_names, as.granges = TRUE)
#'
#' @importFrom readr read_tsv
#' @importFrom GenomicRanges makeGRangesFromDataFrame sort
#' @importFrom GenomeInfoDb `seqlevelsStyle<-` seqlevelsStyle
#'
#' @export
read_pfb <- function(file, col_names = NULL, as.granges = FALSE, style = 'UCSC',
                     ...) {
  # read tsv file
  d <- read_tsv(file, ...)

  # change name when given
  if (!is.null(col_names)) {
    names(d) <- col_names
  }

  if ( as.granges ) {
    # create a GRanges object from signal and sort
    d <- makeGRangesFromDataFrame(
      d,
      start.field = 'pos',
      end.field = 'pos',
      keep.extra.columns = TRUE
    )
    d <- sort(d)
    seqlevelsStyle(d) <- style
  }

  return(d)
}

#' Read cytoBands file
#'
#' Expects a tab separate file with at 5 columns: chrom, start, end, band, stain.
#'
#' @inheritParams read_signal
#' @param as.granges Default TRUE.
#'
#' @examples
#' # locate cytobands file
#' fl <- system.file('extdata/hg38.cytoBand.txt', package = 'cnvr')
#'
#' # read cytobands file
#' read_cytobands(fl)
#'
#' @importFrom readr read_tsv
#' @importFrom GenomicRanges makeGRangesFromDataFrame sort
#' @importFrom GenomeInfoDb `seqlevelsStyle<-` seqlevelsStyle
#'
#' @export
read_cytobands <- function(file, as.granges = TRUE, style = 'UCSC', ...) {
  # read tsv file
  d <- read_tsv(file, col_names = c('chrom', 'start', 'end', 'band', 'stain'))

  if ( as.granges ) {
    # create a GRanges object from signal and sort
    d <- makeGRangesFromDataFrame(d, keep.extra.columns = TRUE)
    d <- sort(d)
    seqlevelsStyle(d) <- style
  }

  return(d)
}

#' Read CNV file
#'
#' Expects tab delimited files similar to PennCNV output. Columns names can be
#' edited using the `col_names` argument.
#'
#' @inheritParams read_signal
#' @param numeric_columns A character vector. Columns to clean as numeric
#' @param character_columns A character vector. Columns to clean as character
#' @param tidy Default TRUE.
#'
#' @examples
#' # locate cnv file
#' fl <- system.file('extdata/offspring.cnv', package = 'cnvr')
#' col_names <- c('region', 'numsnp', 'length', 'cn', 'sample', 'startsnp', 'endsnp', 'conf')
#'
#' # read cnv file
#' read_cnv(fl, col_names)
#'
#' @importFrom readr read_table
#' @importFrom dplyr mutate_at vars
#' @importFrom tidyr separate
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb `seqlevelsStyle<-` seqlevelsStyle
#'
#' @export
read_cnv <- function(file, col_names = FALSE, tidy = TRUE, as.granges = TRUE,
                     numeric_columns = c('numsnp', 'length', 'cn', 'conf'),
                     character_columns = c('startsnp', 'endsnp'),
                     style = 'UCSC', ...) {
  # read tsv file
  d <- readr::read_table(file, col_names = col_names, ...)

  if ( tidy ) {
    d <- mutate_at(d, vars(numeric_columns), ~clean_field(.x))
    d <- mutate_at(d, vars(character_columns), ~clean_field(.x, FALSE))

    d <- separate(d, region, into = c('chrom', 'start', 'end'), remove = FALSE)
    d <- mutate_at(d, vars('start', 'end'), as.integer)
  }

  if ( as.granges ) {
    # create a GRanges object from cnv and sort
    d <- makeGRangesFromDataFrame(
      d,
      keep.extra.columns = TRUE
    )
    d <- sort(d)
    seqlevelsStyle(d) <- style
  }

  return(d)
}

#' Clean CNV fields
#'
#' @examples
#' x <- c('snp=1', 'snp=2')
#' clean_field(x)
#'
#' x <- c('snp=snp1', 'snp=snp2')
#' clean_field(x, FALSE)
#'
#' @importFrom stringr str_remove_all
#'
#' @export
clean_field <- function(x, integer = TRUE) {
  res <- str_remove_all(x, '.*=|,')
  if ( integer ) res <- as.integer(res)

  return(res)
}

