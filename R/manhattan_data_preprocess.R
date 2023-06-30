#' Preprocess GWAS Result
#'
#' Preprocesses a result from Genome Wide Association Study
#' before making a manhattan plot.
#' It accepts a \code{data.frame}, which at bare minimum should
#' contain a chromosome, position, and p-value.
#' Additional options, such as chromosome color, label colum names,
#' and colors for specific variants, are provided here.
#'
#' @param x a data frame or any other extension of data frame (e.g. a tibble).
#'   At bare minimum, it should contain chromosome, position, and p-value.
#' @param ... Additional arguments for manhattan_data_preprocess.
#'
##' @details
#' \code{manhattan_data_preprocess} gathers information needed to plot a manhattan plot
#' and organizes the information as \code{MPdata} S3 object.
#'
#' New positions for each points are calculated, and stored in the data.frame as
#' \code{"new_pos"}. By default, all chromosomes will have the same width, with each
#' point being equally spaced. This behavior is changed when \code{preserve.position = TRUE}.
#' The width of each chromosome will scale to the number of points and the points will
#' reflect the original positions.
#'
#' \code{chr.col} and \code{highlight.col}, maps the data values to colors. If they are
#' an unnamed vector, then the function will try its best to match the values of
#' \code{chr.colname} or \code{highlight.colname} to the colors. If they are a named vector,
#' then they are expected to map all values to a color. If \code{highlight.colname} is
#' supplied, then \code{chr.col} is ignored.
#'
#' While feeding a \code{data.frame} directly into \code{manhattan_plot}
#' does preprocessing & plotting in one step. If you plan on making multiple plots
#' with different graphic options, you have the choice to preprocess separately and
#' then generate plots.
#'
#' @return a MPdata object. This object contains all the necessary info for constructing
#'   a manhattan plot.
#'
#' @examples
#'
#' gwasdat <- data.frame(
#'   "chromosome" = rep(1:5, each = 30),
#'   "position" = c(replicate(5, sample(1:300, 30))),
#'   "pvalue" = rbeta(150, 1, 1)^5
#' )
#'
#'   manhattan_data_preprocess(
#'   gwasdat, pval.colname = "pvalue", chr.colname = "chromosome", pos.colname = "position",
#'   chr.order = as.character(1:5)
#' )
#'
#' @export
manhattan_data_preprocess <- function(x, ...) UseMethod("manhattan_data_preprocess")

#' @rdname manhattan_data_preprocess
#' @method manhattan_data_preprocess default
#' @export
manhattan_data_preprocess.default <- function(x, ...) stop("Provide a valid data.frame or GRanges object to preprocess.")

#' @rdname manhattan_data_preprocess
#' @method manhattan_data_preprocess data.frame
#' @param chromosome a character. This is supplied if a manhattan plot of a single chromosome is
#'   desired. If \code{NULL}, then all the chromosomes in the data will be plotted.
#' @param signif a numeric vector. Significant p-value thresholds to be drawn for
#'   manhattan plot. At least one value should be provided. Default value is c(5e-08, 1e-5)
#' @param pval.colname a character. Column name of \code{x} containing p.value.
#' @param chr.colname a character. Column name of \code{x} containing chromosome number.
#' @param pos.colname a character. Column name of \code{x} containing position.
#' @param pval.log a logical. If \code{FALSE}, p-values will assumed to already be log10-transformed, otherwise they will be converted. Default is \code{TRUE}.
#' @param chr.order a character vector. Order of chromosomes presented in manhattan plot.
#' @param signif.col a character vector of equal length as \code{signif}. It contains
#'   colors for the lines drawn at \code{signif}. If \code{NULL}, the smallest value is colored
#'   black while others are grey.
#' @param chr.col a character vector of equal length as chr.order. It contains colors
#'   for the chromosomes. Name of the vector should match \code{chr.order}. If \code{NULL}, default
#'   colors are applied using \code{RColorBrewer}.
#' @param highlight.colname a character. If you desire to color certain points
#'   (e.g. significant variants) rather than color by chromosome, you can specify the
#'   category in this column, and provide the color mapping in \code{highlight.col}.
#'   Ignored if \code{NULL}.
#' @param highlight.col a character vector. It contains color mapping for the values from
#'   \code{highlight.colname}.
#' @param preserve.position a logical. If \code{TRUE}, the width of each chromosome reflect the
#'   number of variants and the position of each variant is correctly scaled? If \code{FALSE}, the
#'   width of each chromosome is equal and the variants are equally spaced.
#' @param thin a logical. If \code{TRUE}, \code{thinPoints} will be applied. Defaults to \code{TRUE} if
#'   \code{chromosome} is \code{NULL}. Defaults to \code{FALSE} if \code{chromosome} is supplied.
#' @param thin.n an integer. Number of max points per horizontal partitions of the plot.
#'   Defaults to 1000.
#' @importFrom ggplot2 waiver
#' @export
manhattan_data_preprocess.data.frame <- function(
  x, chromosome = NULL, signif = c(5e-8, 1e-5), pval.colname = "pval", pval.log = TRUE,
  chr.colname = "chr", pos.colname = "pos", highlight.colname = NULL, chr.order = NULL,
  signif.col = NULL, chr.col = NULL, highlight.col = NULL, preserve.position = FALSE, thin = NULL,
  thin.n = 1000
) {

  # what manhattan preprocess does:
  # orders the chromosome levels (X axis)
  # calculate x plot position for each variant based on chromosome width, gap.
  # Remove any results with missing chromosome of position
  # also, preprocess data for chromosome color, significance linetype, color, etc.

  # run checks on several arguments
  preprocess_arg_check_out <- preprocess_arg_check(
    x = x, chromosome = chromosome, signif = signif, signif.col = signif.col,
    pval.colname = pval.colname, chr.colname = chr.colname,
    pos.colname = pos.colname, pval.log = pval.log, preserve.position = preserve.position)

  thin <- set_thin_logical(thin, chromosome)

  # remove any results with missing chr, pos, or pval
  x <- remove_na(x, chr.colname, pos.colname, pval.colname)
  x[[pval.colname]] <- replace_0_pval(x[[pval.colname]], pval.log)
  signif.col <- preprocess_arg_check_out$signif.col

  # subset by chromosome if chromosome is specified
  if (!is.null(chromosome)) {
    x <- x[x[[chr.colname]] == chromosome,]
  }

  # factorize chromosome column to set order of chromosomes for the plot
  if (!is.null(chr.order)) {
    x[[chr.colname]] <- factor(x[[chr.colname]], levels = chr.order)
  } else if (!is.factor(x[[chr.colname]])) {
    x[[chr.colname]] <- factor(x[[chr.colname]])
    chr.order <- levels(x[[chr.colname]])
  } else {
    chr.order <- levels(x[[chr.colname]])
  }
  nchr <- length(chr.order)

  # chromosome / highlight color mapping
  chr.col <- set_chr_col(chr.col, nchr, chr.order)
  if (!is.null(highlight.colname)) {
    highlight.col <- set_highlight_col(x, highlight.colname, highlight.col)
  }

  # map each position in the chromosome to new positions
  x <- x[order(x[[chr.colname]], x[[pos.colname]]), ]

  if (preserve.position) {
    # scale the width of chromosome proportional to number of points in chromosome
    # keep original positioning
    chr_width <- table(x[[chr.colname]])
    chr_width <- as.numeric(chr_width); names(chr_width) <- chr.order
    chr_width <- chr_width / sum(chr_width) * nchr

    new_pos <-
      unlist(
        tapply(x[[pos.colname]], x[[chr.colname]], FUN = function(y) sequence_along_chr_scaled(y), simplify = FALSE),
        use.names = FALSE) * unname(chr_width[as.character(x[[chr.colname]])])

  } else {
    # all chromsomes have equal length & all variants are equally spaced
    chr_width <- rep(1, nchr)
    names(chr_width) <- chr.order

    new_pos <-
      unlist(
        tapply(x[[pos.colname]], x[[chr.colname]], FUN = function(y) sequence_along_chr_unscaled(y), simplify = FALSE),
        use.names = FALSE) * unname(chr_width[as.character(x[[chr.colname]])])
  }

  # fix certain widths for each chromosome, and gap for in between chromosomes
  lg <- 0.15 / 26 * nchr # gap between chromosome (should be robust with different lengths of chromosme)
  start_pos <- c(0, cumsum(chr_width)[-nchr]) + ((1:nchr - 1) * lg)
  names(start_pos) <- chr.order # starting x-coordinate for each chr
  end_pos <- start_pos + chr_width # ending x-coordinate for each chr
  center_pos <- (start_pos + end_pos) / 2 # middle x-coordinate for each chr... used for x axis labelling
  x$new_pos <- new_pos + start_pos[as.character(x[[chr.colname]])]

  # -log10(pvalue)
  x$log10pval <- ifelse(pval.log, -log10(x[[pval.colname]]), x[[pval.colname]])

  # thin data points if it set to true
  if (thin) {
    x <- thinPoints(dat = x, value = "log10pval", n = thin.n, nbins = 200, groupBy = chr.colname)
  }

  # Create MPdata Class
  mpdata <- list(
    data = x,
    start_pos = start_pos,
    center_pos = center_pos,
    end_pos = end_pos,
    signif = signif,
    signif.col = signif.col,
    chr.col = chr.col,
    highlight.colname = highlight.colname,
    highlight.col = highlight.col,
    chr.labels = chr.order,
    chr.colname = chr.colname,
    pos.colname = "new_pos",
    true.pos.colname = pos.colname,
    pval.colname = "log10pval"
  )
  class(mpdata) <- "MPdata"

  return(mpdata)
}

#' @rdname manhattan_data_preprocess
#' @method manhattan_data_preprocess GRanges
#' @export
setMethod(
  "manhattan_data_preprocess", signature = "GRanges",
  function(
    x, chromosome = NULL, signif = c(5e-8, 1e-5), pval.colname = "pval", pval.log = TRUE, highlight.colname = NULL, chr.order = NULL,
    signif.col = NULL, chr.col = NULL, highlight.col = NULL, preserve.position = FALSE, thin = NULL,
    thin.n = 100
  ) {
    grdat <- as.data.frame(x)
    grdat$pos <- (grdat$start + grdat$end) %/% 2

    chr.colname <- "seqnames"
    pos.colname <- "pos"

    manhattan_data_preprocess(
      grdat, chromosome = chromosome, signif = signif, pval.colname = pval.colname,
      chr.colname = chr.colname, pos.colname = pos.colname, pval.log = pval.log, highlight.colname = highlight.colname, chr.order = chr.order,
      signif.col = signif.col, chr.col = chr.col, highlight.col = highlight.col, preserve.position = preserve.position, thin = thin,
      thin.n = thin.n
    )
  }
)
