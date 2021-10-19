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
#' library(dplyr)
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
methods::setGeneric("manhattan_data_preprocess", function(x, ...) standardGeneric("manhattan_data_preprocess"), signature = "x")

#' @rdname manhattan_data_preprocess
#' @export
methods::setMethod(
  "manhattan_data_preprocess", signature = "ANY",
  function(x, ...) stop("Provide a valid data.frame or GRanges object to preprocess.")
)

#' @rdname manhattan_data_preprocess
#'
#' @param signif a numeric vector. Significant p-value thresholds to be drawn for
#'   manhattan plot. At least one value should be provided. Default value is c(5e-08, 1e-5)
#' @param pval.colname a character. Column name of \code{x} containing p.value.
#' @param chr.colname a character. Column name of \code{x} containing chromosome number.
#' @param pos.colname a character. Column name of \code{x} containing position.
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
#' @param thin a logical. Reduce number of data points when they are cluttered?
#' @param thin.n an integer. Number of max points per horizontal partitions of the plot.
#'   Defaults to 500.
#' @importFrom ggplot2 waiver
#' @export
methods::setMethod(
  "manhattan_data_preprocess", signature = "data.frame",
  function(
    x, signif = c(5e-8, 1e-5), pval.colname = "pval",
    chr.colname = "chr", pos.colname = "pos", highlight.colname = NULL, chr.order = NULL,
    signif.col = NULL, chr.col = NULL, highlight.col = NULL, preserve.position = FALSE, thin = TRUE,
    thin.n = 500
  ) {

    # what manhattan preprocess does:
    # orders the chromosome levels (X axis)
    # calculate x plot position for each variant based on chromosome width, gap.
    # Remove any results with missing chromosome of position
    # also, preprocess data for chromosome color, significance linetype, color, etc.

    # run checks on several arguments
    preprocess_arg_check_out <- preprocess_arg_check(
      x =x, signif = signif, signif.col = signif.col,
      pval.colname = pval.colname, chr.colname = chr.colname,
      pos.colname = pos.colname, preserve.position = preserve.position)

    # remove any results with missing chr, pos, or pval
    x <- remove_na(x, chr.colname, pos.colname, pval.colname)
    x[[pval.colname]] <- replace_0_pval(x[[pval.colname]])
    signif.col <- preprocess_arg_check_out$signif.col

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
    x$log10pval <- -log10(x[[pval.colname]])

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
)

#' @rdname manhattan_data_preprocess
#' @export
methods::setMethod(
  "manhattan_data_preprocess", signature = "GRanges",
  function(
    x, signif = c(5e-8, 1e-5), pval.colname = "pval", highlight.colname = NULL, chr.order = NULL,
    signif.col = NULL, chr.col = NULL, highlight.col = NULL, preserve.position = FALSE, thin = TRUE,
    thin.n = 500
  ) {
    grdat <- as.data.frame(x)
    grdat$pos <- (grdat$start + grdat$end) %/% 2

    chr.colname <- "seqnames"
    pos.colname <- "pos"

    manhattan_data_preprocess(
      grdat, signif = signif, pval.colname = pval.colname,
      chr.colname = chr.colname, pos.colname = pos.colname, highlight.colname = highlight.colname, chr.order = chr.order,
      signif.col = signif.col, chr.col = chr.col, highlight.col = highlight.col, preserve.position = preserve.position, thin = thin,
      thin.n = thin.n
    )
  }
)

#' Manhattan Plotting
#'
#' A generic function for manhattan plot.
#'
#' @param x a \code{data.frame}, an extension of \code{data.frame} object (e.g.
#'   \code{tibble}), or an \code{MPdata} object.
#' @param ... Additional arguments for manhattan plot.
#' @inheritParams manhattan_data_preprocess
#'
#' @details
#' This generic function accepts a result of a GWAS in the form of \code{data.frame}
#' or a \code{MPdata} object produced by \code{manhattan_data_preprocess}. The
#' function will throw an error if another type of object is passed.
#'
#' Having \code{rescale = TRUE} is useful when there are points with very
#' high -log10(p.value). In this case, the function attempts to split
#' the plot into two different scales, with the split happening near the strictest
#' significance threshold. More precisely, the plot is rescaled when
#' \deqn{-log10(pvalue) / (strictest significance threshold) \ge rescale.ratio.threshold}
#'
#' If you wish to add annotation to the points, provide the name of the column to
#' \code{label.colname}. The labels are added with \code{\link{ggrepel}}.
#'
#' Be careful though: if the annotation column contains
#' a large number of variants, then the plotting could take a long time, and the
#' labels will clutter up the plot. For those points with no annotation, you have the
#' choice to set them as \code{NA} or \code{""}.
#'
#' @return \code{gg} object if \code{is.null(outfn)}, \code{NULL} if \code{!is.null(outf)}
#'
#' @examples
#' library(dplyr)
#'
#' gwasdat <- data.frame(
#'   "chromosome" = rep(1:5, each = 30),
#'   "position" = c(replicate(5, sample(1:300, 30))),
#'   "pvalue" = rbeta(150, 1, 1)^5
#' )
#'
#' manhattan_plot(
#'   gwasdat, pval.colname = "pvalue", chr.colname = "chromosome", pos.colname = "position",
#'   chr.order = as.character(1:5)
#' )
#'
#' mpdata <- manhattan_data_preprocess(
#'   gwasdat, pval.colname = "pvalue", chr.colname = "chromosome", pos.colname = "position",
#'   chr.order = as.character(1:5)
#' )
#'
#' manhattan_plot(mpdata)
#'
#' @rdname manhattan_plot
#' @export
methods::setGeneric("manhattan_plot", function(x, ...) standardGeneric("manhattan_plot"), signature = "x")

#' @rdname manhattan_plot
#' @export
methods::setMethod(
  "manhattan_plot", signature = "ANY",
  function(x, ...) stop("Provide a data.frame to preprocess & plot or a preprocessed MPdata object.")
)

#' @rdname manhattan_plot
#'
#' @export
methods::setMethod(
  "manhattan_plot", signature = "data.frame",
  function(
    x, outfn = NULL, signif = c(5e-8, 1e-5), pval.colname = "pval", chr.colname = "chr",
    pos.colname = "pos", label.colname = NULL, highlight.colname = NULL, chr.order = NULL,
    signif.col = NULL, chr.col = NULL,  highlight.col = NULL,
    rescale = TRUE, rescale.ratio.threshold = 5, signif.rel.pos = 0.4, color.by.highlight = FALSE,
    preserve.position = FALSE, thin = TRUE, thin.n = 500,
    plot.title = ggplot2::waiver(), plot.subtitle = ggplot2::waiver(), plot.width = 10, plot.height = 5,
    point.size = 0.75, label.font.size = 2, max.overlaps = 20,
    x.label = "Chromosome", y.label = expression(-log[10](p)), ...
  ) {

    # preprocess manhattan plot data
    mpdata <- manhattan_data_preprocess(
      x, signif = signif, pval.colname = pval.colname,
      chr.colname = chr.colname, pos.colname = pos.colname, chr.order = chr.order,
      signif.col = signif.col, chr.col = chr.col, highlight.colname = highlight.colname,
      highlight.col = highlight.col, preserve.position = preserve.position, thin = thin
    )

    # manhattan plot
    manhattan_plot(
      x = mpdata, outfn = outfn, rescale = rescale, rescale.ratio.threshold = rescale.ratio.threshold,
      signif.rel.pos = signif.rel.pos, color.by.highlight = color.by.highlight,
      label.colname = label.colname, x.label = x.label, y.label = y.label,
      point.size = point.size, label.font.size = label.font.size,
      max.overlaps = max.overlaps, plot.title = plot.title, plot.subtitle = plot.subtitle,
      plot.width = plot.width, plot.height = plot.height, ...
    )

  }
)

#' @rdname manhattan_plot
#' @method manhattan_plot MPdata
#'
#' @param outfn a character. File name to save the Manhattan Plot. If \code{outfn}
#'  is supplied (i.e. \code{!is.null(outfn)}), then the plot is not drawn in
#'  the graphics window.
#' @param rescale a logical. If \code{TRUE}, the plot will rescale itself depending
#' on the data. More on this in details.
#' @param rescale.ratio.threshold a numeric. Threshold of that triggers the rescale.
#' @param signif.rel.pos a numeric between 0.1 and 0.9. If the plot is rescaled,
#' where should the significance threshold be positioned?
#' @param color.by.highlight a logical. Should the points be colored based on a highlight column?
#' @param label.colname a character. Name of the column in \code{MPdata$data}
#'   to be used for labelling.
#' @param x.label a character. x-axis label
#' @param y.label a character. y-axis label
#' @param point.size a numeric. Size of the points.
#' @param label.font.size a numeric. Size of the labels.
#' @param max.overlaps an integer. Exclude text labels that overlaps too many things.
#' @param plot.title a character. Plot title
#' @param plot.subtitle a character. Plot subtitle
#' @param plot.width a numeric. Plot width in inches.
#' @param plot.height a numeric. Plot height in inches.
#' @param ... additional arguments to be passed onto \code{geom_label_repel}
#'
#' @export
methods::setMethod(
  "manhattan_plot", signature = "MPdata",
  function(
    x, outfn = NULL,
    rescale = TRUE, rescale.ratio.threshold = 5, signif.rel.pos = 0.4, color.by.highlight = FALSE,
    label.colname = NULL, x.label = "Chromosome", y.label = expression(-log[10](p)),
    point.size = 0.75, label.font.size = 2, max.overlaps = 20,
    plot.title = ggplot2::waiver(), plot.subtitle = ggplot2::waiver(),
    plot.width = 10, plot.height = 5, ...
  ) {

    if (all(!is.null(label.colname), !is.na(label.colname), na.rm = TRUE)) {
      if (!(label.colname %in% colnames(x$data))) stop("label.colname not a valid column name for the data.")
      if ((sum(!is.na(x$data[[label.colname]]) & !(x$data[[label.colname]] == ""))) > 20) warning("The plot will generate > 20 labels. The plot may be cluttered & may take a longer to generate.")
    }

    # create transformation object; if rescaling is required, create the according transformation
    trans <- list("trans" = "identity", "breaks" = ggplot2::waiver())
    if (rescale) {
      jump <- get_transform_jump(-log10(x$signif))
      if ((ceiling(max(x$data[[x$pval.colname]])/5)*5)/jump > rescale.ratio.threshold) {
        trans <- get_transform(x$data, jump, x$pval.colname, jump.rel.pos = signif.rel.pos)
      }
    }

    ylimit <- c(0, ifelse(identical(trans$trans, "identity"), NA, max(trans$breaks)))

    # choose whether to use highlight.colname or chr.colname
    if (!is.null(x$highlight.colname) && !is.null(x$highlight.col) && color.by.highlight) {
      point.color <- x$highlight.colname
      point.color.map <- x$highlight.col
    } else {
      point.color <- x$chr.colname
      point.color.map <- x$chr.col
    }

    # manhattan plot without labels
    p <- ggplot2::ggplot(x$data, ggplot2::aes_string(x = x$pos.colname, y = x$pval.colname, color = point.color)) +
      ggplot2::geom_point(size = point.size, pch = 16) +
      ggplot2::scale_discrete_manual(aesthetics = "color", values = point.color.map) +
      ggplot2::scale_y_continuous(
        trans = trans$trans,
        breaks = trans$breaks,
        expand = c(0.02, 0.01),
        limits = ylimit
      ) +
      ggplot2::scale_x_continuous(
        name = x.label,
        breaks = x$center_pos,
        labels = x$chr.labels,
        expand = c(0.01, 0.01)
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(x$signif),
        linetype = 'dashed',
        color = x$signif.col
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      ggplot2::ylab(y.label) +
      ggplot2::ggtitle(label = plot.title, subtitle = plot.subtitle)

    if (all(!is.null(label.colname), !is.na(label.colname), na.rm = TRUE)) {
      p <- p + ggrepel::geom_label_repel(
        ggplot2::aes_string(x = x$pos.colname, y = x$pval.colname, label = label.colname),
        size = label.font.size,
        label.padding = 0.15,
        segment.size = 0.2,
        min.segment.length = 0,
        max.overlaps = max.overlaps,
        ...
      )
    }

    if (rescale & !identical(trans$trans, "identity")) {
      # if the plot is rescaled, change the tick at the "jump" to double line
      jump.tick.size <- 3.5

      p <- p +
        ggplot2::theme(
          axis.ticks.y = ggplot2::element_line(linetype = trans$y_axis_linetype)) +
        ggplot2::annotate(geom = "point", shape = "=", x = -Inf, y = trans$jump, size = jump.tick.size) +
        ggplot2::coord_cartesian(clip = "off")
    }

    if (!is.null(outfn)) {
      ggplot2::ggsave(outfn, plot=p, width=plot.width, height=plot.height, units = "in")
      invisible()
    } else {
      return(p)
    }

  }
)

#' @rdname manhattan_plot
#'
#' @export
methods::setMethod(
  "manhattan_plot", signature = "GRanges",
  function(
    x, outfn = NULL, signif = c(5e-8, 1e-5), pval.colname = "pval", label.colname = NULL,
    highlight.colname = NULL, chr.order = NULL,
    signif.col = NULL, chr.col = NULL,  highlight.col = NULL,
    rescale = TRUE, rescale.ratio.threshold = 5, signif.rel.pos = 0.4, color.by.highlight = FALSE,
    preserve.position = FALSE, thin = TRUE, thin.n = 500,
    plot.title = ggplot2::waiver(), plot.subtitle = ggplot2::waiver(), plot.width = 10, plot.height = 5,
    point.size = 0.75, label.font.size = 2, max.overlaps = 20,
    x.label = "Chromosome", y.label = expression(-log[10](p)), ...
  ) {

    # preprocess manhattan plot data
    mpdata <- manhattan_data_preprocess(
      x, signif = signif, pval.colname = pval.colname, chr.order = chr.order,
      signif.col = signif.col, chr.col = chr.col, highlight.colname = highlight.colname,
      highlight.col = highlight.col, preserve.position = preserve.position, thin = thin
    )

    # manhattan plot
    manhattan_plot(
      x = mpdata, outfn = outfn, rescale = rescale, rescale.ratio.threshold = rescale.ratio.threshold,
      signif.rel.pos = signif.rel.pos, color.by.highlight = color.by.highlight,
      label.colname = label.colname, x.label = x.label, y.label = y.label,
      point.size = point.size, label.font.size = label.font.size,
      max.overlaps = max.overlaps, plot.title = plot.title, plot.subtitle = plot.subtitle,
      plot.width = plot.width, plot.height = plot.height, ...
    )

  }
)

#' Chromosome Manhattan Plotting
#'
#' Manhattan Plot of a chosen chromosome
#'
#' @param x a \code{data.frame}, an extension of \code{data.frame} object (e.g.
#'   \code{tibble}), or an \code{MPdata} object.
#' @param chromosome a character. a specific chromosome in \code{x} to be plotted
#' @param ... Additional arguments for chromosome manhattan plot
#'
#' @inheritParams manhattan_plot
#' @inheritParams manhattan_data_preprocess
#'
#' @rdname manhattan_chromosome
#' @export
#'
manhattan_chromosome <- function(x, chromosome, ...) {
  UseMethod("manhattan_chromosome", x)
}

#' @rdname manhattan_chromosome
#' @method manhattan_chromosome data.frame
#'
#' @param chromosome a character. A specific chromsome in \code{x} to be plotted
#' @return \code{gg} object if \code{is.null(outfn)}, \code{NULL} if \code{!is.null(outf)}
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' gwasdat <- data.frame(
#'   "chromosome" = rep(1:5, each = 30),
#'   "position" = c(replicate(5, sample(1:300, 30))),
#'   "pvalue" = rbeta(150, 1, 1)^5
#' )
#'
#' manhattan_chromosome(
#'   gwasdat, pval.colname = "pvalue", chromosome = "3",
#'   chr.colname = "chromosome", pos.colname = "position",
#'   chr.order = as.character(1:5)
#' )
#'
manhattan_chromosome.data.frame <- function(
  x, chromosome, outfn = NULL, signif = c(5e-8, 1e-5),
  pval.colname = "pval", chr.colname = "chr", pos.colname = "pos",
  chr.order = NULL, signif.col = NULL, chr.col = NULL, highlight.colname = NULL,
  preserve.position = FALSE, thin = FALSE, thin.n = 3000,
  rescale = TRUE, rescale.ratio.threshold = 5,
  signif.rel.pos = 0.4, label.colname = NULL,
  point.size = 0.75, x.label = "Chromosome",
  y.label = expression(-log[10](p)), label.font.size = 2, max.overlaps = 20,
  plot.title = ggplot2::waiver(),
  plot.subtitle = ggplot2::waiver(), plot.width = 10, plot.height = 5, ...
) {

  x <- x[x[[chr.colname]] == chromosome,]

  if (nrow(x) < 1) {
    stop("No matches at the chromosome.")
  }

  mpdata <- manhattan_data_preprocess(
    x, signif = signif, pval.colname = pval.colname,
    chr.colname = chr.colname, pos.colname = pos.colname, chr.order = chromosome,
    signif.col = signif.col, chr.col = chr.col, highlight.colname = highlight.colname,
    preserve.position = preserve.position, thin = thin
  )

  manhattan_chromosome(
    mpdata, chromosome = chromosome, outfn = outfn,
    rescale = rescale, rescale.ratio.threshold = rescale.ratio.threshold, signif.rel.pos = signif.rel.pos,
    label.colname = label.colname, point.size = point.size, label.font.size = label.font.size,
    max.overlaps = max.overlaps, x.label = x.label, y.label = y.label,
    plot.title = plot.title, plot.subtitle = plot.subtitle,
    plot.width = plot.width, plot.height = plot.height, ...
  )
}

#' @rdname manhattan_chromosome
#' @method manhattan_chromosome MPdata
#'
#' @param x MPdata object.
#'
#' @return \code{gg} object if \code{is.null(outfn)}, \code{NULL} if \code{!is.null(outf)}
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' gwasdat <- data.frame(
#'   "chromosome" = rep(1:5, each = 30),
#'   "position" = c(replicate(5, sample(1:300, 30))),
#'   "pvalue" = rbeta(150, 1, 1)^5
#' )
#'
#' mpdata <- manhattan_data_preprocess(
#'   gwasdat, pval.colname = "pvalue", chr.colname = "chromosome", pos.colname = "position",
#'   chr.order = as.character(1:5)
#' )
#'
#' manhattan_chromosome(mpdata, chromosome = "3")
#'
manhattan_chromosome.MPdata <- function(
  x, chromosome, outfn = NULL, label.colname = NULL,
  rescale = TRUE, rescale.ratio.threshold = 5, signif.rel.pos = 0.4, color.by.highlight = FALSE,
  point.size = 0.75, label.font.size = 2, max.overlaps = 20,
  x.label = "Chromosome", y.label = expression(-log[10](p)),
  plot.subtitle = ggplot2::waiver(), plot.title = ggplot2::waiver(),
  plot.width = 10, plot.height = 5, ...
) {

  if (!is.null(label.colname)) {
    if ((sum(!is.na(x$data[[label.colname]]) & !(x$data[[label.colname]] == ""))) > 20) warning("The plot will generate > 20 labels. The plot may be cluttered & may take a longer to generate.")
  }

  if (length(chromosome) > 1) {
    stop("Only one chromosome can be selected.")
  }

  chromosome <- as.character(chromosome)
  # subset data to the chromosome
  x$data <- x$data[x$data[[x$chr.colname]] == chromosome,]
  if (nrow(x$data) == 0) {
    stop("No matching chromosome.")
  }

  # create transformation object; if rescaling is required, create the according transformation
  trans <- list("trans" = "identity", "breaks" = ggplot2::waiver())
  if (rescale) {
    jump <- get_transform_jump(-log10(x$signif))
    if ((ceiling(max(x$data[[x$pval.colname]])/5)*5)/jump > rescale.ratio.threshold) {
      trans <- get_transform(x$data, jump, x$pval.colname, jump.rel.pos = signif.rel.pos)
    }
  }

  ylimit <- c(0, ifelse(identical(trans$trans, "identity"), NA, max(trans$breaks)))

  pos <- x$true.pos.colname
  x.break <- ggplot2::waiver()
  x.break.label <- ggplot2::waiver()

  point.color <- if (!is.null(x$highlight.colname)) x$highlight else x$chr.colname

  # choose whether to use highlight.colname or chr.colname
  if (!is.null(x$highlight.colname) && !is.null(x$highlight.col) && color.by.highlight) {
    point.color <- x$highlight.colname
    point.color.map <- x$highlight.col
  } else {
    point.color <- x$chr.colname
    point.color.map <- x$chr.col
  }

  p <- ggplot2::ggplot(x$data, ggplot2::aes_string(
    x = pos,
    y = x$pval.colname,
    color = point.color
  )) +
    ggplot2::geom_point(size = point.size, pch = 16) +
    ggplot2::scale_discrete_manual(aesthetics = "color", values = point.color.map) +
    ggplot2::scale_y_continuous(trans = trans$trans,
                       breaks = trans$breaks,
                       expand = c(0.02, 0.01),
                       limits = ylimit) +
    ggplot2::scale_x_continuous(
      name = x.label,
      breaks = x.break,
      labels = x.break.label,
      expand = c(0.01, 0.01)
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(x$signif),
      linetype = 'dashed',
      color = x$signif.col
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    ggplot2::ylab(y.label) +
    ggplot2::ggtitle(label = plot.title, subtitle = plot.subtitle)

  if (!is.null(label.colname)) {
    p <- p + ggrepel::geom_label_repel(
      ggplot2::aes_string(x = x$pos.colname, y = x$pval.colname, label = label.colname),
      size = label.font.size,
      label.padding = 0.15,
      segment.size = 0.2,
      min.segment.length = 0,
      max.overlaps = max.overlaps,
      ...
    )
  }

  if (rescale & !identical(trans$trans, "identity")) {
    # if the plot is rescaled, change the tick at the "jump" to double line
    jump.tick.size <- 3.5

    p <- p +
      ggplot2::theme(
        axis.ticks.y = ggplot2::element_line(linetype = trans$y_axis_linetype)) +
      ggplot2::annotate(geom = "point", shape = "=", x = -Inf, y = trans$jump, size = jump.tick.size) +
      ggplot2::coord_cartesian(clip = "off")
  }

  if (!is.null(outfn)) {
    ggplot2::ggsave(outfn, plot=p, width=plot.width, height=plot.height, units = "in")
    invisible()
  } else {
    return(p)
  }
}
