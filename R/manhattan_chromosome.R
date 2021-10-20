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
manhattan_chromosome <- function(x, ...) UseMethod("manhattan_chromosome")

#' @rdname manhattan_chromosome
#' @export
manhattan_chromosome.default <- function(x, ...) stop("Provide a data.frame to preprocess & plot or a preprocessed MPdata object.")

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
    preserve.position = preserve.position, thin = thin, thin.n = thin.n
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

#' @rdname manhattan_chromosome
#' @method manhattan_chomosome GRanges
#' @export
setMethod(
  "manhattan_chromosome", signature = "GRanges",
  function(
    x, chromosome, outfn = NULL, signif = c(5e-8, 1e-5),
    pval.colname = "pval",
    chr.order = NULL, signif.col = NULL, chr.col = NULL, highlight.colname = NULL,
    preserve.position = FALSE, thin = FALSE, thin.n = 3000,
    rescale = TRUE, rescale.ratio.threshold = 5,
    signif.rel.pos = 0.4, label.colname = NULL,
    point.size = 0.75, x.label = "Chromosome",
    y.label = expression(-log[10](p)), label.font.size = 2, max.overlaps = 20,
    plot.title = ggplot2::waiver(),
    plot.subtitle = ggplot2::waiver(), plot.width = 10, plot.height = 5, ...
  ) {

    if (length(chromosome) > 1) {
      stop("Only one chromosome can be selected.")
    }

    if (!(chromosome %in% as.character(GenomicRanges::runValue(GenomicRanges::seqnames(x))))) stop("Chromosome does not exist in x.")

    x <- x[GenomicRanges::seqnames(x) == chromosome]

    mpdata <- manhattan_data_preprocess(
      x, signif = signif, pval.colname = pval.colname,
      chr.order = chromosome, signif.col = signif.col, chr.col = chr.col,
      highlight.colname = highlight.colname, preserve.position = preserve.position, thin = thin,
      thin.n = thin.n
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
)
