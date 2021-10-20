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
manhattan_plot <- function(x, ...) UseMethod("manhattan_plot")

#' @rdname manhattan_plot
#' @export
manhattan_plot.default <- function(x, ...) stop("Provide a data.frame to preprocess & plot or a preprocessed MPdata object.")

#' @rdname manhattan_plot
#' @method manhattan_plot data.frame
#' @export
manhattan_plot.data.frame <- function(
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
    x, signif = signif, pval.colname = pval.colname, chr.colname = chr.colname, pos.colname = pos.colname,
    chr.order = chr.order, signif.col = signif.col, chr.col = chr.col, highlight.colname = highlight.colname,
    highlight.col = highlight.col, preserve.position = preserve.position, thin = thin,
    thin.n = thin.n
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
manhattan_plot.MPdata <- function(
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

#' @rdname manhattan_plot
#' @method manhattan_plot GRanges
#' @export
setMethod(
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
      highlight.col = highlight.col, preserve.position = preserve.position, thin = thin,
      thin.n
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
