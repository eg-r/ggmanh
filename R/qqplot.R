#' Plot Quantile-Quantile Plot of p-values against uniform distribution.
#'
#' @param x a numeric vector of p-values. All values should be between 0 and 1.
#' @param outfn a character. File name to save the QQ Plot. If \code{outfn}
#'  is supplied (i.e. \code{!is.null(outfn)}), then the plot is not drawn in
#'  the graphics window.
#' @param conf.int a numeric between 0 and 1. Confidence band to draw around reference line. Set to \code{NA} to leave it out.
#' @param plot.width a numeric. Plot width in inches.
#' @param plot.height a numeric. Plot height in inches.
#' @param thin a logical. Reduce number of data points when they are cluttered?
#' @param thin.n an integer. Number of max points per horizontal partitions of the plot.
#'   Defaults to 500.
#' @param zero.pval a character. Determine how to treat 0 pvals.
#' "replace" will replace the p-value of zero with the non-zero minimum.
#' "remove" will remove the p-value of zero.
#' @param pval.log a logical. If \code{FALSE}, p-values will assumed to already be log10-transformed, otherwise they will be converted. Default is \code{TRUE}.
#'
#' @return a ggplot object
#'
#' @export
#' @examples
#' x <- rbeta(1000, 1, 1)
#' qqunif(x)
qqunif <- function(
  x, outfn = NULL, conf.int = 0.95,
  plot.width = 5, plot.height = 5, thin = TRUE, thin.n = 500,
  zero.pval = "replace", pval.log = TRUE
) {

  if (!is.atomic(x) || !is.numeric(x)) stop("x should be a numeric vector.")
  if (any((x < 0) | (x > 1))) stop("x should be between 0 and 1")
  if (any(is.na(x))) {
    warning("Removing NAs")
    x <- x[!is.na(x)]
  }
  if (length(conf.int) != 1 || !is.numeric(conf.int)) stop("conf.int should be a numerical value.")
  if (conf.int <= 0 || conf.int >= 1) stop("conf.int should be between 0 and 1.")
  if ((!is.logical(thin) || (length(thin) != 1)) || is.na(thin)) stop("thin should be a logical of length 1.")
  if (!(is.numeric(thin.n) && length(thin.n) == 1)) stop("thin.n should be an integer of length 1.")

  if (zero.pval == "replace") {
    x <- replace_0_pval(x, pval.log = pval.log)
  } else if (zero.pval == "remove") {
    x <- remove_0_pval(x, pval.log = pval.log)
  } else {
    stop("zero.pval should be \"replace\" or \"remove\"")
  }

  x <- sort(x)
  N <- length(x)
  q <- 1:N / N
  qqdf <- data.frame(obs = ifelse(pval.log, -log10(x), x), exp = -log10(q))
  if (thin) {
    qqdf <- thinPoints(qqdf, value = "obs", n = thin.n)
  }

  ## QQ plot using ggplot
  qqtitle <- "QQ Plot"
  p <- ggplot2::ggplot(qqdf)

  # calculate confidence bands
  if (!is.na(conf.int) && !is.null(conf.int)) {
    conf.int.list <- data.frame(
      K <- c(seq(0, -log10(1/N), length.out = round(-log10(1/N)) * 100), -log10(0.9/N))
    )
    ranks <- (10^-K) * N
    conf.int.list$p_025 <- -log10(stats::qbeta((1 - conf.int)/2, ranks, N + 1 - ranks))
    conf.int.list$p_975 <- -log10(stats::qbeta(1/2 + conf.int/2, ranks, N - ranks))

    p <- p +
      ggplot2::geom_line(ggplot2::aes_string(x = "K", y = "p_025"), data = conf.int.list, lty = "dashed", lwd = 0.5) +
      ggplot2::geom_line(ggplot2::aes_string(x = "K", y = "p_975"), data = conf.int.list, lty = "dashed", lwd = 0.5)

    qqtitle <- sprintf("%s with %d%% CI", qqtitle, conf.int * 100)
  }

  p <- p +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "red") +
    ggplot2::geom_point(ggplot2::aes_string(x = "exp", y = "obs"), size = 0.9) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    ggplot2::xlab(expression(paste(Expected, " ", -Log[10], "(pval)"))) +
    ggplot2::ylab(expression(paste(Observed, " ", -Log[10], "(pval)"))) +
    ggplot2::ggtitle(qqtitle)

  if (!is.null(outfn)) {
    ggplot2::ggsave(outfn, plot=p, width=plot.width, height=plot.height, units = "in")
    invisible()
  } else {
    return(p)
  }

}
