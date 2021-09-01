#' Thin Data Points
#'
#' Reduce the number of cluttered data points.
#'
#' @param dat a data frame
#' @param value column name of \code{dat} to be used for partitioning (see details)
#' @param n number of points to sample for each partition
#' @param nbins number of partitions
#' @param groupBy column name of \code{dat} to group by before partitioning
#'
#' @details The result of Genome Wide Association Study can be very large, with the majority
#'   of points being being clustered below significance threshold. This unnecessarily increases the time
#'   to plot while making almost no difference. This function reduces the number of points by partitioning the points
#'   by a numberic column \code{value} into \code{nbins} and sampling \code{n} points.
#'
#' @return a \code{data.frame}
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' library(dplyr)
#' dat <- data.frame(
#'    A1 = c(1:20, 20, 20),
#'    A2 = c(rep(1, 12), rep(1,5), rep(20, 3), 20, 20) ,
#'    B = rep(c("a", "b", "c", "d"), times = c(5, 7, 8, 2))
#' )
#' # partition "A1" into 2 bins and then sample 6 data points
#' thinPoints(dat, value = "A1", n = 6, nbins = 2)
#' # partition "A2" into 2 bins and then sample 6 data points
#' thinPoints(dat, value = "A2", n = 6, nbins = 2)
#' # group by "B", partition "A2" into 2 bins and then sample 3 data points
#' thinPoints(dat, value = "A2", n = 3, nbins = 2, groupBy = "B")
#'
thinPoints <- function(dat, value, n=3000, nbins=200, groupBy=NULL)
{
  if (!(value %in% colnames(dat)) || !is.numeric(dat[[value]])) stop("value should be name of a numeric column.")
  if (!is.numeric(n) || !is.numeric(nbins)) stop("n and nbins should numbers.")
  if (!is.character(value) || (length(value) != 1)) stop("value should be a string.")

  if (!is.null(groupBy))
  {
    if (!is.character(groupBy) || (length(groupBy) != 1)) stop("groupBy should be a string.")
    dat <- dplyr::group_by(dat, .data[[groupBy]])
  }

  dat %>%
    dplyr::mutate(bin = cut(.data[[value]], breaks = nbins, labels = FALSE)) %>%
    dplyr::group_by(.data[["bin"]], .add=TRUE) %>%
    dplyr::slice_sample(n = n, replace = FALSE) %>%
    dplyr::ungroup() %>% # undo all grouping
    dplyr::select(-.data[["bin"]])

}
