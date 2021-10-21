#' ggmanh: A package for visualization of GWAS results.
#'
#' ggmanh provides flexible tools for visualizing GWAS result for downstream analysis.
#'
#' @details
#' Manhattan plot is commonly used to display significant Single Nucleotide Polymorphisms (SNPs)
#' in Genome Wide Association Study (GWAS)
#' This package comes with features useful for manhattan plot creation, including annotation with \code{\link{ggrepel}},
#' truncating data for faster plot generation, and manual rescaling of the y-axis.
#' The manhattan plot is generated in two steps: data preprocessing and plotting. This allows the user to iteratively
#' customize the plot without having the process the GWAS summary data over and over again.
#' Currently, \code{data.frame} and \code{GRanges} from \code{GenomicRanges} are supported.
#'
#' A vignette detailing the usage of the package is accessible by \code{vignette("ggmanh")}
#'
#' @docType package
#' @import methods
#' @name ggmanh
NULL
