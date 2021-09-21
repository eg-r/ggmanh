#' gnomAD Variant Annotation in SeqArray Format
#'
#' \code{ggmanh} provides a GDS file whose path is accessible by \code{default_gds_path}.
#' The original annotation file is from the gnomAD browser v2.1.1 release, available in this link: \url{https://gnomad.broadinstitute.org/downloads}.
#' This gds file contains variants in the exome with the global minor allele frequency \eqn{\ge} 0.0002, and has been manually curated
#' to fit the file size requirement for R Bioconductor packages.
#'
#' @format A GDS file with 1015430 variants with chromosome, position, allele, gene symbol, Ensembl VEP Consequence, and predicted LoF.
#' @name ggmanh_annotation_gds
NULL
