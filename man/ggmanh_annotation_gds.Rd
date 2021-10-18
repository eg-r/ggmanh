% man page for inst/extdata/gnomad.exomes.vep.hg19.v5.gds
\name{ggmanh_annotation_gds}
\docType{data}
\alias{ggmanh_annotation_gds}
\title{gnomAD Variant Annotation in SeqArray Format}
\format{
A GDS file with 1015430 variants with chromosome, position, allele, gene symbol, Ensembl VEP Consequence, and predicted LoF.
}
\description{
\code{ggmanh} provides a GDS file whose path is accessible by \code{default_gds_path}.
The original annotation file is from the gnomAD browser v2.1.1 release, available in this link: \url{https://gnomad.broadinstitute.org/downloads}.
This gds file contains variants in the exome with the global minor allele frequency \eqn{\ge} 0.0002, and has been manually curated
to fit the file size requirement for R Bioconductor packages.
}
