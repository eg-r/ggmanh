#' Path to Default GDS File
#'
#' Find path to the default gds file.
#'
#' @return A character vector.
#'
#' @examples
#' default_gds_path()
#'
#' @export
default_gds_path <- function() {
  message("Using default gds annotation file from ggmanh package. Run '?ggmanh_annotation_gds' for more information.")
  system.file("extdata", "gnomad.exomes.vep.hg19.v5.gds", mustWork = TRUE, package = "ggmanh")
}

# retrieve annotation with position
gds_retrieve_annot_pos <- function(
  gds, chr, pos, ref, alt,
  annotation_names = c("annotation/info/symbol", "annotation/info/consequence", "annotation/info/LoF"),
  concat_char = "/", verbose = FALSE
) {

  # check that the lengths of all parameters are equal
  if (length(chr) != length(pos) | length(chr) != length(ref) | length(chr) != length(alt)) {
    stop("The lengths of chr, pos, ref, and alt should be equal")
  }

  # check that the length of input is at least 1
  if (length(chr) == 0 | length(pos) == 0 | length(ref) == 0 | length(alt) == 0) {
    stop("The length of chr, pos, ref, and alt should be at least 1.")
  }

  # if chr is a factor, convert to character
  if (is.factor(chr)) chr <- as.character(chr)

  # check that annotation_names exist in gds
  if (!gds_node_exists(gds, annotation_names)) stop("Some of annotation_names do not exist in the gds file.")

  # filter the gds file to match chromosome, position, ref, and alt allele
  ret.idx <- SeqArray::seqSetFilterPos(gds, chr = chr, pos = pos, ref = ref, alt = alt, ret.idx = TRUE, verbose = verbose, multi.pos = TRUE)

  # retrieve the annotation
  # the order reflects the order of the filters provided
  annot_list <- SeqArray::seqGetData(gds, annotation_names)
  annot_list <- lapply(annot_list, function(x) x[ret.idx])

  if (is.null(concat_char)) {
    return(as.data.frame(annot_list))
  } else {
    return(concat_list(annot_list, concat_char))
  }

}

# retrieve annotation with rsid
gds_retrieve_annot_rsid <- function(
  gds, rs.id,
  annotation_names = c("annotation/info/symbol", "annotation/info/consequence", "annotation/info/LoF"),
  concat_char = "/", verbose = FALSE
) {
  # check that the lengths of all parameters are equal
  if (length(rs.id) < 1) {
    stop("The length of rs.id should be at least 1.")
  }

  # check that annotation_names exist in gds
  if (!gds_node_exists(gds, annotation_names)) stop("Some of annotation_names do not exist in the gds file.")

  # filter the gds file to match chromosome, position, ref, and alt allele
  ret.idx <- SeqArray::seqSetFilterAnnotID(gds, id = rs.id, ret.idx = TRUE, verbose = verbose)

  if (length(ret.idx) != length(SeqArray::seqGetData(gds, "annotation/id"))) {
    warning("Some variants are multiallelic; only the first allele available in the GDS file is selected. For more precise annotation, try annotating by position with annot.method = \"position\".")
  }

  # retrieve the annotation
  # the order reflects the order of the filters provided
  annot_list <- SeqArray::seqGetData(gds, annotation_names)
  annot_list <- lapply(annot_list, function(x) x[ret.idx])

  if (is.null(concat_char)) {
    return(as.data.frame(annot_list))
  } else {
    return(concat_list(annot_list, concat_char))
  }

}

#' Annotation with GDS File
#'
#' Retrieve variant annotation stored in a GDS file with chromosome location or rs.id.
#'
#' @param x a \code{data.frame} object to be annotated.
#' @param gdsfile a character for GDS filename. If \code{NULL}, the default GDS file included with the package is used.
#' @param annot.method a method for searching variants. "position" requires \code{chr}, \code{pos}, \code{ref}, and \code{alt}. "rs.id" requires \code{rs.id}.
#' @param chr,pos,ref,alt,rs.id column names of \code{x} that contain chromosome, position, reference allele, alternate allele, and rs.id, respectively.
#' @param concat_char a character used to separate multiple annotations returned from the gds file.
#' @param annotation_names a character vector of nodes of the \code{gdsfile} that are to be extracted.
#' @param verbose output messages.
#'
#' @return A character vector the length of \code{nrow(x)} if \code{concat_char} is a character.
#' A data frame with \code{nrow(x)} rows and \code{length(annotation_names)} if \code{concat_char} is null.
#'
#' @examples
#' vardata <- data.frame(
#'   chr = c(11,20,14),
#'   pos = c(12261002, 10033792, 23875025),
#'   ref = c("G", "G", "CG"),
#'   alt = c("A", "A", "C")
#' )
#'
#' annotations <- gds_annotate(
#'   x = vardata, annot.method = "position",
#'   chr = "chr", pos = "pos", ref = "ref", alt = "alt"
#' )
#'
#' print(annotations)
#'
#' @export
gds_annotate <- function(
  x, gdsfile = NULL, annot.method = "position",
  chr = NULL, pos = NULL, ref = NULL, alt = NULL, rs.id = NULL,
  concat_char = "/", verbose = TRUE,
  annotation_names = c("annotation/info/symbol", "annotation/info/consequence", "annotation/info/LoF")
) {

  # check that x is data.frame
  if (!is.data.frame(x)) stop("x should be a data.frame.")

  # annot.method is either "position" or "rs.id"
  if (!(annot.method %in% c("position", "rs.id"))) stop("annot.method should be \"position\" or \"rs.id\".")

  # check & open gdsfile
  if (is.null(gdsfile)) {
    gdsfile <- default_gds_path()
  }

  stopifnot(file.exists(gdsfile))
  gds <- SeqArray::seqOpen(gdsfile, readonly = TRUE)

  # close gds on exit
  on.exit(SeqArray::seqClose(gds))

  if (annot.method == "position") {
    if (any(!(as.character(c(chr, pos, ref, alt)) %in% colnames(x)))) stop("One of chr, position, ref, alt is not a column name of x.")
    annot <-
      gds_retrieve_annot_pos(
        gds, chr = x[[chr]], pos = x[[pos]], ref = x[[ref]], alt = x[[alt]],
        annotation_names = annotation_names, concat_char = concat_char,
        verbose = verbose
      )
  } else if (annot.method == "rs.id") {
    if (is.null(gdsfile)) stop("Default gds file does not have rs.id.")
    if (any(!(as.character(rs.id) %in% colnames(x))))
    annot <- gds_retrieve_annot_rsid(
      gds, rs.id = x[[rs.id]],
      annotation_names = annotation_names, concat_char = concat_char,
      verbose = verbose
    )
  } else {
    stop("annot.method should be either \"position\" or \"rs.id\".")
  }

  return(annot)
}
