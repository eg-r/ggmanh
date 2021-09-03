# check preprocessing arguments
preprocess_arg_check <- function(
  x, signif, signif.col, pval.colname, chr.colname, pos.colname, scale.chr.width
) {
  preprocess_checklist <- list(signif.col = NULL)
  # check significance cutoff exists
  if (length(signif) < 1) stop("At least one significance threshold should be provided.")

  # check that significance cutoff is numeric
  if (!is.numeric(signif)) stop("signif should be a numeric vector.")

  # check signif.col
  if (is.null(signif.col)) {
    preprocess_checklist$signif.col <- rep("grey", length(signif))
    preprocess_checklist$signif.col[which.max(-log10(signif))] <- "black"
  } else if (!all(valid_colors)) {
    warning("invalid signif.col colors... using default colors")
    preprocess_checklist$signif.col <- rep("grey", length(signif))
    preprocess_checklist$signif.col[which.max(-log10(signif))] <- "black"
  }

  # check that the supplied column names exist
  if (!all(c(is.character(pval.colname), is.character(chr.colname), is.character(pos.colname)))) {
    stop("Column names should be characters")
  }
  if (!all(c(pval.colname, chr.colname, pos.colname) %in% colnames(x))) {
    stop("Column name(s) not in data.frame.")
  }
  if (pval.colname == "log10pval") {
    stop("Choose a different name for pvalue column name.")
  }

  if (any(x[[pval.colname]] < 0, na.rm = TRUE) | any(x[[pval.colname]] > 1, na.rm = TRUE)) stop("p.value is a probability between 0 and 1.")


  # check that column names are valid
  if (!is.numeric(x[[pval.colname]])) stop(pval.colname, " should be a numeric column.")
  if (!is.numeric(x[[pos.colname]])) stop(pos.colname, " should be a numeric column.")

  # check that values in p value column are correct
  if (any(x[[pval.colname]] < 0, na.rm = TRUE) | any(x[[pval.colname]] > 1, na.rm = TRUE)) stop("p.value is a probability between 0 and 1.")

  if (length(scale.chr.width) != 1 | !is.logical(scale.chr.width)) {
    stop("scale.chr.width should be TRUE or FALSE.")
  }

  return(preprocess_checklist)
}

# remove entries where position, chromosome, or pvalue is missing
remove_na <- function(x, chr.colname, pos.colname, pval.colname) {
  na_remove <- which(is.na(x[[chr.colname]]) | is.na(x[[pos.colname]]) | is.na(x[[pval.colname]]))
  if (length(na_remove) > 0) {
    warning("Removed ", length(na_remove), " rows due to missing chromosome/position/pvalue.\n")
    x <- x[-na_remove,]
  }
  if (nrow(x) < 1) {
    stop("Empty rows after omitting missing chromosome/position/pvalue.\n")
  }
  return(x)
}

# TEMPORARY: remove entries where p-value is zero
remove_0_pval <- function(x, pval.colname) {
  zero_pval <- which(x[[pval.colname]] == 0)
  if (length(zero_pval) > 0) {
    warning("Removing observations with p.value = 0.")
    x <- x[-zero_pval,]
  }
  return(x)
}

set_chr_col <- function(chr.col, nchr, chr.order) {
  if (is.null(chr.col)) {
    chr.col <- stats::setNames(rep_len(RColorBrewer::brewer.pal(8, "Dark2"), nchr), chr.order)
  } else {
    if (!all(valid_colors(chr.col))) {
      warning("Invalid chr.col colors. Using default colors")
      chr.col <- stats::setNames(rep_len(RColorBrewer::brewer.pal(8, "Dark2"), nchr), chr.order)
    }
    if (!is.null(names(chr.col))) {
      if (!all(chr.order %in% names(chr.col))) {
        stop("names(chr.col) is missing values from chr.colname.")
      }
    } else {
      if (nchr > length(chr.col)) {
        warning("chr.col is recycled to match chr.order length")
        chr.col <- rep(chr.col, length.out = nchr)
        names(chr.col) <- chr.order
      } else {
        warning(paste0("Using first ", nchr, " colors for chr.col."))
        chr.col <- chr.col[1:nchr]
        names(chr.col) <- chr.order
      }
    }
  }
  return(chr.col)
}

set_highlight_col <- function(x, highlight.colname, highlight.col) {
  if (!(highlight.colname %in% colnames(x))) stop(paste0(highlight.colname, " not in data."))
  highlight.levels <- unique(x[[highlight.colname]])
  if (is.null(highlight.col)) {
    highlight.col <- RColorBrewer::brewer.pal(length(highlight.levels), "Dark2")
    names(highlight.col) <- as.character(highlight.levels)
  } else {
    if (!all(valid_colors(highlight.col))) stop("Please provide valid colors.")
    if (!is.null(names(highlight.col))) {
      if (!all(highlight.levels %in% names(highlight.col))) {
        stop("names(highlight.col) is missing values from column ", highlight.colname, ".")
      }
    } else {
      if (length(highlight.levels) > length(highlight.col)) {
        warning("highlight.col is recycled to match unique values of ", highlight.colname, ".")
        highlight.col <- rep(highlight.col, length.out = length(highlight.levels))
        names(highlight.col) <- highlight.levels
      } else if (length(highlight.levels) < length(highlight.col)) {
        warning("Using first ", length(highlight.levels), " colors.")
        highlight.col <- highlight.col[1:length(highlight.levels)]
        names(highlight.col) <- highlight.levels
      } else {
        names(highlight.col) <- highlight.levels
      }
    }
  }

  return(highlight.col)
}

# check that a character string is a valid color
valid_colors_ <- function(clr) {
  tryCatch(is.matrix(grDevices::col2rgb(clr)), error = function(clr) FALSE)
}

# vectorized version of valid_colors_
valid_colors <- function(clr) {
  vapply(clr, valid_colors_, logical(1))
}

# create spaced points of length(..)
sequence_along_chr_scaled <- function(pos) {
  pos <- pos - min(pos)
  if (max(pos) != 0) {
    pos <- pos / max(pos)
    return(pos)
  } else {
    return(pos)
  }
}

# create an equally spaced points of length(..)
sequence_along_chr_unscaled <- function(pos) {
  return(seq(from = 0, to = 1, length.out = length(pos) + 2)[-c(1,length(pos) + 2)])
}

# concatenate elements across the list
concat_list <- function(dflist, concat_char = "/") {

  if (!all(unlist(lapply(dflist, is.vector)) | unlist(lapply(dflist, is.factor)))) {
    stop("All elements in the list should be a vector.")
  }

  check_lengths <- lapply(dflist, length)
  if (length(unique(unlist(check_lengths))) != 1) {
    stop("Length of all list elements should be equal.")
  }

  if (!is.character(concat_char)) stop("concat_char should be of character type.")

  if (length(concat_char) > 1) {
    warning("concat_char should be a character vector of length 1. Using first element.")
    concat_char <- concat_char[1]
  } else if (length(concat_char) < 1) {
    warning("concat_char should be a character of length 1. Using \"/\" by default.")
  }

  if (nchar(concat_char) < 1) {
    warning("concat_char should be a character of length 1. Using \"/\" by default.")
  }

  dflist <- lapply(dflist, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    return(x)
  })
  dflist <- unname(dflist)
  dflist$sep <- concat_char
  dflist <- do.call(paste, dflist)
  dflist <- gsub(paste0("(", concat_char, ")", "+$"), "", dflist)
  dflist <- gsub(paste0("^", "(", concat_char, ")", "+"), "", dflist)
  return(dflist)
}

# concatenate columns of data.frame and produce a character vector
concat_df_cols <- function(df, concat_char = "/") {
  if (!is.data.frame(df)) stop("df should be a data.frame.")
  if (nrow(df) == 0) {
    return("")
  }
  if (length(concat_char) > 1) {
    warning("concat_char should be a character vector of length 1. Using first element.")
    concat_char <- concat_char[1]
  } else if (length(concat_char) < 1) {
    warning("concat_char should be a character of length 1. Using \"/\" by default.")
  }
  if (nchar(concat_char) < 1) {
    warning("concat_char should be a character of length 1. Using \"/\" by default.")
  }

  return(concat_list(as.list(df), concat_char))

}

# check that gds node exists
gds_node_exists <- function(gds, nodes) {
  all(nodes %in% gdsfmt::ls.gdsn(gds, recursive = TRUE))
}
