# Shared internal helpers ------------------------------------------------------

# Wrap column names in backticks so that any name works inside a formula.
backtick <- function(x) {
  if (length(x) == 0L) {
    return(character())
  }
  paste0("`", gsub("`", "\\\\`", x), "`")
}

# Observed levels of one stressor as text, with the control first. Factors keep
# their level order, numbers are sorted numerically, text alphabetically.
stressor_levels <- function(x, control) {
  if (is.factor(x)) {
    lv <- levels(droplevels(x))
  } else if (is.numeric(x)) {
    lv <- as.character(sort(unique(x)))
  } else {
    lv <- sort(unique(as.character(x)))
  }
  c(control, setdiff(lv, control))
}

# Label of a treatment combination: "temperature=25 | nitrate=high".
# `levels` is a data frame or a named list with one element per stressor.
make_cell_label <- function(levels, stressors) {
  parts <- lapply(stressors, function(s) paste0(s, "=", levels[[s]]))
  do.call(paste, c(parts, sep = " | "))
}

# Short label with only the stressors away from control:
# "herbicide=diuron + nutrients=added".
compact_levels <- function(labels, control) {
  at_control <- paste0(names(control), "=", control)
  vapply(strsplit(labels, " | ", fixed = TRUE), function(parts) {
    paste(parts[!parts %in% at_control], collapse = " + ")
  }, character(1))
}

# Labels for messages: only the stressors away from control ("control" if
# none), preceded by the `by` group when there is one.
short_cells <- function(cell, control, context = "") {
  short <- compact_levels(cell, control)
  short[!nzchar(short)] <- "control"
  context <- rep_len(context, length(short))
  with_context <- nzchar(context)
  short[with_context] <- paste0(context[with_context], ": ", short[with_context])
  short
}

# Label of the context defined by the `by` columns ("" when there are none).
make_context_label <- function(data, by) {
  if (length(by) == 0L) {
    return(rep("", nrow(data)))
  }
  parts <- lapply(by, function(b) paste0(b, "=", as.character(data[[b]])))
  do.call(paste, c(parts, sep = " | "))
}

# Key that never collides with the content of the labels.
paste_key <- function(...) {
  paste(..., sep = "\r")
}

# Number of distinct values of `x` within each group of `g` (named vector).
n_distinct_by <- function(x, g) {
  vapply(split(x, g), function(v) length(unique(v)), integer(1))
}

# Format a range of integers as "3" or "2-5".
format_range <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return("none")
  }
  if (min(x) == max(x)) as.character(min(x)) else paste0(min(x), "-", max(x))
}

# All subsets of the positions 1..m, ordered by size (empty set first).
all_subsets <- function(m) {
  out <- list(integer(0))
  for (size in seq_len(m)) {
    out <- c(out, utils::combn(m, size, simplify = FALSE))
  }
  out
}

# Bind a list of data frames that may be empty.
bind_rows_base <- function(x, template = NULL) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0L) {
    return(template)
  }
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}
