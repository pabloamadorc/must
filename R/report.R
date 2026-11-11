# Methods paragraph and table of decisions ------------------------------------
#
# Everything in the report comes from the settings of the objects, so that
# the text describes the analysis that was run (roadmap, line 4: the paragraph
# must agree with the table of decisions and with `settings`).

#' Report the methods and the decisions of an analysis
#'
#' @description
#' Writes, from the output of [estimate_interactions()], the methods of the
#' statistical analysis as paragraphs of text and a table with every decision
#' behind it: the design, the model of each response, the null models and why
#' they were chosen, the inference and, when given, the classification of the
#' interactions and the agreement among null models. The text cites the
#' methods used and lists their references. Every value comes from the
#' settings of the objects, so the description matches the analysis that was
#' run.
#'
#' @param x The output of [estimate_interactions()].
#' @param classes Optional: the output of [classify_interactions()] for `x`,
#'   to describe the classification of the interactions.
#' @param agreement Optional: the output of [compare_null_models()] for `x`,
#'   to describe the agreement among the null models.
#' @param justification Optional text, one or more sentences, that justifies
#'   the choice of null models; it is added to the methods as given. When the
#'   null models came from [suggest_null_models()], the criteria behind them
#'   are described without it.
#'
#' @details
#' The methods have up to five paragraphs: the design (stressors, units,
#' times and the interactions the design identifies); the software and the
#' model of each response (for a model fitted with [fit_responses()], the
#' candidates for its error distribution, how one was chosen and the random
#' effects; for a model fitted outside must, its class and family); the null
#' models and why they were chosen (the criteria of [suggest_null_models()],
#' the `justification` given or, for the default, why it is the default);
#' how the interactions were estimated and tested; and, when `classes` or
#' `agreement` are given, how the interactions were classified and compared.
#' They are written in English, in the past tense, for the methods section of
#' a manuscript, and are meant to be edited: the description of the
#' experiment itself (organisms, exposure, measurements) is left to the
#' authors.
#'
#' The table of decisions has one row per decision, with the function
#' (`step`) and argument that control it, the value used, the default and
#' its `source`: `"default"`, `"set"` (given by the user), `"chosen"` (by
#' must: the distribution selected by [fit_responses()], a seed drawn at
#' random), `"estimated"` (criteria of [suggest_null_models()] evaluated from
#' the effects estimated), `"suggest_null_models()"` (the null models
#' suggested) or `"recorded"` (versions of the software). The versions are
#' those recorded when the models were fitted and the interactions estimated;
#' objects created by earlier versions of must lack them and the summary of
#' the models, and the report then says less.
#'
#' @return An object of class `must_report`, a list with `methods` (the
#'   paragraphs, a character vector), `decisions` (a data frame with columns
#'   `step`, `decision`, `argument`, `value`, `default` and `source`) and
#'   `references` (the works cited, formatted, a character vector). `format()`
#'   gives the whole report as lines of plain text or, with `format =
#'   "markdown"`, of markdown for an R Markdown or Quarto document (in a chunk
#'   with `results = "asis"`); `print()` shows the plain text.
#'
#' @seealso [estimate_interactions()], [suggest_null_models()],
#'   [classify_interactions()], [compare_null_models()].
#'
#' @examples
#' jars <- expand.grid(
#'   salinity = c("ambient", "high"),
#'   warming = c("no", "yes"),
#'   replicate = 1:5,
#'   stringsAsFactors = FALSE
#' )
#' jars$jar <- seq_len(nrow(jars))
#' mu <- 40 * ifelse(jars$salinity == "high", 0.5, 1) *
#'   ifelse(jars$warming == "yes", 0.5, 1)
#' jars$algae <- stats::qpois(c(0.2, 0.35, 0.5, 0.65, 0.8)[jars$replicate], mu)
#' design <- describe_design(jars, c("salinity", "warming"),
#'   control = c(salinity = "ambient", warming = "no"), unit = "jar")
#' fits <- fit_responses(jars, "algae", design, residual_checks = FALSE,
#'   seed = 1, quiet = TRUE)
#' est <- estimate_interactions(fits, null_models = c("AD", "MU"), n_sim = 2000,
#'   seed = 1)
#' report <- report_interactions(est, classes = classify_interactions(est))
#' report
#' report$decisions[, c("decision", "value", "source")]
#' # Markdown, for an R Markdown document (in a chunk with results = "asis")
#' cat(format(report, format = "markdown"), sep = "\n")
#' @export
report_interactions <- function(x, classes = NULL, agreement = NULL, justification = NULL) {
  if (!inherits(x, "must_interactions")) {
    abort_input(sprintf(
      "`x` must be the output of estimate_interactions(), not %s.", describe_value(x)
    ))
  }
  responses <- unique(x$estimates$response)
  if (!is.null(classes) && !(is.data.frame(classes) &&
                             all(c("class", "verdict", "null_model") %in% names(classes)) &&
                             all(classes$response %in% responses) &&
                             all(classes$null_model %in% x$settings$null_models))) {
    abort_input(paste(
      "`classes` must be NULL or the output of classify_interactions() for `x`",
      "(the same responses and null models)."
    ))
  }
  if (!is.null(agreement) && !(is.data.frame(agreement) &&
                               all(c("agreement", "direction") %in% names(agreement)) &&
                               all(agreement$response %in% responses) &&
                               all(compared_null_models(agreement) %in% x$settings$null_models))) {
    abort_input(paste(
      "`agreement` must be NULL or the output of compare_null_models() for `x`",
      "(the same responses and null models)."
    ))
  }
  if (!is.null(justification)) {
    if (!is.character(justification) || length(justification) == 0L || anyNA(justification) ||
        !any(nzchar(trimws(justification)))) {
      abort_input("`justification` must be NULL or text (a character string).")
    }
    justification <- paste(trimws(justification), collapse = " ")
  }
  refs <- new_citations()
  null_models <- report_null_models(x, justification, refs)
  estimation <- report_estimation(x, refs)
  # Without reasons for the null models, their sentence opens the paragraph on
  # the estimation.
  tested <- if (length(null_models) == 1L) {
    paste(null_models, estimation)
  } else {
    c(paste(null_models, collapse = " "), estimation)
  }
  methods <- c(
    report_design(x),
    report_models(x, refs),
    tested,
    report_classes(x, classes, agreement, refs)
  )
  structure(
    list(
      methods = methods[nzchar(methods)],
      decisions = report_decisions(x, classes, agreement, justification),
      references = refs$formatted()
    ),
    class = "must_report"
  )
}

#' @export
format.must_report <- function(x, format = "text", width = 80L, ...) {
  check_choice(format, c("text", "markdown"))
  check_number(width, lower = 20)
  d <- x$decisions
  if (identical(format, "markdown")) {
    cell <- function(v) gsub("|", "\\|", v, fixed = TRUE)
    table <- c(
      "| Step | Decision | Argument | Value | Default | Source |",
      "|---|---|---|---|---|---|",
      sprintf("| %s | %s | %s | %s | %s | %s |",
              ifelse(grepl("()", d$step, fixed = TRUE), paste0("`", d$step, "`"), d$step),
              cell(d$decision),
              ifelse(nzchar(d$argument), paste0("`", cell(d$argument), "`"), ""),
              cell(d$value), cell(d$default), cell(d$source))
    )
    return(c(
      unlist(lapply(x$methods, function(p) c(p, ""))),
      table, "",
      "**References**", "",
      unlist(lapply(x$references, function(r) c(r, "")))
    ))
  }
  wrap <- function(s, ...) strwrap(strip_markup(s), width = width, ...)
  decisions <- unlist(lapply(unique(d$step), function(st) {
    part <- d[d$step == st, , drop = FALSE]
    c(st, unlist(lapply(seq_len(nrow(part)), function(i) {
      where <- paste(c(if (nzchar(part$argument[[i]])) part$argument[[i]], part$source[[i]]),
                     collapse = ", ")
      strwrap(sprintf("%s: %s [%s]", part$decision[[i]], part$value[[i]], where),
              width = width, indent = 2L, exdent = 4L)
    })), "")
  }))
  c(
    "Methods", "",
    unlist(lapply(x$methods, function(p) c(wrap(p), ""))),
    "Decisions", "",
    decisions,
    "References", "",
    unlist(lapply(x$references, function(r) c(wrap(r, exdent = 2L), "")))
  )
}

# The null models compared in the output of compare_null_models().
compared_null_models <- function(agreement) {
  compared <- attr(agreement, "null_models")
  if (is.null(compared)) {
    compared <- intersect(c("AD", "MU", "DO", "BL", "SAM", "CA"), names(agreement))
  }
  compared
}

#' @export
print.must_report <- function(x, ...) {
  cat("<must_report>", format(x, format = "text"), sep = "\n")
  invisible(x)
}

# Citations ---------------------------------------------------------------------

# Works cited: cite(key) gives the citation inside parentheses ("Madge
# Pimentel et al. 2025"), citet(key) the narrative one ("Madge Pimentel et al.
# (2025)"), and both record the work for the list of references.
new_citations <- function() {
  env <- new.env(parent = emptyenv())
  env$entries <- list()
  entry <- function(key) {
    if (is.null(env$entries[[key]])) {
      env$entries[[key]] <- report_reference(key)
    }
    env$entries[[key]]
  }
  list(
    cite = function(key) {
      e <- entry(key)
      if (is.na(e$year)) e$who else paste(e$who, e$year)
    },
    citet = function(key) {
      e <- entry(key)
      if (is.na(e$year)) e$who else sprintf("%s (%s)", e$who, e$year)
    },
    formatted = function() {
      full <- vapply(env$entries, `[[`, character(1), "full")
      unname(full[order(tolower(strip_markup(full)), method = "radix")])
    }
  )
}

# Text without the emphasis of markdown (`*...*`, `_..._`).
strip_markup <- function(x) {
  x <- gsub("\\*([^*]+)\\*", "\\1", x)
  gsub("(?<![[:alnum:]])_([^_]+)_(?![[:alnum:]])", "\\1", x, perl = TRUE)
}

report_reference <- function(key) {
  static <- report_static_references()
  if (!is.null(static[[key]])) {
    e <- static[[key]]
    return(list(who = e[[1L]], year = e[[2L]], full = e[[3L]]))
  }
  switch(
    key,
    must = must_reference(),
    R = bib_reference(utils::citation()),
    DHARMa = package_reference("DHARMa", "Hartig", paste(
      "Hartig, F. (n.d.). *DHARMa: Residual diagnostics for hierarchical",
      "(multi-level/mixed) regression models*. R package.",
      "<https://CRAN.R-project.org/package=DHARMa>"
    )),
    emmeans = package_reference("emmeans", "Lenth", paste(
      "Lenth, R. V. (n.d.). *emmeans: Estimated marginal means, aka least-squares",
      "means*. R package. <https://CRAN.R-project.org/package=emmeans>"
    )),
    abort_internal(sprintf("no reference for '%s'.", key))
  )
}

# Who, year and full reference of a bibentry (its first entry), the full one
# in the style of the other references for manuals.
bib_reference <- function(bib) {
  b <- bib[[1L]]
  family <- vapply(b$author, function(p) {
    paste(if (length(p$family) > 0L) p$family else p$given, collapse = " ")
  }, character(1))
  who <- if (length(family) == 1L) {
    family
  } else if (length(family) == 2L) {
    paste(family, collapse = " & ")
  } else {
    paste(family[[1L]], "et al.")
  }
  year <- if (is.null(b$year)) "n.d." else as.character(b$year)
  list(who = who, year = year, full = bib_full(b))
}

bib_full <- function(b) {
  if (!tolower(attr(unclass(b)[[1L]], "bibtype")) %in% c("manual", "misc")) {
    return(gsub("\\s+", " ", paste(format(b, style = "text"), collapse = " ")))
  }
  authors <- vapply(b$author, function(p) {
    if (length(p$family) == 0L) {
      return(paste(p$given, collapse = " "))
    }
    initials <- paste0(substr(p$given, 1L, 1L), ".", collapse = " ")
    family <- paste(p$family, collapse = " ")
    if (length(p$given) == 0L) family else sprintf("%s, %s", family, initials)
  }, character(1))
  authors <- if (length(authors) == 1L) {
    authors
  } else {
    paste0(paste(authors[-length(authors)], collapse = ", "), ", & ", authors[[length(authors)]])
  }
  publisher <- c(b$organization, b$address)
  pieces <- c(
    sprintf("%s (%s).", authors, if (is.null(b$year)) "n.d." else b$year),
    sprintf("*%s*.", sub("\\.$", "", b$title)),
    if (!is.null(b$note)) sprintf("%s.", sub("\\.$", "", b$note)),
    if (length(publisher) > 0L) sprintf("%s.", paste(publisher, collapse = ", ")),
    if (!is.null(b$url)) sprintf("<%s>", b$url)
  )
  paste(pieces, collapse = " ")
}

# The citation of an installed package, or a fallback without a year.
package_reference <- function(pkg, who, fallback) {
  bib <- if (is_installed(pkg)) {
    tryCatch(utils::citation(pkg), error = function(e) NULL)
  } else {
    NULL
  }
  if (is.null(bib)) list(who = who, year = "n.d.", full = fallback) else bib_reference(bib)
}

# The citation of must, from its CITATION file (its first entry).
must_reference <- function() {
  file <- system.file("CITATION", package = "must")
  bib <- if (nzchar(file)) {
    tryCatch(
      utils::readCitationFile(file, meta = list(Version = must_version(), Encoding = "UTF-8")),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  if (is.null(bib)) {
    return(list(who = "Amador Crespo", year = "n.d.", full = sprintf(paste(
      "Amador Crespo, P. (n.d.). *must: Multiple-Stressor Interactions Under",
      "Alternative Null Models*. R package version %s.",
      "<https://github.com/pabloamadorc/must>"
    ), must_version())))
  }
  bib_reference(bib)
}

must_version <- function() {
  package_version_text("must")
}

package_version_text <- function(pkg) {
  tryCatch(as.character(getNamespaceVersion(pkg)), error = function(e) {
    tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA_character_)
  })
}

# Versions recorded with the objects: must, R and the packages given.
analysis_versions <- function(packages = character()) {
  c(must = must_version(), R = paste(R.version$major, R.version$minor, sep = "."),
    vapply(packages, package_version_text, character(1)))
}

# What the report needs from the models of fit_responses(), kept in the
# settings of estimate_interactions() (the models themselves are not).
fit_summary <- function(fit) {
  s <- fit$settings
  responses <- lapply(fit$fits, function(f) {
    ok <- identical(f$status, "ok")
    candidates <- f$selection
    if (!is.null(candidates)) {
      candidates <- candidates[, intersect(c("family", "link", "converged", "selected"),
                                           names(candidates)), drop = FALSE]
      given <- s$candidate_families[[f$response]]
      candidates$distribution <- if (is.null(given)) {
        candidates$family
      } else {
        unname(given[candidates$family])
      }
    }
    list(
      response = f$response, status = f$status, reason = f$reason, family = f$family,
      family_name = if (ok) model_family(f$model, f$family) else NA_character_,
      link = f$link, response_type = f$response_type, type_source = f$type_source,
      kind = f$kind, formula = f$formula, candidates = candidates,
      n_boundary = length(f$boundary), checks_passed = f$checks_passed,
      reml = !is.null(f$model_reml)
    )
  })
  c(s, list(responses = unname(responses)))
}

model_family <- function(model, fallback) {
  out <- tryCatch(stats::family(model)$family, error = function(e) NULL)
  if (is.character(out) && length(out) == 1L) out else fallback
}

report_static_references <- function() {
  list(
    anscombe1948 = list("Anscombe", "1948", paste(
      "Anscombe, F. J. (1948). The transformation of Poisson, binomial and",
      "negative-binomial data. *Biometrika*, 35(3\u20134), 246\u2013254.",
      "<https://doi.org/10.1093/biomet/35.3-4.246>"
    )),
    benjamini1995 = list("Benjamini & Hochberg", "1995", paste(
      "Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate:",
      "A practical and powerful approach to multiple testing. *Journal of the",
      "Royal Statistical Society: Series B*, 57(1), 289\u2013300.",
      "<https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>"
    )),
    benjamini2001 = list("Benjamini & Yekutieli", "2001", paste(
      "Benjamini, Y., & Yekutieli, D. (2001). The control of the false discovery",
      "rate in multiple testing under dependency. *The Annals of Statistics*,",
      "29(4), 1165\u20131188. <https://doi.org/10.1214/aos/1013699998>"
    )),
    beppler2016 = list("Beppler et al.", "2016", paste(
      "Beppler, C., Tekin, E., Mao, Z., White, C., McDiarmid, C., Vargas, E.,",
      "Miller, J. H., Savage, V. M., & Yeh, P. J. (2016). Uncovering emergent",
      "interactions in three-way combinations of stressors. *Journal of the Royal",
      "Society Interface*, 13, 20160800. <https://doi.org/10.1098/rsif.2016.0800>"
    )),
    bliss1939 = list("Bliss", "1939", paste(
      "Bliss, C. I. (1939). The toxicity of poisons applied jointly. *Annals of",
      "Applied Biology*, 26(3), 585\u2013615.",
      "<https://doi.org/10.1111/j.1744-7348.1939.tb06990.x>"
    )),
    brooks2017 = list("Brooks et al.", "2017", paste(
      "Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg, C.",
      "W., Nielsen, A., Skaug, H. J., M\u00e4chler, M., & Bolker, B. M. (2017).",
      "glmmTMB balances speed and flexibility among packages for zero-inflated",
      "generalized linear mixed modeling. *The R Journal*, 9(2), 378\u2013400.",
      "<https://doi.org/10.32614/RJ-2017-066>"
    )),
    hochberg1988 = list("Hochberg", "1988", paste(
      "Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of",
      "significance. *Biometrika*, 75(4), 800\u2013802.",
      "<https://doi.org/10.1093/biomet/75.4.800>"
    )),
    holm1979 = list("Holm", "1979", paste(
      "Holm, S. (1979). A simple sequentially rejective multiple test procedure.",
      "*Scandinavian Journal of Statistics*, 6, 65\u201370.",
      "<https://www.jstor.org/stable/4615733>"
    )),
    hommel1988 = list("Hommel", "1988", paste(
      "Hommel, G. (1988). A stagewise rejective multiple test procedure based on a",
      "modified Bonferroni test. *Biometrika*, 75(2), 383\u2013386.",
      "<https://doi.org/10.1093/biomet/75.2.383>"
    )),
    krinsky1986 = list("Krinsky & Robb", "1986", paste(
      "Krinsky, I., & Robb, A. L. (1986). On approximating the statistical",
      "properties of elasticities. *The Review of Economics and Statistics*,",
      "68(4), 715\u2013719. <https://doi.org/10.2307/1924536>"
    )),
    liess2016 = list("Liess et al.", "2016", paste(
      "Liess, M., Foit, K., Knillmann, S., Sch\u00e4fer, R. B., & Liess, H.-D. (2016).",
      "Predicting the synergy of multiple stress effects. *Scientific Reports*,",
      "6, 32965. <https://doi.org/10.1038/srep32965>"
    )),
    loewe1953 = list("Loewe", "1953", paste(
      "Loewe, S. (1953). The problem of synergism and antagonism of combined",
      "drugs. *Arzneimittelforschung*, 3, 285\u2013290."
    )),
    madge2025 = list("Madge Pimentel et al.", "2025", paste(
      "Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S.",
      "J., Matthaei, C. D., Orr, J. A., Piggott, J. J., & Sch\u00e4fer, R. B. (2025).",
      "Hypothesis-driven research on multiple stressors: An analytical framework",
      "for stressor interactions. *Ecology and Evolution*, 15, e71959.",
      "<https://doi.org/10.1002/ece3.71959>"
    )),
    piggott2015 = list("Piggott et al.", "2015", paste(
      "Piggott, J. J., Townsend, C. R., & Matthaei, C. D. (2015).",
      "Reconceptualizing synergism and antagonism among multiple stressors.",
      "*Ecology and Evolution*, 5(7), 1538\u20131547.",
      "<https://doi.org/10.1002/ece3.1465>"
    )),
    satterthwaite1946 = list("Satterthwaite", "1946", paste(
      "Satterthwaite, F. E. (1946). An approximate distribution of estimates of",
      "variance components. *Biometrics Bulletin*, 2(6), 110\u2013114.",
      "<https://doi.org/10.2307/3002019>"
    )),
    schafer2018 = list("Sch\u00e4fer & Piggott", "2018", paste(
      "Sch\u00e4fer, R. B., & Piggott, J. J. (2018). Advancing understanding and",
      "prediction in multiple stressor research through a mechanistic basis for",
      "null models. *Global Change Biology*, 24(5), 1817\u20131826.",
      "<https://doi.org/10.1111/gcb.14073>"
    )),
    schuirmann1987 = list("Schuirmann", "1987", paste(
      "Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure",
      "and the power approach for assessing the equivalence of average",
      "bioavailability. *Journal of Pharmacokinetics and Biopharmaceutics*, 15(6),",
      "657\u2013680. <https://doi.org/10.1007/BF01068419>"
    ))
  )
}

# Pieces of text -------------------------------------------------------------------

# Values in the table of decisions: NULL as "none", vectors separated by commas.
show_value <- function(v) {
  if (is.list(v)) {
    v <- unlist(v)
  }
  if (is.null(v) || length(v) == 0L) {
    return("none")
  }
  if (is.numeric(v)) {
    return(paste(format(v, scientific = FALSE, trim = TRUE, drop0trailing = TRUE),
                 collapse = ", "))
  }
  paste(as.character(v), collapse = ", ")
}

# Default of an argument of a function, as text ("required" if it has none).
default_text <- function(fun, arg) {
  f <- formals(fun)
  # An argument without a default holds the empty symbol, which cannot be
  # assigned to a variable.
  if (is.symbol(f[[arg]]) && !nzchar(as.character(f[[arg]]))) {
    return("required")
  }
  show_value(eval(f[[arg]], envir = baseenv()))
}

# Readable names of the families of glmmTMB and of the candidates of
# fit_responses().
family_name <- function(family) {
  labels <- c(gaussian = "Gaussian", gaussian_log = "Gaussian", Gamma = "Gamma",
              gamma_log = "Gamma", lognormal = "lognormal", tweedie = "Tweedie",
              poisson = "Poisson", nbinom2 = "negative binomial (NB2)",
              nbinom1 = "negative binomial (NB1)", binomial = "binomial",
              betabinomial = "beta-binomial", beta = "beta", beta_family = "beta")
  out <- unname(labels[family])
  ifelse(is.na(out), family, out)
}

# "Gamma (log link)", "negative binomial (NB2, log link)".
family_label <- function(family, link) {
  name <- family_name(family)
  vapply(seq_along(name), function(i) {
    if (grepl(")", name[[i]], fixed = TRUE)) {
      sub(")$", sprintf(", %s link)", link[[i]]), name[[i]])
    } else {
      sprintf("%s (%s link)", name[[i]], link[[i]])
    }
  }, character(1))
}

# A family recorded as "poisson (log)" in the settings, as a label.
settings_family_label <- function(x) {
  family <- sub(" \\(.*$", "", x)
  link <- sub("^[^(]*\\(([^)]*)\\).*$", "\\1", x)
  family_label(family, link)
}

# A list of phrases that may themselves contain commas or "and".
phrase_list <- function(x) {
  if (length(x) > 1L && any(grepl(",| and ", x))) paste(x, collapse = "; ") else and_list(x)
}

# "Gamma (log link) for a and b; Tweedie (log link) for c".
grouped_labels <- function(labels, responses) {
  u <- unique(labels)
  phrase_list(vapply(u, function(l) {
    sprintf("%s for %s", l, and_list(responses[labels == l]))
  }, character(1), USE.NAMES = FALSE))
}

# Times in numeric order when they are numbers.
sort_times <- function(times) {
  times <- unique(as.character(times))
  num <- suppressWarnings(as.numeric(times))
  if (anyNA(num)) times else times[order(num)]
}

plural <- function(n, one, many) {
  if (n == 1L) one else many
}

# "The response was", "The responses a and b were".
responses_subject <- function(who, all) {
  n <- length(who)
  if (setequal(who, all)) {
    plural(n, "The response was", "The responses were")
  } else {
    sprintf("The %s %s %s", plural(n, "response", "responses"), and_list(who),
            plural(n, "was", "were"))
  }
}

range_text <- function(x) {
  sub("-", " to ", x, fixed = TRUE)
}

version_suffix <- function(versions, pkg) {
  if (!is.null(versions) && pkg %in% names(versions) && !is.na(versions[[pkg]])) {
    paste0(" ", versions[[pkg]])
  } else {
    ""
  }
}

# The design ------------------------------------------------------------------------

report_design <- function(x) {
  d <- x$design
  spec <- d$spec
  stressors <- spec$stressors
  stressor_text <- vapply(stressors, function(s) {
    sprintf("%s (%s, with %s as control)", s, and_list(spec$levels[[s]]), spec$control[[s]])
  }, character(1), USE.NAMES = FALSE)
  combos <- prod(lengths(spec$levels[stressors]))
  observed <- nrow(unique(d$cells[, stressors, drop = FALSE]))
  n_missing <- if (is.null(d$missing_cells)) 0L else nrow(d$missing_cells)
  layout <- if (n_missing == 0L && observed == combos) {
    sprintf("a full factorial design (%d treatment combinations%s)", combos,
            if (length(spec$by) > 0L) sprintf(" in each group of %s", and_list(spec$by)) else "")
  } else if (length(spec$by) == 0L) {
    sprintf("%d of the %d combinations of their levels", observed, combos)
  } else {
    sprintf(paste(
      "an incomplete factorial design (%d combinations of their levels, %d of them",
      "missing within the groups of %s)"
    ), combos, n_missing, and_list(spec$by))
  }
  per_cell <- range_text(format_range(d$cells$n_units))
  replication <- if (is.null(spec$unit)) {
    sprintf("%s replicates per combination", per_cell)
  } else {
    sprintf("%s experimental units (%s) per combination", per_cell, spec$unit)
  }
  text <- sprintf("The experiment combined %d stressors, %s, in %s, with %s",
                  length(stressors), and_list(stressor_text), layout, replication)
  baseline <- baseline_times(spec)
  if (!is.null(spec$time)) {
    times <- x$estimates$time
    tested <- times[!is.na(times) & times != "overall"]
    times <- sort_times(c(unlist(baseline), tested))
    if (length(times) > 0L) {
      text <- paste0(text, sprintf(", sampled at %d %s (%s %s)", length(times),
                                   plural(length(times), "time", "times"), spec$time,
                                   and_list(times)))
    }
  }
  text <- paste0(text, ".")
  if (length(spec$by) > 0L) {
    n_groups <- length(unique(d$sets$context))
    text <- paste(text, sprintf(
      "Interactions were estimated within %s of %s%s.",
      plural(n_groups, "the only group", sprintf("each of the %d groups", n_groups)),
      and_list(spec$by),
      if (length(spec$shared) > 0L) {
        sprintf(", with the rows shared by every group (%s) used in each of them",
                paste(sprintf("%s = %s", names(spec$shared), unlist(spec$shared)),
                      collapse = ", "))
      } else {
        ""
      }
    ))
  }
  before <- Filter(length, baseline)
  if (length(before) > 0L) {
    text <- paste(text, if (common_baseline(before, stressors)) {
      sprintf("The %s before the stressors were applied (%s %s) %s not tested.",
              plural(length(before[[1L]]), "sampling", "samplings"), spec$time,
              and_list(sort_times(before[[1L]])), plural(length(before[[1L]]), "was", "were"))
    } else {
      sprintf(
        paste(
          "Interactions that involve a stressor were not tested at the samplings",
          "before it was applied (%s)."
        ),
        paste(vapply(names(before), function(s) {
          sprintf("%s: %s %s", s, spec$time, and_list(sort_times(before[[s]])))
        }, character(1)), collapse = "; ")
      )
    })
  }
  text <- paste(c(text, group_sentences(d$groups)), collapse = " ")
  paste(text, set_sentences(d, spec))
}

# Whether the samplings before the stressors were applied are the same for
# every stressor.
common_baseline <- function(before, stressors) {
  length(before) == length(stressors) && length(unique(lapply(before, sort_times))) == 1L
}

group_sentences <- function(groups) {
  if (is.null(groups) || nrow(groups) == 0L) {
    return(character())
  }
  out <- character()
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, , drop = FALSE]
    between <- strsplit(g$between, ", ", fixed = TRUE)[[1L]]
    n <- length(between)
    grouped <- sprintf("The units were grouped by %s (%d groups)", g$group, g$n_groups)
    out <- c(out, switch(
      g$structure,
      blocks = paste0(grouped, ", within which the treatment combinations were compared."),
      whole_units = sprintf(
        paste(
          "%s, and %s %s applied to whole groups, with %s groups per %s, so the groups",
          "are the replicates for the effects that involve %s."
        ),
        grouped, and_list(between), plural(n, "was", "were"),
        range_text(g$groups_per_combination),
        plural(n, "level", "combination of their levels"), plural(n, "it", "them")
      ),
      confounded = sprintf(
        paste(
          "%s, and %s %s applied to whole groups, with a single group for some %s, so",
          "the effects that involve %s cannot be separated from differences among groups."
        ),
        grouped, and_list(between), plural(n, "was", "were"),
        plural(n, "levels", "combinations of their levels"), plural(n, "it", "them")
      ),
      NULL
    ))
  }
  out
}

set_sentences <- function(d, spec) {
  sets <- d$sets[d$sets$net_identifiable, , drop = FALSE]
  if (nrow(sets) == 0L) {
    return("The design identified no interaction.")
  }
  by_order <- table(sets$order)
  order_text <- vapply(names(by_order), function(o) {
    n <- by_order[[o]]
    what <- switch(o, "2" = "pairwise", "3" = "three-way", "4" = "four-way", paste0(o, "-way"))
    sprintf("%d %s %s", n, what, plural(n, "interaction", "interactions"))
  }, character(1), USE.NAMES = FALSE)
  where <- if (identical(spec$background, "all")) {
    ", in the absence and in the presence of the other stressors"
  } else if (length(spec$stressors) > 2L && any(sets$order < length(spec$stressors))) {
    ", with the stressors outside each set at their control"
  } else {
    ""
  }
  text <- if (nrow(sets) == 1L) {
    sprintf("The design identified one contrast set, a %s interaction%s.",
            sub("^1 ", "", sub(" interaction$", "", order_text)), where)
  } else {
    sprintf("The design identified %d contrast sets: %s%s.", nrow(sets), and_list(order_text),
            where)
  }
  no_pure <- sum(sets$order >= 3L & sets$pure_identifiable %in% FALSE)
  if (no_pure > 0L) {
    text <- paste(text, sprintf(
      paste(
        "For %d of the sets of three or more stressors the design lacks lower-order",
        "combinations, so only their net deviation could be estimated."
      ),
      no_pure
    ))
  }
  not_identified <- sum(!d$sets$net_identifiable)
  if (not_identified > 0L) {
    text <- paste(text, sprintf(
      "%d more %s could not be estimated because the design lacks combinations %s.",
      not_identified, plural(not_identified, "interaction", "interactions"),
      plural(not_identified, "it needs", "they need")
    ))
  }
  confounded <- sum(sets$confounded %in% TRUE)
  if (confounded > 0L) {
    text <- paste(text, sprintf(
      "%d %s confounded with differences among groups of units.",
      confounded, plural(confounded, "set was", "sets were")
    ))
  }
  text
}

# Software and models ------------------------------------------------------------------

report_models <- function(x, refs) {
  s <- x$settings
  software <- sprintf("The analysis was run with the R package must%s (%s) in R%s (%s).",
                      version_suffix(s$versions, "must"), refs$cite("must"),
                      version_suffix(s$versions, "R"), refs$cite("R"))
  if (!is.null(s$model)) {
    return(paste(software, report_external_model(x, refs)))
  }
  if (is.null(s$fit)) {
    # Interactions estimated by an earlier version of must, which kept only the
    # family of each model.
    fams <- s$families[!is.na(s$families)]
    return(paste(software, sprintf(
      paste(
        "Each response was modelled with a generalized linear model with one mean",
        "per treatment combination, fitted with glmmTMB (%s): %s."
      ),
      refs$cite("brooks2017"), grouped_labels(settings_family_label(fams), names(fams))
    )))
  }
  paste(software, report_fit_model(x, s$fit, refs))
}

report_fit_model <- function(x, fit, refs) {
  spec <- x$design$spec
  ok <- Filter(function(r) identical(r$status, "ok"), fit$responses)
  failed <- Filter(function(r) !identical(r$status, "ok"), fit$responses)
  failed_text <- if (length(failed) > 0L) {
    sprintf("No model could be fitted to %s (%s).",
            and_list(vapply(failed, `[[`, character(1), "response")),
            paste(vapply(failed, `[[`, character(1), "reason"), collapse = " "))
  } else {
    character()
  }
  if (length(ok) == 0L) {
    return(failed_text)
  }
  responses <- vapply(ok, `[[`, character(1), "response")
  by_text <- if (length(spec$by) > 0L) sprintf(" and group of %s", and_list(spec$by)) else ""
  means <- if (is.null(spec$time)) {
    sprintf("one mean per treatment combination%s (a cell-means model)", by_text)
  } else if (identical(fit$time_structure, "additive")) {
    sprintf(paste(
      "one mean per treatment combination%s (a cell-means model) plus an effect of %s",
      "common to all of them"
    ), by_text, spec$time)
  } else {
    sprintf("one mean per treatment combination%s and sampling time (a cell-means model)",
            by_text)
  }
  random <- vapply(ok, function(r) random_description(r$formula, spec, fit$random_requested),
                   character(1))
  has_random <- random != "no random effects"
  random_text <- if (!any(has_random)) {
    ""
  } else if (length(unique(random)) == 1L) {
    paste(" and", random[[1L]])
  } else {
    paste0(", with ", grouped_labels(random, responses))
  }
  dispersion <- switch(
    fit$dispersion,
    by = sprintf(", with one dispersion parameter per group of %s", and_list(spec$by)),
    cell = ", with one dispersion parameter per treatment combination",
    ""
  )
  text <- sprintf(
    "%s (%s) was modelled with a %s with %s%s%s, fitted by maximum likelihood with glmmTMB%s (%s).",
    if (length(ok) == 1L) "The response" else "Each response",
    and_list(responses),
    if (any(has_random)) "generalized linear mixed model" else "generalized linear model",
    means, random_text, dispersion, version_suffix(fit$versions, "glmmTMB"),
    refs$cite("brooks2017")
  )
  if (!is.null(spec$time) && identical(fit$time_structure, "additive")) {
    text <- paste(text, paste(
      "Treatment effects were therefore assumed constant over the sampling times on the",
      "scale of the link function."
    ))
  }
  dropped <- dropped_random(ok, fit, spec)
  if (length(dropped) > 0L) {
    text <- paste(text, sprintf(
      paste(
        "Random intercepts whose variance was estimated at zero were removed, with",
        "the same maximum likelihood: %s."
      ),
      phrase_list(dropped)
    ))
  }
  trials <- unlist(fit$trials)
  if (length(trials) > 0L) {
    text <- paste(text, sprintf(
      "%s modelled as successes out of a number of trials (%s).",
      responses_subject(names(trials), responses),
      paste(sprintf("%s out of %s", names(trials), trials), collapse = ", ")
    ))
  }
  offset <- unlist(fit$offset)
  if (length(offset) > 0L) {
    text <- paste(text, sprintf(
      "%s modelled per unit of sampling effort, with the logarithm of the effort as an offset (%s).",
      responses_subject(names(offset), responses),
      paste(sprintf("%s per %s", names(offset), offset), collapse = ", ")
    ))
  }
  text <- paste(text, report_selection(ok, fit, refs))
  boundary <- vapply(ok, function(r) as.integer(r$n_boundary), integer(1))
  if (any(boundary > 0L)) {
    text <- paste(text, sprintf(
      paste(
        "%s with only zeros%s%s %s fixed at the boundary of the link function; a",
        "contrast that compares such a combination with its expectation was tested",
        "from the probability of observing only zeros under that expectation, and",
        "contrasts that need such a combination to compute it were not estimated."
      ),
      if (length(ok) == 1L && boundary[[1L]] == 1L) {
        "The mean of the treatment combination"
      } else if (length(ok) == 1L) {
        sprintf("The means of the %d treatment combinations", boundary[[1L]])
      } else {
        "The means of the treatment combinations"
      },
      if (any(vapply(ok, function(r) identical(r$kind, "binomial"), logical(1)))) {
        " or only successes"
      } else {
        ""
      },
      if (length(ok) > 1L) {
        sprintf(" (%s)", paste(sprintf("%d for %s", boundary[boundary > 0L],
                                       responses[boundary > 0L]), collapse = ", "))
      } else {
        ""
      },
      plural(sum(boundary), "was", "were")
    ))
  }
  if (identical(x$settings$correction, "small_sample")) {
    no_reml <- responses[!vapply(ok, function(r) isTRUE(r$reml), logical(1))]
    text <- paste(text, sprintf(
      "For the inference, the selected %s refitted by restricted maximum likelihood (REML)%s.",
      plural(length(ok), "model was", "models were"),
      if (length(no_reml) > 0L) {
        sprintf("; the refit failed for %s, whose maximum-likelihood fit was used", and_list(no_reml))
      } else {
        ""
      }
    ))
  }
  paste(c(text, failed_text), collapse = " ")
}

# How the error distribution of each response was chosen.
report_selection <- function(ok, fit, refs) {
  responses <- vapply(ok, `[[`, character(1), "response")
  one <- length(ok) == 1L
  given <- vapply(ok, function(r) identical(r$type_source, "families"), logical(1))
  candidates <- vapply(ok, function(r) {
    cd <- r$candidates
    if (is.null(cd) || nrow(cd) == 0L) "" else and_list(family_label(cd$distribution, cd$link))
  }, character(1))
  key <- paste(given, candidates, sep = "\r")
  mixed <- any(given) && !all(given)
  groups <- vapply(unique(key), function(k) {
    i <- which(key == k)
    sprintf("%s%s%s", candidates[[i[[1L]]]],
            if (mixed && given[[i[[1L]]]]) ", given by the authors" else "",
            if (one) "" else sprintf(", for %s", and_list(responses[i])))
  }, character(1), USE.NAMES = FALSE)
  counts <- !given & vapply(ok, function(r) identical(r$response_type, "count"), logical(1))
  text <- sprintf(
    paste(
      "The error distribution was chosen among %s (%s): the candidate that converged",
      "with a positive-definite Hessian%s and had the lowest AIC was selected%s."
    ),
    if (all(given)) {
      "the candidates given by the authors"
    } else if (mixed) {
      "candidates suited to the values of each response or given by the authors"
    } else {
      sprintf("candidates suited to the values of %s", if (one) "the response" else "each response")
    },
    paste(groups, collapse = "; "),
    if (isTRUE(fit$residual_checks)) {
      sprintf(paste(
        ", passed three checks of its simulated residuals (uniformity, dispersion",
        "and zero inflation, p >= 0.05, with %d simulations; %s)"
      ), as.integer(fit$n_sim_residuals), refs$cite("DHARMa"))
    } else {
      ""
    },
    if (all(counts)) {
      ", except that the Poisson was used only if the negative binomial did not converge"
    } else if (any(counts)) {
      sprintf(paste(
        ", except that for the counts (%s) the Poisson was used only if the",
        "negative binomial did not converge"
      ), and_list(responses[counts]))
    } else {
      ""
    }
  )
  failed_checks <- vapply(ok, function(r) isFALSE(r$checks_passed), logical(1))
  if (isTRUE(fit$residual_checks) && any(failed_checks)) {
    text <- paste(text, sprintf(
      "No candidate passed every check for %s, and the one that passed the most was used.",
      and_list(responses[failed_checks])
    ))
  }
  selected <- vapply(ok, function(r) family_label(r$family_name, r$link), character(1))
  paste(text, if (one) {
    sprintf("The model selected was %s.", selected[[1L]])
  } else {
    sprintf("The models selected were %s.", grouped_labels(selected, responses))
  })
}

# Random terms of a formula: "(1 | .unit)" and the like.
random_terms <- function(formula) {
  if (is.null(formula) || is.na(formula)) {
    return(character())
  }
  regmatches(formula, gregexpr("\\([^()]*\\|[^()]*\\)", formula))[[1L]]
}

normalise_term <- function(x) {
  gsub("[` ]", "", x)
}

random_description <- function(formula, spec, requested = "automatic") {
  if (!is.null(requested) && !requested %in% c("automatic", "none")) {
    return(sprintf("the random part %s given by the authors", requested))
  }
  terms <- random_terms(formula)
  if (length(terms) == 0L) {
    return("no random effects")
  }
  and_list(vapply(terms, function(t) {
    group <- normalise_term(sub("\\)$", "", sub("^.*\\|", "", t)))
    who <- if (identical(group, ".unit")) {
      if (is.null(spec$unit)) "unit" else spec$unit
    } else {
      sub("^\\.g_", "", group)
    }
    sprintf("a random intercept for each %s", who)
  }, character(1), USE.NAMES = FALSE))
}

# Automatic random intercepts removed from the selected model of a response.
dropped_random <- function(ok, fit, spec) {
  if (!identical(fit$random_requested, "automatic") || identical(fit$random, "none")) {
    return(character())
  }
  wanted <- normalise_term(random_terms(fit$random))
  out <- vapply(ok, function(r) {
    lost <- setdiff(wanted, normalise_term(random_terms(r$formula)))
    if (length(lost) == 0L) {
      return("")
    }
    sprintf("%s for %s", random_description(paste(lost, collapse = " + "), spec), r$response)
  }, character(1))
  out[nzchar(out)]
}

report_external_model <- function(x, refs) {
  s <- x$settings
  label <- s$families[[1L]]
  pattern <- "^(\\S+) \\(([^)]*)\\)(; response transformed with (.*))?$"
  parts <- regmatches(label, regexec(pattern, label))[[1L]]
  family <- if (length(parts) > 0L) {
    transformed <- parts[[5L]]
    if (nzchar(transformed)) {
      transformed <- sprintf(", on a response transformed with %s", transformed)
    }
    sprintf("%s family and %s link%s", family_name(parts[[2L]]), parts[[3L]], transformed)
  } else {
    label
  }
  sprintf(
    paste(
      "The response (%s) was modelled outside must, with a model of class %s (%s),",
      "and the means of the treatment combinations were its estimated marginal",
      "means from emmeans%s (%s)%s."
    ),
    names(s$families)[[1L]], s$model, family, version_suffix(s$versions, "emmeans"),
    refs$cite("emmeans"),
    if (length(s$averaged_over) > 0L) sprintf(", averaged over %s", and_list(s$averaged_over)) else ""
  )
}

# Null models and inference ----------------------------------------------------------------

null_model_phrase <- function(code, s, refs) {
  switch(
    code,
    AD = "simple addition (AD)",
    MU = "the multiplicative model (MU)",
    DO = "dominance (DO)",
    BL = sprintf("Bliss independence (BL; %s)", refs$cite("bliss1939")),
    SAM = sprintf("the stress addition model (SAM; %s)", refs$cite("liess2016")),
    CA = sprintf("concentration addition (CA; %s)", refs$cite("loewe1953"))
  )
}

# What the null models that need more than the means assume.
null_model_details <- function(nm, s) {
  c(
    if ("SAM" %in% nm) {
      sprintf(paste(
        "Under SAM, the stress capacities of the individuals followed a symmetric beta",
        "distribution with both shape parameters equal to %s."
      ), format(s$sam_shape))
    },
    if ("CA" %in% nm) {
      sprintf(paste(
        "Under CA, two-parameter log-logistic curves were fitted to the responses at",
        "the doses of each stressor tested alone%s."
      ), if (identical(s$inference, "simulation")) ", and refitted in each simulation draw" else "")
    }
  )
}

report_null_models <- function(x, justification, refs) {
  s <- x$settings
  nm <- s$null_models
  text <- sprintf("Each interaction was tested against %s.", and_list(vapply(
    nm, null_model_phrase, character(1), s = s, refs = refs, USE.NAMES = FALSE
  )))
  why <- c(
    if (!is.null(s$null_model_choice)) null_model_reason(s$null_model_choice, refs),
    justification
  )
  if (length(why) == 0L && identical(nm, "AD")) {
    why <- paste(
      "No hypothesis about how the stressors act was stated, so the default null",
      "model of must was used: in simulated experiments, among the null models",
      "defined in two-level factorials, it flagged the fewest false interactions",
      "on average over five mechanisms of no interaction."
    )
  }
  # The sentences on the null models, and those on why they were chosen.
  c(paste(c(text, null_model_details(nm, s)), collapse = " "), if (length(why) > 0L) {
    paste(why, collapse = " ")
  })
}

# How the interactions were estimated and tested.
report_estimation <- function(x, refs) {
  s <- x$settings
  nm <- s$null_models
  scales <- c(
    if (any(nm %in% c("AD", "DO"))) {
      sprintf("as a difference on the scale of the response (%s)",
              and_list(intersect(c("AD", "DO"), nm)))
    },
    if ("MU" %in% nm) "as a difference on the log scale (MU)",
    if (any(nm %in% relative_null_models())) {
      sprintf("on the responses relative to the reference combination (%s)",
              and_list(intersect(relative_null_models(), nm)))
    }
  )
  n_models <- length(unique(x$estimates$response))
  means <- if (!is.null(s$model)) {
    "the estimated marginal means of the model"
  } else if (fit_has_random(s$fit)) {
    sprintf("the means predicted by the %s with the random effects at zero",
            plural(n_models, "model", "models"))
  } else {
    sprintf("the means predicted by the %s", plural(n_models, "model", "models"))
  }
  text <- sprintf(
    paste(
      "Following %s, each interaction was computed from %s (post-estimation",
      "inference), as the deviation of the combination of stressors from the",
      "expectation of the null model, %s."
    ),
    refs$citet("madge2025"), means, and_list(scales)
  )
  sets <- x$design$sets
  if (any(sets$order >= 3L & sets$net_identifiable)) {
    pure <- intersect(nm, c("AD", "MU", "BL"))
    text <- paste(text, if (length(pure) > 0L) {
      sprintf(
        paste(
          "For sets of three or more stressors, both the net deviation from the single",
          "effects and the pure highest-order interaction (under %s%s) were estimated."
        ),
        and_list(pure),
        if ("BL" %in% pure) {
          sprintf("; under BL, the emergent interaction of %s, for three stressors",
                  refs$citet("beppler2016"))
        } else {
          ""
        }
      )
    } else {
      "For sets of three or more stressors, the net deviation from the single effects was estimated."
    })
    if ("DO" %in% nm) {
      text <- paste(text, "Dominance is defined for two stressors only.")
    }
  }
  paste(text, report_inference(x, refs))
}

# Whether some model of fit_responses() has random effects.
fit_has_random <- function(fit) {
  if (is.null(fit)) {
    return(FALSE)
  }
  any(vapply(fit$responses, function(r) {
    identical(r$status, "ok") && length(random_terms(r$formula)) > 0L
  }, logical(1)))
}

# Why these null models: the criteria of suggest_null_models().
null_model_reason <- function(choice, refs) {
  cr <- choice$criteria
  value <- stats::setNames(cr$value, cr$criterion)
  stated <- list(
    similar_action = c(yes = "the stressors act similarly", no = "the stressors do not act similarly"),
    specific_action = c(yes = "one of them acts specifically on physiology",
                        no = "none of them acts specifically on physiology"),
    mortality = c(yes = "the effect is mortality", no = "the effect is not mortality"),
    stressor_effect = c(yes = "the stressor-effect relationships are known",
                        no = "the stressor-effect relationships are not known"),
    sensitivity_distribution = c(yes = "the distribution of sensitivities is known",
                                 no = "the distribution of sensitivities is not known"),
    correlation = c(positive = "the sensitivities are highly positively correlated",
                    negative = "the sensitivities are highly negatively correlated",
                    weak = "the sensitivities are weakly correlated")
  )
  open_question <- c(similar_action = "whether the stressors act similarly",
                     specific_action = "whether one of them acts specifically on physiology",
                     mortality = "whether the effect is mortality",
                     stressor_effect = "whether the stressor-effect relationships are known",
                     sensitivity_distribution = "whether the distribution of sensitivities is known",
                     correlation = "how the sensitivities are correlated")
  user <- names(open_question)
  given <- user[value[user] != "unknown"]
  not_known <- user[value[user] == "unknown"]
  text <- sprintf(
    paste(
      "The null models were chosen with the guidelines of %s, in which must adds two",
      "conditions: concentration addition requires that the stressors alone act in",
      "the same direction, and stress addition that none of them increases the",
      "response."
    ),
    refs$citet("schafer2018")
  )
  if (length(given) > 0L) {
    text <- paste(text, sprintf("The authors stated that %s.", and_list(vapply(given, function(n) {
      unname(stated[[n]][[value[[n]]]])
    }, character(1), USE.NAMES = FALSE))))
  }
  if (length(not_known) > 0L) {
    text <- paste(text, sprintf("They did not state %s.", or_list(unname(open_question[not_known]))))
  }
  criteria <- sprintf(paste(
    "whether one effect reaches %s%% of the effect limit, whether their sum exceeds",
    "the limit, whether a stressor alone increases the response and whether two act",
    "in opposite directions"
  ), format(100 * choice$near_limit))
  text <- paste(text, switch(
    choice$effects,
    estimated = sprintf(
      paste(
        "The criteria that depend on the effects of the stressors alone (%s) were",
        "evaluated from the effects estimated by the models in each contrast set%s, %s."
      ),
      criteria,
      if (any(grepl("times", cr$detail, fixed = TRUE))) " and sampling time" else "",
      if (is.null(choice$level)) {
        "from their point estimates"
      } else {
        sprintf(paste(
          "with %s%% intervals; a criterion was not known when the interval spanned",
          "its threshold"
        ), format(100 * choice$level))
      }
    ),
    given = sprintf(paste(
      "The criteria that depend on the effects of the stressors alone (%s) were",
      "evaluated from effects given by the authors."
    ), criteria),
    sprintf(paste(
      "The effects of the stressors alone were not used, so the criteria that depend",
      "on them (%s) were not known."
    ), criteria)
  ))
  if (any(cr$value %in% c("unknown", "uncertain", "varies")) && length(choice$null_models) > 1L) {
    text <- paste(text, sprintf(
      paste(
        "Each criterion that was not known%s took each of its values, and every null",
        "model reached was kept."
      ),
      if (any(cr$value == "varies")) " (or that differed among contrast sets)" else ""
    ))
  }
  if (length(choice$complement) > 0L) {
    text <- paste(text, sprintf(
      paste(
        "For %s, which the stressors alone increased, the guidelines were applied to",
        "the complement of the proportion (1 - p)."
      ),
      and_list(choice$complement)
    ))
  }
  text
}

or_list <- function(x) {
  n <- length(x)
  if (n <= 1L) {
    return(paste(x, collapse = ""))
  }
  paste(paste(x[-n], collapse = ", "), "or", x[[n]])
}

report_inference <- function(x, refs) {
  s <- x$settings
  est <- x$estimates
  tested <- est[est$component == "interaction" & est$status %in% c("ok", "boundary"), ,
                drop = FALSE]
  finite <- round(tested$df[is.finite(tested$df)], 1)
  df_text <- if (length(finite) == 0L) {
    ""
  } else if (min(finite) == max(finite)) {
    format(min(finite))
  } else {
    sprintf("%s to %s", format(min(finite)), format(max(finite)))
  }
  text <- if (identical(s$inference, "simulation")) {
    sprintf(paste(
      "Uncertainty was obtained from %s draws of the fixed effects from their",
      "estimated distribution (%s), which gave percentile intervals and two-sided",
      "p-values."
    ), format(s$n_sim, big.mark = ",", scientific = FALSE), refs$cite("krinsky1986"))
  } else {
    "Uncertainty was obtained with the delta method (Wald intervals and tests)."
  }
  normal <- unique(tested$response[is.infinite(tested$df)])
  correction <- if (!identical(s$correction, "small_sample")) {
    sprintf(paste(
      "No correction for small samples was applied: the covariance was that of the",
      "maximum-likelihood fit and the distributions were normal (as in %s)."
    ), refs$cite("madge2025"))
  } else if (!is.null(s$model)) {
    if (nzchar(df_text)) {
      sprintf("To correct for small samples, the distributions were t, with %s (%s).",
              if (is.numeric(s$df)) {
                "the degrees of freedom given by the authors"
              } else if (identical(s$df, "residual")) {
                "the residual degrees of freedom of the model"
              } else {
                "the degrees of freedom given by emmeans"
              }, df_text)
    } else {
      "The distributions were normal, because emmeans gave no finite degrees of freedom."
    }
  } else if (!nzchar(df_text)) {
    sprintf(paste(
      "To correct for small samples, the covariance came from the REML refit; the",
      "distributions stayed normal, because the %s no dispersion parameter or",
      "random effects."
    ), plural(length(unique(est$response)), "model has", "models have"))
  } else {
    whole <- !is.null(x$design$groups) && any(x$design$groups$structure == "whole_units")
    sprintf(paste(
      "To correct for small samples, the covariance came from the REML refit and the",
      "distributions were t, with the degrees of freedom given by the replication of",
      "the experimental units%s (%s)%s."
    ), if (whole) " and of the groups to which stressors were applied whole" else "",
    df_text, if (length(normal) > 0L) {
      sprintf("; for %s, without a dispersion parameter or random effects, they stayed normal",
              and_list(normal))
    } else {
      ""
    })
  }
  text <- paste(text, correction)
  if (identical(s$fit$dispersion, "cell") && identical(s$correction, "small_sample")) {
    text <- paste(text, sprintf(
      paste(
        "With one dispersion parameter per combination, the degrees of freedom of each",
        "contrast followed the Welch-Satterthwaite approximation (%s)."
      ),
      refs$cite("satterthwaite1946")
    ))
  }
  has_time <- any(!is.na(est$time) & est$time != "overall")
  adjust <- if (identical(s$p_adjust, "none")) {
    "p-values were not adjusted for multiple testing"
  } else {
    sprintf(
      "p-values were adjusted %s within each response, null model and estimate, across contrast sets%s",
      p_adjust_text(s$p_adjust, refs), if (has_time) " and sampling times" else ""
    )
  }
  text <- paste(text, sprintf(
    "Intervals were computed at the %s%% level, and %s%s.",
    format(100 * s$level), adjust,
    if (isTRUE(s$overall) && any(est$time %in% "overall")) {
      "; effects averaged over the sampling times were tested and adjusted separately"
    } else {
      ""
    }
  ))
  sensitivity <- as.character(s$sensitivity)
  if (length(sensitivity) > 0L) {
    complement <- intersect(sensitivity, est$response[grepl("complement", est$note, fixed = TRUE)])
    other <- setdiff(sensitivity, complement)
    if (length(complement) > 0L) {
      text <- paste(text, sprintf(
        paste(
          "For %s, a proportion in which an increase means a stronger effect, the null",
          "models were applied to its complement (1 - p)."
        ),
        and_list(complement)
      ))
    }
    if (length(other) > 0L) {
      text <- paste(text, sprintf(
        paste(
          "For %s, in which an increase means a stronger effect, only simple addition",
          "was tested, because the other null models assume that an increase means",
          "tolerance (%s)."
        ),
        and_list(other), refs$cite("madge2025")
      ))
    }
  }
  if (!is.null(s$margin)) {
    text <- paste(text, sprintf(
      paste(
        "Each interaction was also tested for equivalence with two one-sided tests",
        "(%s) of whether its deviation relative to the response of the reference",
        "combination, (observed - expected) / reference, lay between -%s and %s; each",
        "test measured the signed distance of the means from the closest means at the",
        "bound, on a scale that stabilizes their variance (%s)."
      ),
      refs$cite("schuirmann1987"), format(s$margin), format(s$margin),
      refs$cite("anscombe1948")
    ))
  }
  if (identical(s$inference, "simulation") || any(est$status %in% "boundary")) {
    text <- paste(text, sprintf("Random numbers were generated with seed %s.", format(s$seed)))
  }
  text
}

p_adjust_text <- function(method, refs) {
  switch(
    method,
    BH = , fdr = sprintf("for the false discovery rate (%s)", refs$cite("benjamini1995")),
    BY = sprintf("for the false discovery rate under dependence (%s)", refs$cite("benjamini2001")),
    holm = sprintf("with the method of %s", refs$citet("holm1979")),
    hochberg = sprintf("with the method of %s", refs$citet("hochberg1988")),
    hommel = sprintf("with the method of %s", refs$citet("hommel1988")),
    bonferroni = "with the Bonferroni correction",
    sprintf("with the method \"%s\" of p.adjust()", method)
  )
}

# Classification and agreement ------------------------------------------------------

report_classes <- function(x, classes, agreement, refs) {
  parts <- character()
  if (!is.null(classes)) {
    alpha <- attr(classes, "alpha")
    text <- sprintf(
      paste(
        "Interactions were classified with the directional scheme of %s: each single,",
        "combined and interaction effect was coded as positive, negative or not",
        "significant from its adjusted p-value%s, and each net interaction was classed",
        "as additive (no deviation detected), as a positive or negative synergism or as",
        "a positive or negative antagonism (their Table 1)."
      ),
      refs$citet("piggott2015"),
      if (is.null(alpha)) "" else sprintf(" (alpha = %s)", format(alpha))
    )
    if ("DO" %in% x$settings$null_models) {
      text <- paste(text, paste(
        "Dominance, which does not predict the combination from the single effects,",
        "was not classified."
      ))
    }
    if (any(x$design$sets$order >= 3L & x$design$sets$net_identifiable)) {
      text <- paste(text, paste(
        "For three or more stressors, the same rules were applied to the net deviation",
        "and the single effects of all the stressors of the set; pure interactions were",
        "not classified."
      ))
    }
    if (!is.null(x$settings$margin)) {
      text <- paste(text, paste(
        "With the equivalence margin, each net interaction was also a deviation, a",
        "negligible deviation (detected but within the margin), equivalent (no",
        "deviation detected and within the margin) or inconclusive."
      ))
    }
    parts <- c(parts, text)
  }
  if (!is.null(agreement)) {
    alpha <- attr(agreement, "alpha")
    compared <- compared_null_models(agreement)
    parts <- c(parts, sprintf(
      paste(
        "Across the null models compared (%s), an interaction was robust when every",
        "null model tested flagged a deviation with the same sign (adjusted p <",
        "%s), dependent on the null model when only some did, and conflicting when",
        "they flagged deviations of opposite sign."
      ),
      and_list(compared), if (is.null(alpha)) "alpha" else format(alpha)
    ))
  }
  paste(parts, collapse = " ")
}

# The table of decisions ----------------------------------------------------------------

# One row of the table of decisions; without `origin`, "default" when the
# value is the default and "set" otherwise.
decision_row <- function(step, decision, argument, value, default, origin = NULL) {
  if (is.null(origin)) {
    origin <- if (identical(value, default)) "default" else "set"
  }
  data.frame(step = step, decision = decision, argument = argument, value = value,
             default = default, source = origin, stringsAsFactors = FALSE)
}

report_decisions <- function(x, classes, agreement, justification) {
  rows <- c(decisions_design(x), decisions_fit(x), decisions_suggestion(x),
            decisions_estimate(x), decisions_report(justification),
            decisions_classes(x, classes, agreement), decisions_software(x))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

decisions_design <- function(x) {
  spec <- x$design$spec
  step <- "describe_design()"
  row <- function(decision, arg, value) {
    decision_row(step, decision, arg, value, default_text(describe_design, arg))
  }
  k <- length(spec$stressors)
  rows <- list(
    decision_row(step, "Stressors and their control levels", "stressors, control",
                 paste(sprintf("%s (control %s)", spec$stressors,
                               unlist(spec$control[spec$stressors])), collapse = ", "),
                 "required", "set"),
    row("Experimental unit", "unit", show_value(spec$unit)),
    row("Sampling time", "time", show_value(spec$time))
  )
  if (!is.null(spec$time)) {
    before <- Filter(length, baseline_times(spec))
    rows <- c(rows, list(row(
      "Samplings before the stressors were applied", "baseline",
      if (length(before) == 0L) {
        "none"
      } else if (common_baseline(before, spec$stressors)) {
        show_value(sort_times(before[[1L]]))
      } else {
        paste(sprintf("%s: %s", names(before),
                      vapply(before, function(b) show_value(sort_times(b)), character(1))),
              collapse = "; ")
      }
    )))
  }
  rows <- c(rows, list(
    row("Groups of units (blocks, rooms, whole plots)", "groups", show_value(spec$groups)),
    row("Factors within which interactions are estimated", "by", show_value(spec$by))
  ))
  if (length(spec$by) > 0L) {
    rows <- c(rows, list(row(
      "Rows shared by every group", "shared",
      if (length(spec$shared) == 0L) {
        "none"
      } else {
        paste(sprintf("%s = %s", names(spec$shared), unlist(spec$shared)), collapse = ", ")
      }
    )))
  }
  c(rows, list(
    row("Level of the stressors outside each set", "background", spec$background),
    decision_row(step, "Highest number of stressors in a set", "max_order",
                 if (spec$max_order == k) sprintf("%d (all stressors)", k) else format(spec$max_order),
                 "all stressors", if (spec$max_order == k) "default" else "set")
  ))
}

decisions_fit <- function(x) {
  fit <- x$settings$fit
  if (is.null(fit)) {
    return(list())
  }
  step <- "fit_responses()"
  row <- function(decision, arg, value, origin = NULL) {
    decision_row(step, decision, arg, value, default_text(fit_responses, arg), origin)
  }
  resp <- fit$responses
  ok <- Filter(function(r) identical(r$status, "ok"), resp)
  rows <- list()
  if (length(ok) > 0L) {
    given <- vapply(ok, function(r) identical(r$type_source, "families"), logical(1))
    candidates <- vapply(ok, function(r) {
      cd <- r$candidates
      text <- if (is.null(cd) || nrow(cd) == 0L) {
        "none"
      } else {
        and_list(family_label(cd$distribution, cd$link))
      }
      if (r$type_source %in% c("detected", "declared")) {
        text <- sprintf("%s (%s, %s)", text, r$response_type, r$type_source)
      }
      text
    }, character(1))
    names_ok <- vapply(ok, `[[`, character(1), "response")
    key <- paste(given, candidates, sep = "\r")
    for (k in unique(key)) {
      i <- which(key == k)
      rows <- c(rows, list(decision_row(
        step, sprintf("Candidate distributions (%s)", and_list(names_ok[i])), "families",
        candidates[[i[[1L]]]], "chosen from the values", if (given[[i[[1L]]]]) "set" else "default"
      )))
    }
    declared <- vapply(ok, function(r) identical(r$type_source, "declared"), logical(1))
    if (any(declared)) {
      rows <- c(rows, list(decision_row(
        step, "Response types", "response_types",
        paste(sprintf("%s: %s", names_ok[declared],
                      vapply(ok[declared], `[[`, character(1), "response_type")), collapse = ", "),
        "detected from the values", "set"
      )))
    }
  }
  if (length(unlist(fit$trials)) > 0L) {
    rows <- c(rows, list(row("Responses out of a number of trials", "trials",
                             paste(sprintf("%s: %s", names(fit$trials), unlist(fit$trials)),
                                   collapse = ", "))))
  }
  if (length(unlist(fit$offset)) > 0L) {
    rows <- c(rows, list(row("Sampling effort (offset)", "offset",
                             paste(sprintf("%s: %s", names(fit$offset), unlist(fit$offset)),
                                   collapse = ", "))))
  }
  if (!is.null(x$design$spec$time)) {
    rows <- c(rows, list(row("Treatment effects over time", "time_structure",
                             fit$time_structure)))
  }
  automatic <- identical(fit$random_requested, "automatic")
  rows <- c(rows, list(
    row("Observations that share a dispersion parameter", "dispersion", fit$dispersion),
    decision_row(step, "Random effects", "random",
                 if (automatic) sprintf("automatic: %s", fit$random) else fit$random_requested,
                 "automatic", if (automatic) "default" else "set"),
    decision_row(step, "Limits of the optimiser", "control", fit$control,
                 "scaled iteration limits"),
    row("Residual checks with DHARMa", "residual_checks", show_value(fit$residual_checks))
  ))
  if (isTRUE(fit$residual_checks)) {
    rows <- c(rows, list(
      row("Simulations for the residual checks", "n_sim_residuals",
          show_value(fit$n_sim_residuals)),
      decision_row(step, "Seed of the residual checks", "seed", show_value(fit$seed),
                   "drawn at random", seed_source(fit$seed_given))
    ))
  }
  for (r in resp) {
    rows <- c(rows, list(decision_row(
      step, sprintf("Distribution selected (%s)", r$response), "",
      if (identical(r$status, "ok")) {
        paste0(family_label(r$family_name, r$link),
               if (isFALSE(r$checks_passed)) " (no candidate passed every residual check)" else "")
      } else {
        "none: the model could not be fitted"
      },
      "", "chosen"
    )))
  }
  rows
}

seed_source <- function(given) {
  if (is.null(given)) "recorded" else if (isTRUE(given)) "set" else "chosen"
}

decisions_suggestion <- function(x) {
  choice <- x$settings$null_model_choice
  if (is.null(choice)) {
    return(list())
  }
  step <- "suggest_null_models()"
  cr <- choice$criteria
  user <- c("similar_action", "specific_action", "mortality", "stressor_effect",
            "sensitivity_distribution", "correlation")
  rows <- lapply(which(cr$criterion %in% user), function(i) {
    decision_row(step, cr$question[[i]], cr$criterion[[i]], cr$value[[i]], "NA (unknown)",
                 if (identical(cr$source[[i]], "given")) "set" else "default")
  })
  rows <- c(rows, list(decision_row(
    step, "Effects of the stressors alone", "effects",
    switch(choice$effects, estimated = "estimated from the models", given = "given",
           "not known"),
    "estimated from the models, or not known without them",
    if (identical(choice$effects, "given")) "set" else "default"
  )))
  if (identical(choice$effects, "estimated")) {
    rows <- c(rows, list(decision_row(
      step, "Level of the intervals of the effects", "level",
      if (is.null(choice$level)) "none (point estimates)" else show_value(choice$level),
      default_text(suggest_null_models, "level")
    )))
  }
  rows <- c(rows, list(decision_row(
    step, "Share of the effect limit that counts as close to it", "near_limit",
    show_value(choice$near_limit), default_text(suggest_null_models, "near_limit")
  )))
  effect_source <- c(estimated = "estimated", given = "set", unknown = "default")
  rows <- c(rows, lapply(which(!cr$criterion %in% user), function(i) {
    decision_row(step, cr$question[[i]], "",
                 if (nzchar(cr$detail[[i]])) {
                   sprintf("%s (%s)", cr$value[[i]], cr$detail[[i]])
                 } else {
                   cr$value[[i]]
                 },
                 "", unname(effect_source[cr$source[[i]]]))
  }))
  if (length(choice$complement) > 0L) {
    rows <- c(rows, list(decision_row(
      step, "Guidelines applied to the complement (1 - p)", "", and_list(choice$complement),
      "", "chosen"
    )))
  }
  rows
}

decisions_estimate <- function(x) {
  s <- x$settings
  step <- "estimate_interactions()"
  external <- !is.null(s$model)
  fun <- if (external) estimate_interactions.default else estimate_interactions.must_fit
  row <- function(decision, arg, value, origin = NULL) {
    decision_row(step, decision, arg, value, default_text(fun, arg), origin)
  }
  rows <- list()
  if (external) {
    rows <- c(rows, list(decision_row(step, "Model fitted outside must", "fit",
                                      sprintf("%s: %s", s$model, s$families[[1L]]),
                                      "required", "set")))
  }
  rows <- c(rows, list(
    row("Null models", "null_models", show_value(s$null_models),
        if (!is.null(s$null_model_choice)) "suggest_null_models()"),
    row("Inference", "inference", s$inference)
  ))
  if (identical(s$inference, "simulation") || any(x$estimates$status %in% "boundary")) {
    rows <- c(rows, list(row("Number of simulation draws", "n_sim", show_value(s$n_sim))))
  }
  rows <- c(rows, list(row("Small-sample correction", "correction", s$correction)))
  if (external && !identical(s$correction, "none")) {
    rows <- c(rows, list(decision_row(
      step, "Degrees of freedom", "df",
      if (is.numeric(s$df)) {
        show_value(s$df)
      } else if (identical(s$df, "residual")) {
        "residual degrees of freedom of the model"
      } else {
        "given by emmeans"
      },
      "given by emmeans", if (is.numeric(s$df)) "set" else "default"
    )))
  }
  rows <- c(rows, list(
    row("Confidence level", "level", show_value(s$level)),
    row("Adjustment of p-values", "p_adjust", s$p_adjust)
  ))
  if (!is.null(x$design$spec$time)) {
    rows <- c(rows, list(row("Averages over the sampling times", "overall",
                             show_value(s$overall))))
  }
  sensitivity <- as.character(s$sensitivity)
  rows <- c(rows, list(decision_row(
    step, "Responses in which an increase means sensitivity", "sensitivity",
    show_value(sensitivity), "none", if (length(sensitivity) == 0L) "default" else "set"
  )))
  if ("SAM" %in% s$null_models) {
    rows <- c(rows, list(row("Shape of the distribution of stress capacities (SAM)",
                             "sam_shape", show_value(s$sam_shape))))
  }
  rows <- c(rows, list(row("Equivalence margin", "margin", show_value(s$margin))))
  if (identical(s$inference, "simulation") || any(x$estimates$status %in% "boundary")) {
    rows <- c(rows, list(decision_row(step, "Seed of the simulations", "seed",
                                      show_value(s$seed), "drawn at random",
                                      seed_source(s$seed_given))))
  }
  rows
}

decisions_report <- function(justification) {
  list(decision_row("report_interactions()", "Justification of the null models",
                    "justification", if (is.null(justification)) "none" else "given",
                    "none"))
}

decisions_classes <- function(x, classes, agreement) {
  rows <- list()
  if (!is.null(classes)) {
    alpha <- attr(classes, "alpha")
    rows <- c(rows, list(decision_row(
      "classify_interactions()", "Significance level of the classes", "alpha",
      if (is.null(alpha)) "not recorded" else show_value(alpha),
      default_text(classify_interactions, "alpha"), if (is.null(alpha)) "recorded" else NULL
    )))
  }
  if (!is.null(agreement)) {
    alpha <- attr(agreement, "alpha")
    compared <- compared_null_models(agreement)
    rows <- c(rows, list(
      decision_row("compare_null_models()", "Null models compared", "null_models",
                   show_value(compared), "all those estimated",
                   if (setequal(compared, x$settings$null_models)) "default" else "set"),
      decision_row("compare_null_models()", "Significance level of the agreement", "alpha",
                   if (is.null(alpha)) "not recorded" else show_value(alpha),
                   default_text(compare_null_models, "alpha"),
                   if (is.null(alpha)) "recorded" else NULL)
    ))
  }
  rows
}

decisions_software <- function(x) {
  v <- c(x$settings$versions, x$settings$fit$versions)
  v <- v[!duplicated(names(v)) & !is.na(v)]
  lapply(names(v), function(pkg) {
    decision_row("software", sprintf("Version of %s", pkg), "", v[[pkg]], "", "recorded")
  })
}
