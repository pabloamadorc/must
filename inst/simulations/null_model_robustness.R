# Robustness of the null models of must across mechanisms of no interaction
#
# The default null model of estimate_interactions() was chosen with this
# simulation: the null model that flags the fewest false interactions when the
# data are generated without interaction by different mechanisms (simple
# addition, multiplication or Bliss independence, dominance, stress addition
# and concentration addition), without losing the power to detect
# interactions that no mechanism explains.
#
# What it does:
#   1. Simulates experiments in a 2 x 2 factorial and in a 3 x 3 dose
#      factorial, four units per combination, with survival (alive out of 20
#      animals per unit, binomial) or overdispersed counts (Poisson with
#      log-normal noise), and moderate or weak single effects.
#   2. The true response of each combination follows one of five mechanisms
#      without interaction, or an interaction that no mechanism explains:
#      synergy (half the response that the most harmful mechanism predicts)
#      or antagonism (halfway between the most affected single stressor and
#      the control). The predictions of the mechanisms are computed here,
#      independently of the package (uniroot() for concentration addition).
#   3. Fits a cell-means glm (binomial or quasi-Poisson) and estimates the
#      interactions with estimate_interactions(model, design) under AD, MU,
#      DO, BL, SAM and CA (CA in the dose factorial only).
#   4. Summarises, for each null model, the rate of false interactions
#      (p < 0.05 without adjustment, in each combination) under each
#      mechanism and the power against the interactions that no mechanism
#      explains, and the same for a consensus criterion (an interaction only
#      if every null model flags it with the same sign).
#
# Rule fixed before the results were seen: among the null models defined in
# two-level factorials and with an average power of at least 80% of the
# highest, the one with the lowest average rate of false interactions (mean
# over the five mechanisms, each with the same weight); if two are within one
# percentage point, the more powerful one.
#
# Usage:
#   source(system.file("simulations", "null_model_robustness.R", package = "must"))
# or, from a shell, Rscript null_model_robustness.R [replicates] [cores].
# Settings can be given before sourcing the script, without editing it:
#   REPLICATES <- 20; CORES <- 4; OUTPUT_DIR <- "robustness"; source(...)
# Run from the root of the package source, it uses the working version of
# must (through pkgload); otherwise, the installed one. emmeans is needed.
#
# Results go to OUTPUT_DIR (by default, a folder in the temporary directory
# of the session): null_model_robustness.csv, with the rates of every
# scenario and null model, and null_model_robustness.md, the report.
# source(...)$value also returns them. The files installed next to this
# script come from 100 replicates per scenario. Each replicate has its own
# seed, so the results do not depend on CORES (forked processes, on
# Unix-alikes only). SAVE_RAW keeps the result of every test in an .rds file;
# FROM_RAW rebuilds the summary and the report from such a file.

local({

  ## Settings -----------------------------------------------------------------

  option <- function(name, value) {
    if (exists(name, envir = globalenv(), inherits = FALSE)) {
      get(name, envir = globalenv())
    } else {
      value
    }
  }
  arguments <- commandArgs(trailingOnly = TRUE)
  REPLICATES <- option("REPLICATES",
                       if (length(arguments) >= 1L) as.integer(arguments[[1L]]) else 100L)
  CORES <- option("CORES", if (length(arguments) >= 2L) as.integer(arguments[[2L]]) else 1L)
  N_SIM <- option("N_SIM", 1000L)
  ALPHA <- 0.05
  SEED <- option("SEED", 2026L)
  OUTPUT_DIR <- option("OUTPUT_DIR", file.path(tempdir(), "null_model_robustness"))
  SAVE_RAW <- option("SAVE_RAW", NULL)
  FROM_RAW <- option("FROM_RAW", NULL)

  in_source <- file.exists("DESCRIPTION") &&
    identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]), "must")
  if (in_source && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(must)
  }
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("The simulation needs emmeans: install.packages(\"emmeans\").", call. = FALSE)
  }

  ## Dose-response curves and mechanisms (independent of the package) --------

  # Log-logistic curve relative to the control through w(1) = w1 and
  # w(2) = w2: 1 / w - 1 = (x / ec50)^b.
  loglogistic <- function(w1, w2) {
    b <- log((1 / w2 - 1) / (1 / w1 - 1)) / log(2)
    ec50 <- 1 / (1 / w1 - 1)^(1 / b)
    list(b = b, ec50 = ec50, w = function(x) 1 / (1 + (x / ec50)^b))
  }
  EFFECTS <- list(
    moderate = list(a = loglogistic(0.80, 0.50), b = loglogistic(0.85, 0.60)),
    weak = list(a = loglogistic(0.92, 0.75), b = loglogistic(0.95, 0.80))
  )

  SAM_SHAPE <- 3.2
  predict_sam <- function(w) {
    stress <- sum(stats::qbeta(pmin(pmax(1 - w, 0), 1), SAM_SHAPE, SAM_SHAPE))
    1 - stats::pbeta(min(stress, 1), SAM_SHAPE, SAM_SHAPE)
  }
  predict_ca <- function(x, curves) {
    active <- x > 0
    if (sum(active) == 1L) {
      return(curves[[which(active)]]$w(x[active]))
    }
    ec <- function(cv, w) cv$ec50 * ((1 - w) / w)^(1 / cv$b)
    f <- function(w) sum(vapply(which(active), function(i) x[[i]] / ec(curves[[i]], w), 0)) - 1
    stats::uniroot(f, c(1e-12, 1 - 1e-12), tol = 1e-12)$root
  }
  # Relative response of the combination with doses x under each mechanism.
  prediction <- function(mechanism, x, curves) {
    w <- c(curves$a$w(x[[1L]]), curves$b$w(x[[2L]]))
    switch(
      mechanism,
      AD = max(1 + sum(w - 1), 0.02),
      MU = prod(w),
      DO = min(w),
      SAM = predict_sam(w),
      CA = predict_ca(x, list(curves$a, curves$b))
    )
  }
  MECHANISMS <- c("AD", "MU", "DO", "SAM", "CA")
  TRUTHS <- c(MECHANISMS, "synergy", "antagonism")

  # True relative response of a cell (doses x) under a truth.
  true_response <- function(truth, x, curves) {
    if (any(x == 0)) {
      return(if (x[[1L]] > 0) curves$a$w(x[[1L]]) else if (x[[2L]] > 0) curves$b$w(x[[2L]]) else 1)
    }
    if (truth %in% MECHANISMS) {
      return(prediction(truth, x, curves))
    }
    predictions <- vapply(MECHANISMS, prediction, numeric(1), x = x, curves = curves)
    w <- c(curves$a$w(x[[1L]]), curves$b$w(x[[2L]]))
    if (identical(truth, "synergy")) {
      # Half the response of the most harmful mechanism. If that one already
      # predicts zero, the combination has only zeros and cannot be tested: it
      # does not count.
      0.5 * min(predictions)
    } else {
      min(w) + 0.5 * (1 - min(w))
    }
  }

  ## Designs and responses ----------------------------------------------------

  DESIGNS <- list(
    factorial_2x2 = expand.grid(a = c(0, 1), b = c(0, 1)),
    dose_3x3 = expand.grid(a = c(0, 1, 2), b = c(0, 1, 2))
  )
  UNITS <- 4L
  RESPONSES <- c("survival", "counts")
  NULL_MODELS <- c("AD", "MU", "DO", "BL", "SAM", "CA")

  layout <- function(design) {
    cells <- DESIGNS[[design]]
    d <- cells[rep(seq_len(nrow(cells)), each = UNITS), ]
    d$unit <- seq_len(nrow(d))
    rownames(d) <- NULL
    d
  }
  MUST_DESIGNS <- lapply(names(DESIGNS), function(name) {
    describe_design(layout(name), c("a", "b"), control = c(a = "0", b = "0"), unit = "unit")
  })
  names(MUST_DESIGNS) <- names(DESIGNS)

  simulate_data <- function(d, response, w) {
    if (identical(response, "survival")) {
      d$n <- 20L
      d$alive <- stats::rbinom(nrow(d), 20L, 0.9 * w)
    } else {
      d$y <- stats::rpois(nrow(d), 50 * w * exp(stats::rnorm(nrow(d), 0, 0.25)))
    }
    d
  }
  fit_model <- function(d, response) {
    if (identical(response, "survival")) {
      stats::glm(cbind(alive, n - alive) ~ factor(a) * factor(b), family = stats::binomial(),
                 data = d)
    } else {
      stats::glm(y ~ factor(a) * factor(b), family = stats::quasipoisson(), data = d)
    }
  }

  ## Simulation ---------------------------------------------------------------

  scenarios <- expand.grid(design = names(DESIGNS), response = RESPONSES,
                           effects = names(EFFECTS), truth = TRUTHS,
                           stringsAsFactors = FALSE)
  tasks <- merge(scenarios, data.frame(replicate = seq_len(REPLICATES)))
  tasks$seed <- SEED + seq_len(nrow(tasks))

  run_task <- function(k) {
    task <- tasks[k, ]
    curves <- EFFECTS[[task$effects]]
    d <- layout(task$design)
    w <- vapply(seq_len(nrow(d)), function(i) {
      true_response(task$truth, c(d$a[[i]], d$b[[i]]), curves)
    }, numeric(1))
    set.seed(task$seed)
    d <- simulate_data(d, task$response, w)
    model <- suppressWarnings(fit_model(d, task$response))
    null_models <- if (identical(task$design, "dose_3x3")) NULL_MODELS else setdiff(NULL_MODELS, "CA")
    est <- tryCatch(
      suppressMessages(estimate_interactions(
        model, MUST_DESIGNS[[task$design]], null_models = null_models, n_sim = N_SIM,
        seed = task$seed, overall = FALSE, p_adjust = "none"
      )),
      error = function(e) NULL
    )
    if (is.null(est)) {
      return(NULL)
    }
    e <- est$estimates
    e <- e[e$component == "interaction" & e$estimand %in% "net", , drop = FALSE]
    data.frame(
      task[rep(1L, nrow(e)), c("design", "response", "effects", "truth", "replicate")],
      combination = e$levels, null_model = e$null_model, estimate = e$estimate,
      p_value = e$p_value, status = e$status, stringsAsFactors = FALSE, row.names = NULL
    )
  }

  if (is.null(FROM_RAW)) {
    start <- Sys.time()
    # A failure in one replicate (an impossible fit, say) does not stop the rest.
    safe <- function(k) tryCatch(run_task(k), error = function(e) NULL)
    parts <- if (CORES > 1L && .Platform$OS.type == "unix") {
      parallel::mclapply(seq_len(nrow(tasks)), safe, mc.cores = CORES, mc.preschedule = TRUE)
    } else {
      lapply(seq_len(nrow(tasks)), safe)
    }
    valid <- vapply(parts, is.data.frame, logical(1))
    if (any(!valid)) {
      message(sum(!valid), " replicates without results.")
    }
    raw <- do.call(rbind, parts[valid])
    attr(raw, "minutes") <- as.numeric(difftime(Sys.time(), start, units = "mins"))
    attr(raw, "n_sim") <- N_SIM
    attr(raw, "date") <- format(Sys.Date())
    attr(raw, "r_version") <- R.version.string
    attr(raw, "must_version") <- as.character(utils::packageVersion("must"))
    if (!is.null(SAVE_RAW)) {
      saveRDS(raw, SAVE_RAW)
    }
  } else {
    raw <- readRDS(FROM_RAW)
    REPLICATES <- max(raw$replicate)
    if (!is.null(attr(raw, "n_sim"))) N_SIM <- attr(raw, "n_sim")
  }
  run <- lapply(c(date = "date", must_version = "must_version", r_version = "r_version"),
                function(a) if (is.null(attr(raw, a))) "unknown" else attr(raw, a))
  run$minutes <- attr(raw, "minutes")

  ## Summary ------------------------------------------------------------------

  tested <- raw$status %in% c("ok", "boundary") & is.finite(raw$p_value)
  raw$signal <- ifelse(tested & raw$p_value < ALPHA, sign(raw$estimate), 0)
  raw$tested <- tested

  # Consensus: an interaction only if every null model tested in that
  # combination flags it with the same sign.
  key <- paste(raw$design, raw$response, raw$effects, raw$truth, raw$replicate,
               raw$combination, sep = "|")
  consensus <- do.call(rbind, lapply(split(seq_len(nrow(raw)), key), function(i) {
    r <- raw[i, ]
    all_tested <- all(r$tested)
    s <- if (all_tested && all(r$signal == 1)) 1 else if (all_tested && all(r$signal == -1)) -1 else 0
    data.frame(r[1L, c("design", "response", "effects", "truth", "replicate", "combination")],
               null_model = "consensus", estimate = NA_real_, p_value = NA_real_,
               status = if (all_tested) "ok" else "not_estimable", signal = s,
               tested = all_tested, stringsAsFactors = FALSE, row.names = NULL)
  }))
  rows <- rbind(as.data.frame(raw), consensus)

  expected <- ifelse(rows$truth == "synergy", -1, ifelse(rows$truth == "antagonism", 1, 0))
  rows$false_interaction <- expected == 0 & rows$signal != 0
  rows$detected <- expected != 0 & rows$signal == expected

  # One row per design, response, effects, truth and null model. `rate` is the
  # share of tested combinations with a false interaction (under a mechanism)
  # or with the interaction detected (under synergy or antagonism);
  # `share_negative`, the share of flagged interactions with a negative
  # estimate (synergy, for these responses).
  by_scenario <- do.call(rbind, lapply(
    split(seq_len(nrow(rows)), paste(rows$design, rows$response, rows$effects, rows$truth,
                                     rows$null_model, sep = "|")),
    function(i) {
      x <- rows[i, ]
      flagged <- x$tested & x$signal != 0
      data.frame(
        x[1L, c("design", "response", "effects", "truth", "null_model")],
        replicates = REPLICATES, n_sim = N_SIM,
        combinations = nrow(x),
        tested = mean(x$tested),
        rate = if (x$truth[[1L]] %in% MECHANISMS) {
          mean(x$false_interaction[x$tested])
        } else {
          mean(x$detected[x$tested])
        },
        share_negative = if (any(flagged)) mean(x$signal[flagged] < 0) else NA_real_,
        stringsAsFactors = FALSE, row.names = NULL
      )
    }
  ))
  by_scenario <- by_scenario[order(by_scenario$design, by_scenario$response,
                                   by_scenario$effects, by_scenario$truth,
                                   by_scenario$null_model, method = "radix"), ]
  rownames(by_scenario) <- NULL

  # Means: each mechanism with the same weight; for each null model, over the
  # designs, responses and effect sizes in which it is tested.
  COLUMNS <- c(NULL_MODELS, "consensus")
  mean_rate <- function(r) {
    m <- tapply(r$rate, r$null_model, mean, na.rm = TRUE)
    stats::setNames(as.numeric(m[COLUMNS]), COLUMNS)
  }
  false_rates <- do.call(rbind, lapply(MECHANISMS, function(m) {
    mean_rate(by_scenario[by_scenario$truth == m, ])
  }))
  rownames(false_rates) <- MECHANISMS
  power <- rbind(
    synergy = mean_rate(by_scenario[by_scenario$truth == "synergy", ]),
    antagonism = mean_rate(by_scenario[by_scenario$truth == "antagonism", ])
  )
  # Direction of the false interactions, weighted by the number flagged.
  negative <- do.call(rbind, lapply(MECHANISMS, function(m) {
    r <- by_scenario[by_scenario$truth == m, ]
    flagged <- r$rate * r$tested * r$combinations
    vapply(COLUMNS, function(nm) {
      j <- which(r$null_model == nm & flagged > 0)
      if (length(j) > 0L) sum(r$share_negative[j] * flagged[j]) / sum(flagged[j]) else NA_real_
    }, numeric(1))
  }))
  rownames(negative) <- MECHANISMS
  summary_table <- data.frame(
    null_model = COLUMNS,
    false_mean = colMeans(false_rates, na.rm = TRUE),
    false_max = apply(false_rates, 2L, max, na.rm = TRUE),
    power_mean = colMeans(power, na.rm = TRUE),
    in_two_level_factorials = COLUMNS != "CA",
    stringsAsFactors = FALSE, row.names = NULL
  )

  # The rule. The consensus is not a null model: it is reported, but not part
  # of the rule.
  models <- summary_table[summary_table$null_model != "consensus", ]
  candidates <- models[models$in_two_level_factorials &
                         models$power_mean >= 0.8 * max(models$power_mean), ]
  candidates <- candidates[order(candidates$false_mean, -candidates$power_mean), ]
  best <- candidates$null_model[[1L]]
  tied <- candidates$null_model[abs(candidates$false_mean - candidates$false_mean[[1L]]) < 0.01]
  if (length(tied) > 1L) {
    best <- tied[which.max(candidates$power_mean[match(tied, candidates$null_model)])]
  }

  ## Report -------------------------------------------------------------------

  percent <- function(x) ifelse(is.na(x), "-", sprintf("%.1f", 100 * x))
  table_md <- function(m) {
    head <- paste0("| ", paste(c("", colnames(m)), collapse = " | "), " |")
    rule <- paste0("|", paste(rep("---", ncol(m) + 1L), collapse = "|"), "|")
    body <- vapply(seq_len(nrow(m)), function(i) {
      paste0("| ", rownames(m)[[i]], " | ", paste(percent(m[i, ]), collapse = " | "), " |")
    }, character(1))
    c(head, rule, body)
  }
  # One row per design, response and effects; one column per null model with
  # "false interactions / power".
  breakdown_md <- function(bs) {
    groups <- unique(bs[, c("design", "response", "effects")])
    cells <- vapply(seq_len(nrow(groups)), function(i) {
      g <- bs[bs$design == groups$design[[i]] & bs$response == groups$response[[i]] &
                bs$effects == groups$effects[[i]], ]
      values <- vapply(COLUMNS, function(m) {
        x <- g[g$null_model == m, ]
        if (nrow(x) == 0L) {
          return("-")
        }
        paste(percent(mean(x$rate[x$truth %in% MECHANISMS], na.rm = TRUE)),
              percent(mean(x$rate[!x$truth %in% MECHANISMS], na.rm = TRUE)), sep = " / ")
      }, character(1))
      paste0("| ", paste(groups$design[[i]], groups$response[[i]], groups$effects[[i]]),
             " | ", paste(values, collapse = " | "), " |")
    }, character(1))
    c(paste0("| Scenario | ", paste(COLUMNS, collapse = " | "), " |"),
      paste0("|", paste(rep("---", length(COLUMNS) + 1L), collapse = "|"), "|"),
      cells)
  }

  report <- c(
    "# Robustness of the null models across mechanisms",
    "",
    sprintf(paste0("Simulated with null_model_robustness.R on %s, with must %s and %s: ",
                   "%d replicate%s per scenario and %d Krinsky-Robb draws per test%s."),
            run$date, run$must_version, run$r_version, REPLICATES,
            if (REPLICATES == 1L) "" else "s", N_SIM,
            if (is.null(run$minutes)) {
              ""
            } else if (run$minutes < 1) {
              " (under a minute)"
            } else {
              sprintf(" (%.0f minutes)", run$minutes)
            }),
    "",
    paste("Rates in %, over the combinations tested: 2 x 2 factorial and 3 x 3 dose designs,",
          "survival (binomial) and counts (quasi-Poisson), moderate and weak effects;",
          "p < 0.05 without adjustment in each combination. CA is tested in the dose design",
          "only. Consensus: an interaction only if every null model flags it with the same",
          "sign."),
    "",
    "## False interactions (no interaction, by the mechanism that generated the data)",
    "",
    table_md(false_rates),
    "",
    "## Power against interactions that no mechanism explains",
    "",
    table_md(power),
    "",
    "## Direction of the false interactions",
    "",
    "Share of the false interactions with a negative estimate (synergy, for these responses), in %.",
    "",
    table_md(negative),
    "",
    "## Rule",
    "",
    paste("Among the null models defined in two-level factorials and with an average power of",
          "at least 80% of the highest, the one with the lowest average rate of false",
          "interactions (the five mechanisms with the same weight); within one percentage",
          "point, the more powerful one."),
    "",
    sprintf("Result: **%s**.", best),
    "",
    "## Summary by null model",
    "",
    "| Null model | False interactions, mean | False interactions, maximum | Power, mean |",
    "|---|---|---|---|",
    sprintf("| %s | %s | %s | %s |", summary_table$null_model, percent(summary_table$false_mean),
            percent(summary_table$false_max), percent(summary_table$power_mean)),
    "",
    "## By design, response and size of the effects",
    "",
    "Mean rate of false interactions (five mechanisms) / mean power (synergy and antagonism), in %.",
    "",
    breakdown_md(by_scenario),
    ""
  )

  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(by_scenario, file.path(OUTPUT_DIR, "null_model_robustness.csv"),
                   row.names = FALSE)
  writeLines(report, file.path(OUTPUT_DIR, "null_model_robustness.md"))
  cat(report, sep = "\n")
  message("Results written to ", normalizePath(OUTPUT_DIR))
  invisible(list(by_scenario = by_scenario, summary = summary_table, best = best))
})
