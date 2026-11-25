# Resampling tests of the interactions of each response (internal) ------------
#
# With `inference = "resampling"`, estimate_interactions() keeps the
# post-estimation estimates and intervals and replaces the p-values of the
# net interactions under AD and MU with resampling tests of the null model
# (R/resampling.R): a GLM of the family of the fitted model, with the null
# model as a constraint on the means of the set at each time, against the
# cell means, resampling whole units with the PIT-trap. Rows that cannot be
# tested this way keep the post-estimation test and say why.

# Resampling family of a fitted response, or the reason why there is none.
resampling_family_of <- function(f) {
  if (!identical(f$kind, "plain") && !identical(f$kind, "offset")) {
    return(list(problem = paste(
      "proportions cannot be tested by resampling: the identity and log links of AD",
      "and MU leave the interval from 0 to 1 (Madge Pimentel et al. 2025)."
    )))
  }
  model <- f$model
  name <- stats::family(model)$family
  dispersion_by_cell <- identical(f$dispersion, "cell")
  fam <- switch(
    name,
    poisson = rs_family("poisson"),
    nbinom2 = rs_family("nbinom2", theta = if (dispersion_by_cell) NA_real_ else {
      unname(stats::sigma(model))
    }),
    nbinom1 = rs_family("nbinom1"),
    gaussian = rs_family("gaussian"),
    Gamma = rs_family("gamma"),
    lognormal = rs_family("lognormal"),
    tweedie = rs_family("tweedie", power = unname(glmmTMB::family_params(model)[[1L]])),
    NULL
  )
  if (is.null(fam)) {
    return(list(problem = sprintf(
      paste(
        "the %s distribution is not available for resampling (Poisson, negative",
        "binomial, Gaussian, gamma, lognormal and Tweedie are)."
      ),
      name
    )))
  }
  list(family = fam)
}

# The samples of a fitted response: value, cell, time, unit, effort and the
# levels of the stressors of its cell.
response_samples <- function(f, design) {
  frame <- f$model$frame
  grid <- f$grid
  by_time <- ".cell_time" %in% names(frame)
  key <- as.character(if (by_time) frame$.cell_time else frame$.group)
  at <- match(key, as.character(if (by_time) grid$.cell_time else grid$.group))
  obs <- observation_table(f)
  out <- data.frame(
    y = as.numeric(frame[[f$response]]), .context = as.character(grid$.context[at]),
    .cell = as.character(grid$.cell[at]), .time = as.character(grid$.time[at]),
    .unit = obs$unit, stringsAsFactors = FALSE
  )
  offset <- grep("^offset\\(", names(frame), value = TRUE)
  out$.effort <- if (length(offset) > 0L) exp(frame[[offset[[1L]]]]) else NA_real_
  cells <- design$cells
  pos <- match(paste_key(out$.context, out$.cell), paste_key(cells$context, cells$cell))
  for (st in design$spec$stressors) {
    out[[st]] <- as.character(cells[[st]][pos])
  }
  out
}

# Resampling p-values for the net interactions under AD and MU of one
# response (`rows`, its estimates).
resample_interactions <- function(rows, f, design, opts) {
  if (is.null(rows) || nrow(rows) == 0L) {
    return(rows)
  }
  target <- which(rows$component == "interaction" & rows$estimand %in% "net" &
                    rows$null_model %in% c("AD", "MU") & rows$status %in% c("ok", "boundary"))
  if (length(target) == 0L) {
    return(rows)
  }
  fam <- resampling_family_of(f)
  if (!is.null(fam$problem)) {
    return(keep_post_estimation(rows, target, fam$problem))
  }
  data <- response_samples(f, design)
  effort <- if (all(is.na(data$.effort))) NULL else data$.effort
  keys <- unique(rows[target, c("set_id", "time"), drop = FALSE])
  for (k in seq_len(nrow(keys))) {
    set_id <- keys$set_id[[k]]
    tt <- keys$time[[k]]
    here <- target[rows$set_id[target] == set_id & rows$time[target] %in% tt]
    s <- design$sets[design$sets$set_id == set_id, , drop = FALSE]
    sc <- design$set_cells[[set_id]]
    contexts <- set_cell_contexts(sc, s)
    use <- c(which(sc$size == 0L), which(sc$size == 1L), which(sc$size == s$order))
    in_set <- strsplit(s$stressors, " x ", fixed = TRUE)[[1L]]
    in_cells <- paste_key(data$.context, data$.cell) %in% paste_key(contexts[use], sc$cell[use])
    all_times <- identical(tt, "overall")
    times <- if (all_times) {
      setdiff(unique(rows$time[rows$set_id == set_id & rows$component == "interaction"]),
              c("overall", NA))
    } else {
      tt
    }
    sel <- in_cells & (if (is.na(tt)) TRUE else data$.time %in% times)
    whole <- whole_group_problem(design, in_set)
    if (nzchar(whole)) {
      rows <- keep_post_estimation(rows, here, whole)
      next
    }
    smp <- rs_samples(data[sel, , drop = FALSE], in_set, design$spec$control, all_times)
    if (nzchar(smp$problem)) {
      rows <- keep_post_estimation(rows, here, smp$problem)
      next
    }
    fam_k <- fam$family
    if (identical(fam_k$name, "nbinom2") && is.na(fam_k$theta)) {
      fam_k$theta <- nb2_size(data$y[sel], if (is.null(effort)) rep(1, sum(sel)) else effort[sel],
                              smp$layout$group)
    }
    Y <- matrix(data$y[sel], ncol = 1L, dimnames = list(NULL, f$response))
    E <- if (is.null(effort)) NULL else matrix(effort[sel], ncol = 1L)
    for (i in here) {
      res <- tryCatch(rs_test(Y, E, smp$layout, smp$unit, list(fam_k), rows$null_model[[i]],
                              opts$n_boot),
                      error = function(e) e)
      if (inherits(res, "error")) {
        rows <- keep_post_estimation(rows, i, conditionMessage(res))
        next
      }
      rows$p_value[[i]] <- res$p_value
      rows$df[[i]] <- NA_real_
      rows$mc_se[[i]] <- mc_error(res$p_value, opts$n_boot)
      rows$inference[[i]] <- "resampling"
      if (identical(rows$status[[i]], "boundary")) {
        rows$status[[i]] <- "ok"
        rows$reason[[i]] <- ""
      }
      if (all_times) {
        rows$note[[i]] <- add_note(rows$note[[i]], "joint_over_times")
      }
    }
  }
  rows
}

keep_post_estimation <- function(rows, idx, why) {
  for (i in idx) {
    rows$note[[i]] <- add_note(rows$note[[i]], "resampling_not_available")
    rows$reason[[i]] <- trimws(paste(rows$reason[[i]], sprintf(
      "Not tested by resampling: %s", sub("\\.?$", ".", why)
    )))
  }
  rows
}
