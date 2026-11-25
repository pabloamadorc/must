# estimate_community(): the null models of Thompson et al. (2018) on the
# examples of their Fig. 1, the point estimates against the means of the
# fitted models, the correlation of the taxa, the statuses and the methods
# shared with estimate_interactions().

community_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      d <- community_counts()
      des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank",
                             time = "day")
      cache <<- fit_responses(d, c("sp1", "sp2", "sp3", "sp4", "sp5", "mass"), des,
                              residual_checks = FALSE, seed = 1, quiet = TRUE)
    }
    cache
  }
})

# Statistics of a community of deterministic means (one row per cell:
# reference, a, b and their combination; one column per taxon). With
# Poisson taxa of mean 1000 or 0, a taxon is present when its mean is
# positive, as in the simulated communities of Thompson et al. (2018).
thompson_statistics <- function(means, property = "richness") {
  n_taxa <- ncol(means)
  mu <- array(means, dim = c(1L, 4L, n_taxa))
  info <- lapply(seq_len(n_taxa), function(i) {
    list(family = "poisson", power = NA_real_, disp = rep(NA_real_, 4L))
  })
  names(info) <- paste0("sp", seq_len(n_taxa))
  layout <- list(m = 2L, pos = matrix(1:4, ncol = 1L), cols = 1:4, complete = TRUE)
  res <- community_statistics(mu, layout, info, property, community_null_models(),
                              point = TRUE)[[property]][[1L]]
  out <- res[1L, ]
  names(out)[1:2] <- c("a", "b")
  as.list(c(out, reference = attr(res, "reference"), observed = attr(res, "observed")))
}

test_that("the null models reproduce the three cases of Thompson et al. (2018, Fig. 1)", {
  full <- rep(1000, 12)
  # (a) Redundant effects: both stressors eliminate the same 7 of 12 taxa.
  a <- c(rep(0, 7), rep(1000, 5))
  res <- thompson_statistics(rbind(full, a, a, a))
  expect_equal(c(res$a, res$b), rep(-7 / 12, 2))
  expect_equal(res$combined, -7 / 12)
  # Additive on the property: -14/12 (116 % of loss in the figure), limited to -1.
  expect_equal(res$AD, -1)
  expect_equal(res$MU, (5 / 12)^2 - 1)
  expect_equal(res$AD_taxa, -7 / 12)
  expect_equal(res$MU_taxa, -7 / 12)
  # (b) Effects sub-lethal alone for four taxa, two of which cannot persist
  # with both stressors (60 % of loss with each).
  a <- c(rep(0, 4), rep(1000, 4), 400, 400, 700, 700)
  b <- c(rep(1000, 4), rep(0, 4), 400, 400, 700, 700)
  ab <- c(rep(0, 10), 400, 400)
  res <- thompson_statistics(rbind(full, a, b, ab))
  expect_equal(res$combined, -10 / 12)
  expect_equal(res$AD, -8 / 12)
  expect_equal(res$MU, (8 / 12)^2 - 1)
  expect_equal(res$AD_taxa, -10 / 12)
  expect_equal(res$MU_taxa, -8 / 12)
  # (c) Opposite effects: b offsets the loss of two of the six taxa that a
  # eliminates.
  a <- c(rep(0, 6), rep(1000, 6))
  b <- c(2500, 2500, rep(1000, 10))
  ab <- c(1500, 1500, rep(0, 4), rep(1000, 6))
  res <- thompson_statistics(rbind(full, a, b, ab))
  expect_equal(c(res$a, res$b), c(-0.5, 0))
  expect_equal(res$combined, -4 / 12)
  expect_equal(res$AD, -0.5)
  expect_equal(res$MU, -0.5)
  expect_equal(res$AD_taxa, -4 / 12)
  expect_equal(res$MU_taxa, -0.5)
  # Total abundance of case (c): additive on each taxon and on the total
  # agree when no taxon is expected below zero.
  res <- thompson_statistics(rbind(full, a, b, ab), "abundance")
  expect_equal(res$AD_taxa, 9000 / 12000 - 1)
  expect_equal(res$AD, res$AD_taxa)
  expect_equal(res$reference, 12000)
  expect_equal(res$observed, 9000)
})

test_that("the multiplicative model of each taxon needs the taxon in the reference", {
  x0 <- matrix(c(0, 0, 10), nrow = 1L)
  xs <- list(matrix(c(0, 5, 5), nrow = 1L), matrix(c(0, 2, 20), nrow = 1L))
  out <- expected_taxa_means(x0, xs, "MU_taxa")
  expect_equal(out[1L, 1L], 0)
  expect_true(is.nan(out[1L, 2L]))
  expect_equal(out[1L, 3L], 10 * (5 / 10) * (20 / 10))
  expect_equal(expected_taxa_means(x0, xs, "AD_taxa")[1L, ], c(0, 7, 15))
})

test_that("taxa whose persistence is uncertain are those eliminated in some draws", {
  taxa <- c("a", "b", "c", "d")
  expect_equal(uncertain_persistence(c(0, 300, 980, 1000), 1000, taxa), "b")
  expect_equal(uncertain_persistence(c(200, 800, 201, 799), 1000, taxa), c("c", "d"))
  expect_equal(uncertain_persistence(NULL, 1000, taxa), character())
  values <- persistence_values(c("b", "c"), 0.95)
  expect_equal(values$status, "not_estimable")
  expect_true(is.na(values$estimate))
  expect_match(values$reason, "'b', 'c'")
})

test_that("the probability of absence of a NB1 count matches dnbinom", {
  expect_equal(prob_boundary_obs("nbinom1", 3, NA, 2, NA, "zero"),
               stats::dnbinom(0, size = 3 / 2, mu = 3))
  expect_equal(presence_probability(c(0, 2), "nbinom2", 1.5, NA),
               1 - stats::dnbinom(0, size = 1.5, mu = c(0, 2)))
  p <- 1.4
  expect_equal(presence_probability(2, "tweedie", 0.8, p),
               1 - exp(-2^(2 - p) / (0.8 * (2 - p))))
})

test_that("the shrinkage of the correlation follows Schäfer and Strimmer (2005)", {
  set.seed(3)
  n <- 30
  z <- stats::rnorm(n)
  D <- cbind(z + stats::rnorm(n), z + stats::rnorm(n), stats::rnorm(n))
  D <- sweep(D, 2L, colMeans(D))
  # One observation per unit: their estimator (eq. 10 and Appendix A).
  out <- shrink_correlation(D, unit = seq_len(n))
  x <- scale(D)
  r <- stats::cor(D)
  var_r <- matrix(0, 3, 3)
  for (i in 1:3) for (j in 1:3) {
    w <- x[, i] * x[, j]
    var_r[i, j] <- n / (n - 1)^3 * sum((w - mean(w))^2)
  }
  off <- upper.tri(r)
  lambda <- min(1, max(0, sum(var_r[off]) / sum(r[off]^2)))
  expect_equal(out$lambda, lambda)
  expect_equal(out$raw, r, ignore_attr = TRUE)
  expect_equal(out$matrix[off], (1 - lambda) * r[off])
  expect_equal(unname(diag(out$matrix)), rep(1, 3))
  # Repeated observations of the same units do not count as independent.
  twice <- rbind(D, D)
  clustered <- shrink_correlation(twice, unit = rep(seq_len(n), 2))
  naive <- shrink_correlation(twice, unit = seq_len(2 * n))
  expect_equal(clustered$raw, naive$raw)
  expect_gt(clustered$lambda, naive$lambda)
  # Without correlation there is nothing to keep.
  expect_equal(shrink_correlation(cbind(c(1, -1, 0, 0), c(0, 0, 1, -1)), 1:4)$lambda, 1)
})

# Means of each taxon (columns) in the cells of the first set at day 7 (rows:
# reference, a, b and a + b), from predict() on the maximum likelihood fits.
hand_means <- function(fit, taxa) {
  sc <- fit$design$set_cells[[fit$design$sets$set_id[[1L]]]]
  grid <- fit$fits[[taxa[[1L]]]]$grid
  ct <- as.character(grid$.cell_time[match(paste(sc$cell, "7"), paste(grid$.cell, grid$.time))])
  x <- vapply(fit$fits[taxa], function(f) {
    p <- stats::predict(f$model, re.form = NA, type = "response")
    p[match(ct, as.character(f$model$frame$.cell_time))]
  }, numeric(4))
  x[x < 1e-6] <- 0
  x
}

# Probability that each taxon is found in a sample, from its model.
hand_presence <- function(fit, taxa, means) {
  vapply(seq_along(taxa), function(i) {
    m <- fit$fits[[taxa[[i]]]]$model
    if (stats::family(m)$family == "nbinom2") {
      1 - stats::dnbinom(0, size = stats::sigma(m), mu = means[[i]])
    } else {
      1 - stats::dpois(0, means[[i]])
    }
  }, numeric(1))
}

# The rows of day 7 of a property against the hand calculation, for the null
# models in `models`: R is the property in the four cells, `xhat` the means
# that AD_taxa and MU_taxa expect, `value` the property of a set of means.
expect_hand <- function(e, prop, R, xhat, value, models) {
  row <- function(model, component, stressor = NA) {
    e[e$response == prop & e$time == "7" & e$null_model == model & e$component == component &
        (is.na(stressor) | e$stressor %in% stressor), ]
  }
  single <- R[2:3] / R[1] - 1
  combined <- R[[4L]] / R[[1L]] - 1
  expected <- c(AD_taxa = value(xhat$AD_taxa) / R[[1L]] - 1,
                MU_taxa = value(xhat$MU_taxa) / R[[1L]] - 1,
                AD = max(-1, sum(single)), MU = prod(1 + single) - 1)
  for (nm in models) {
    expect_equal(row(nm, "single", "a")$estimate, single[[1L]], tolerance = 1e-6)
    expect_equal(row(nm, "single", "b")$estimate, single[[2L]], tolerance = 1e-6)
    expect_equal(row(nm, "combined")$estimate, combined, tolerance = 1e-6)
    expect_equal(row(nm, "expected")$estimate, expected[[nm]], tolerance = 1e-6)
    int <- row(nm, "interaction")
    expect_equal(int$estimate, combined - expected[[nm]], tolerance = 1e-6)
    expect_equal(int$observed, R[[4L]], tolerance = 1e-6)
    expect_equal(int$expected, R[[1L]] * (1 + expected[[nm]]), tolerance = 1e-6)
  }
}

expected_means <- function(x) {
  list(AD_taxa = pmax(x[1L, ] + (x[2L, ] - x[1L, ]) + (x[3L, ] - x[1L, ]), 0),
       MU_taxa = x[2L, ] * x[3L, ] / x[1L, ])
}

test_that("the effects on the properties follow the means of the fitted models", {
  fit <- community_fit()
  taxa <- c("sp1", "sp2", "sp3", "sp4")
  # Maximum likelihood fits (no REML refit), whose means predict() gives.
  cm <- estimate_community(fit, taxa = taxa, correction = "none", n_sim = 1000, seed = 1)
  expect_s3_class(cm, c("must_community", "must_interactions"))
  e <- cm$estimates
  expect_setequal(unique(e$response), c("abundance", "richness"))
  expect_setequal(unique(e$null_model), community_null_models())
  expect_true(all(e$scale == "response / control"))
  x <- hand_means(fit, taxa)
  expect_equal(unname(x[4L, "sp3"]), 0)
  xhat <- expected_means(x)
  expect_hand(e, "abundance", rowSums(x), xhat, sum, community_null_models())
  R <- vapply(1:4, function(k) sum(hand_presence(fit, taxa, x[k, ])), numeric(1))
  expect_hand(e, "richness", R, xhat, function(v) sum(hand_presence(fit, taxa, v)),
              c("MU_taxa", "AD", "MU"))
  # sp1 is expected at 7.5 with both stressors, far less than its uncertainty:
  # its persistence under AD_taxa is too uncertain for the richness.
  flagged <- e[e$response == "richness" & e$null_model == "AD_taxa" &
                 e$component %in% c("expected", "interaction"), ]
  expect_true(all(flagged$status == "not_estimable"))
  expect_match(flagged$reason[flagged$time == "7"], "'sp1'")
  expect_true(all(e$status[e$response == "abundance"] == "ok"))
  # Intervals, tests and the average over times.
  ok <- e[e$status == "ok", ]
  expect_true(all(ok$lower <= ok$estimate & ok$estimate <= ok$upper))
  expect_true(all(is.na(ok$p_value[ok$component == "expected"])))
  expect_true(all(ok$p_value[ok$component != "expected"] > 0 &
                    ok$p_value[ok$component != "expected"] <= 1))
  expect_true(all(is.finite(ok$p_adjusted[ok$component != "expected"])))
  by_time <- e[e$response == "abundance" & e$null_model == "AD" & e$component == "interaction", ]
  expect_equal(by_time$estimate[by_time$time == "overall"],
               mean(by_time$estimate[by_time$time != "overall"]))
  expect_equal(unique(ok$df), Inf)
  # Taxa whose persistence is clear: the richness under AD_taxa is estimated.
  taxa <- c("sp2", "sp4", "sp5")
  cm <- estimate_community(fit, taxa = taxa, null_models = c("AD_taxa", "AD"),
                           correction = "none", n_sim = 1000, seed = 1)
  e <- cm$estimates
  expect_true(all(e$status == "ok"))
  x <- hand_means(fit, taxa)
  xhat <- expected_means(x)
  expect_hand(e, "abundance", rowSums(x), xhat, sum, c("AD_taxa", "AD"))
  R <- vapply(1:4, function(k) sum(hand_presence(fit, taxa, x[k, ])), numeric(1))
  expect_hand(e, "richness", R, xhat, function(v) sum(hand_presence(fit, taxa, v)),
              c("AD_taxa", "AD"))
})

test_that("the correlation of the taxa links their draws but not their estimates", {
  fit <- community_fit()
  taxa <- c("sp1", "sp2", "sp3", "sp4")
  linked <- estimate_community(fit, taxa = taxa, n_sim = 500, seed = 4)
  apart <- estimate_community(fit, taxa = taxa, n_sim = 500, seed = 4,
                              taxa_correlation = "none")
  again <- estimate_community(fit, taxa = taxa, n_sim = 500, seed = 4)
  expect_identical(linked$estimates, again$estimates)
  tested <- linked$estimates[linked$estimates$status == "ok", ]
  expect_equal(unique(tested$df), min(vapply(fit$fits[taxa], function(f) {
    set_df(f$df_info, fit$design$sets[1L, ], "small_sample")
  }, numeric(1))))
  expect_equal(linked$estimates$estimate, apart$estimates$estimate)
  expect_true(is.na(apart$settings$shrinkage))
  expect_equal(unname(apart$diagnostics$correlation), diag(4))
  corr <- linked$diagnostics$correlation
  expect_equal(dimnames(corr), list(taxa, taxa))
  expect_equal(unname(diag(corr)), rep(1, 4))
  expect_true(linked$settings$shrinkage >= 0 && linked$settings$shrinkage <= 1)
  expect_equal(linked$diagnostics$n_units, 16L)
  expect_equal(linked$diagnostics$n_obs, 32L)
  expect_equal(nrow(fit$fits$sp1$observations), fit$fits$sp1$n_obs)
  expect_equal(fit$fits$sp1$observations$unit, community_counts()$tank)
  # Fits made before `observations` was kept give the same correlation.
  old <- fit
  for (taxon in taxa) old$fits[[taxon]]$observations <- NULL
  expect_equal(residual_correlation(old$fits[taxa]), residual_correlation(fit$fits[taxa]))
})

test_that("a single sampling gives no averages over times", {
  d <- community_counts()
  d <- d[d$day == 7, ]
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank")
  fit <- fit_responses(d, c("sp1", "sp2", "sp4"), des, residual_checks = FALSE, seed = 1,
                       quiet = TRUE)
  cm <- estimate_community(fit, n_sim = 500, seed = 1)
  e <- cm$estimates
  expect_equal(nrow(e), 2L * 4L * 5L)
  expect_true(all(is.na(e$time)))
  compositional_richness <- e$response == "richness" & e$null_model == "AD_taxa"
  expect_true(all(e$status[!compositional_richness] == "ok"))
  expect_equal(cm$diagnostics$n_units, 16L)
  expect_equal(fit$fits$sp1$observations$cell_time,
               as.character(fit$fits$sp1$model$frame$.group))
})

test_that("statuses explain what the null models and properties cannot give", {
  fit <- community_fit()
  # A taxon only found with `a` alone: no multiplicative prediction.
  cm <- estimate_community(fit, taxa = c("sp1", "sp5"), properties = "abundance",
                           n_sim = 500, seed = 1)
  e <- cm$estimates
  mu_taxa <- e[e$null_model == "MU_taxa" & e$component %in% c("expected", "interaction"), ]
  expect_true(all(mu_taxa$status == "not_defined"))
  expect_match(mu_taxa$reason[[1L]], "'sp5'")
  expect_true(all(e$status[e$null_model != "MU_taxa"] == "ok"))
  expect_true(all(e$status[e$component %in% c("single", "combined")] == "ok"))
  # A biomass without zeros has no probability of absence.
  cm <- estimate_community(fit, taxa = c("sp1", "mass"), null_models = c("AD_taxa", "AD"),
                           n_sim = 500, seed = 1)
  e <- cm$estimates
  expect_true(all(e$status[e$response == "richness"] == "not_defined"))
  expect_match(e$reason[e$response == "richness"][[1L]], "'mass'")
  expect_true(all(e$status[e$response == "abundance"] == "ok"))
  expect_match(cm$diagnostics$richness, "tweedie")
  expect_output(print(cm), "Richness not defined")
})

test_that("community estimates work with the methods of estimate_interactions()", {
  fit <- community_fit()
  cm <- estimate_community(fit, taxa = c("sp1", "sp2", "sp3", "sp4"), n_sim = 500, seed = 1)
  expect_output(print(cm), "<must_community> 4 taxa")
  expect_output(print(cm), "additive on each taxon")
  cls <- classify_interactions(cm)
  expect_setequal(unique(cls$null_model), community_null_models())
  # The richness under AD_taxa is not estimable (persistence of sp1).
  classified <- !(cls$response == "richness" & cls$null_model == "AD_taxa")
  expect_true(all(cls$status[classified] == "ok"))
  expect_true(all(cls$status[!classified] == "not_estimable"))
  expect_true(all(cls$class[classified] %in% c("AD", "+S", "-S", "+A", "-A")))
  expect_match(class_meaning("+S", "AD_taxa"), "adding the effects on each taxon")
  cmp <- compare_null_models(cm, null_models = c("AD_taxa", "AD"))
  expect_named(cmp, c("response", "context", "time", "set_id", "order", "stressors", "levels",
                      "background", "estimand", "AD_taxa", "AD", "n_tested", "n_deviating",
                      "agreement", "direction"))
  expect_error(compare_null_models(cm, null_models = c("AD", "BL")), class = "must_error_input")
  s <- summary(cm)
  expect_setequal(unique(s$interactions$response), c("abundance", "richness"))
  grDevices::pdf(NULL)
  device <- grDevices::dev.cur()
  withr::defer(grDevices::dev.off(device))
  expect_invisible(plot(cm, response = "richness", null_model = "AD_taxa"))
  expect_invisible(plot(cm, type = "models"))
  expect_invisible(plot(cm, type = "time"))
  expect_error(report_interactions(cm), "estimate_community", class = "must_error_input")
})

test_that("estimate_community() checks its input", {
  fit <- community_fit()
  expect_error(estimate_community(list()), "fit_responses", class = "must_error_input")
  expect_error(estimate_community(fit, taxa = "sp1"), "two or more", class = "must_error_input")
  expect_error(estimate_community(fit, taxa = c("sp1", "other")), class = "must_error_input")
  expect_error(estimate_community(fit, properties = "evenness"), class = "must_error_input")
  expect_error(estimate_community(fit, null_models = c("AD", "AD")), class = "must_error_input")
  expect_error(estimate_community(fit, null_models = "BL"), class = "must_error_input")
  expect_error(estimate_community(fit, taxa_correlation = "full"), class = "must_error_input")
  expect_error(estimate_community(fit, n_sim = 10), class = "must_error_input")
  broken <- fit
  broken$fits$sp2$status <- "failed"
  expect_error(estimate_community(broken, taxa = c("sp1", "sp2")), "'sp2' has no fitted model",
               class = "must_error_input")
  broken <- fit
  broken$fits$sp2$kind <- "binomial"
  expect_error(estimate_community(broken, taxa = c("sp1", "sp2")), "proportions",
               class = "must_error_input")
  broken <- fit
  broken$fits$sp2$non_negative <- FALSE
  expect_error(estimate_community(broken, taxa = c("sp1", "sp2")), "negative",
               class = "must_error_input")
  broken <- fit
  broken$fits$sp2$grid <- broken$fits$sp2$grid[-1L, ]
  expect_error(estimate_community(broken, taxa = c("sp1", "sp2")), "same treatment",
               class = "must_error_input")
})
