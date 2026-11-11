fit_counts <- function() {
  d <- two_by_two_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                residual_checks = FALSE, quiet = TRUE)
}

# Set the adjusted p-value, estimate and status of the interaction of each
# null model, to test the agreement rules on known verdicts.
set_verdicts <- function(est, p, estimate, status = rep("ok", length(p))) {
  e <- est$estimates
  models <- est$settings$null_models
  for (k in seq_along(models)) {
    i <- e$component == "interaction" & e$null_model == models[[k]]
    e$p_adjusted[i] <- p[[k]]
    e$estimate[i] <- estimate[[k]]
    e$status[i] <- status[[k]]
  }
  est$estimates <- e
  est
}

test_that("agreement follows the verdicts of the null models", {
  est <- estimate_interactions(fit_counts(), null_models = c("AD", "MU", "DO"), n_sim = 500,
                               seed = 1)
  agree <- function(...) compare_null_models(set_verdicts(est, ...))

  r <- agree(p = c(0.01, 0.02, 0.001), estimate = c(-1, -0.3, -5))
  expect_equal(nrow(r), 1L)
  expect_equal(unlist(r[c("AD", "MU", "DO")]), c(AD = "-", MU = "-", DO = "-"))
  expect_equal(r$agreement, "robust")
  expect_equal(r$direction, "negative")
  expect_equal(c(r$n_tested, r$n_deviating), c(3L, 3L))

  r <- agree(p = c(0.01, 0.2, 0.001), estimate = c(2, 0.3, 5))
  expect_equal(r$MU, "0")
  expect_equal(r$agreement, "model_dependent")
  expect_equal(r$direction, "positive")

  r <- agree(p = c(0.01, 0.01, 0.5), estimate = c(2, -0.3, 1))
  expect_equal(r$agreement, "conflicting")
  expect_true(is.na(r$direction))

  r <- agree(p = c(0.3, 0.2, 0.5), estimate = c(2, -0.3, 1))
  expect_equal(r$agreement, "none")
  expect_true(is.na(r$direction))

  # Null models that could not be tested do not count
  r <- agree(p = c(0.01, NA, 0.01), estimate = c(-2, NA, -1),
             status = c("ok", "not_estimable", "boundary"))
  expect_true(is.na(r$MU))
  expect_equal(r$n_tested, 2L)
  expect_equal(r$agreement, "robust")
  r <- agree(p = c(0.01, NA, NA), estimate = c(2, NA, NA),
             status = c("ok", "not_estimable", "not_applicable"))
  expect_equal(r$n_tested, 1L)
  expect_true(is.na(r$agreement))

  # A subset of the null models, at another level
  r <- compare_null_models(set_verdicts(est, p = c(0.04, 0.02, 0.2), estimate = c(-1, -1, -1)),
                           null_models = c("MU", "AD"), alpha = 0.03)
  expect_true(all(c("AD", "MU") %in% names(r)))
  expect_false("DO" %in% names(r))
  expect_equal(r$agreement, "model_dependent")
})

test_that("compare_null_models() checks its input", {
  fit <- fit_counts()
  one <- estimate_interactions(fit, n_sim = 200, seed = 1)
  expect_error(compare_null_models(one), "two or more null models", class = "must_error_input")
  two <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
  expect_error(compare_null_models(two, null_models = c("AD", "BL")), "not estimated under BL",
               class = "must_error_input")
  expect_error(compare_null_models(two, null_models = "XX"), class = "must_error_input")
  expect_error(compare_null_models(list()), "estimate_interactions", class = "must_error_input")
  expect_error(compare_null_models(two, alpha = 1), class = "must_error_input")
})

test_that("net and pure interactions of three stressors are compared separately", {
  d <- three_way_gaussian()
  des <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                         unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("AD", "MU", "DO", "SAM"), n_sim = 500,
                               seed = 3)
  cmp <- compare_null_models(est)
  int <- est$estimates[est$estimates$component == "interaction", ]
  expect_equal(nrow(cmp), nrow(unique(int[c("set_id", "estimand")])))
  pure <- cmp[cmp$order == 3L & cmp$estimand == "pure", ]
  # DO and SAM do not define the pure interaction
  expect_true(is.na(pure$DO) && is.na(pure$SAM))
  expect_equal(pure$n_tested, 2L)
  # Each verdict is that of the adjusted p-value of its null model
  for (nm in c("AD", "MU", "DO", "SAM")) {
    rows <- int[int$null_model == nm, ]
    key <- match(paste(cmp$set_id, cmp$estimand), paste(rows$set_id, rows$estimand))
    tested <- rows$status[key] %in% c("ok", "boundary")
    expect_identical(is.na(cmp[[nm]]), !tested)
    flagged <- tested & rows$p_adjusted[key] < 0.05
    expect_identical(cmp[[nm]] %in% c("+", "-"), flagged)
  }
})

test_that("summary() reports the agreement between null models", {
  fit <- fit_counts()
  est <- estimate_interactions(fit, null_models = c("AD", "MU", "DO"), n_sim = 500, seed = 1)
  sm <- summary(est)
  cmp <- compare_null_models(est)
  expect_equal(sum(unlist(sm$agreement[c("n_robust", "n_model_dependent", "n_conflicting",
                                         "n_none", "n_not_compared")])), nrow(cmp))
  expect_output(print(sm), "Agreement between null models \\(AD, MU, DO\\)")
  # Not with a single null model
  single <- summary(estimate_interactions(fit, n_sim = 200, seed = 1))
  expect_null(single$agreement)
  expect_false(any(grepl("Agreement", utils::capture.output(print(single)))))
})
