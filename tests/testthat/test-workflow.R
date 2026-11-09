# Integration test: the documented workflow (design, fit, estimate, classify)
# on the example data.

test_that("the documented workflow runs end to end on polazzo_zooplankton", {
  # Day -5 precedes the pesticides but not the nutrients: it is left out
  zoo <- polazzo_zooplankton[polazzo_zooplankton$day > 0, ]
  design <- describe_design(
    zoo,
    stressors = c("insecticide", "herbicide", "nutrients"),
    control = c(insecticide = "absent", herbicide = "absent", nutrients = "ambient"),
    unit = "mesocosm",
    time = "day"
  )
  expect_equal(nrow(design$sets), 4L)
  expect_true(all(design$sets$net_identifiable))
  expect_true(design$sets$pure_identifiable[design$sets$order == 3L])
  expect_false(any(design$checks$severity == "problem"))

  fits <- fit_responses(zoo, c("total_density", "daphnia"), design,
                        residual_checks = FALSE, seed = 1, quiet = TRUE)
  expect_equal(vapply(fits$fits, `[[`, character(1), "status"),
               c(total_density = "ok", daphnia = "ok"))
  expect_match(fits$fits$total_density$formula, "(1 | .unit)", fixed = TRUE)

  est <- estimate_interactions(fits, n_sim = 500, seed = 1)
  e <- est$estimates
  # 4 sets x (4 days + their average) for each response
  int <- e[e$component == "interaction" & e$null_model == "AD" & e$estimand == "net", ]
  expect_equal(nrow(int), 2L * 4L * 5L)
  # 24 mesocosms in 8 treatment combinations
  expect_true(all(int$df[int$response == "total_density"] == 16))
  # Net = pure + pairwise, at every day and on average, for total density under AD
  td <- e[e$response == "total_density" & e$null_model == "AD" & e$component == "interaction", ]
  for (d in unique(td$time)) {
    net <- td$estimate[td$time == d & td$order == 3L & td$estimand == "net"]
    pure <- td$estimate[td$time == d & td$order == 3L & td$estimand == "pure"]
    pairs <- sum(td$estimate[td$time == d & td$order == 2L])
    expect_equal(net, pure + pairs, tolerance = 1e-6)
  }
  # Daphnia: the herbicide alone eliminates it, which is tested from the
  # probability of zeros, and the interactions that need that cell are not
  # estimable
  dap <- e[e$response == "daphnia", ]
  expect_true(any(dap$status == "boundary"))
  expect_true(any(dap$status == "not_estimable"))

  cls <- classify_interactions(est)
  expect_true(any(cls$status[cls$order == 3L] == "ok"))
  expect_true(all(cls$class[cls$status == "ok"] %in% c("AD", "+S", "-S", "+A", "-A")))
  sm <- summary(est)
  expect_s3_class(sm, "summary.must_interactions")
  expect_output(print(sm), "Classes of the net deviations")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(est, response = "total_density", null_model = "MU"))
})

test_that("summary() and plot() check their input", {
  d <- two_by_two_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
  sm <- summary(est)
  expect_equal(sm$interactions$n_sets, c(1L, 1L))
  expect_equal(sum(sm$classes[, c("AD", "+S", "-S", "+A", "-A")]), 2L)
  expect_error(plot(est, null_model = "DO"), class = "must_error_input")
  expect_error(plot(est, response = "other"), class = "must_error_input")
  expect_error(summary(est, alpha = 2), class = "must_error_input")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par("mar")
  plot(est, null_model = "MU", main = "MU")
  expect_equal(graphics::par("mar"), old)
})
