design_2x2 <- function(d, unit = "id") {
  describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = unit)
}

test_that("Gaussian cell means equal the sample means", {
  d <- two_by_two_gaussian()
  fit <- fit_responses(d, "y", design_2x2(d), families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  f <- fit$fits$y
  expect_s3_class(fit, "must_fit")
  expect_equal(f$status, "ok")
  beta <- glmmTMB::fixef(f$model)$cond
  expect_equal(unname(beta[paste0(".group", c("a=0 | b=0", "a=1 | b=0", "a=0 | b=1", "a=1 | b=1"))]),
               c(10, 8, 7, 3), tolerance = 1e-6)
  # One row per unit: no random effect is added
  expect_equal(f$formula, "`y` ~ 0 + .group")
  expect_equal(f$n_obs, 20L)
})

test_that("default candidates follow the values of the response", {
  expect_named(default_families(c(0, 3, 5))$families, c("poisson", "nbinom2"))
  expect_named(default_families(c(0.5, 3, 5))$families,
               c("gaussian", "gaussian_log", "gamma_log", "lognormal"))
  expect_named(default_families(c(0, 0.5, 3))$families, c("gaussian", "gaussian_log", "tweedie"))
  expect_named(default_families(c(-1, 0.5, 3))$families, "gaussian")
  # Whole numbers that are not counts can be declared continuous
  expect_equal(default_families(c(14, 15, 20), "continuous")$type, "continuous")
  expect_equal(default_families(c(14, 15, 20), "continuous")$source, "declared")
  expect_equal(default_families(c(14, 15, 20))$source, "detected")
  expect_named(default_families(c(1, 5, 3), kind = "binomial")$families,
               c("binomial", "betabinomial"))
  # With an offset only log links are kept
  fams <- default_families(c(0.5, 3, 5), kind = "offset")$families
  expect_true(all(vapply(fams, function(f) f$link, character(1)) == "log"))
})

test_that("missing responses are dropped for that response and counted", {
  d <- two_by_two_gaussian()
  d$y[c(2, 9)] <- NA
  fit <- fit_responses(d, "y", design_2x2(d), families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  expect_equal(fit$fits$y$n_missing, 2L)
  expect_equal(fit$fits$y$n_obs, 18L)
})

test_that("repeated measures get a random intercept per unit", {
  d <- repeated_counts()
  # Clear differences between tanks, so that the variance is not zero
  d$count <- round(d$count * c(0.5, 1.6, 0.8, 1.3)[as.integer(factor(d$tank)) %% 4 + 1])
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank", time = "day")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  expect_equal(fit$settings$random, "(1 | .unit)")
  expect_match(fit$fits$count$formula, "0 + .cell_time + (1 | .unit)", fixed = TRUE)
  additive <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                            time_structure = "additive", residual_checks = FALSE,
                            quiet = TRUE)
  expect_match(additive$fits$count$formula, "0 + .group + .time", fixed = TRUE)
  none <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                        random = "none", residual_checks = FALSE, quiet = TRUE)
  expect_equal(none$settings$random, "none")
})

test_that("residual checks need DHARMa and say how to skip them", {
  d <- two_by_two_gaussian()
  local_mocked_bindings(is_installed = function(pkg) FALSE)
  expect_error(fit_responses(d, "y", design_2x2(d), quiet = TRUE),
               "residual_checks = FALSE", class = "must_error_missing_package")
})

test_that("candidates that pass the residual checks are preferred", {
  d <- two_by_two_gaussian()
  fams <- list(gaussian = stats::gaussian(), gaussian_log = stats::gaussian(link = "log"))
  local_mocked_bindings(
    is_installed = function(pkg) TRUE,
    residual_tests = function(model, n_sim) {
      ok <- identical(stats::family(model)$link, "log")
      list(p_uniformity = if (ok) 0.5 else 0.01, p_dispersion = 0.5,
           p_zero_inflation = 0.5, ok = ok, message = "")
    }
  )
  fit <- fit_responses(d, "y", design_2x2(d), families = fams, quiet = TRUE)
  expect_equal(fit$fits$y$family, "gaussian_log")
  expect_true(fit$fits$y$checks_passed)
})

test_that("if no candidate passes, the lowest AIC is kept and flagged", {
  d <- two_by_two_gaussian()
  fams <- list(gaussian = stats::gaussian(), gaussian_log = stats::gaussian(link = "log"))
  local_mocked_bindings(
    is_installed = function(pkg) TRUE,
    residual_tests = function(model, n_sim) {
      list(p_uniformity = 0.01, p_dispersion = 0.5, p_zero_inflation = 0.5,
           ok = FALSE, message = "")
    }
  )
  fit <- fit_responses(d, "y", design_2x2(d), families = fams, quiet = TRUE)
  expect_false(fit$fits$y$checks_passed)
  sel <- fit$fits$y$selection
  expect_equal(sel$aic[sel$selected], min(sel$aic))
  expect_match(fit$fits$y$reason, "lowest AIC")
})

test_that("combinations with only zeros are fixed at the boundary and reported", {
  d <- two_by_two_counts()
  d$count[d$a == "1" & d$b == "1"] <- 0
  fit <- fit_responses(d, "count", design_2x2(d),
                       families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  f <- fit$fits$count
  expect_equal(f$status, "ok")
  expect_match(f$reason, "only zeros")
  expect_equal(names(f$boundary), "a=1 | b=1")
  beta <- glmmTMB::fixef(f$model)$cond
  expect_equal(unname(beta[[".groupa=1 | b=1"]]), -20)
  # The fixed parameter is counted in the AIC
  expect_equal(f$selection$aic, stats::AIC(f$model) + 2)
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
  uses_ab <- est$estimates$component %in% c("combined", "interaction")
  expect_true(all(est$estimates$status[uses_ab] == "boundary"))
  expect_true(all(est$estimates$status[!uses_ab] == "ok"))
})

test_that("boundary columns are found whatever the characters of the labels", {
  d <- zero_combination_counts()
  d$site <- paste0("Schleim", "\u00fc", "nde")
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id",
                         by = "site")
  withr::with_locale(c(LC_CTYPE = "C"), {
    fit <- fit_responses(d, "y", des, families = list(poisson = stats::poisson()),
                         residual_checks = FALSE, quiet = TRUE)
  })
  expect_equal(fit$fits$y$status, "ok")
  expect_length(fit$fits$y$boundary, 1L)
})

test_that("with additive time, only combinations with zeros at every time are fixed", {
  d <- repeated_counts()
  d$count[d$a == "1" & d$b == "1" & d$day == 7] <- 0
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank", time = "day")
  frame <- design_model_frame(d, des)
  boundary <- boundary_cells(frame, "count", "plain")
  expect_length(boundary, 1L)
  expect_length(effective_boundary(boundary, frame, "0 + .group + .time", "log"), 0L)
  expect_length(effective_boundary(boundary, frame, "0 + .cell_time", "log"), 1L)
  expect_length(effective_boundary(boundary, frame, "0 + .cell_time", "identity"), 0L)
  d$count[d$a == "1" & d$b == "1"] <- 0
  frame <- design_model_frame(d, des)
  boundary <- boundary_cells(frame, "count", "plain")
  expect_length(effective_boundary(boundary, frame, "0 + .group + .time", "log"), 3L)
})

test_that("binomial responses need whole numbers of successes and trials", {
  d <- survival_jars()
  des <- design_2x2(d)
  fit <- fit_responses(d, "alive", des, trials = c(alive = "n"), residual_checks = FALSE,
                       quiet = TRUE)
  f <- fit$fits$alive
  expect_equal(f$status, "ok")
  expect_true(f$family %in% c("binomial", "betabinomial"))
  expect_equal(f$type_source, "trials")
  expect_match(f$formula, "cbind(`alive`, .failures)", fixed = TRUE)
  p <- d
  p$alive <- p$alive / p$n
  expect_error(fit_responses(p, "alive", design_2x2(p), trials = c(alive = "n"),
                             residual_checks = FALSE), "number of successes",
               class = "must_error_input")
  expect_error(fit_responses(p, "alive", design_2x2(p),
                             families = list(binomial = stats::binomial()),
                             residual_checks = FALSE), "trials",
               class = "must_error_input")
  expect_error(fit_responses(d, "alive", des, trials = c(alive = "n"),
                             families = list(poisson = stats::poisson()),
                             residual_checks = FALSE), "binomial",
               class = "must_error_input")
  expect_error(fit_responses(d, "alive", des, trials = c(other = "n"),
                             residual_checks = FALSE), class = "must_error_input")
})

test_that("an offset turns counts into rates", {
  d <- two_by_two_counts()
  d$volume <- rep(c(1, 2), 10)
  d$count <- d$count * d$volume
  fit <- fit_responses(d, "count", design_2x2(d), offset = c(count = "volume"),
                       families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  f <- fit$fits$count
  expect_match(f$formula, "offset(log(.effort))", fixed = TRUE)
  beta <- glmmTMB::fixef(f$model)$cond
  beta <- beta[paste0(".group", c("a=0 | b=0", "a=1 | b=0", "a=0 | b=1", "a=1 | b=1"))]
  rate <- tapply(d$count, paste(d$a, d$b), sum) / tapply(d$volume, paste(d$a, d$b), sum)
  expect_equal(unname(exp(beta)), as.numeric(rate[c("0 0", "1 0", "0 1", "1 1")]), tolerance = 1e-5)
  expect_error(fit_responses(d, "count", design_2x2(d), offset = c(count = "volume"),
                             families = list(gaussian = stats::gaussian()),
                             residual_checks = FALSE), "log", class = "must_error_input")
})

test_that("equal AIC between links is broken towards the log link and reported", {
  d <- two_by_two_gaussian()
  fams <- list(gaussian = stats::gaussian(), gaussian_log = stats::gaussian(link = "log"))
  fit <- fit_responses(d, "y", design_2x2(d), families = fams, residual_checks = FALSE,
                       quiet = TRUE)
  expect_equal(fit$fits$y$family, "gaussian_log")
  expect_match(fit$fits$y$reason, "fit equally well")
  z <- two_by_two_gaussian()
  z$y[z$a == "1" & z$b == "1"] <- 0
  fit_z <- fit_responses(z, "y", design_2x2(z), families = fams, residual_checks = FALSE,
                         quiet = TRUE)
  expect_equal(fit_z$fits$y$family, "gaussian")
  expect_match(fit_z$fits$y$reason, "only zeros")
})

test_that("the dispersion can vary between `by` groups or cells", {
  d <- rbind(transform(two_by_two_gaussian(), site = "A"),
             transform(two_by_two_gaussian(dev = c(-3, -1.5, 0, 1.5, 3)), site = "B"))
  d$id <- seq_len(nrow(d))
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id",
                         by = "site")
  fams <- list(gaussian = stats::gaussian())
  common <- fit_responses(d, "y", des, families = fams, residual_checks = FALSE, quiet = TRUE)
  by_site <- fit_responses(d, "y", des, families = fams, dispersion = "by",
                           residual_checks = FALSE, quiet = TRUE)
  expect_equal(by_site$settings$dispersion, "by")
  expect_lt(by_site$fits$y$selection$aic, common$fits$y$selection$aic)
  est <- estimate_interactions(by_site, null_models = "AD", inference = "delta")
  ad <- est$estimates[est$estimates$component == "interaction", ]
  # One variance per site, estimated from its 20 observations and 4 means
  expect_equal(ad$df, c(16, 16))
  expect_equal(ad$se[ad$site == "B"] / ad$se[ad$site == "A"], 3, tolerance = 1e-4)
  expect_error(fit_responses(two_by_two_gaussian(), "y", design_2x2(two_by_two_gaussian()),
                             dispersion = "by", residual_checks = FALSE),
               "by", class = "must_error_input")
})

test_that("every value of a response with an offset needs an effort", {
  d <- two_by_two_counts()
  d$volume <- 1
  d$volume[3] <- NA
  expect_error(fit_responses(d, "count", design_2x2(d), offset = c(count = "volume"),
                             residual_checks = FALSE), "sampling effort",
               class = "must_error_input")
})

test_that("one dispersion per cell skips only the candidates that have one", {
  d <- two_by_two_counts()
  d <- d[-c(17:20), ]
  des <- design_2x2(d)
  pois <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                        dispersion = "cell", residual_checks = FALSE, quiet = TRUE)
  expect_equal(pois$fits$count$status, "ok")
  both <- fit_responses(d, "count", des, dispersion = "cell", residual_checks = FALSE,
                        quiet = TRUE)
  sel <- both$fits$count$selection
  expect_equal(both$fits$count$family, "poisson")
  expect_match(sel$message[sel$family == "nbinom2"], "two or more units")
})

test_that("the seed is drawn when not given and stored", {
  d <- two_by_two_gaussian()
  fit <- fit_responses(d, "y", design_2x2(d), families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  expect_true(is.integer(fit$settings$seed))
  expect_output(print(fit), "seed")
  expect_error(fit_responses(d, "y", design_2x2(d), seed = 1.5, residual_checks = FALSE),
               class = "must_error_input")
})

test_that("the REML refit is kept for the small-sample correction", {
  d <- two_by_two_gaussian()
  fit <- fit_responses(d, "y", design_2x2(d), families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  f <- fit$fits$y
  expect_false(is.null(f$model_reml))
  # REML variance: RSS / (N - 4)
  expect_equal(stats::sigma(f$model_reml)^2, 4 * 2.5 / 16, tolerance = 1e-5)
  expect_equal(stats::sigma(f$model)^2, 4 * 2.5 / 20, tolerance = 1e-5)
})

test_that("invalid arguments give classed errors", {
  d <- two_by_two_gaussian()
  des <- design_2x2(d)
  expect_error(fit_responses(d, "zz", des, residual_checks = FALSE), class = "must_error_input")
  expect_error(fit_responses(d, "a", des, residual_checks = FALSE), class = "must_error_input")
  expect_error(fit_responses(d, "y", list(), residual_checks = FALSE), class = "must_error_input")
  expect_error(fit_responses(d, "y", des, time_structure = "add", residual_checks = FALSE),
               class = "must_error_input")
  expect_error(fit_responses(d, "y", des, families = list(gaussian()), residual_checks = FALSE),
               class = "must_error_input")
  other <- d[d$a == "0", ]
  expect_error(fit_responses(other, "y", des, residual_checks = FALSE),
               "describe_design", class = "must_error_input")
})

test_that("print methods summarise fits and interactions", {
  d <- two_by_two_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  expect_output(print(fit), "count: poisson \\(log link; user, families\\)")
  expect_output(print(fit), "Residual checks: not run")
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
  expect_output(print(est), "AD = differences")
})

test_that("failed responses are listed by print()", {
  d <- two_by_two_counts()
  d$bad <- NA_real_
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, c("count", "bad"), des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  expect_equal(fit$fits$bad$status, "failed")
  expect_output(print(fit), "bad: FAILED")
  est <- estimate_interactions(fit, null_models = "MU", inference = "delta")
  expect_equal(est$diagnostics$failed_responses$response, "bad")
  expect_output(print(est), "Responses without a fitted model: bad")
})

test_that("the optimiser limits grow with the number of fixed effects", {
  expect_equal(optimizer_control(10)$optCtrl, list(iter.max = 300L, eval.max = 400L))
  expect_equal(optimizer_control(200)$optCtrl, list(iter.max = 4000L, eval.max = 6000L))
  # 153 cell means: more than glmmTMB's default 300 iterations are needed
  d <- many_single_cells()
  s <- paste0("s", 1:12)
  des <- describe_design(d, s, control = stats::setNames(rep("0", 12), s), unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  expect_equal(fit$fits$y$status, "ok")
  expect_equal(fit$settings$control, "scaled iteration limits")
  frame <- design_model_frame(d, des)
  means <- tapply(frame$y, frame$.group, mean)
  beta <- glmmTMB::fixef(fit$fits$y$model)$cond
  expect_equal(unname(beta[paste0(".group", names(means))]), unname(as.numeric(means)),
               tolerance = 1e-5)
})

test_that("a user `control` is passed to glmmTMB", {
  d <- two_by_two_gaussian()
  ctrl <- glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1, eval.max = 1))
  fit <- fit_responses(d, "y", design_2x2(d), families = list(gaussian = stats::gaussian()),
                       control = ctrl, residual_checks = FALSE, quiet = TRUE)
  expect_false(fit$fits$y$selection$converged)
  expect_equal(fit$settings$control, "user")
  expect_error(fit_responses(d, "y", design_2x2(d), control = list(a = 1),
                             residual_checks = FALSE),
               "glmmTMBControl", class = "must_error_input")
})

test_that("groups get a random intercept only when it can be estimated", {
  spec <- list(stressors = c("temperature", "nutrients"),
               control = c(temperature = "15", nutrients = "no"))
  plots <- jars_in_rooms(rooms_per_level = 2L)
  des <- describe_design(plots, spec$stressors, control = spec$control, unit = "id",
                         groups = "room")
  fit <- fit_responses(plots, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  expect_equal(fit$fits$y$formula, "`y` ~ 0 + .group + (1 | `.g_room`)")
  expect_equal(fit$fits$y$status, "ok")
  single <- jars_in_rooms(rooms_per_level = 1L)
  des1 <- describe_design(single, spec$stressors, control = spec$control, unit = "id",
                          groups = "room")
  fit1 <- fit_responses(single, "y", des1, families = list(gaussian = stats::gaussian()),
                        residual_checks = FALSE, quiet = TRUE)
  expect_equal(fit1$fits$y$formula, "`y` ~ 0 + .group")
})

test_that("an autocorrelation structure can replace the random intercept", {
  d <- autocorrelated_gaussian()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank", time = "day")
  fams <- list(gaussian = stats::gaussian())
  ar1 <- fit_responses(d, "y", des, families = fams, random = "ar1(.time + 0 | .unit)",
                       residual_checks = FALSE, quiet = TRUE)
  intercept <- fit_responses(d, "y", des, families = fams, residual_checks = FALSE,
                             quiet = TRUE)
  expect_equal(ar1$fits$y$formula, "`y` ~ 0 + .cell_time + ar1(.time + 0 | .unit)")
  expect_equal(ar1$fits$y$status, "ok")
  # The autocorrelated noise is better described by the AR1 structure
  expect_lt(ar1$fits$y$selection$aic, intercept$fits$y$selection$aic)
  est <- estimate_interactions(ar1, null_models = "AD", inference = "delta")
  # Six days and their average
  expect_equal(nrow(est$estimates[est$estimates$component == "interaction", ]), 7L)
})

test_that("an automatic random intercept with variance at zero is removed", {
  tanks <- expand.grid(a = c("0", "1"), b = c("0", "1"), rep = 1:4, stringsAsFactors = FALSE)
  tanks$tank <- paste0("T", seq_len(nrow(tanks)))
  d <- tanks[rep(seq_len(nrow(tanks)), times = 2), ]
  d$day <- rep(c(1, 2), each = nrow(tanks))
  # Deviations of opposite sign on the two days: no variation between tanks
  e <- rep(c(-0.3, 0.1, 0.4, -0.2), each = 4)[match(d$tank, tanks$tank)]
  d$y <- 10 - 2 * (d$a == "1") - (d$b == "1") + ifelse(d$day == 1, e, -e)
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank", time = "day")
  fams <- list(gaussian = stats::gaussian())
  fit <- fit_responses(d, "y", des, families = fams, residual_checks = FALSE, quiet = TRUE)
  f <- fit$fits$y
  expect_equal(f$status, "ok")
  expect_equal(f$selection$random, "none")
  expect_equal(f$formula, "`y` ~ 0 + .cell_time")
  expect_match(f$reason, "variance estimated at zero")
  expect_equal(fit$settings$random, "(1 | .unit)")
  # Terms given by the user are not removed
  given <- fit_responses(d, "y", des, families = fams, random = "(1 | .unit)",
                         residual_checks = FALSE, quiet = TRUE)
  expect_equal(given$fits$y$selection$random, "(1 | .unit)")
})

test_that("residual checks run with DHARMa when it is installed", {
  skip_if_not_installed("DHARMa")
  d <- two_by_two_counts()
  fit <- fit_responses(d, "count", design_2x2(d), seed = 1, n_sim_residuals = 250, quiet = TRUE)
  sel <- fit$fits$count$selection
  expect_setequal(sel$family, c("poisson", "nbinom2"))
  p <- unlist(sel[sel$converged, c("p_uniformity", "p_dispersion", "p_zero_inflation")])
  expect_true(all(is.na(p) | (p >= 0 & p <= 1)))
  expect_true(all(sel$n_tests_passed[sel$converged] %in% 0:3))
  expect_false(is.na(fit$fits$count$checks_passed))
  again <- fit_responses(d, "count", design_2x2(d), seed = 1, n_sim_residuals = 250, quiet = TRUE)
  expect_identical(again$fits$count$selection, sel)
})

test_that("the default Poisson for counts is selected only if the negative binomial fails", {
  sel <- data.frame(family = c("poisson", "nbinom2"), link = "log", converged = TRUE,
                    aic = c(100, 101), residuals_ok = NA, n_tests_passed = NA_integer_,
                    stringsAsFactors = FALSE)
  # Lower AIC, but only a fallback
  expect_equal(select_candidate(sel, FALSE, FALSE, fallback = "poisson")$best, 2L)
  # Without the rule (candidates given in `families`), the AIC decides
  expect_equal(select_candidate(sel, FALSE, FALSE)$best, 1L)
  # The negative binomial did not converge: the Poisson, with a note
  failed <- sel
  failed$converged[2] <- FALSE
  failed$aic[2] <- NA
  chosen <- select_candidate(failed, FALSE, FALSE, fallback = "poisson")
  expect_equal(chosen$best, 1L)
  expect_match(chosen$notes, "'nbinom2' did not converge")
  # A negative binomial that fails the residual checks still comes first
  checked <- sel
  checked$residuals_ok <- c(TRUE, FALSE)
  checked$n_tests_passed <- c(3L, 2L)
  chosen <- select_candidate(checked, TRUE, FALSE, fallback = "poisson")
  expect_equal(chosen$best, 2L)
  expect_false(chosen$checks_passed)
})

test_that("fit_responses() says when it falls back to the Poisson", {
  d <- two_by_two_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  f <- fit_responses(d, "count", des, residual_checks = FALSE, quiet = TRUE)$fits$count
  # These counts are less variable than a Poisson: the negative binomial
  # cannot estimate its overdispersion
  if (identical(f$family, "poisson")) {
    expect_match(f$reason, "'nbinom2' did not converge")
  } else {
    expect_equal(f$family, "nbinom2")
  }
  # With the candidates given, no fallback rule
  user <- fit_responses(d, "count", des, residual_checks = FALSE, quiet = TRUE,
                        families = list(poisson = stats::poisson(), nbinom2 = glmmTMB::nbinom2()))
  expect_false(grepl("did not converge", user$fits$count$reason))
})
