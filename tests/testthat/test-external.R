# Models fitted outside must, through emmeans.

design_2x2 <- function(d, ...) {
  describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), ...)
}

interactions_of <- function(est, null_model) {
  e <- est$estimates
  e[e$component == "interaction" & e$null_model == null_model, , drop = FALSE]
}

test_that("lm(): the AD interaction is the interaction coefficient, with its t test", {
  skip_if_not_installed("emmeans")
  d <- two_by_two_gaussian(dev = c(-1.3, -0.4, 0, 0.6, 1.1))
  model <- stats::lm(y ~ a * b, data = d)
  est <- estimate_interactions(model, design_2x2(d), null_models = c("AD", "MU"),
                               inference = "delta")
  ad <- interactions_of(est, "AD")
  coefs <- summary(model)$coefficients
  expect_equal(ad$estimate, coefs["a1:b1", "Estimate"])
  expect_equal(ad$se, coefs["a1:b1", "Std. Error"])
  expect_equal(ad$df, 16)
  expect_equal(ad$p_value, coefs["a1:b1", "Pr(>|t|)"])
  expect_equal(c(ad$lower, ad$upper), unname(stats::confint(model)["a1:b1", ]))
  expect_equal(ad$response, "y")
  single <- est$estimates[est$estimates$component == "single" & est$estimates$null_model == "AD", ]
  expect_equal(single$estimate, unname(coefs[c("a1", "b1"), "Estimate"]))
  mu <- interactions_of(est, "MU")
  expect_equal(mu$estimate, log(3) - log(8) - log(7) + log(10), tolerance = 1e-6)
  expect_equal(mu$expected, 8 * 7 / 10, tolerance = 1e-6)
  expect_equal(est$settings$model, "lm")
  expect_equal(unname(est$settings$families), "gaussian (identity)")
  expect_equal(nrow(est$diagnostics$means), 4L)
  expect_equal(est$diagnostics$means$mean, c(10, 8, 7, 3), tolerance = 1e-6)
  none <- estimate_interactions(model, design_2x2(d), null_models = "AD", inference = "delta",
                                correction = "none")
  expect_equal(interactions_of(none, "AD")$df, Inf)
  fixed <- estimate_interactions(model, design_2x2(d), null_models = "AD", inference = "delta",
                                 df = 4)
  expect_equal(interactions_of(fixed, "AD")$df, 4)
  expect_equal(interactions_of(fixed, "AD")$p_value,
               2 * stats::pt(-abs(ad$estimate / ad$se), 4))
})

test_that("glm() with an estimated dispersion uses its residual degrees of freedom", {
  skip_if_not_installed("emmeans")
  d <- two_by_two_counts()
  des <- design_2x2(d, unit = "id")
  qp <- stats::glm(count ~ a * b, family = stats::quasipoisson(), data = d)
  est <- estimate_interactions(qp, des, null_models = "MU", inference = "delta")
  mu <- interactions_of(est, "MU")
  coefs <- summary(qp)$coefficients
  # The t test of summary.glm(), which emmeans would replace by a z test
  expect_equal(mu$df, qp$df.residual)
  expect_equal(mu$estimate, coefs["a1:b1", "Estimate"], tolerance = 1e-8)
  expect_equal(mu$p_value, coefs["a1:b1", "Pr(>|t|)"], tolerance = 1e-8)
  expect_identical(est$settings$df, "residual")
  expect_output(print(est), "residual degrees of freedom of the glm")
  # emmeans' degrees of freedom are still reported for the means
  expect_equal(unique(est$diagnostics$means$df), Inf)
  # Poisson, as in summary.glm(): normal
  po <- stats::glm(count ~ a * b, family = stats::poisson(), data = d)
  est_po <- estimate_interactions(po, des, null_models = "MU", inference = "delta")
  expect_equal(interactions_of(est_po, "MU")$df, Inf)
  expect_identical(est_po$settings$df, "emmeans")
  # Without the correction, or with `df`, the residual ones are not used
  none <- estimate_interactions(qp, des, null_models = "MU", inference = "delta",
                                correction = "none")
  expect_equal(interactions_of(none, "MU")$df, Inf)
  expect_identical(none$settings$df, "emmeans")
  fixed <- estimate_interactions(qp, des, null_models = "MU", inference = "delta", df = 4)
  expect_equal(interactions_of(fixed, "MU")$df, 4)
})

test_that("glm(): the MU interaction of a Poisson model is its interaction coefficient", {
  skip_if_not_installed("emmeans")
  d <- two_by_two_counts()
  model <- stats::glm(count ~ a * b, family = stats::poisson(), data = d)
  est <- estimate_interactions(model, design_2x2(d), null_models = c("AD", "MU", "BL"),
                               inference = "delta")
  mu <- interactions_of(est, "MU")
  coefs <- summary(model)$coefficients
  expect_equal(mu$estimate, coefs["a1:b1", "Estimate"], tolerance = 1e-8)
  expect_equal(mu$se, coefs["a1:b1", "Std. Error"], tolerance = 1e-6)
  expect_equal(mu$p_value, coefs["a1:b1", "Pr(>|z|)"], tolerance = 1e-6)
  expect_equal(mu$df, Inf)
  bl <- interactions_of(est, "BL")
  expect_equal(bl$estimate, 0.25 - 0.5 * 0.7, tolerance = 1e-8)
  expect_equal(bl$rescaled, -0.1 / 0.35, tolerance = 1e-6)
  # The same model fitted by fit_responses() gives the same estimates
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  ref <- estimate_interactions(fit, null_models = c("AD", "MU", "BL"), inference = "delta")
  expect_equal(est$estimates$estimate, ref$estimates$estimate, tolerance = 1e-5)
  expect_equal(est$estimates$se, ref$estimates$se, tolerance = 1e-4)
  sim <- estimate_interactions(model, design_2x2(d), null_models = c("MU", "DO"),
                               n_sim = 20000, seed = 4)
  expect_equal(interactions_of(sim, "MU")$se, mu$se, tolerance = 0.03)
  expect_equal(interactions_of(sim, "DO")$status, "ok")
  expect_output(print(est), "Model: glm, poisson \\(log\\); means from emmeans")
  expect_output(print(est), "degrees of freedom from emmeans")
  expect_output(print(summary(est)), "degrees of freedom from emmeans")
  fixed <- estimate_interactions(model, design_2x2(d), null_models = "MU", inference = "delta",
                                 df = 15)
  expect_output(print(fixed), "t with 15 degrees of freedom")
})

test_that("a quantity fixed at zero by the model is not tested", {
  skip_if_not_installed("emmeans")
  d <- two_by_two_counts()
  model <- stats::glm(count ~ a + b, family = stats::poisson(), data = d)
  est <- estimate_interactions(model, design_2x2(d), null_models = c("AD", "MU", "BL"),
                               inference = "delta")
  int <- est$estimates[est$estimates$component == "interaction", ]
  expect_equal(int$status[int$null_model %in% c("MU", "BL")], c("not_applicable", "not_applicable"))
  expect_match(int$reason[int$null_model == "MU"], "fixes this quantity at zero")
  # Simple addition is not constrained by a log-link model without interaction
  expect_equal(int$status[int$null_model == "AD"], "ok")
  expect_true(all(est$estimates$status[est$estimates$component != "interaction"] == "ok"))
})

test_that("means at the boundary of the link or not estimable are named", {
  skip_if_not_installed("emmeans")
  d <- two_by_two_counts()
  d$count[d$a == "1" & d$b == "1"] <- 0
  model <- suppressWarnings(stats::glm(count ~ a * b, family = stats::poisson(), data = d))
  est <- estimate_interactions(model, design_2x2(d), null_models = c("AD", "MU"),
                               inference = "delta")
  e <- est$estimates
  expect_true(all(e$status[e$component %in% c("combined", "interaction")] == "not_estimable"))
  expect_match(e$reason[e$component == "interaction"][1], "'a=1 \\+ b=1' is at the boundary")
  expect_true(all(e$status[e$component == "single"] == "ok"))
  expect_equal(est$diagnostics$means$problem, c("", "", "", "boundary"))
  # A combination missing at one time in a model with all interactions
  r <- repeated_counts()
  r <- r[!(r$a == "1" & r$b == "1" & r$day == 21), ]
  des <- describe_design(r, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank",
                         time = "day")
  model <- stats::glm(count ~ a * b * factor(day), family = stats::poisson(), data = r)
  est <- estimate_interactions(model, des, null_models = "MU", inference = "delta")
  int <- interactions_of(est, "MU")
  expect_equal(int$time, c("7", "14", "21", "overall"))
  expect_equal(int$status, c("ok", "ok", "not_estimable", "not_estimable"))
  expect_match(int$reason[3], "cannot estimate the mean of 'a=1 \\+ b=1'")
})

test_that("numeric stressors are evaluated at each observed value", {
  skip_if_not_installed("emmeans")
  d <- expand.grid(temp = c(15, 20, 25), nit = c(0, 2), rep = 1:3)
  d$y <- round(exp(3 - 0.05 * (d$temp - 15) - 0.2 * d$nit - 0.02 * (d$temp - 15) * d$nit) *
                 c(0.9, 1, 1.1)[d$rep])
  des <- describe_design(d, c("temp", "nit"), control = c(temp = "15", nit = "0"))
  model <- stats::glm(y ~ temp * nit, family = stats::poisson(), data = d)
  est <- estimate_interactions(model, des, null_models = "MU", inference = "delta")
  int <- interactions_of(est, "MU")
  expect_equal(int$levels, c("temp=20 | nit=2", "temp=25 | nit=2"))
  expect_equal(int$estimate, unname(stats::coef(model)["temp:nit"]) * c(5, 10) * 2,
               tolerance = 1e-8)
})

test_that("sampling times, baselines and averages over times follow the design", {
  skip_if_not_installed("emmeans")
  r <- repeated_counts()
  r <- rbind(transform(r[r$day == 7, ], day = 0), r)
  des <- describe_design(r, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank",
                         time = "day", baseline = 0)
  model <- stats::glm(count ~ a * b * factor(day), family = stats::poisson(), data = r)
  est <- estimate_interactions(model, des, null_models = "MU", inference = "delta")
  int <- interactions_of(est, "MU")
  expect_equal(int$time, c("7", "14", "21", "overall"))
  by_time <- vapply(c("7", "14", "21"), function(tt) {
    sub <- r[r$day == as.numeric(tt), ]
    unname(stats::coef(stats::glm(count ~ a * b, family = stats::poisson(), data = sub))["a1:b1"])
  }, numeric(1))
  expect_equal(int$estimate[1:3], unname(by_time), tolerance = 1e-6)
  expect_equal(int$estimate[4], mean(by_time), tolerance = 1e-6)
  # Time as a number in the model: the same means
  model2 <- stats::glm(count ~ a * b * poly(day, 3), family = stats::poisson(), data = r)
  est2 <- estimate_interactions(model2, des, null_models = "MU", inference = "delta")
  expect_equal(interactions_of(est2, "MU")$estimate, int$estimate, tolerance = 1e-6)
})

test_that("`by` groups and other predictors are handled as emmeans does", {
  skip_if_not_installed("emmeans")
  d <- expand.grid(a = c("0", "1"), b = c("0", "1"), sp = c("x", "y"), rep = 1:4,
                   stringsAsFactors = FALSE)
  d$y <- round(ifelse(d$sp == "x", 20, 40) * ifelse(d$a == "1", 0.5, 1) *
                 ifelse(d$b == "1", 0.6, 1) * c(0.8, 0.9, 1.1, 1.2)[d$rep])
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), by = "sp")
  model <- stats::glm(y ~ a * b * sp, family = stats::poisson(), data = d)
  est <- estimate_interactions(model, des, null_models = "MU", inference = "delta")
  int <- interactions_of(est, "MU")
  expect_equal(int$sp, c("x", "y"))
  cf <- stats::coef(model)
  expect_equal(int$estimate, unname(c(cf["a1:b1"], cf["a1:b1"] + cf["a1:b1:spy"])),
               tolerance = 1e-8)
  # Without `by` in the design, the model averages over sp
  x <- d[d$sp == "x", ]
  x$block <- rep(c("B1", "B2"), each = 8)
  blocked <- stats::glm(y ~ a * b + block, family = stats::poisson(), data = x)
  est2 <- estimate_interactions(blocked, design_2x2(x), null_models = "MU",
                                inference = "delta")
  expect_equal(est2$settings$averaged_over, "block")
  expect_equal(interactions_of(est2, "MU")$estimate,
               unname(stats::coef(blocked)["a1:b1"]), tolerance = 1e-8)
})

test_that("proportions are analysed as their complement when higher means sensitivity", {
  skip_if_not_installed("emmeans")
  d <- expand.grid(a = c("0", "1"), b = c("0", "1"), rep = 1:4, stringsAsFactors = FALSE)
  dead <- c("00" = 3, "10" = 7, "01" = 5, "11" = 11)
  d$dead <- unname(dead[paste0(d$a, d$b)]) + c(-1, 0, 0, 1)[d$rep]
  d$alive <- 20 - d$dead
  model <- stats::glm(cbind(dead, alive) ~ a * b, family = stats::binomial(), data = d)
  est <- estimate_interactions(model, design_2x2(d), null_models = c("MU", "BL"),
                               inference = "delta", sensitivity = TRUE)
  p <- tapply(d$dead, list(d$a, d$b), sum) / tapply(d$dead + d$alive, list(d$a, d$b), sum)
  s <- 1 - p
  mu <- interactions_of(est, "MU")
  expect_equal(mu$estimate, log(s["1", "1"]) - log(s["1", "0"]) - log(s["0", "1"]) +
                 log(s["0", "0"]), tolerance = 1e-6)
  expect_equal(mu$note, "complement")
  expect_equal(est$settings$sensitivity, "cbind(dead, alive)")
  # A response that is not a proportion: MU is not applicable
  linear <- stats::lm(dead ~ a * b, data = d)
  est2 <- estimate_interactions(linear, design_2x2(d), null_models = c("AD", "MU"),
                                inference = "delta", sensitivity = TRUE)
  expect_equal(interactions_of(est2, "MU")$status, "not_applicable")
  expect_equal(interactions_of(est2, "AD")$status, "ok")
})

test_that("three stressors: pure interactions are the three-way coefficients", {
  skip_if_not_installed("emmeans")
  d <- three_way_gaussian()
  des <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"))
  linear <- stats::lm(y ~ n * s * f, data = d)
  est <- estimate_interactions(linear, des, null_models = c("AD", "MU", "BL"),
                               inference = "delta")
  rows <- est$estimates[est$estimates$order == 3L & est$estimates$component == "interaction", ]
  coefs <- summary(linear)$coefficients
  ad <- rows[rows$null_model == "AD" & rows$estimand == "pure", ]
  expect_equal(ad$estimate, coefs["n1:s1:f1", "Estimate"])
  expect_equal(ad$se, coefs["n1:s1:f1", "Std. Error"])
  mu <- tapply(d$y, paste0(d$n, d$s, d$f), mean)
  w <- mu / mu[["000"]]
  e3 <- w[["111"]] - w[["100"]] * w[["011"]] - w[["010"]] * w[["101"]] -
    w[["001"]] * w[["110"]] + 2 * w[["100"]] * w[["010"]] * w[["001"]]
  expect_equal(rows$estimate[rows$null_model == "BL" & rows$estimand == "pure"], e3,
               tolerance = 1e-8)
  counts <- d
  counts$y <- round(counts$y)
  pois <- stats::glm(y ~ n * s * f, family = stats::poisson(), data = counts)
  est2 <- estimate_interactions(pois, des, null_models = "MU", inference = "delta")
  pure <- est2$estimates[est2$estimates$order == 3L & est2$estimates$estimand %in% "pure", ]
  expect_equal(pure$estimate, unname(stats::coef(pois)["n1:s1:f1"]), tolerance = 1e-8)
})

test_that("transformed responses are back-transformed as emmeans does", {
  skip_if_not_installed("emmeans")
  d <- two_by_two_counts()
  model <- stats::lm(log(count + 1) ~ a * b, data = d)
  est <- estimate_interactions(model, design_2x2(d), null_models = "AD", inference = "delta")
  means <- est$diagnostics$means
  expect_equal(means$mean, exp(means$link_estimate) - 1)
  expect_match(unname(est$settings$families), "response transformed with log\\(y \\+ 1\\)")
  expect_equal(est$estimates$response[[1L]], "log(count + 1)")
})

test_that("mixed models are analysed through emmeans", {
  skip_if_not_installed("emmeans")
  r <- repeated_counts()
  des <- describe_design(r, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank",
                         time = "day")
  tmb <- glmmTMB::glmmTMB(count ~ a * b * factor(day) + (1 | tank),
                          family = stats::poisson(), data = r)
  est <- estimate_interactions(tmb, des, null_models = "MU", inference = "delta")
  int <- interactions_of(est, "MU")
  cf <- glmmTMB::fixef(tmb)$cond
  expect_equal(int$estimate[1], unname(cf["a1:b1"]), tolerance = 1e-6)
  expect_equal(est$settings$model, "glmmTMB")
  skip_if_not_installed("lme4")
  glmer <- suppressMessages(lme4::glmer(count ~ a * b * factor(day) + (1 | tank),
                                        family = stats::poisson(), data = r))
  est2 <- estimate_interactions(glmer, des, null_models = "MU", inference = "delta")
  expect_equal(interactions_of(est2, "MU")$estimate[1],
               unname(lme4::fixef(glmer)["a1:b1"]), tolerance = 1e-6)
})

test_that("mgcv::gam(): parametric and smooth terms are analysed through emmeans", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("mgcv")
  d <- two_by_two_gaussian(dev = c(-1.3, -0.4, 0, 0.6, 1.1))
  d$temp <- rep(c(12, 15, 18, 14, 16), 4)
  model <- mgcv::gam(y ~ a * b + s(temp, k = 4), data = d)
  est <- estimate_interactions(model, design_2x2(d), inference = "delta")
  ad <- interactions_of(est, "AD")
  coefs <- summary(model)$p.table
  expect_equal(ad$estimate, coefs["a1:b1", "Estimate"], tolerance = 1e-8)
  expect_equal(ad$se, coefs["a1:b1", "Std. Error"], tolerance = 1e-6)
  # emmeans gives the residual degrees of freedom of the GAM
  expect_equal(ad$df, model$df.residual, tolerance = 1e-6)
  expect_equal(ad$p_value, coefs["a1:b1", "Pr(>|t|)"], tolerance = 1e-6)
  # A smooth over time for each treatment combination
  r <- repeated_counts()
  r$ab <- interaction(r$a, r$b)
  des <- describe_design(r, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank",
                         time = "day")
  smooth_time <- mgcv::gam(count ~ a * b + s(day, by = ab, k = 3), family = stats::poisson(),
                           data = r)
  est_time <- suppressMessages(estimate_interactions(smooth_time, des, null_models = "MU",
                                                     inference = "delta", overall = FALSE))
  mu <- interactions_of(est_time, "MU")
  grid <- expand.grid(a = c("0", "1"), b = c("0", "1"), day = c(7, 14, 21),
                      stringsAsFactors = FALSE)
  grid$ab <- interaction(grid$a, grid$b)
  eta <- stats::predict(smooth_time, newdata = grid, type = "link")
  # log yAB - log yA - log yB + log y0 at each day
  hand <- tapply(eta * ifelse(grid$a == grid$b, 1, -1), grid$day, sum)
  expect_equal(mu$time, c("7", "14", "21"))
  expect_equal(mu$estimate, as.numeric(hand), tolerance = 1e-6)
})

test_that("models fitted outside must are checked against the design", {
  skip_if_not_installed("emmeans")
  d <- two_by_two_counts()
  des <- design_2x2(d)
  model <- stats::glm(count ~ a * b, family = stats::poisson(), data = d)
  expect_error(estimate_interactions(model), "describe_design", class = "must_error_input")
  expect_error(estimate_interactions(list(1)), "not the output of fit_responses",
               class = "must_error_input")
  expect_error(estimate_interactions(stats::glm(count ~ a, family = stats::poisson(), data = d),
                                     des),
               "'b' of `design` is not a predictor", class = "must_error_input")
  other <- d
  other$b[other$b == "1"] <- "2"
  expect_error(estimate_interactions(model, design_2x2(other)), "no level '2' of 'b'",
               class = "must_error_input")
  expect_error(estimate_interactions(model, des, at = list(a = "1")), "`at` cannot set 'a'",
               class = "must_error_input")
  expect_error(estimate_interactions(model, des, cov.keep = "a"), "`cov.keep`",
               class = "must_error_input")
  expect_error(estimate_interactions(model, des, regrid = "response"), "`regrid`",
               class = "must_error_input")
  expect_error(estimate_interactions(model, des, df = 3, correction = "none"),
               class = "must_error_input")
  expect_error(estimate_interactions(model, des, df = -1), class = "must_error_input")
  expect_error(estimate_interactions(model, des, sensitivity = "count"),
               class = "must_error_input")
  # A stressor nested in another factor gives several means per combination
  nested <- d
  nested$block <- ifelse(nested$a == "1", "B2", "B1")
  nested_model <- stats::glm(count ~ a * b + block, family = stats::poisson(), data = nested)
  expect_error(suppressMessages(estimate_interactions(nested_model, design_2x2(nested))),
               "several means for the same treatment combination", class = "must_error_input")
  err <- expect_error(estimate_interactions(model, des, null_models = "XX"),
                      class = "must_error_input")
  expect_identical(conditionCall(err)[[1L]], as.name("estimate_interactions"))
})

test_that("emmeans is required for models fitted outside must", {
  d <- two_by_two_counts()
  model <- stats::glm(count ~ a * b, family = stats::poisson(), data = d)
  local_mocked_bindings(is_installed = function(pkg) FALSE)
  expect_error(estimate_interactions(model, design_2x2(d)), "emmeans",
               class = "must_error_missing_package")
})

test_that("links and transformations of emmeans are inverted", {
  expect_equal(external_link(list(tran = "log"), NULL)$linkinv(log(5)), 5)
  expect_equal(external_link(list(tran = "logit"), NULL)$linkinv(0), 0.5)
  expect_equal(external_link(list(), NULL)$name, "identity")
  shifted <- external_link(list(tran = "log", tran.offset = 1), NULL)
  expect_equal(shifted$linkinv(log(6)), 5)
  expect_equal(shifted$mu.eta(log(6)), 6)
  expect_equal(external_link(list(tran = "log10"), NULL)$linkinv(2), 100)
  expect_equal(external_link(list(tran = "sqrt"), NULL)$linkinv(c(-1, 3)), c(0, 9))
  expect_error(external_link(list(tran = "unknown"), NULL), class = "must_error_input")
  expect_error(external_link(list(tran = "log", tran2 = "logit"), NULL),
               class = "must_error_input")
  expect_equal(external_problems(c(1, -20, 25, NA, 0), c(0.1, 500, 500, 1, NA), "logit"),
               c("", "boundary", "boundary", "nonestimable", "no_se"))
  expect_equal(external_problems(c(25, -20), c(500, 500), "log"), c("", "boundary"))
  expect_equal(external_problems(-20, 500, "identity"), "")
  expect_equal(match_numeric(c(0.1, 0.3, 2), c(2, 0.1)), c(2L, NA, 1L))
})

test_that("concentration addition and stress addition work with models fitted outside must", {
  skip_if_not_installed("emmeans")
  d <- ca_factorial()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"))
  model <- stats::lm(y ~ factor(a) * factor(b), data = d)
  est <- estimate_interactions(model, des, null_models = c("SAM", "CA"), inference = "delta")
  ca <- est$estimates[est$estimates$component == "interaction" & est$estimates$null_model == "CA", ]
  expect_true(all(abs(ca$estimate) < 1e-6))
  expect_true(all(ca$status == "ok"))
  sam <- est$estimates[est$estimates$component == "interaction" &
                         est$estimates$null_model == "SAM", ]
  expect_true(all(sam$status == "ok"))
})

test_that("lm(): a control shared by two `by` groups enters the sets of both", {
  skip_if_not_installed("emmeans")
  d <- shared_regimes()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id",
                         by = "regime", shared = c(regime = "both"))
  model <- stats::lm(y ~ regime * a * b, data = d)
  ad <- interactions_of(estimate_interactions(model, des, data = d, null_models = "AD",
                                              inference = "delta"), "AD")
  mu <- tapply(d$y, paste(d$regime, d$a, d$b), mean)
  expected <- vapply(ad$regime, function(r) {
    mu[[paste(r, "1 1")]] - mu[[paste(r, "1 0")]] - mu[[paste(r, "0 1")]] + mu[["both 0 0"]]
  }, numeric(1))
  expect_equal(ad$estimate, unname(expected), tolerance = 1e-8)
  expect_equal(ad$df, c(21, 21))
  expect_true(all(ad$status == "ok"))
})
