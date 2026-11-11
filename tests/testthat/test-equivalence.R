fit_2x2_eq <- function(d, response, families) {
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit_responses(d, response, des, families = families, residual_checks = FALSE, quiet = TRUE)
}

net_rows <- function(est, null_model) {
  e <- est$estimates
  e[e$component == "interaction" & e$null_model == null_model & e$estimand %in% "net", ]
}

# Two one-sided tests from the signed distances to the lower and upper bounds.
tost_from_distances <- function(lower, upper, df) {
  max(stats::pt(lower, df, lower.tail = FALSE), stats::pt(upper, df))
}

# Signed distance, on the scale `h` (variances `v` of the cell means on that
# scale, independent cells), from the estimated means `y` (control, A, B, AB)
# to the means with multiplicative relative deviation y_AB / y_0 - y_A y_B /
# y_0^2 equal to `bound`: y_AB is eliminated and the other three means are
# optimised, independently of the package's search.
mu_distance_by_elimination <- function(y, v, bound, h = identity, h_inv = identity) {
  z0 <- h(y)
  obj <- function(u) {
    m <- h_inv(u)
    yk <- m[2] * m[3] / m[1] + bound * m[1]
    if (!is.finite(yk) || yk < 0) return(1e10)
    sum((c(u, h(yk)) - z0)^2 / v)
  }
  o <- stats::optim(z0[1:3], obj, method = "BFGS", control = list(reltol = 1e-15, maxit = 1000))
  o <- stats::optim(o$par, obj, method = "Nelder-Mead", control = list(reltol = 1e-15, maxit = 5000))
  d0 <- y[4] / y[1] - y[2] * y[3] / y[1]^2
  sign(d0 - bound) * sqrt(o$value)
}

test_that("relative deviations and equivalence tests match hand calculations", {
  d <- two_by_two_gaussian()
  fit <- fit_2x2_eq(d, "y", list(gaussian = stats::gaussian()))
  est <- estimate_interactions(fit, null_models = c("AD", "MU", "BL"), inference = "delta",
                               margin = 0.1)
  expect_equal(est$settings$margin, 0.1)
  # REML variance RSS / (N - p) = 10 / 16; variance of each cell mean, over 5
  v <- 10 / 16 / 5
  y <- c(10, 8, 7, 3)
  crit <- stats::qt(0.975, 16)

  # AD: relative deviation = (yAB - yA - yB + y0) / y0, Wald interval
  ad <- net_rows(est, "AD")
  theta <- y[4] - y[2] - y[3] + y[1]
  rel <- theta / y[1]
  g <- c(1 / y[1] - theta / y[1]^2, -1 / y[1], -1 / y[1], 1 / y[1])
  se <- sqrt(sum(g^2) * v)
  expect_equal(ad$relative, rel, tolerance = 1e-8)
  expect_equal(c(ad$relative_lower, ad$relative_upper), rel + c(-1, 1) * crit * se,
               tolerance = 1e-6)
  # Equivalence: the bounds are linear in the means, so the tests are the
  # Fieller-type t tests of yAB - yA - yB + (1 +/- margin) y0
  fieller <- function(m) {
    cv <- c(1 + m, -1, -1, 1)
    sum(cv * y) / sqrt(sum(cv^2) * v)
  }
  expect_equal(ad$p_equivalence, tost_from_distances(fieller(0.1), fieller(-0.1), 16),
               tolerance = 1e-6)

  # MU and BL: wAB - wA wB, the same relative deviation and the same test
  bl <- net_rows(est, "BL")
  mu <- net_rows(est, "MU")
  rel <- y[4] / y[1] - y[2] * y[3] / y[1]^2
  expect_equal(bl$relative, bl$estimate)
  expect_equal(mu$relative, rel, tolerance = 1e-8)
  expect_equal(bl$relative, rel, tolerance = 1e-8)
  lower <- mu_distance_by_elimination(y, rep(v, 4), -0.1)
  upper <- mu_distance_by_elimination(y, rep(v, 4), 0.1)
  expect_equal(mu$p_equivalence, tost_from_distances(lower, upper, 16), tolerance = 1e-5)
  expect_equal(bl$p_equivalence, mu$p_equivalence, tolerance = 1e-8)

  # Without a margin there is no equivalence test, but the relative deviation is given
  plain <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  expect_null(plain$settings$margin)
  expect_equal(net_rows(plain, "AD")$relative, 0.2 * -1, tolerance = 1e-8)
  expect_true(is.na(net_rows(plain, "AD")$p_equivalence))
  # Single and combined effects have no relative deviation
  others <- est$estimates[est$estimates$component != "interaction", ]
  expect_true(all(is.na(others$relative) & is.na(others$p_equivalence)))
})

test_that("counts are tested on the square-root scale", {
  d <- two_by_two_counts()
  fit <- fit_2x2_eq(d, "count", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = "MU", inference = "delta", margin = 0.25)
  mu <- net_rows(est, "MU")
  # Poisson cell means of five units: the variance of the square root of a
  # mean is 1 / (4 x 5), whatever the mean
  y <- c(100, 50, 70, 25) / 5
  lower <- mu_distance_by_elimination(y, rep(1 / 20, 4), -0.25, sqrt, function(z) pmax(z, 0)^2)
  upper <- mu_distance_by_elimination(y, rep(1 / 20, 4), 0.25, sqrt, function(z) pmax(z, 0)^2)
  expect_equal(mu$df, Inf)
  # (up to the precision of the fitted means)
  expect_equal(mu$p_equivalence,
               max(stats::pnorm(lower, lower.tail = FALSE), stats::pnorm(upper)),
               tolerance = 1e-4)
})

test_that("models of a transformed response are tested on the scale of the model", {
  skip_if_not_installed("emmeans")
  d <- two_by_two_gaussian(means = c(10, 8, 7, 4), dev = c(-0.4, -0.2, 0, 0.2, 0.4))
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- stats::lm(log(y) ~ a * b, data = d)
  est <- estimate_interactions(fit, des, null_models = "MU", inference = "delta", margin = 0.25)
  mu <- net_rows(est, "MU")
  # Means back-transformed with exp(); the variance of each mean is constant
  # on the log scale, RSS / 16 / 5
  cell <- interaction(d$a, d$b)
  y <- exp(tapply(log(d$y), cell, mean))[c("0.0", "1.0", "0.1", "1.1")]
  v <- sum(stats::residuals(fit)^2) / 16 / 5
  lower <- mu_distance_by_elimination(unname(y), rep(v, 4), -0.25, log, exp)
  upper <- mu_distance_by_elimination(unname(y), rep(v, 4), 0.25, log, exp)
  expect_equal(mu$p_equivalence, tost_from_distances(lower, upper, 16), tolerance = 1e-5)
})

test_that("the equivalence test does not depend on the inference", {
  d <- two_by_two_gaussian(means = c(10, 8, 7, 5.4))
  fit <- fit_2x2_eq(d, "y", list(gaussian = stats::gaussian()))
  delta <- estimate_interactions(fit, null_models = c("AD", "BL"), inference = "delta",
                                 margin = 0.15)
  sim <- estimate_interactions(fit, null_models = c("AD", "BL"), n_sim = 40000, seed = 4,
                               margin = 0.15)
  for (nm in c("AD", "BL")) {
    expect_equal(net_rows(sim, nm)$relative, net_rows(delta, nm)$relative)
    expect_equal(net_rows(sim, nm)$p_equivalence, net_rows(delta, nm)$p_equivalence,
                 tolerance = 1e-10)
    # Percentile and Wald intervals of a ratio differ a little
    expect_lt(abs(net_rows(sim, nm)$relative_upper - net_rows(delta, nm)$relative_upper), 0.02)
  }
  # Adjusted within the same families as the tests of deviation
  e <- sim$estimates
  ok <- is.finite(e$p_equivalence)
  expect_true(all(e$p_equivalence_adjusted[ok] >= e$p_equivalence[ok]))
  expect_identical(match("p_equivalence_adjusted", names(e)), match("p_equivalence", names(e)) + 1L)
})

test_that("the scale of the tests follows the error distribution", {
  near <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-12))
  expect_true(near(stabilizing_scale(stats::gaussian(link = "log"))$h(4), 4))
  expect_true(near(stabilizing_scale(stats::quasipoisson())$h(4), 2))
  expect_true(near(stabilizing_scale(stats::binomial())$h(0.25), asin(0.5)))
  expect_true(near(stabilizing_scale(stats::Gamma(link = "log"))$h(4), log(4)))
  nb <- stabilizing_scale(structure(list(family = "Negative Binomial(2.5)", link = "log"),
                                    class = "family"))
  expect_true(near(nb$h(5), asinh(sqrt(2))))
  expect_true(near(nb$h_inv(nb$h(7)), 7))
  # derivative of h
  for (fam in list(stats::poisson(), stats::binomial(), stats::Gamma(link = "log"))) {
    st <- stabilizing_scale(fam)
    x <- 0.3
    expect_equal(st$dh(x), (st$h(x + 1e-6) - st$h(x - 1e-6)) / 2e-6, tolerance = 1e-6)
  }
  expect_equal(nb$dh(5), (nb$h(5 + 1e-6) - nb$h(5 - 1e-6)) / 2e-6, tolerance = 1e-6)
  # Unknown families and transformed responses use the scale of the model
  expect_null(stabilizing_scale(stats::inverse.gaussian()))
  expect_null(stabilizing_scale(stats::gaussian(), transformed = TRUE))
  expect_null(stabilizing_scale(NULL))
})

test_that("signed distances to linear and curved bounds", {
  V <- diag(c(0.5, 0.25))
  # linear bound: the t statistic of the linear combination
  g <- function(Z) Z[, 1] + 2 * Z[, 2] - 1
  z0 <- c(1, 1)
  expect_equal(signed_distance(z0, V, g), 2 / sqrt(0.5 + 4 * 0.25), tolerance = 1e-8)
  expect_equal(signed_distance(c(0, 0), V, g), -1 / sqrt(0.5 + 4 * 0.25), tolerance = 1e-8)
  expect_identical(signed_distance(c(1, 0), V, g), 0)
  # circle of radius 1 around the origin, from (3, 0) with V = I: distance 2
  circle <- function(Z) Z[, 1]^2 + Z[, 2]^2 - 1
  expect_equal(signed_distance(c(3, 0.001), diag(2), circle), 2, tolerance = 1e-6)
  # a correlated metric needs no inverse
  W <- matrix(c(1, 0.9, 0.9, 1), 2)
  lin <- function(Z) Z[, 1] - Z[, 2] - 1
  expect_equal(signed_distance(c(0, 0), W, lin), -1 / sqrt(2 - 2 * 0.9), tolerance = 1e-8)
})

test_that("verdicts combine the tests of deviation and of equivalence", {
  d <- two_by_two_counts()
  fit <- fit_2x2_eq(d, "count", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), n_sim = 500, seed = 1,
                               margin = 0.1)
  verdict_of <- function(p, p_eq) {
    e <- est$estimates
    i <- e$component == "interaction" & e$null_model == "AD"
    e$p_adjusted[i] <- p
    e$p_equivalence_adjusted[i] <- p_eq
    x <- est
    x$estimates <- e
    cls <- classify_interactions(x)
    cls$verdict[cls$null_model == "AD"]
  }
  expect_equal(verdict_of(0.01, 0.6), "deviation")
  expect_equal(verdict_of(0.01, 0.01), "negligible")
  expect_equal(verdict_of(0.4, 0.01), "equivalent")
  expect_equal(verdict_of(0.4, 0.6), "inconclusive")
  # The level of classify_interactions() applies to both tests
  e <- est$estimates
  i <- e$component == "interaction" & e$null_model == "AD"
  e$p_adjusted[i] <- 0.03
  e$p_equivalence_adjusted[i] <- 0.6
  est$estimates <- e
  expect_equal(classify_interactions(est, alpha = 0.01)$verdict[1L], "inconclusive")

  sm <- summary(est)
  expect_equal(sum(unlist(sm$verdicts[c("n_deviation", "n_negligible", "n_equivalent",
                                        "n_inconclusive")])), 2L)
  expect_output(print(sm), "equivalence margin of \\+/-0.1")
  expect_output(print(est), "Equivalence margin: \\+/-0.1")

  # Without a margin there is no verdict
  plain <- estimate_interactions(fit, null_models = "AD", n_sim = 500, seed = 1)
  expect_true(is.na(classify_interactions(plain)$verdict))
  expect_null(summary(plain)$verdicts)
})

test_that("combinations with only zeros have a relative deviation but no equivalence test", {
  d <- zero_combination_counts()
  fit <- fit_2x2_eq(d, "y", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), n_sim = 2000, seed = 1,
                               margin = 0.1)
  ad <- net_rows(est, "AD")
  mu <- net_rows(est, "MU")
  expect_equal(ad$status, "boundary")
  expect_equal(ad$relative, (0 - 15 - 12 + 20) / 20, tolerance = 1e-4)
  # MU: (observed - expected) / reference, finite although the log ratio is -Inf
  expect_equal(mu$relative, (0 - 15 * 12 / 20) / 20, tolerance = 1e-4)
  expect_true(is.na(ad$p_equivalence) && is.na(ad$relative_lower))
  cls <- classify_interactions(est)
  expect_equal(cls$verdict, c("deviation", "deviation"))
})

test_that("averages over times carry the relative deviation", {
  d <- repeated_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank", time = "day")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  for (inference in c("delta", "simulation")) {
    est <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = inference,
                                 n_sim = 2000, seed = 2, margin = 0.2)
    ad <- net_rows(est, "AD")
    by_time <- ad[ad$time != "overall", ]
    overall <- ad[ad$time == "overall", ]
    expect_equal(overall$relative, mean(by_time$relative), tolerance = 1e-8)
    expect_true(is.finite(overall$p_equivalence))
    expect_true(overall$relative_lower < overall$relative &&
                  overall$relative < overall$relative_upper)
  }
})

test_that("the relative deviation needs a positive reference", {
  d <- two_by_two_gaussian(means = c(-1, 2, 3, 4))
  fit <- fit_2x2_eq(d, "y", list(gaussian = stats::gaussian()))
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta", margin = 0.1)
  ad <- net_rows(est, "AD")
  expect_equal(ad$status, "ok")
  expect_true(is.na(ad$relative) && is.na(ad$p_equivalence))
  expect_match(ad$note, "reference_not_positive")
  expect_true(classify_interactions(est)$verdict %in% c("deviation", "inconclusive"))
})

test_that("the margin is checked", {
  d <- two_by_two_counts()
  fit <- fit_2x2_eq(d, "count", list(poisson = stats::poisson()))
  for (bad in list(-0.1, 0, "0.1", c(0.1, 0.2), Inf, NA_real_)) {
    expect_error(estimate_interactions(fit, inference = "delta", margin = bad),
                 class = "must_error_input")
  }
})

test_that("models fitted outside must get the same equivalence tests", {
  skip_if_not_installed("emmeans")
  d <- two_by_two_gaussian()
  fit <- fit_2x2_eq(d, "y", list(gaussian = stats::gaussian()))
  inside <- estimate_interactions(fit, null_models = c("AD", "BL"), inference = "delta",
                                  margin = 0.1)
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  outside <- estimate_interactions(stats::lm(y ~ a * b, data = d), des,
                                   null_models = c("AD", "BL"), inference = "delta",
                                   margin = 0.1)
  for (nm in c("AD", "BL")) {
    expect_equal(net_rows(outside, nm)$relative, net_rows(inside, nm)$relative,
                 tolerance = 1e-6)
    expect_equal(net_rows(outside, nm)$p_equivalence, net_rows(inside, nm)$p_equivalence,
                 tolerance = 1e-6)
  }
})
