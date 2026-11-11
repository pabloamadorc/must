test_that("two-stressor weights are the classical interaction contrast", {
  w <- contrast_weights(2)
  # Subset order: control, A, B, AB
  expect_equal(w$net, c(1, -1, -1, 1))
  expect_equal(w$pure, c(1, -1, -1, 1))
  expect_equal(w$single[[1]], c(-1, 1, 0, 0))
  expect_equal(w$combined, c(-1, 0, 0, 1))
})

test_that("three-stressor net deviation equals pure plus the pairwise interactions", {
  w <- contrast_weights(3)
  subsets <- all_subsets(3)
  y <- c(20, 17, 16, 18, 12.5, 14, 14.5, 7.5)
  pos <- function(t) which(vapply(subsets, function(s) identical(s, t), logical(1)))
  # Pairwise interaction of stressors i and j, third stressor at control
  pair <- function(v, i, j) v[pos(c(i, j))] - v[pos(i)] - v[pos(j)] + v[pos(integer(0))]
  for (f in list(identity, log)) {
    fy <- f(y)
    pairwise <- pair(fy, 1L, 2L) + pair(fy, 1L, 3L) + pair(fy, 2L, 3L)
    expect_equal(sum(w$net * fy), sum(w$pure * fy) + pairwise)
  }
  # By hand: yABC - yA - yB - yC + 2 y0 = 7.5 - 17 - 16 - 18 + 40
  expect_equal(sum(w$net * y), -3.5)
})

test_that("dominance interaction is min(yAB, y0) - min(yA, yB)", {
  M <- rbind(c(10, 8, 6, 3), c(10, 12, 11, 15), c(10, 12, 6, 9))
  # 3 - 6 = -3; min(15, 10) - 11 = -1; 9 - 6 = 3
  expect_equal(dominance_interaction(M), c(-3, -1, 3))
})

test_that("expected combined responses follow each null model", {
  mu <- c(10, 8, 6, 3)
  expect_equal(expected_combined(mu, 2, "AD"), 4)
  expect_equal(expected_combined(mu, 2, "MU"), 4.8)
  expect_equal(expected_combined(mu, 2, "DO"), 6)
  # Dominance has no point expectation when no stressor reduces the response
  expect_true(is.na(expected_combined(c(10, 12, 11, 15), 2, "DO")))
  # Multiplicative needs positive means
  expect_true(is.na(expected_combined(c(10, -1, 6, 3), 2, "MU")))
})

test_that("simulation summaries give percentile intervals and bounded p-values", {
  draws <- c(seq(0.5, 9.5, by = 1))
  s <- summarise_draws(draws, level = 0.8)
  expect_equal(s$lower, unname(stats::quantile(draws, 0.1)))
  expect_equal(s$upper, unname(stats::quantile(draws, 0.9)))
  # No draw crosses zero: p is bounded by 1 / n
  expect_equal(s$p_value, 1 / length(draws))
  centred <- summarise_draws(c(-2, -1, 1, 2), level = 0.5)
  expect_equal(centred$p_value, 1)
})

test_that("delta gradients and Wald summaries use the link scale and the t distribution", {
  # Identity link: gradient = weights
  V <- diag(c(1, 2, 3, 4))
  g <- delta_gradient(c(1, -1, -1, 1), mu = c(5, 4, 3, 2), dmu = rep(1, 4), null_model = "AD")
  expect_equal(sqrt(drop(t(g) %*% V %*% g)), sqrt(10))
  # Log link with MU: d log(mu) / d eta = 1
  mu <- c(5, 4, 3, 2)
  g2 <- delta_gradient(c(1, -1, -1, 1), mu = mu, dmu = mu, null_model = "MU")
  expect_equal(g2, c(1, -1, -1, 1))
  normal <- summarise_wald(2, 1, Inf, 0.95)
  expect_equal(normal$p_value, 2 * stats::pnorm(-2))
  expect_equal(normal$upper, 2 + stats::qnorm(0.975))
  t5 <- summarise_wald(2, 1, 5, 0.95)
  expect_equal(t5$p_value, 2 * stats::pt(-2, 5))
  expect_equal(t5$lower, 2 - stats::qt(0.975, 5))
  expect_true(is.na(summarise_wald(1, 0, 5, 0.95)$p_value))
})

test_that("degrees of freedom follow the strata of the design", {
  info <- list(has_dispersion = TRUE, has_random = FALSE, repeated = FALSE,
               dispersion = "common", rank_total = 4L,
               contexts = data.frame(context = "", n_obs = 20L, rank = 4L, n_units = 20L,
                                     n_unit_groups = 4L, block_df = 0),
               whole = list())
  s <- data.frame(context = "", stressors = "a x b")
  expect_equal(set_df(info, s, "small_sample"), 16)
  expect_equal(set_df(info, s, "none"), Inf)
  info$repeated <- TRUE
  info$contexts$n_units <- 8L
  expect_equal(set_df(info, s, "small_sample"), 4)
  info$contexts$block_df <- 1
  expect_equal(set_df(info, s, "small_sample"), 3)
  info$whole <- list(list(group = "room", between = "a",
                          contexts = data.frame(context = "", n_groups = 4L, n_combos = 2L,
                                                block_df = 0)))
  expect_equal(set_df(info, s, "small_sample"), 2)
  # Sets without the stressor applied to whole groups lose the whole-group
  # stratum (4 groups - 2 combinations)
  expect_equal(set_df(info, data.frame(context = "", stressors = "b x c"), "small_sample"), 1)
  info$dispersion <- "cell"
  expect_true(is.na(set_df(info, s, "small_sample")))
  poisson <- list(has_dispersion = FALSE, has_random = FALSE)
  expect_equal(set_df(poisson, s, "small_sample"), Inf)
  # Without a dispersion parameter, one dispersion per `by` group changes
  # nothing
  two <- list(has_dispersion = FALSE, has_random = TRUE, repeated = TRUE, dispersion = "by",
              rank_total = 8L,
              contexts = data.frame(context = c("x", "y"), n_obs = 32L, rank = 4L,
                                    n_units = 16L, n_unit_groups = 4L, block_df = 0),
              whole = list())
  expect_equal(set_df(two, data.frame(context = "x", stressors = "a x b"), "small_sample"), 24)
})

test_that("Welch-Satterthwaite degrees of freedom match the two-sample formula", {
  v <- c(4 / 5, 9 / 7)
  expect_equal(welch_df(v, c(5, 7)), sum(v)^2 / (v[1]^2 / 4 + v[2]^2 / 6))
  expect_equal(welch_df(c(1, 1, 1, 1), rep(5, 4)), 16)
  expect_equal(welch_df(c(1, 0), c(1, 5)), 0)
  expect_equal(welch_df(c(0, 0), c(5, 5)), Inf)
})

test_that("Bliss contrasts follow Beppler et al. (2016)", {
  # Two stressors: w = (1, 0.5, 0.7, 0.25)
  M <- matrix(c(20, 10, 14, 5), nrow = 1)
  w <- contrast_weights(2)
  specs <- component_specs("BL", 2, c("a", "b"), w)
  values <- vapply(specs, function(sp) bliss_contrast(M, 2, sp), numeric(1))
  # Single and combined effects w - 1; net interaction wAB - wA wB (eq. 2.1)
  expect_equal(values, c(-0.5, -0.3, -0.75, 0.25 - 0.5 * 0.7))
  # Three stressors: net DA (eq. 2.2) and emergent E3 (eq. 2.4)
  y <- c(20, 16, 12, 18, 9, 13.5, 10, 6)
  W <- y / y[[1L]]
  M3 <- matrix(y, nrow = 1)
  w3 <- contrast_weights(3)
  net <- list(component = "interaction", estimand = "net", weights = w3$net)
  pure <- list(component = "interaction", estimand = "pure", weights = w3$pure)
  expect_equal(bliss_contrast(M3, 3, net), W[[8]] - W[[2]] * W[[3]] * W[[4]])
  e3 <- W[[8]] - W[[2]] * W[[7]] - W[[3]] * W[[6]] - W[[4]] * W[[5]] +
    2 * W[[2]] * W[[3]] * W[[4]]
  expect_equal(bliss_contrast(M3, 3, pure), e3)
  # E3 = DA minus the pairwise interactions weighted by the third stressor
  # (eq. 2.3)
  pair <- function(xy, x, yy) W[[xy]] - W[[x]] * W[[yy]]
  expect_equal(e3, bliss_contrast(M3, 3, net) - W[[2]] * pair(7, 3, 4) -
                 W[[3]] * pair(6, 2, 4) - W[[4]] * pair(5, 2, 3))
  # With Bliss-additive pairs, E3 and DA coincide
  y2 <- c(10, 8, 5, 4, 4, 3.2, 2, 1)
  expect_equal(bliss_contrast(matrix(y2, nrow = 1), 3, pure),
               bliss_contrast(matrix(y2, nrow = 1), 3, net))
  expect_error(bliss_contrast(matrix(1, nrow = 1, ncol = 16), 4, pure),
               class = "must_error_internal")
})

test_that("rescaled Bliss interactions reach -1 at lethal synergy and 1 at full buffering", {
  # w = (1, 0.5, 0.7, wAB): lethal synergy (wAB = 0), buffering (wAB = 0.5)
  M <- function(ab) matrix(c(20, 10, 14, 20 * ab), nrow = 1)
  value <- function(ab) ab - 0.5 * 0.7
  expect_equal(bliss_rescaled(value(0), M(0), 2), -1)
  expect_equal(bliss_rescaled(value(0.5), M(0.5), 2), 1)
  expect_equal(bliss_rescaled(value(0.25), M(0.25), 2), (0.25 - 0.35) / 0.35)
  expect_equal(bliss_rescaled(value(0.4), M(0.4), 2), (0.4 - 0.35) / (0.5 - 0.35))
  expect_equal(bliss_rescaled(0, M(0.35), 2), 0)
  expect_equal(bliss_class(c(-0.8, -0.5, 0.2, 0.51, NA)),
               c("synergy", "additive", "additive", "antagonism", NA))
  # Three stressors: the denominator replaces wXYZ by 0 or by min(wX, wY, wZ)
  y <- c(20, 16, 12, 18, 9, 13.5, 10, 6)
  W <- y / y[[1L]]
  e3 <- W[[8]] - W[[2]] * W[[7]] - W[[3]] * W[[6]] - W[[4]] * W[[5]] +
    2 * W[[2]] * W[[3]] * W[[4]]
  expect_equal(bliss_rescaled(e3, matrix(y, nrow = 1), 3), e3 / abs(e3 - W[[8]]))
})

test_that("Bliss expectations multiply the relative single effects", {
  expect_equal(expected_combined(c(20, 10, 14, 5), 2, "BL"), 20 * 0.5 * 0.7)
  expect_true(is.na(expected_combined(c(0, 1, 2, 3), 2, "BL")))
  expect_equal(null_model_scale("BL"), "response / control")
})

test_that("the stress addition model adds general stress (Liess et al. 2016)", {
  # A single stressor is predicted by itself: the beta quantile and the beta
  # distribution function are inverse
  expect_equal(sam_expected(matrix(0.7, nrow = 1), 3.2), 0.7)
  # No stress from a stressor at or above the reference
  expect_equal(sam_expected(matrix(c(1, 0.8), nrow = 1), 3.2), 0.8)
  expect_equal(sam_expected(matrix(c(1.2, 0.8), nrow = 1), 3.2), 0.8)
  # Two stressors killing half each exert a total stress of 1: no survivors
  expect_equal(sam_expected(matrix(c(0.5, 0.5), nrow = 1), 3.2), 0)
  # By hand, eqs. 2-5: S = qbeta(0.1) + qbeta(0.2), survival 1 - pbeta(S)
  s <- stats::qbeta(0.1, 3.2, 3.2) + stats::qbeta(0.2, 3.2, 3.2)
  expect_equal(sam_expected(matrix(c(0.9, 0.8), nrow = 1), 3.2),
               1 - stats::pbeta(s, 3.2, 3.2))
  # More synergistic than Bliss independence for these effects
  expect_lt(sam_expected(matrix(c(0.9, 0.8), nrow = 1), 3.2), 0.9 * 0.8)
  # Rows are evaluated separately
  expect_equal(sam_expected(rbind(c(0.7, 1), c(0.5, 0.5)), 3.2), c(0.7, 0))
  sp <- list(component = "interaction", estimand = "net", shape = 3.2)
  M <- matrix(c(20, 18, 16, 10), nrow = 1)
  expect_equal(sam_contrast(M, 2, sp), 0.5 - sam_expected(matrix(c(0.9, 0.8), 1), 3.2))
})

test_that("log-logistic curves are fitted by least squares", {
  dose <- c(1, 4, 16)
  w <- 1 / (1 + (dose / 4)^2)
  fit <- fit_loglogistic(matrix(w, nrow = 1), log(dose))
  expect_equal(unname(fit[1, ]), c(2 * log(4), -2), tolerance = 1e-8)
  noisy <- w + c(0.03, -0.04, 0.02)
  fit2 <- fit_loglogistic(matrix(noisy, nrow = 1), log(dose))
  ref <- stats::nls(noisy ~ stats::plogis(a + b * log(dose)), start = list(a = 2, b = -1.5))
  expect_equal(unname(fit2[1, ]), unname(stats::coef(ref)), tolerance = 1e-5)
  # Responses at or above the reference are kept in the fit
  flat <- fit_loglogistic(matrix(c(1.02, 0.6, 0.1), nrow = 1), log(dose))
  expect_true(flat[1, 2] < 0)
})

test_that("concentration addition solves sum(x_i / EC_i(w)) = 1", {
  a <- fit_loglogistic(matrix(1 / (1 + (c(1, 4, 16) / 4)^2), nrow = 1), log(c(1, 4, 16)))
  b <- fit_loglogistic(matrix(1 / (1 + (c(2, 10, 50) / 10)), nrow = 1), log(c(2, 10, 50)))
  equation <- function(w) 4 / (4 * ((1 - w) / w)^(1 / 2)) + 10 / (10 * (1 - w) / w) - 1
  expect_equal(ca_expected(list(a, b), c(4, 10)),
               stats::uniroot(equation, c(1e-9, 1 - 1e-9), tol = 1e-12)$root, tolerance = 1e-8)
  # A sham combination of a substance with itself is its curve at the summed dose
  expect_equal(ca_expected(list(a, a), c(4, 4)), 1 / (1 + (8 / 4)^2), tolerance = 1e-8)
  # A curve that does not decrease adds no dose equivalent
  flat <- cbind(alpha = 0, beta = 0.5)
  expect_equal(ca_expected(list(a, flat), c(4, 10)), 0.5, tolerance = 1e-8)
  expect_equal(ca_expected(list(flat, flat), c(4, 10)), 1)
})
