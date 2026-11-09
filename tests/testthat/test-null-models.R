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
