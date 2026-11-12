fit_2x2 <- function(d, response, families) {
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit_responses(d, response, des, families = families, residual_checks = FALSE, quiet = TRUE)
}

interaction_rows <- function(est, null_model, estimand = "net") {
  e <- est$estimates
  e[e$component == "interaction" & e$null_model == null_model & e$estimand %in% estimand, ]
}

test_that("Gaussian AD and MU estimates and delta SEs match hand calculations", {
  d <- two_by_two_gaussian()
  fit <- fit_2x2(d, "y", list(gaussian = stats::gaussian()))
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta",
                               correction = "none")
  # ML residual variance: RSS / N = 4 * 2.5 / 20
  sigma2 <- 4 * 2.5 / 20
  ad <- interaction_rows(est, "AD")
  expect_equal(ad$estimate, 3 - 8 - 7 + 10, tolerance = 1e-6)
  expect_equal(ad$se, sqrt(sigma2 * 4 / 5), tolerance = 1e-4)
  expect_equal(ad$observed, 3, tolerance = 1e-6)
  expect_equal(ad$expected, 8 + 7 - 10, tolerance = 1e-6)
  mu <- interaction_rows(est, "MU")
  expect_equal(mu$estimate, log(3) - log(8) - log(7) + log(10), tolerance = 1e-6)
  # Identity link: d log(mu) / d eta = 1 / mu
  expect_equal(mu$se, sqrt(sigma2 / 5 * sum(1 / c(10, 8, 7, 3)^2)), tolerance = 1e-4)
  expect_equal(mu$expected, 8 * 7 / 10, tolerance = 1e-6)
  singles <- est$estimates[est$estimates$component == "single" & est$estimates$null_model == "AD", ]
  expect_equal(singles$estimate, c(8 - 10, 7 - 10), tolerance = 1e-6)
})

test_that("the small-sample correction reproduces the t intervals of lm()", {
  d <- two_by_two_gaussian(dev = c(-1.3, -0.4, 0, 0.6, 1.1))
  fit <- fit_2x2(d, "y", list(gaussian = stats::gaussian()))
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  ad <- interaction_rows(est, "AD")
  ref <- stats::lm(y ~ a * b, data = d)
  expect_equal(ad$df, 16)
  expect_equal(ad$se, unname(summary(ref)$coefficients[4L, 2L]), tolerance = 1e-5)
  expect_equal(ad$p_value, unname(summary(ref)$coefficients[4L, 4L]), tolerance = 1e-5)
  expect_equal(c(ad$lower, ad$upper), unname(stats::confint(ref)[4L, ]), tolerance = 1e-5)
  sim <- estimate_interactions(fit, null_models = "AD", n_sim = 40000, seed = 5)
  ks <- interaction_rows(sim, "AD")
  expect_equal(c(ks$lower, ks$upper), unname(stats::confint(ref)[4L, ]), tolerance = 0.03)
  expect_true(is.finite(ks$mc_se) && ks$mc_se > 0)
  expect_true(is.na(ad$mc_se))
  none <- estimate_interactions(fit, null_models = "AD", inference = "delta", correction = "none")
  expect_equal(interaction_rows(none, "AD")$df, Inf)
  expect_lt(interaction_rows(none, "AD")$se, ad$se)
  expect_equal(est$settings$correction, "small_sample")
  expect_output(print(none), "no small-sample correction")
})

test_that("families without dispersion keep the normal distribution", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  a <- estimate_interactions(fit, null_models = "MU", inference = "delta")
  b <- estimate_interactions(fit, null_models = "MU", inference = "delta", correction = "none")
  expect_equal(a$estimates, b$estimates)
  expect_true(all(a$estimates$df == Inf))
})

test_that("degrees of freedom follow the units with repeated measures and whole plots", {
  d <- repeated_counts()
  d$y <- d$count + rep(c(0.3, -0.2, 0.1, -0.4), length.out = nrow(d))
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank", time = "day")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  # 16 tanks in 4 treatment combinations
  expect_true(all(est$estimates$df == 12))
  rooms <- jars_in_rooms(rooms_per_level = 2L)
  des2 <- describe_design(rooms, c("temperature", "nutrients"),
                          control = c(temperature = "15", nutrients = "no"),
                          unit = "id", groups = "room")
  fit2 <- fit_responses(rooms, "y", des2, families = list(gaussian = stats::gaussian()),
                        residual_checks = FALSE, quiet = TRUE)
  est2 <- estimate_interactions(fit2, null_models = "AD", inference = "delta")
  # Temperature is applied to 4 rooms (2 per level): 4 - 2 degrees of freedom
  expect_true(all(est2$estimates$df == 2))
  # The standard error stays finite with 2 degrees of freedom
  sim2 <- estimate_interactions(fit2, null_models = "AD", n_sim = 4000, seed = 1)
  expect_equal(sim2$estimates$se, est2$estimates$se, tolerance = 0.05)
})

test_that("split plots: within-group sets lose the whole-group stratum", {
  # Temperature applied to 4 rooms (2 per level); nutrients x salinity within
  # each room, 2 jars per combination
  d <- expand.grid(room_rep = 1:2, temperature = c("15", "20"), nutrients = c("no", "yes"),
                   salinity = c("low", "high"), jar = 1:2, stringsAsFactors = FALSE)
  d$room <- paste0("R", d$temperature, "-", d$room_rep)
  d$id <- seq_len(nrow(d))
  room_effect <- c(0.4, -0.3, 0.2, -0.1)[match(d$room, sort(unique(d$room)))]
  d$y <- 10 - 2 * (d$temperature == "20") - (d$nutrients == "yes") +
    0.5 * (d$salinity == "high") + room_effect + c(-0.2, 0.2)[d$jar] +
    0.1 * (d$nutrients == "yes") * c(-1, 1)[d$jar]
  des <- describe_design(d, c("temperature", "nutrients", "salinity"),
                         control = c(temperature = "15", nutrients = "no", salinity = "low"),
                         unit = "id", groups = "room")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  e <- est$estimates[est$estimates$component == "interaction" & est$estimates$estimand == "net", ]
  # 32 jars in 8 combinations; rooms: 4 - 2 temperature levels
  expect_equal(e$df[e$stressors == "nutrients x salinity"], 32 - 8 - 2)
  expect_true(all(e$df[grepl("temperature", e$stressors)] == 2))
})

test_that("units sampled once each are nested in combination x time", {
  # Three new jars per combination on each of three days, two subsamples each
  jars <- expand.grid(a = c("0", "1"), b = c("0", "1"), day = c(1, 2, 3), rep = 1:3,
                      stringsAsFactors = FALSE)
  jars$jar <- paste0("J", seq_len(nrow(jars)))
  d <- jars[rep(seq_len(nrow(jars)), each = 2), ]
  d$y <- 10 - (d$a == "1") - 2 * (d$b == "1") + 0.3 * d$day +
    rep(c(-0.4, 0.1, 0.3, -0.2, 0.5, -0.3), length.out = nrow(d))
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "jar", time = "day")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  # 36 jars in 12 combinations x days
  expect_true(all(est$estimates$df == 24))
})

test_that("partly confounded groups with a random intercept bound the degrees of freedom", {
  # Temperature 15 in two rooms, 20 in a single room
  d <- expand.grid(room = c("R15a", "R15b", "R20"), nutrients = c("no", "yes"), jar = 1:3,
                   stringsAsFactors = FALSE)
  d$temperature <- ifelse(d$room == "R20", "20", "15")
  d$id <- seq_len(nrow(d))
  d$y <- 10 - 2 * (d$temperature == "20") - (d$nutrients == "yes") +
    c(R15a = 0.3, R15b = -0.3, R20 = 0.1)[d$room] + c(-0.2, 0, 0.2)[d$jar]
  des <- describe_design(d, c("temperature", "nutrients"),
                         control = c(temperature = "15", nutrients = "no"),
                         unit = "id", groups = "room")
  expect_equal(des$groups$structure, "confounded")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  expect_match(fit$fits$y$formula, ".g_room", fixed = TRUE)
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  # 3 rooms - 2 temperatures
  expect_true(all(est$estimates$df == 1))
})

test_that("blocks given in `random` count like automatic ones", {
  d <- expand.grid(a = c("0", "1"), b = c("0", "1"), block = paste0("B", 1:4),
                   stringsAsFactors = FALSE)
  d$id <- seq_len(nrow(d))
  d$y <- 5 + (d$a == "1") - (d$b == "1") + c(B1 = 0.5, B2 = -0.2, B3 = 0.1, B4 = -0.4)[d$block] +
    c(0.1, -0.1, 0.05, -0.05)[seq_len(nrow(d)) %% 4 + 1]
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id",
                         groups = "block")
  fams <- list(gaussian = stats::gaussian())
  auto <- fit_responses(d, "y", des, families = fams, residual_checks = FALSE, quiet = TRUE)
  given <- fit_responses(d, "y", des, families = fams, random = "(1 | block)",
                         residual_checks = FALSE, quiet = TRUE)
  df_auto <- estimate_interactions(auto, null_models = "AD", inference = "delta")$estimates$df
  df_given <- estimate_interactions(given, null_models = "AD", inference = "delta")$estimates$df
  # 16 plots - 4 combinations - 3 block degrees of freedom
  expect_true(all(df_auto == 9))
  expect_equal(df_given, df_auto)
})

test_that("averages over times with one dispersion per cell follow their t distribution", {
  d <- repeated_counts()
  d$y <- d$count * c(1, 1.3, 0.8, 1.1)[as.integer(factor(d$tank)) %% 4 + 1] +
    0.5 * (d$a == "1")
  d$tank_day <- paste(d$tank, d$day)
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank_day", time = "day")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       dispersion = "cell", residual_checks = FALSE, quiet = TRUE)
  delta <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  sim <- estimate_interactions(fit, null_models = "AD", n_sim = 20000, seed = 4)
  pick <- function(x) {
    x$estimates[x$estimates$component == "interaction" & x$estimates$time == "overall", ]
  }
  width <- function(r) r$upper - r$lower
  expect_equal(width(pick(sim)), width(pick(delta)), tolerance = 0.04)
  expect_equal(pick(sim)$df, pick(delta)$df)
})

test_that("one dispersion per cell gives Welch-Satterthwaite degrees of freedom", {
  d <- two_by_two_gaussian(dev = c(-1.3, -0.4, 0, 0.6, 1.1))
  d$y <- d$mean + c(1, 1, 3, 0.5)[match(paste(d$a, d$b), c("0 0", "1 0", "0 1", "1 1"))] *
    (d$y - d$mean)
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       dispersion = "cell", residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  ad <- interaction_rows(est, "AD")
  v <- tapply(d$y, paste(d$a, d$b), stats::var) / 5
  expect_equal(ad$se, sqrt(sum(v)), tolerance = 1e-4)
  expect_equal(ad$df, sum(v)^2 / sum(v^2 / 4), tolerance = 1e-4)
})

test_that("Poisson MU and AD standard errors match hand calculations", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
  # Cell means 20, 10, 14, 5 (totals 100, 50, 70, 25 over 5 units)
  mu <- interaction_rows(est, "MU")
  expect_equal(mu$estimate, log(5) - log(10) - log(14) + log(20), tolerance = 1e-4)
  # Var(log mean) = 1 / cell total
  expect_equal(mu$se, sqrt(1 / 100 + 1 / 50 + 1 / 70 + 1 / 25), tolerance = 1e-4)
  ad <- interaction_rows(est, "AD")
  expect_equal(ad$estimate, 5 - 10 - 14 + 20, tolerance = 1e-4)
  # Var(mean) = mean / n
  expect_equal(ad$se, sqrt((20 + 10 + 14 + 5) / 5), tolerance = 1e-4)
})

test_that("Krinsky-Robb simulation agrees with the delta method and gives DO", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  delta <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
  sim <- estimate_interactions(fit, null_models = c("AD", "MU", "DO"), n_sim = 20000, seed = 42)
  expect_equal(interaction_rows(sim, "MU")$se, interaction_rows(delta, "MU")$se,
               tolerance = 0.03)
  do <- interaction_rows(sim, "DO")
  # min(yAB, y0) - min(yA, yB) = 5 - 10
  expect_equal(do$estimate, -5, tolerance = 1e-4)
  expect_equal(do$expected, 10, tolerance = 1e-4)
  expect_true(do$lower < do$estimate && do$estimate < do$upper)
})

test_that("a seed makes the simulation reproducible and restores the RNG state", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  set.seed(1)
  before <- .Random.seed
  a <- estimate_interactions(fit, n_sim = 500, seed = 7)
  expect_identical(.Random.seed, before)
  b <- estimate_interactions(fit, n_sim = 500, seed = 7)
  expect_identical(a$estimates, b$estimates)
})

test_that("AD is the default and each null model ignores the others requested", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  default <- estimate_interactions(fit, n_sim = 500, seed = 7)
  expect_identical(default$settings$null_models, "AD")
  expect_identical(unique(default$estimates$null_model), "AD")
  several <- estimate_interactions(fit, null_models = c("SAM", "MU", "AD", "DO"), n_sim = 500,
                                   seed = 7)
  # Returned in the canonical order, whatever the order requested
  expect_identical(several$settings$null_models, c("AD", "MU", "DO", "SAM"))
  ad <- several$estimates[several$estimates$null_model == "AD", ]
  rownames(ad) <- NULL
  expect_identical(ad, default$estimates)
  mu <- estimate_interactions(fit, null_models = "MU", n_sim = 500, seed = 7)$estimates
  mu_several <- several$estimates[several$estimates$null_model == "MU", ]
  rownames(mu) <- rownames(mu_several) <- NULL
  expect_identical(mu_several, mu)
})

test_that("the delta method is refused for the dominance model", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  expect_error(estimate_interactions(fit, null_models = c("AD", "DO"), inference = "delta"),
               "dominance", class = "must_error_input")
  expect_error(estimate_interactions(fit, null_models = "XX"), class = "must_error_input")
})

test_that("MU is not estimable when an identity-link mean is not positive", {
  d <- two_by_two_gaussian(means = c(2, 1, 1, -1))
  fit <- fit_2x2(d, "y", list(gaussian = stats::gaussian()))
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
  expect_equal(interaction_rows(est, "MU")$status, "not_estimable")
  expect_equal(interaction_rows(est, "AD")$status, "ok")
})

test_that("simple addition below zero is flagged for non-negative responses", {
  d <- two_by_two_counts()
  d$count <- c(18, 20, 22, 19, 21, 3, 5, 4, 4, 4, 4, 4, 3, 5, 4, 1, 1, 0, 2, 1)
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
  ad <- interaction_rows(est, "AD")
  expect_equal(ad$expected, 4 + 4 - 20, tolerance = 1e-4)
  expect_equal(ad$note, "expected_below_zero")
})

test_that("three stressors: net and pure interactions match hand calculations", {
  d <- three_way_gaussian()
  des <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                         unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  triple <- est$estimates[est$estimates$order == 3L & est$estimates$component == "interaction", ]
  # Cell means of the generating model (deviations sum to zero in each cell)
  expect_equal(triple$estimate[triple$estimand == "net"], 9.5 - 17 - 16 - 18 + 2 * 20,
               tolerance = 1e-6)
  expect_equal(triple$estimate[triple$estimand == "pure"], -2.5, tolerance = 1e-6)
  pairs <- est$estimates[est$estimates$order == 2L & est$estimates$component == "interaction", ]
  expect_equal(sort(pairs$estimate), sort(c(1.5, -1, 0.5)), tolerance = 1e-6)
  # Net = pure + sum of the pairwise interactions
  expect_equal(triple$estimate[triple$estimand == "net"],
               triple$estimate[triple$estimand == "pure"] + sum(pairs$estimate),
               tolerance = 1e-6)
  sim <- estimate_interactions(fit, null_models = "DO", n_sim = 200, seed = 1)
  do_triple <- sim$estimates[sim$estimates$order == 3L, ]
  expect_equal(do_triple$status, "not_defined")
})

test_that("pure interactions are not identifiable without the pairwise cells", {
  d <- three_way_gaussian()
  d <- d[!(rowSums(d[c("n", "s", "f")] == "1") == 2L), ]
  des <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                         unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  rows <- est$estimates[est$estimates$component == "interaction", ]
  expect_equal(rows$status[rows$estimand == "pure"], "not_identifiable")
  expect_equal(rows$status[rows$estimand == "net"], "ok")
})

test_that("repeated measures give one estimate per sampling time", {
  d <- repeated_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank", time = "day")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       time_structure = "additive", residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
  mu <- interaction_rows(est, "MU")
  expect_equal(mu$time, c("7", "14", "21", "overall"))
  # Additive structure on a log link: the MU interaction is the same every day
  expect_equal(diff(mu$estimate), c(0, 0, 0), tolerance = 1e-8)
  expect_true(all(est$estimates$p_adjusted >= est$estimates$p_value, na.rm = TRUE))
})

test_that("averages over times combine the estimates of every time", {
  d <- repeated_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank", time = "day")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
  ad <- interaction_rows(est, "AD")
  by_time <- ad[ad$time != "overall", ]
  overall <- ad[ad$time == "overall", ]
  expect_equal(overall$estimate, mean(by_time$estimate), tolerance = 1e-8)
  expect_lt(overall$se, max(by_time$se))
  # Times and averages are adjusted separately
  expect_equal(overall$p_adjusted, overall$p_value)
  sim <- estimate_interactions(fit, null_models = "AD", n_sim = 4000, seed = 3)
  s_ad <- interaction_rows(sim, "AD")
  expect_equal(s_ad$se[s_ad$time == "overall"], overall$se, tolerance = 0.05)
  none <- estimate_interactions(fit, null_models = "AD", inference = "delta", overall = FALSE)
  expect_false("overall" %in% none$estimates$time)
})

test_that("baseline times are modelled but not estimated", {
  d <- repeated_counts()
  base <- d[d$day == 7, ]
  base$day <- 0
  base$count <- round(30 * c(0.9, 1.1, 1, 1.05)[as.integer(factor(base$tank)) %% 4 + 1])
  d <- rbind(base, d)
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank", time = "day", baseline = 0)
  expect_output(print(des), "baseline: 0")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  expect_false("0" %in% est$estimates$time)
  expect_setequal(unique(est$estimates$time), c("7", "14", "21", "overall"))
  expect_error(describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                               unit = "tank", time = "day", baseline = 3),
               "not in column", class = "must_error_input")
  expect_error(describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                               unit = "tank", baseline = 0),
               "needs `time`", class = "must_error_input")
  hint <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                          unit = "tank", time = "day")
  expect_true("baseline" %in% hint$checks$check)
})

test_that("interactions can be estimated with other stressors away from control", {
  d <- three_way_gaussian()
  des <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                         unit = "id", background = "all")
  expect_equal(nrow(des$sets), 7L)
  expect_equal(sum(nzchar(des$sets$background)), 3L)
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  ns_f <- est$estimates[est$estimates$stressors == "n x s" & est$estimates$background == "f=1" &
                          est$estimates$component == "interaction", ]
  # n x s at f = 1: y(1,1,1) - y(1,0,1) - y(0,1,1) + y(0,0,1)
  m <- function(n, s, f) {
    20 - 3 * n - 4 * s - 2 * f + 1.5 * n * s - 1 * n * f + 0.5 * s * f - 2.5 * n * s * f
  }
  expect_equal(ns_f$estimate, m(1, 1, 1) - m(1, 0, 1) - m(0, 1, 1) + m(0, 0, 1),
               tolerance = 1e-6)
  single <- est$estimates[est$estimates$stressors == "n x s" & est$estimates$background == "f=1" &
                            est$estimates$component == "single" & est$estimates$stressor == "n", ]
  expect_equal(single$estimate, m(1, 0, 1) - m(0, 0, 1), tolerance = 1e-6)
  cls <- classify_interactions(est)
  expect_equal(sum(nzchar(cls$background)), 3L)
  expect_output(print(est), "background")
})

test_that("a combination with only zeros is tested from the probability of zeros", {
  d <- zero_combination_counts()
  fit <- fit_2x2(d, "y", list(poisson = stats::poisson()))
  expect_equal(unname(fit$fits$y$boundary), "zero")
  est <- estimate_interactions(fit, null_models = c("AD", "MU", "DO"), n_sim = 2000, seed = 1)
  ad <- interaction_rows(est, "AD")
  expect_equal(ad$status, "boundary")
  expect_equal(ad$inference, "boundary test")
  expect_equal(ad$estimate, 0 - 15 - 12 + 20, tolerance = 1e-4)
  expect_equal(ad$expected, 15 + 12 - 20, tolerance = 1e-4)
  expect_true(is.na(ad$lower))
  # P(five zeros | Poisson mean 7) = exp(-35), doubled
  expect_equal(ad$p_value, 2 * exp(-35), tolerance = 0.5)
  mu <- interaction_rows(est, "MU")
  expect_equal(mu$estimate, -Inf)
  expect_equal(mu$expected, 15 * 12 / 20, tolerance = 1e-4)
  expect_lt(mu$p_value, 1e-10)
  singles <- est$estimates[est$estimates$component == "single", ]
  expect_true(all(singles$status == "ok"))
  combined <- est$estimates[est$estimates$component == "combined", ]
  expect_true(all(combined$status == "boundary"))
  expect_true(all(est$estimates$p_adjusted >= est$estimates$p_value, na.rm = TRUE))
})

test_that("the covariance of the fixed effects is rebuilt when vcov() of glmmTMB fails", {
  beta <- c(a = 1, b = -20, c = 3, d = -20)
  est <- matrix(c(1, 0.5, 0.5, 2), 2, dimnames = list(c("a", "c"), c("a", "c")))
  # As glmmTMB 1.1.11 to 1.1.15 with two or more parameters fixed via `map`
  failing_default <- function(model, include_nonest = TRUE) {
    if (include_nonest) stop("number of items to replace is not a multiple of replacement length")
    list(cond = est)
  }
  V <- conditional_covariance(NULL, beta, vcov_fun = failing_default)
  expect_identical(dimnames(V), list(names(beta), names(beta)))
  expect_identical(V[c("a", "c"), c("a", "c")], est)
  expect_true(all(is.na(V[c("b", "d"), ])) && all(is.na(V[, c("b", "d")])))
  # A vcov() that works is used as it is
  full <- matrix(NA_real_, 4, 4, dimnames = list(names(beta), names(beta)))
  full[c("a", "c"), c("a", "c")] <- est
  expect_identical(conditional_covariance(NULL, beta, function(model, ...) list(cond = full)),
                   full)
  # When both fail, the error says what glmmTMB said
  expect_error(conditional_covariance(NULL, beta, function(model, ...) stop("boom")),
               "boom). Updating glmmTMB", fixed = TRUE, class = "must_error")
})

test_that("the probability of zeros allows for overdispersion and random effects", {
  expect_equal(prob_boundary_obs("poisson", 2, NA, NA, NA, "zero"), exp(-2))
  expect_equal(prob_boundary_obs("nbinom2", 2, NA, 1.5, NA, "zero"), stats::dnbinom(0, size = 1.5, mu = 2))
  expect_equal(prob_boundary_obs("binomial", 0.2, 10, NA, NA, "zero"), 0.8^10)
  expect_equal(prob_boundary_obs("binomial", 0.9, 10, NA, NA, "one"), 0.9^10)
  rows <- data.frame(.unit = c("u1", "u1", "u2"))
  # Without random effects: product over observations
  expect_equal(prob_all_boundary(2, rows, "poisson", "log", rep(NA, 3), NA, 0, "zero"),
               exp(-6))
  # A shared random intercept makes zeros in the same unit more likely together
  with_re <- prob_all_boundary(2, rows, "poisson", "log", rep(NA, 3), NA, 1, "zero")
  expect_gt(with_re, exp(-6))
  b <- seq(-8, 8, length.out = 4001)
  w <- stats::dnorm(b) * diff(b)[1]
  ref <- sum(w * exp(-2 * 2 * exp(b))) * sum(w * exp(-2 * exp(b)))
  expect_equal(with_re, ref, tolerance = 1e-4)
  expect_equal(prob_all_boundary(c(0, -1), rows, "poisson", "log", rep(NA, 3), NA, 0, "zero"),
               c(1, 1))
})

test_that("observations that share a random-effect level share their intercept", {
  rows <- data.frame(.unit = c("j1", "j2", "j3", "j4"), .g_room = c("R1", "R1", "R2", "R2"))
  expect_equal(random_components(rows, ".unit"), 1:4)
  expect_equal(random_components(rows, ".g_room"), c(1L, 1L, 2L, 2L))
  expect_equal(random_components(rows, c(".unit", ".g_room")), c(1L, 1L, 2L, 2L))
  crossed <- data.frame(.unit = 1:4, a = c(1, 1, 2, 2), b = c(1, 2, 1, 2))
  expect_equal(random_components(crossed, c("a", "b")), rep(1L, 4))
  expect_equal(random_components(rows, "missing"), rep(1L, 4))
  # Jars in two rooms with a room intercept: zeros in the same room go together
  by_room <- prob_all_boundary(1.5, rows, "poisson", "log", rep(NA, 4), NA, 0.8, "zero",
                               factors = ".g_room")
  by_jar <- prob_all_boundary(1.5, rows, "poisson", "log", rep(NA, 4), NA, 0.8, "zero",
                              factors = ".unit")
  expect_gt(by_room, by_jar)
  b <- seq(-8, 8, length.out = 4001)
  w <- stats::dnorm(b) * diff(b)[1]
  room <- sum(w * exp(-2 * 1.5 * exp(0.8 * b)))
  expect_equal(by_room, room^2, tolerance = 1e-4)
})

test_that("a single stressor that eliminates the response is tested, its interaction is not", {
  d <- zero_combination_counts()
  d$y[d$a == "1" & d$b == "0"] <- 0
  fit <- fit_2x2(d, "y", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = "AD", n_sim = 1000, seed = 1)
  single_a <- est$estimates[est$estimates$component == "single" & est$estimates$stressor == "a", ]
  expect_equal(single_a$status, "boundary")
  expect_lt(single_a$p_value, 1e-6)
  expect_equal(interaction_rows(est, "AD")$status, "not_estimable")
  expect_match(interaction_rows(est, "AD")$reason, "a=1")
})

test_that("binomial responses with trials and the complement for sensitivity", {
  d <- survival_jars()
  d$dead <- d$n - d$alive
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, c("alive", "dead"), des, trials = c(alive = "n", dead = "n"),
                       families = list(binomial = stats::binomial()),
                       residual_checks = FALSE, quiet = TRUE)
  survival <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
  mortality <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta",
                                     sensitivity = "dead")
  s <- survival$estimates[survival$estimates$response == "alive", ]
  m <- mortality$estimates[mortality$estimates$response == "dead", ]
  expect_true(all(grepl("complement", m$note)))
  expect_equal(m$estimate, s$estimate, tolerance = 1e-5)
  expect_equal(m$p_value, s$p_value, tolerance = 1e-5)
  raw <- mortality$estimates[mortality$estimates$response == "alive", ]
  expect_false(any(grepl("complement", raw$note)))
  sim <- estimate_interactions(fit, null_models = c("AD", "MU", "DO"), n_sim = 2000, seed = 2,
                               sensitivity = "dead")
  do <- sim$estimates[sim$estimates$response == "dead" & sim$estimates$null_model == "DO", ]
  expect_equal(do$status, "ok")
})

test_that("sensitivity of a non-binomial response makes MU and DO not applicable", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = c("AD", "MU", "DO"), n_sim = 500, seed = 1,
                               sensitivity = "count")
  e <- est$estimates
  expect_true(all(e$status[e$null_model %in% c("MU", "DO")] == "not_applicable"))
  expect_true(all(e$status[e$null_model == "AD"] == "ok"))
  expect_error(estimate_interactions(fit, sensitivity = "other"), class = "must_error_input")
})

test_that("tests of zeros with the delta method record their draws", {
  d <- zero_combination_counts()
  fit <- fit_2x2(d, "y", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta",
                               n_sim = 1000, seed = 1)
  expect_equal(est$settings$n_sim, 1000L)
  expect_true(any(est$estimates$status == "boundary"))
})

test_that("baselines with additive time are coded as untreated", {
  d <- repeated_counts()
  base <- d[d$day == 7, ]
  base$day <- 0
  base$count <- round(30 * c(0.9, 1.1, 1, 1.05)[as.integer(factor(base$tank)) %% 4 + 1])
  with_base <- rbind(base, d)
  des_base <- describe_design(with_base, c("a", "b"), control = c(a = "0", b = "0"),
                              unit = "tank", time = "day", baseline = 0)
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank", time = "day")
  fams <- list(poisson = stats::poisson())
  fit_base <- fit_responses(with_base, "count", des_base, families = fams,
                            time_structure = "additive", residual_checks = FALSE, quiet = TRUE)
  fit <- fit_responses(d, "count", des, families = fams, time_structure = "additive",
                       residual_checks = FALSE, quiet = TRUE)
  pick <- function(f) {
    e <- estimate_interactions(f, null_models = "MU", inference = "delta")$estimates
    e$estimate[e$time == "14"]
  }
  # Single, single, combined, interaction: the baseline does not dilute them
  expect_equal(pick(fit_base), pick(fit), tolerance = 0.02)
  expect_equal(pick(fit)[1:2], log(c(0.6, 0.5)), tolerance = 0.05)
})

test_that("without a seed, the one drawn is stored and repeats the estimates", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  a <- estimate_interactions(fit, n_sim = 500)
  expect_true(is.integer(a$settings$seed))
  b <- estimate_interactions(fit, n_sim = 500, seed = a$settings$seed)
  expect_identical(a$estimates, b$estimates)
})

test_that("print() shows settings and interaction rows", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
  expect_output(print(est), "delta method")
  expect_output(print(est), "Interactions")
})

test_that("estimates of confounded sets carry a note and their classes are flagged", {
  d <- jars_in_rooms(rooms_per_level = 1L)
  des <- describe_design(d, c("temperature", "nutrients"),
                         control = c(temperature = "15", nutrients = "no"),
                         unit = "id", groups = "room")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
  expect_true(all(est$estimates$note == "confounded_with_groups"))
  cls <- classify_interactions(est)
  expect_true(all(cls$confounded))
  expect_output(print(est), "cannot be separated")
})

test_that("combinations without values of the response are named in the reason", {
  d <- two_by_two_gaussian()
  d$y[d$a == "1" & d$b == "1"] <- NA
  fit <- fit_2x2(d, "y", list(gaussian = stats::gaussian()))
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  ab <- est$estimates[est$estimates$component %in% c("combined", "interaction"), ]
  expect_true(all(ab$status == "not_estimable"))
  expect_true(all(ab$reason == "No values of the response in 'a=1 + b=1'."))
})

test_that("a design-level reason comes before a data-level one", {
  # Pure interaction not identifiable (no pairwise cells) and, in addition,
  # no values for one single-stressor cell
  d <- three_way_gaussian()
  d <- d[!(rowSums(d[c("n", "s", "f")] == "1") == 2L), ]
  d$y[d$n == "1" & d$s == "0" & d$f == "0"] <- NA
  des <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                         unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  rows <- est$estimates[est$estimates$component == "interaction", ]
  expect_equal(rows$status[rows$estimand == "pure"], "not_identifiable")
  expect_match(rows$reason[rows$estimand == "pure"], "'n=1 + s=1'", fixed = TRUE)
  expect_equal(rows$status[rows$estimand == "net"], "not_estimable")
  expect_match(rows$reason[rows$estimand == "net"], "No values of the response in 'n=1'",
               fixed = TRUE)
})

test_that("Bliss estimates and delta standard errors match hand calculations", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = c("MU", "BL"), inference = "delta")
  bl <- est$estimates[est$estimates$null_model == "BL", ]
  # Cell means 20, 10, 14 and 5: w = (1, 0.5, 0.7, 0.25)
  w <- c(1, 0.5, 0.7, 0.25)
  expect_equal(bl$estimate, c(-0.5, -0.3, -0.75, 0.25 - 0.35), tolerance = 1e-5)
  expect_equal(bl$scale, rep("response / control", 4))
  int <- bl[bl$component == "interaction", ]
  expect_equal(int$expected, 20 * 0.35, tolerance = 1e-5)
  expect_equal(int$rescaled, -0.1 / 0.35, tolerance = 1e-4)
  expect_equal(int$rescaled_class, "additive")
  # Poisson with a log link: var(eta) = 1 / cell total; gradient of
  # wAB - wA wB with respect to (eta0, etaA, etaB, etaAB)
  g <- c(2 * w[2] * w[3] - w[4], -w[2] * w[3], -w[2] * w[3], w[4])
  expect_equal(int$se, sqrt(sum(g^2 / c(100, 50, 70, 25))), tolerance = 1e-4)
  expect_true(all(is.na(est$estimates$rescaled[est$estimates$null_model == "MU"])))
})

test_that("simulated Bliss and multiplicative p-values coincide for the net interaction", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = c("MU", "BL"), n_sim = 4000, seed = 2)
  e <- est$estimates
  # w_K - prod(w_i) and log(w_K) - sum(log(w_i)) have the same sign in every
  # draw
  expect_equal(e$p_value[e$null_model == "BL"], e$p_value[e$null_model == "MU"])
  expect_output(print(est), "BL = relative to the reference cell")
  expect_output(print(est), "rescaled")
})

test_that("the emergent Bliss interaction of three stressors matches eq. 2.4", {
  d <- three_way_gaussian()
  des <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                         unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "BL", inference = "delta")
  rows <- est$estimates[est$estimates$component == "interaction" & est$estimates$order == 3L, ]
  mu <- tapply(d$y, paste0(d$n, d$s, d$f), mean)
  w <- mu / mu[["000"]]
  e3 <- w[["111"]] - w[["100"]] * w[["011"]] - w[["010"]] * w[["101"]] -
    w[["001"]] * w[["110"]] + 2 * w[["100"]] * w[["010"]] * w[["001"]]
  expect_equal(rows$estimate[rows$estimand == "pure"], e3, tolerance = 1e-6)
  expect_equal(rows$estimate[rows$estimand == "net"],
               w[["111"]] - w[["100"]] * w[["010"]] * w[["001"]], tolerance = 1e-6)
  expect_equal(rows$expected[rows$estimand == "pure"], mu[["111"]] - mu[["000"]] * e3,
               tolerance = 1e-6)
  expect_true(all(is.finite(rows$rescaled)))
})

test_that("the emergent Bliss interaction is not defined for four stressors", {
  grid <- expand.grid(a = c("0", "1"), b = c("0", "1"), c = c("0", "1"), e = c("0", "1"),
                      stringsAsFactors = FALSE)
  d <- grid[rep(seq_len(16), each = 2), ]
  d$y <- 10 + rep(c(-0.5, 0.5), 16) + rowSums(d[c("a", "b", "c", "e")] == "1")
  d$id <- seq_len(nrow(d))
  des <- describe_design(d, c("a", "b", "c", "e"),
                         control = c(a = "0", b = "0", c = "0", e = "0"), unit = "id",
                         max_order = 4)
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("MU", "BL"), inference = "delta")
  four <- est$estimates[est$estimates$order == 4L & est$estimates$component == "interaction", ]
  expect_equal(four$status[four$null_model == "BL" & four$estimand == "pure"], "not_defined")
  expect_equal(four$status[four$null_model == "BL" & four$estimand == "net"], "ok")
  expect_equal(four$status[four$null_model == "MU" & four$estimand == "pure"], "ok")
})

test_that("a combination with only zeros gives a rescaled Bliss interaction of -1", {
  d <- two_by_two_counts()
  d$count[d$a == "1" & d$b == "1"] <- 0
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = c("MU", "BL"), n_sim = 2000, seed = 1)
  int <- est$estimates[est$estimates$component == "interaction", ]
  expect_equal(int$status, c("boundary", "boundary"))
  expect_equal(int$p_value[2], int$p_value[1])
  expect_equal(int$estimate[2], -0.5 * 0.7, tolerance = 1e-5)
  expect_equal(int$rescaled[2], -1)
  expect_equal(int$rescaled_class[2], "synergy")
})

test_that("Bliss independence is not applied to responses where higher means sensitivity", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = c("AD", "BL"), inference = "delta",
                               sensitivity = "count")
  expect_true(all(est$estimates$status[est$estimates$null_model == "BL"] == "not_applicable"))
  expect_true(all(est$estimates$status[est$estimates$null_model == "AD"] == "ok"))
})

test_that("arguments that a must_fit does not use are refused", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  expect_error(estimate_interactions(fit, n_sims = 100), "'n_sims'",
               class = "must_error_input")
  err <- expect_error(estimate_interactions(fit, null_models = "XX"), class = "must_error_input")
  expect_identical(conditionCall(err)[[1L]], as.name("estimate_interactions"))
})

test_that("stress addition estimates match the model of Liess et al. (2016)", {
  d <- two_by_two_counts()
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = c("BL", "SAM"), inference = "delta")
  e <- est$estimates
  sam <- e[e$null_model == "SAM", ]
  expect_equal(sam$estimate[sam$component != "interaction"],
               e$estimate[e$null_model == "BL" & e$component != "interaction"])
  w_sam <- sam_expected(matrix(c(0.5, 0.7), nrow = 1), 3.2)
  int <- sam[sam$component == "interaction", ]
  expect_equal(int$estimate, 0.25 - w_sam, tolerance = 1e-5)
  expect_equal(int$expected, 20 * w_sam, tolerance = 1e-4)
  expect_equal(int$scale, "response / control")
  expect_true(is.finite(int$se) && int$se > 0)
  expect_true(is.na(int$rescaled))
  # Another shape of the beta distribution
  flatter <- estimate_interactions(fit, null_models = "SAM", inference = "delta", sam_shape = 1)
  # With a uniform distribution (shape 1) the general stresses are the mortalities
  expect_equal(flatter$estimates$estimate[flatter$estimates$component == "interaction"],
               0.25 - (1 - 0.5 - 0.3), tolerance = 1e-4)
  expect_equal(flatter$settings$sam_shape, 1)
  expect_error(estimate_interactions(fit, null_models = "SAM", sam_shape = 0),
               class = "must_error_input")
  expect_output(print(est), "SAM = relative to the reference cell")
})

test_that("stress addition needs stressors that reduce the response", {
  d <- two_by_two_counts()
  d$count[d$a == "1" & d$b == "0"] <- d$count[d$a == "1" & d$b == "0"] * 3
  fit <- fit_2x2(d, "count", list(poisson = stats::poisson()))
  est <- estimate_interactions(fit, null_models = "SAM", inference = "delta")
  int <- est$estimates[est$estimates$component == "interaction", ]
  expect_equal(int$status, "not_applicable")
  expect_match(int$reason, "'a=1' increases it")
})

test_that("concentration addition is not defined without dose-response curves", {
  d <- two_by_two_counts()
  d$a <- ifelse(d$a == "1", "5", "0")
  d$b <- ifelse(d$b == "1", "2", "0")
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "CA", inference = "delta")
  e <- est$estimates
  expect_equal(e$status[e$component == "interaction"], "not_defined")
  expect_match(e$reason[e$component == "interaction"], "at least two doses of 'a' alone")
  expect_true(all(e$status[e$component != "interaction"] == "ok"))
  text <- two_by_two_counts()
  text$a <- ifelse(text$a == "1", "high", "none")
  des2 <- describe_design(text, c("a", "b"), control = c(a = "none", b = "0"), unit = "id")
  fit2 <- fit_responses(text, "count", des2, families = list(poisson = stats::poisson()),
                        residual_checks = FALSE, quiet = TRUE)
  est2 <- estimate_interactions(fit2, null_models = "CA", inference = "delta")
  expect_match(est2$estimates$reason[est2$estimates$component == "interaction"],
               "levels of 'a' are not numbers")
})

test_that("concentration addition finds no interaction in data that follow it", {
  d <- ca_factorial()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("BL", "CA"), inference = "delta")
  int <- est$estimates[est$estimates$component == "interaction", ]
  ca <- int[int$null_model == "CA", ]
  expect_equal(nrow(ca), 9L)
  expect_true(all(ca$status == "ok"))
  expect_true(all(abs(ca$estimate) < 1e-5))
  expect_true(all(ca$se > 0))
  # Bliss independence does detect departures from it
  expect_gt(max(abs(int$estimate[int$null_model == "BL"])), 0.05)
  # Expected = y0 * w_CA; the simulation agrees with the delta method
  expect_equal(ca$expected, ca$observed, tolerance = 1e-4)
  sim <- estimate_interactions(fit, null_models = "CA", n_sim = 4000, seed = 2)
  sim_int <- sim$estimates[sim$estimates$component == "interaction", ]
  expect_equal(sim_int$se, ca$se, tolerance = 0.1)
})

test_that("concentration addition uses every dose tested alone, as computed by hand", {
  ray <- data.frame(f = c(0, 2, 20, 200, 0, 0, 0, 2, 20, 200),
                    i = c(0, 0, 0, 0, 2, 20, 200, 2, 20, 200))
  ray <- ray[rep(seq_len(nrow(ray)), each = 3), ]
  ray$y <- round(60 / (1 + (ray$f / 30)^1.2) / (1 + ray$i / 60) * rep(c(0.9, 1, 1.1), 10))
  ray$id <- seq_len(nrow(ray))
  des <- describe_design(ray, c("f", "i"), control = c(f = "0", i = "0"), unit = "id")
  fit <- fit_responses(ray, "y", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "CA", inference = "delta")
  int <- est$estimates[est$estimates$component == "interaction", ]
  expect_equal(int$levels, c("f=2 | i=2", "f=20 | i=20", "f=200 | i=200"))
  mu <- tapply(ray$y, paste(ray$f, ray$i), mean)
  w_f <- mu[c("2 0", "20 0", "200 0")] / mu[["0 0"]]
  w_i <- mu[c("0 2", "0 20", "0 200")] / mu[["0 0"]]
  curves <- list(fit_loglogistic(matrix(w_f, nrow = 1), log(c(2, 20, 200))),
                 fit_loglogistic(matrix(w_i, nrow = 1), log(c(2, 20, 200))))
  hand <- as.numeric(mu[c("2 2", "20 20", "200 200")]) / mu[["0 0"]] -
    vapply(c(2, 20, 200), function(x) ca_expected(curves, c(x, x)), numeric(1))
  expect_equal(int$estimate, hand, tolerance = 1e-4)
})

test_that("a combination with only zeros is tested against concentration addition", {
  d <- ca_factorial(a = c(0, 1, 4), b = c(0, 2, 10))
  d$y <- round(d$y / 4)
  d$y[d$a == 4 & d$b == 10] <- 0
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "y", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("SAM", "CA"), n_sim = 2000, seed = 1)
  int <- est$estimates[est$estimates$component == "interaction" &
                         est$estimates$levels == "a=4 | b=10", ]
  expect_equal(int$status, c("boundary", "boundary"))
  expect_true(all(int$observed == 0))
  expect_true(all(is.finite(int$p_value)))
  expect_lt(int$expected[int$null_model == "CA"], 25 / 4)
})

test_that("stress addition and concentration addition give no pure interaction", {
  d <- three_way_gaussian()
  des <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                         unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("SAM", "CA"), inference = "delta")
  three <- est$estimates[est$estimates$order == 3L & est$estimates$component == "interaction", ]
  expect_equal(three$status[three$estimand == "pure"], c("not_defined", "not_defined"))
  expect_equal(three$status[three$null_model == "SAM" & three$estimand == "net"], "ok")
  mu <- tapply(d$y, paste0(d$n, d$s, d$f), mean)
  w <- mu / mu[["000"]]
  expect_equal(three$estimate[three$null_model == "SAM" & three$estimand == "net"],
               w[["111"]] - sam_expected(matrix(w[c("100", "010", "001")], nrow = 1), 3.2),
               tolerance = 1e-6)
})

test_that("a control shared by two `by` groups enters the sets of both", {
  d <- shared_regimes()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id",
                         by = "regime", shared = c(regime = "both"))
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  ad <- est$estimates[est$estimates$component == "interaction", ]
  expect_setequal(ad$regime, c("pulse", "press"))
  mu <- tapply(d$y, paste(d$regime, d$a, d$b), mean)
  for (r in c("pulse", "press")) {
    expected <- mu[[paste(r, "1 1")]] - mu[[paste(r, "1 0")]] - mu[[paste(r, "0 1")]] +
      mu[["both 0 0"]]
    expect_equal(ad$estimate[ad$regime == r], expected, tolerance = 1e-6)
  }
  # One residual variance: 28 units in 7 cells
  expect_equal(ad$df, c(21, 21))
  # One variance per `by` group: each set uses its group and the shared rows
  by_group <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                            dispersion = "by", residual_checks = FALSE, quiet = TRUE)
  ad_by <- interaction_rows(estimate_interactions(by_group, null_models = "AD",
                                                  inference = "delta"), "AD")
  expect_equal(ad_by$df, c(12, 12))
  expect_equal(ad_by$estimate, ad$estimate, tolerance = 1e-6)
})

test_that("sets that involve a stressor are not estimated before it was applied", {
  d <- three_way_repeated()
  ctrl <- c(n = "0", s = "0", f = "0")
  des <- describe_design(d, c("n", "s", "f"), control = ctrl, unit = "tank", time = "day",
                         baseline = list(f = 0))
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta", overall = FALSE)
  int <- est$estimates[est$estimates$component == "interaction", ]
  expect_setequal(unique(int$stressors[int$time == "0"]), "n x s")
  expect_setequal(unique(int$stressors[int$time == "10"]),
                  c("n x s", "n x f", "s x f", "n x s x f"))
  # With other stressors as background, a set with f away from control waits too
  all_bg <- describe_design(d, c("n", "s", "f"), control = ctrl, unit = "tank", time = "day",
                            baseline = list(f = 0), background = "all")
  fit_bg <- fit_responses(d, "y", all_bg, families = list(gaussian = stats::gaussian()),
                          residual_checks = FALSE, quiet = TRUE)
  bg <- estimate_interactions(fit_bg, null_models = "AD", inference = "delta", overall = FALSE)
  day0 <- bg$estimates[bg$estimates$component == "interaction" & bg$estimates$time == "0", ]
  expect_true(all(day0$stressors == "n x s" & day0$background == ""))
  # With additive time, f is coded at its control on day 0
  frame <- design_model_frame(d, des)
  frame <- baseline_as_control(frame, des)
  expect_true(all(grepl("f=0", frame$.cell[frame$day == 0], fixed = TRUE)))
  expect_identical(frame$.cell[frame$day == 10], design_model_frame(d, des)$.cell[d$day == 10])
  additive <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                            time_structure = "additive", residual_checks = FALSE, quiet = TRUE)
  add_est <- estimate_interactions(additive, null_models = "AD", inference = "delta",
                                   overall = FALSE)
  expect_true(any(add_est$estimates$time == "0"))
})
