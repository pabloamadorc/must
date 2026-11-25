# Resampling tests of AD and MU (R/resampling.R, R/interaction-resampling.R)
# and the resampling inference of estimate_interactions() and
# estimate_community(): deviances, fits and distributions against stats and
# MASS, constrained additive fits against brute force, and the routes.

family_names <- c("poisson", "nbinom2", "nbinom1", "gaussian", "gamma", "lognormal", "tweedie")

make_family <- function(name) {
  rs_family(name, power = 1.5, theta = 3)
}

test_that("unit deviances are those of the GLM families", {
  y <- c(0, 1, 3, 7.5, 12)
  mu <- c(0.5, 2, 3, 6, 15)
  pos <- y > 0
  expect_equal(rs_unit_deviance(y, mu, make_family("poisson")),
               stats::poisson()$dev.resids(y, mu, 1))
  expect_equal(rs_unit_deviance(y, mu, make_family("gaussian")),
               stats::gaussian()$dev.resids(y, mu, 1))
  expect_equal(rs_unit_deviance(y[pos], mu[pos], make_family("gamma")),
               stats::Gamma()$dev.resids(y[pos], mu[pos], 1))
  expect_equal(rs_unit_deviance(y, mu, make_family("tweedie")),
               tweedie_unit_deviance(y, mu, 1.5))
  skip_if_not_installed("MASS")
  expect_equal(rs_unit_deviance(y, mu, make_family("nbinom2")),
               MASS::negative.binomial(3)$dev.resids(y, mu, 1))
})

test_that("rates of the full model are the maximum likelihood cell means", {
  set.seed(3)
  group <- factor(rep(1:3, each = 5))
  e <- stats::runif(15, 0.5, 2)
  y <- stats::rpois(15, 4 * e * c(1, 2, 3)[group])
  for (name in c("poisson", "gaussian", "gamma", "tweedie")) {
    fam <- make_family(name)
    yy <- if (name == "gamma") y + 0.5 else y
    fit <- stats::glm(yy ~ 0 + group + offset(log(e)), family = rs_glm_family(fam),
                      control = stats::glm.control(epsilon = 1e-12, maxit = 100))
    expect_equal(rs_rates(yy, e, group, fam), unname(exp(stats::coef(fit))), tolerance = 1e-6)
  }
  skip_if_not_installed("MASS")
  fit <- stats::glm(y ~ 0 + group + offset(log(e)), family = MASS::negative.binomial(3),
                    control = stats::glm.control(epsilon = 1e-12, maxit = 100))
  expect_equal(rs_rates(y, e, group, make_family("nbinom2")),
               unname(exp(stats::coef(fit))), tolerance = 1e-6)
})

test_that("PIT residuals and quantiles invert each other", {
  set.seed(4)
  mu <- c(0.4, 3, 12, 40)
  for (name in family_names) {
    fam <- make_family(name)
    disp <- switch(name, gaussian = 4, gamma = 0.3, lognormal = 0.3, tweedie = 1.2,
                   nbinom1 = 0.5, NA_real_)
    u <- stats::runif(4, 0.05, 0.95)
    q <- rs_quantile(u, mu, disp, fam)
    u2 <- rs_pit(q, mu, disp, fam)
    if (fam$discrete) {
      expect_equal(q, round(q))
      # The randomised residual of q falls in the step of the distribution
      # function that contains u.
      expect_equal(rs_quantile(u2, mu, disp, fam), q)
    } else if (identical(name, "tweedie")) {
      expect_equal(u2[q > 0], u[q > 0], tolerance = 1e-8)
    } else {
      expect_equal(u2, u, tolerance = 1e-8)
    }
  }
})

test_that("the multiplicative null model is the log-link GLM with main effects", {
  d <- two_by_two_counts()
  env <- data.frame(.s1 = factor(d$a == "1", levels = c(FALSE, TRUE)),
                    .s2 = factor(d$b == "1", levels = c(FALSE, TRUE)))
  layout <- rs_layout(env, c(".s1", ".s2"), all_times = FALSE)
  expect_equal(layout$role, ifelse(d$a == "1", ifelse(d$b == "1", 3L, 1L),
                                   ifelse(d$b == "1", 2L, 0L)))
  fit <- rs_null(d$count, rep(1, 20), layout, make_family("poisson"), "MU")
  glm_fit <- stats::glm(count ~ a + b, family = stats::poisson(), data = d)
  expect_equal(fit$mu, unname(stats::fitted(glm_fit)), tolerance = 1e-6)
  expect_equal(fit$deviance, stats::deviance(glm_fit), tolerance = 1e-6)
})

test_that("the additive null model is the identity-link GLM, or its boundary", {
  d <- two_by_two_counts()
  role <- ifelse(d$a == "1", ifelse(d$b == "1", 3L, 1L), ifelse(d$b == "1", 2L, 0L))
  # Interior: the Poisson GLM with an identity link.
  fit <- rs_ad_time(d$count, rep(1, 20), role, 2L, make_family("poisson"), truncate = FALSE)
  glm_fit <- stats::glm(count ~ a + b, family = stats::poisson(link = "identity"), data = d,
                        start = c(20, -10, -6))
  expect_equal(fit$mu, unname(stats::fitted(glm_fit)), tolerance = 1e-5)
  # Gaussian: least squares, also below zero.
  g <- two_by_two_gaussian(means = c(10, 3, 4, 0))
  role_g <- ifelse(g$a == "1", ifelse(g$b == "1", 3L, 1L), ifelse(g$b == "1", 2L, 0L))
  fit_g <- rs_ad_time(g$y, rep(1, 20), role_g, 2L, make_family("gaussian"), truncate = FALSE)
  expect_equal(fit_g$mu, unname(stats::fitted(stats::lm(y ~ a + b, data = g))))
  # Boundary: AD expects 2 + 3 - 20 < 0 for a combination that is absent.
  y <- c(19, 21, 20, 20, 2, 1, 3, 2, 3, 2, 4, 3, 0, 0, 0, 0)
  r <- rep(0:3, each = 4)
  for (name in c("poisson", "nbinom2", "tweedie")) {
    fam <- make_family(name)
    fit <- rs_ad_time(y, rep(1, 16), r, 2L, fam, truncate = FALSE)
    brute <- stats::optim(c(20, 2, 3), function(p) {
      m <- c(p[[1L]], p[[2L]], p[[3L]], p[[2L]] + p[[3L]] - p[[1L]])
      if (any(m < 0)) {
        return(1e10)
      }
      sum(rs_unit_deviance(y, pmax(m[r + 1L], 1e-300), fam))
    }, control = list(maxit = 5000, reltol = 1e-14))
    expect_lte(sum(rs_unit_deviance(y, fit$mu, fam)), brute$value + 1e-6)
    expect_equal(fit$mu[r == 3L], rep(0, 4))
    # Truncated at zero, the null model fits the absence.
    trunc <- rs_ad_time(y, rep(1, 16), r, 2L, fam, truncate = TRUE)
    expect_equal(trunc$mu, rep(unname(tapply(y, r, mean)), each = 4))
  }
})

test_that("resampling tests are reproducible and have the size of their statistic", {
  d <- two_by_two_counts()
  env <- data.frame(.s1 = factor(d$a == "1", levels = c(FALSE, TRUE)),
                    .s2 = factor(d$b == "1", levels = c(FALSE, TRUE)))
  layout <- rs_layout(env, c(".s1", ".s2"), all_times = FALSE)
  unit <- factor(d$id)
  Y <- matrix(d$count, ncol = 1L, dimnames = list(NULL, "count"))
  one <- withr::with_seed(1, rs_test(Y, NULL, layout, unit, list(make_family("poisson")), "AD",
                                     n_boot = 99))
  again <- withr::with_seed(1, rs_test(Y, NULL, layout, unit, list(make_family("poisson")),
                                       "AD", n_boot = 99))
  expect_identical(one, again)
  expect_equal(one$df, 1)
  expect_true(one$p_value > 0 && one$p_value <= 1)
  # The drop in deviance divided by the dispersion estimated from the
  # deviance of the null model: df0 (D0 - D1) / D0, also for counts.
  cells <- stats::ave(d$count, d$a, d$b)
  null <- rs_null(d$count, rep(1, 20), layout, make_family("poisson"), "AD")
  d1 <- sum(stats::poisson()$dev.resids(d$count, cells, 1))
  expect_equal(one$statistic, 17 * (null$deviance - d1) / null$deviance, tolerance = 1e-6)
  g <- two_by_two_gaussian()
  Yg <- matrix(g$y, ncol = 1L)
  res <- withr::with_seed(1, rs_test(Yg, NULL, layout, factor(g$id),
                                     list(make_family("gaussian")), "AD", n_boot = 99))
  rss0 <- sum(stats::residuals(stats::lm(y ~ a + b, data = g))^2)
  rss1 <- sum(stats::residuals(stats::lm(y ~ a * b, data = g))^2)
  expect_equal(res$statistic, 17 * (rss0 - rss1) / rss0)
})

# estimate_interactions(inference = "resampling") -------------------------------

test_that("resampling replaces the tests of the net interactions under AD and MU", {
  d <- community_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank",
                         time = "day")
  fit <- fit_responses(d, c("sp1", "sp3"), des, residual_checks = FALSE, quiet = TRUE)
  pei <- estimate_interactions(fit, null_models = c("AD", "MU"), n_sim = 1000, seed = 1)
  res <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "resampling",
                               n_sim = 1000, n_boot = 99, seed = 1)
  est <- res$estimates
  net <- est$component == "interaction" & est$status == "ok"
  expect_true(all(est$inference[net] == "resampling"))
  expect_equal(est$inference[!net], pei$estimates$inference[!net])
  expect_equal(est$estimate, pei$estimates$estimate)
  expect_equal(est$lower, pei$estimates$lower)
  expect_true(all(is.na(est$df[net])))
  expect_true(all(grepl("joint_over_times", est$note[net & est$time == "overall"])))
  # Single effects keep the post-estimation tests.
  single <- est$component == "single"
  expect_equal(est$p_value[single], pei$estimates$p_value[single])
  # sp3 has no individuals with both stressors: the post-estimation test is a
  # boundary one; the resampling test fits the zeros.
  sp3 <- est$response == "sp3" & est$component == "interaction" & est$time == "7"
  expect_equal(est$status[sp3], c("ok", "ok"))
  expect_equal(pei$estimates$status[sp3], c("boundary", "boundary"))
  expect_equal(res$settings$inference, "resampling")
  expect_equal(res$settings$n_boot, 99L)
  expect_output(print(res), "tested by resampling the units (99 resamples)", fixed = TRUE)
  again <- estimate_interactions(fit, null_models = c("AD", "MU"), inference = "resampling",
                                 n_sim = 1000, n_boot = 99, seed = 1)
  expect_identical(again$estimates, res$estimates)
  expect_s3_class(classify_interactions(res), "data.frame")
})

test_that("resampling tests use the family of the fitted model", {
  d <- community_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank",
                         time = "day")
  fit <- fit_responses(d, "sp1", des, residual_checks = FALSE, quiet = TRUE,
                       families = list(sp1 = list(nbinom2 = glmmTMB::nbinom2())))
  res <- estimate_interactions(fit, null_models = "MU", inference = "resampling",
                               n_sim = 1000, n_boot = 99, seed = 2, overall = FALSE)
  day7 <- d[d$day == 7, ]
  row <- res$estimates$component == "interaction" & res$estimates$time == "7"
  theta <- unname(stats::sigma(fit$fits$sp1$model))
  env <- data.frame(.s1 = factor(day7$a == "1", levels = c(FALSE, TRUE)),
                    .s2 = factor(day7$b == "1", levels = c(FALSE, TRUE)))
  layout <- rs_layout(env, c(".s1", ".s2"), all_times = FALSE)
  direct <- withr::with_seed(2, {
    rs_test(matrix(day7$sp1, ncol = 1L), NULL, layout, factor(day7$tank),
            list(rs_family("nbinom2", theta = theta)), "MU", n_boot = 99)
  })
  expect_true(is.finite(res$estimates$p_value[row]))
  expect_true(direct$p_value > 0 && direct$p_value <= 1)
})

test_that("resampling says what it cannot test", {
  d <- community_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank",
                         time = "day")
  fit <- fit_responses(d, "sp1", des, residual_checks = FALSE, quiet = TRUE)
  expect_error(estimate_interactions(fit, null_models = c("AD", "DO"), inference = "resampling"),
               "DO", class = "must_error_input")
  expect_error(estimate_interactions(fit, inference = "resampling", n_boot = 10),
               class = "must_error_input")
  additive <- fit_responses(d, "sp1", des, residual_checks = FALSE, quiet = TRUE,
                            time_structure = "additive")
  expect_error(estimate_interactions(additive, inference = "resampling"), "additive",
               class = "must_error_input")
  # Proportions keep the post-estimation test.
  jars <- survival_jars()
  jdes <- describe_design(jars, c("a", "b"), control = c(a = "0", b = "0"))
  jfit <- fit_responses(jars, "alive", jdes, trials = c(alive = "n"),
                        residual_checks = FALSE, quiet = TRUE)
  res <- estimate_interactions(jfit, null_models = "AD", inference = "resampling",
                               n_sim = 1000, n_boot = 99, seed = 1)
  net <- res$estimates[res$estimates$component == "interaction", ]
  expect_equal(net$inference, "simulation")
  expect_match(net$note, "resampling_not_available")
  expect_match(net$reason, "proportions")
  # A stressor applied to whole rooms: the rooms, not the jars, are its
  # replicates.
  rooms <- jars_in_rooms(rooms_per_level = 2L)
  rooms$y2 <- withr::with_seed(5, rooms$y + stats::rnorm(nrow(rooms), sd = 0.5))
  rdes <- describe_design(rooms, c("temperature", "nutrients"),
                          control = c(temperature = "15", nutrients = "no"),
                          unit = "id", groups = "room")
  rfit <- fit_responses(rooms, c("y", "y2"), rdes, families = list(gaussian = stats::gaussian()),
                        residual_checks = FALSE, quiet = TRUE)
  rres <- estimate_interactions(rfit, inference = "resampling", n_sim = 1000, n_boot = 99,
                                seed = 1)
  rnet <- rres$estimates[rres$estimates$component == "interaction", ]
  expect_equal(rnet$inference, c("simulation", "simulation"))
  expect_match(rnet$reason, "temperature' is applied to whole groups of room", fixed = TRUE)
  com <- estimate_community(rfit, properties = "abundance", inference = "resampling",
                            n_boot = 99, seed = 1)
  expect_true(all(com$estimates$status == "not_estimable"))
  expect_match(com$estimates$reason, "whole groups of room")
  # A model fitted outside must.
  model <- stats::glm(count ~ a * b, family = stats::poisson(), data = two_by_two_counts())
  ddes <- describe_design(two_by_two_counts(), c("a", "b"), control = c(a = "0", b = "0"))
  expect_error(estimate_interactions(model, ddes, inference = "resampling"), "fit_responses",
               class = "must_error_input")
})

# estimate_community(inference = "resampling") -----------------------------------

test_that("community resampling uses the means of the samples and resampled units", {
  d <- community_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank",
                         time = "day")
  fit <- fit_responses(d, c("sp1", "sp2", "sp4"), des, residual_checks = FALSE, quiet = TRUE)
  res <- estimate_community(fit, inference = "resampling", n_boot = 199, seed = 1)
  est <- res$estimates
  expect_true(all(est$inference == "resampling"))
  expect_equal(res$settings$inference, "resampling")
  expect_equal(res$settings$taxa_correlation, "resampled")
  day7 <- d[d$day == 7, ]
  cell <- paste(day7$a, day7$b)
  total <- tapply(day7$sp1 + day7$sp2 + day7$sp4, cell, mean)
  row <- est$response == "abundance" & est$time == "7" & est$component == "combined" &
    est$null_model == "AD"
  expect_equal(est$estimate[row], unname(total[["1 1"]] / total[["0 0"]] - 1))
  ok <- est$status == "ok" & est$component != "expected"
  expect_true(all(est$p_value[ok] > 0 & est$p_value[ok] <= 1))
  expect_true(all(est$lower[ok] <= est$estimate[ok] + 1e-12 |
                    est$upper[ok] >= est$estimate[ok] - 1e-12))
  again <- estimate_community(fit, inference = "resampling", n_boot = 199, seed = 1)
  expect_identical(again$estimates, res$estimates)
  expect_output(print(res), "resampling of whole units")
  # Units sampled a different number of times cannot be resampled whole.
  uneven <- d[!(d$tank == "T1" & d$day == 14), ]
  ufit <- fit_responses(uneven, c("sp1", "sp2"), describe_design(
    uneven, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank", time = "day"
  ), residual_checks = FALSE, quiet = TRUE)
  ures <- estimate_community(ufit, properties = "abundance", inference = "resampling",
                             n_boot = 99, seed = 1)
  expect_true(all(ures$estimates$status == "not_estimable"))
  expect_match(ures$estimates$reason[[1L]], "every sampling time")
})
