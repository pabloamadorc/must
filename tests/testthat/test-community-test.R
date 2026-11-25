# test_community(): the test of the multiplicative null model (MU) on every
# taxon of a community, with multivariate GLMs of counts (mvabund) or with
# Tweedie models of densities; against glm(), the cell means and the series
# of the Tweedie distribution; resampling, statuses and errors.

community_design <- function(d) {
  describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank", time = "day")
}

test_that("the community statistic adds up the likelihood ratios of the taxa", {
  skip_if_not_installed("mvabund")
  d <- community_counts()
  taxa <- c("sp1", "sp2", "sp3", "sp4", "sp5")
  out <- test_community(d, taxa, community_design(d), n_boot = 99, seed = 1, quiet = TRUE)
  expect_s3_class(out, "must_community_test")
  tests <- out$tests
  expect_equal(tests$time, c("7", "14", "overall"))
  expect_true(all(tests$status == "ok"))
  expect_equal(tests$null_model, rep("MU", 3))
  expect_equal(tests$df, c(1, 1, 2))
  expect_equal(tests$n_taxa, rep(5L, 3))
  for (tt in tests$time) {
    rows <- out$taxa$time == tt
    expect_equal(tests$statistic[tests$time == tt], sum(out$taxa$statistic[rows]),
                 tolerance = 1e-6)
  }
  expect_true(all(tests$p_value > 0 & tests$p_value <= 1))
  expect_equal(tests$p_adjusted[1:2], stats::p.adjust(tests$p_value[1:2], "BH"))
  expect_equal(tests$p_adjusted[[3L]], tests$p_value[[3L]])
  # Net interactions under MU: log ratios of the cell means.
  day7 <- d[d$day == 7, ]
  means <- function(taxon) tapply(day7[[taxon]], paste(day7$a, day7$b), mean)
  est <- out$taxa[out$taxa$time == "7", ]
  expect_equal(unique(est$null_model), "MU")
  expect_equal(unique(est$scale), "log(response)")
  m <- means("sp1")
  expect_equal(est$estimate[est$taxon == "sp1"],
               unname(log(m[["1 1"]]) - log(m[["1 0"]]) - log(m[["0 1"]]) + log(m[["0 0"]])),
               tolerance = 1e-5)
  # sp3 has no individuals with both stressors and sp5 none in the reference.
  expect_true(is.na(est$estimate[est$taxon == "sp3"]))
  expect_true(is.na(est$estimate[est$taxon == "sp5"]))
  expect_true(all(out$taxa$p_adjusted > 0 & out$taxa$p_adjusted <= 1))
  expect_output(print(out), "multiplicative \\(MU\\)")
  expect_null(out$diagnostics$power)
  again <- test_community(d, taxa, community_design(d), n_boot = 99, seed = 1, quiet = TRUE)
  expect_identical(again$tests, out$tests)
  expect_equal(out$settings$seed, 1L)
  expect_equal(out$settings$null_model, "MU")
})

test_that("with Poisson counts each likelihood ratio is that of glm()", {
  skip_if_not_installed("mvabund")
  d <- community_counts()
  out <- test_community(d, c("sp2", "sp4"), community_design(d), family = "poisson",
                        overall = FALSE, n_boot = 99, seed = 1, quiet = TRUE)
  expect_equal(out$tests$time, c("7", "14"))
  day7 <- d[d$day == 7, ]
  day7$cell <- paste(day7$a, day7$b)
  for (taxon in c("sp2", "sp4")) {
    y <- day7[[taxon]]
    null <- stats::glm(y ~ a + b, family = stats::poisson(), data = day7)
    full <- stats::glm(y ~ cell, family = stats::poisson(), data = day7)
    lr <- stats::deviance(null) - stats::deviance(full)
    row <- out$taxa$time == "7" & out$taxa$taxon == taxon
    expect_equal(out$taxa$statistic[row], lr, tolerance = 1e-4)
  }
})

test_that("the sampling effort enters as an offset and gives densities", {
  skip_if_not_installed("mvabund")
  d <- community_counts()
  d$volume <- rep(c(1, 2, 0.5, 1.5), length.out = nrow(d))
  des <- community_design(d)
  out <- test_community(d, c("sp1", "sp2", "sp4"), des, effort = "volume", overall = FALSE,
                        n_boot = 99, seed = 1, quiet = TRUE)
  day7 <- d[d$day == 7, ]
  day7$cell <- factor(paste(day7$a, day7$b), levels = c("0 0", "1 0", "0 1", "1 1"))
  p_out <- test_community(d, c("sp1", "sp2"), des, family = "poisson", effort = "volume",
                          overall = FALSE, n_boot = 99, seed = 1, quiet = TRUE)
  for (taxon in c("sp1", "sp2")) {
    y <- day7[[taxon]]
    null <- stats::glm(y ~ a + b + offset(log(volume)), family = stats::poisson(), data = day7)
    full <- stats::glm(y ~ 0 + cell + offset(log(volume)), family = stats::poisson(),
                       data = day7)
    row <- p_out$taxa$time == "7" & p_out$taxa$taxon == taxon
    expect_equal(p_out$taxa$statistic[row], stats::deviance(null) - stats::deviance(full),
                 tolerance = 1e-4)
    rate <- exp(stats::coef(full))
    expect_equal(p_out$taxa$estimate[row],
                 unname(log(rate[["cell1 1"]]) - log(rate[["cell1 0"]]) -
                          log(rate[["cell0 1"]]) + log(rate[["cell0 0"]])),
                 tolerance = 1e-5)
  }
  expect_output(print(out), "offset log\\(effort\\)")
  expect_equal(out$settings$effort, rep("volume", 3))
  # A constant effort changes nothing.
  d$same <- 2
  plain <- test_community(d, c("sp1", "sp2"), des, overall = FALSE, n_boot = 99, seed = 1,
                          quiet = TRUE)
  constant <- test_community(d, c("sp1", "sp2"), des, effort = "same", overall = FALSE,
                             n_boot = 99, seed = 1, quiet = TRUE)
  expect_equal(constant$tests$statistic, plain$tests$statistic, tolerance = 1e-6)
  expect_equal(constant$taxa$estimate, plain$taxa$estimate, tolerance = 1e-6)
  # One effort per taxon.
  d$other <- 1
  per_taxon <- test_community(d, c("sp1", "sp2"), des, effort = c("volume", "other"),
                              overall = FALSE, n_boot = 99, seed = 1, quiet = TRUE)
  alone <- test_community(d, c("sp1", "sp4"), des, effort = c("volume", "volume"),
                          overall = FALSE, n_boot = 99, seed = 1, quiet = TRUE)
  expect_equal(per_taxon$taxa$statistic[per_taxon$taxa$taxon == "sp1"],
               alone$taxa$statistic[alone$taxa$taxon == "sp1"], tolerance = 1e-6)
  expect_equal(per_taxon$taxa$statistic[per_taxon$taxa$taxon == "sp2"],
               plain$taxa$statistic[plain$taxa$taxon == "sp2"], tolerance = 1e-6)
})

test_that("tests that cannot be done say why", {
  skip_if_not_installed("mvabund")
  d <- community_counts()
  d$empty <- 0
  out <- test_community(d, c("sp1", "sp2", "empty"), community_design(d), n_boot = 99,
                        seed = 1, quiet = TRUE)
  expect_equal(out$tests$n_taxa, rep(2L, 3))
  expect_setequal(out$diagnostics$excluded$taxon, "empty")
  expect_match(out$diagnostics$excluded$reason[[1L]], "Not found")
  # A tank sampled once: the units cannot be resampled whole over the times.
  uneven <- d[!(d$tank == "T1" & d$day == 14), ]
  out <- test_community(uneven, c("sp1", "sp2"), community_design(uneven), n_boot = 99,
                        seed = 1, quiet = TRUE)
  overall <- out$tests[out$tests$time == "overall", ]
  expect_equal(overall$status, "not_estimable")
  expect_match(overall$reason, "every unit sampled at every time")
  expect_true(all(out$tests$status[out$tests$time != "overall"] == "ok"))
  # A single tank in a combination at a time.
  few <- d[!(d$a == "1" & d$b == "1" & d$day == 7 & d$tank != "T4"), ]
  out <- test_community(few, c("sp1", "sp2"), community_design(few), overall = FALSE,
                        n_boot = 99, seed = 1, quiet = TRUE)
  expect_equal(out$tests$status, c("not_estimable", "ok"))
  expect_match(out$tests$reason[[1L]], "fewer than two units")
  # A stressor applied to whole racks of tanks: the racks are its replicates.
  racks <- d
  racks$rack <- paste0("R", racks$a, "-", ifelse(racks$rep <= 2L, 1L, 2L))
  rdes <- describe_design(racks, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank",
                          time = "day", groups = "rack")
  out <- test_community(racks, c("sp1", "sp2"), rdes, overall = FALSE, n_boot = 99, seed = 1,
                        quiet = TRUE)
  expect_equal(out$tests$status, rep("not_estimable", 2))
  expect_match(out$tests$reason, "whole groups of rack", fixed = TRUE)
})

# Tweedie ----------------------------------------------------------------------

test_that("the Tweedie series agree with the tweedie package", {
  skip_if_not_installed("tweedie")
  cases <- list(c(5, 2, 1.5), c(0.3, 1, 1.3), c(200, 0.5, 1.6), c(50, 0.2, 1.1), c(2, 5, 1.9))
  for (cs in cases) {
    mu <- cs[[1L]]
    phi <- cs[[2L]]
    power <- cs[[3L]]
    y <- c(0, stats::qgamma(c(0.05, 0.3, 0.6, 0.9, 0.99), shape = 2, scale = mu / 2))
    expect_equal(ptweedie_series(y, mu, phi, power),
                 tweedie::ptweedie(y, mu = mu, phi = phi, power = power), tolerance = 1e-8)
    u <- c(0.02, 0.3, 0.5, 0.8, 0.97)
    expect_equal(qtweedie_series(u, mu, phi, power),
                 tweedie::qtweedie(u, mu = mu, phi = phi, power = power), tolerance = 1e-6)
  }
})

test_that("the Tweedie quantile inverts the distribution function", {
  u <- c(0.001, 0.2, 0.5, 0.9, 0.999)
  for (power in c(1.05, 1.4, 1.95)) {
    q <- qtweedie_series(u, mu = 3, phi = 1.5, power = power)
    p0 <- exp(-3^(2 - power) / (1.5 * (2 - power)))
    expect_true(all(q[u <= p0] == 0))
    expect_equal(ptweedie_series(q[u > p0], 3, 1.5, power), u[u > p0], tolerance = 1e-8)
  }
  # Beyond 1000 expected amounts, the gamma with the same mean and variance.
  q <- qtweedie_series(u, mu = 5000, phi = 0.01, power = 1.5)
  expect_equal(q, stats::qgamma(u, shape = 5000^0.5 / 0.01, scale = 0.01 * 5000^0.5))
  expect_equal(ptweedie_series(q, 5000, 0.01, 1.5), u)
  expect_equal(tweedie_unit_deviance(c(0, 2), c(1, 2), 1.5), c(2 / 0.5, 0))
})

test_that("the Tweedie family of glm.fit() is that of statmod", {
  skip_if_not_installed("statmod")
  d <- community_counts()
  d <- d[d$day == 7, ]
  ours <- stats::glm(mass ~ a + b, family = tweedie_glm_family(1.4), data = d)
  theirs <- stats::glm(mass ~ a + b, family = statmod::tweedie(var.power = 1.4, link.power = 0),
                       data = d)
  expect_equal(stats::coef(ours), stats::coef(theirs))
  expect_equal(stats::deviance(ours), stats::deviance(theirs))
})

test_that("Tweedie tests scale the drop in deviance of each taxon", {
  d <- community_counts()
  d$sp1 <- d$sp1 / 5
  d$sp2 <- d$sp2 / 5
  d$sp3 <- d$sp3 / 5
  taxa <- c("sp1", "sp2", "sp3", "mass")
  out <- test_community(d, taxa, community_design(d), family = "tweedie", n_boot = 99,
                        seed = 1, quiet = TRUE)
  tests <- out$tests
  expect_true(all(tests$status == "ok"))
  expect_equal(tests$df, c(1, 1, 2))
  power <- out$diagnostics$power
  expect_equal(power$taxon, taxa)
  expect_true(all(power$power > 1 & power$power < 2))
  day7 <- d[d$day == 7, ]
  day7$cell <- paste(day7$a, day7$b)
  for (taxon in taxa) {
    p <- power$power[power$taxon == taxon]
    y <- day7[[taxon]]
    null <- stats::glm(y ~ a + b, family = tweedie_glm_family(p), data = day7,
                       control = stats::glm.control(maxit = 100))
    # The full model: the cell means (zero where the taxon is absent).
    d1 <- sum(tweedie_glm_family(p)$dev.resids(y, stats::ave(y, day7$cell), 1))
    d0 <- stats::deviance(null)
    expected <- stats::df.residual(null) * (d0 - d1) / d0
    row <- out$taxa$time == "7" & out$taxa$taxon == taxon
    expect_equal(out$taxa$statistic[row], expected, tolerance = 1e-5)
  }
  # A taxon that one stressor eliminates alone follows MU (zero with the
  # combination): its means at the boundary make no deviance.
  d$gone <- ifelse(d$a == "1", 0, d$sp1)
  gone <- test_community(d, c("gone", "mass"), community_design(d), family = "tweedie",
                         overall = FALSE, n_boot = 99, seed = 1, quiet = TRUE)
  expect_equal(gone$taxa$statistic[gone$taxa$taxon == "gone"], c(0, 0), tolerance = 1e-6)
  for (tt in tests$time) {
    expect_equal(tests$statistic[tests$time == tt],
                 sum(out$taxa$statistic[out$taxa$time == tt]))
  }
  est <- out$taxa[out$taxa$time == "7" & out$taxa$taxon == "mass", ]
  m <- tapply(day7$mass, day7$cell, mean)
  expect_equal(est$estimate,
               unname(log(m[["1 1"]]) - log(m[["1 0"]]) - log(m[["0 1"]]) + log(m[["0 0"]])))
  expect_true(all(tests$p_value > 0 & tests$p_value <= 1))
  expect_true(all(out$taxa$p_adjusted > 0 & out$taxa$p_adjusted <= 1))
  expect_output(print(out), "tweedie GLM")
  again <- test_community(d, taxa, community_design(d), family = "tweedie", n_boot = 99,
                          seed = 1, quiet = TRUE)
  expect_identical(again$tests, out$tests)
  expect_equal(out$settings$versions[["glmmTMB"]], as.character(utils::packageVersion("glmmTMB")))
})

test_that("Tweedie tests take the effort and leave out taxa they cannot model", {
  d <- community_counts()
  d$volume <- rep(c(1, 2, 0.5, 1.5), length.out = nrow(d))
  d$dens <- d$sp1 / d$volume
  d$rare <- 0
  d$rare[[1L]] <- 0.4
  d$empty <- 0
  out <- test_community(d, c("dens", "mass", "rare", "empty"), community_design(d),
                        family = "tweedie", effort = c("volume", "volume", "volume", "volume"),
                        overall = FALSE, n_boot = 99, seed = 2, quiet = TRUE)
  expect_equal(out$tests$n_taxa, c(2L, 2L))
  ex <- out$diagnostics$excluded
  expect_setequal(ex$taxon[ex$time == "14"], c("rare", "empty"))
  expect_match(ex$reason[ex$taxon == "rare" & ex$time == "7"], "power")
  expect_match(ex$reason[ex$taxon == "empty" & ex$time == "7"], "Not found")
  expect_true(is.na(out$diagnostics$power$power[out$diagnostics$power$taxon == "rare"]))
  # With the effort, the full model is the glm with an offset.
  p <- out$diagnostics$power$power[out$diagnostics$power$taxon == "dens"]
  day7 <- d[d$day == 7, ]
  day7$cell <- paste(day7$a, day7$b)
  full <- rs_full(day7$dens, day7$volume, factor(day7$cell), rs_family("tweedie", power = p))
  fit <- stats::glm(dens ~ 0 + factor(cell) + offset(log(volume)),
                    family = tweedie_glm_family(p), data = day7,
                    control = stats::glm.control(epsilon = 1e-12, maxit = 100))
  expect_equal(full$deviance, stats::deviance(fit), tolerance = 1e-6)
  # Tweedie tests do not need mvabund.
  local_mocked_bindings(is_installed = function(pkg) FALSE)
  expect_no_error(test_community(d, c("dens", "mass"), community_design(d), family = "tweedie",
                                 overall = FALSE, n_boot = 99, seed = 1, quiet = TRUE))
})

# AD_taxa -----------------------------------------------------------------------

test_that("AD_taxa compares each taxon with the identity-link GLM of AD", {
  d <- community_counts()
  out <- test_community(d, c("sp1", "sp4"), community_design(d), family = "poisson",
                        null_models = "AD_taxa", overall = FALSE, n_boot = 99, seed = 1,
                        quiet = TRUE)
  expect_equal(out$tests$null_model, c("AD_taxa", "AD_taxa"))
  expect_true(all(out$tests$status == "ok"))
  expect_equal(out$tests$df, c(1, 1))
  day7 <- d[d$day == 7, ]
  day7$cell <- paste(day7$a, day7$b)
  for (taxon in c("sp1", "sp4")) {
    y <- day7[[taxon]]
    # With an interior solution, the additive null model is the Poisson GLM
    # with an identity link.
    null <- stats::glm(y ~ a + b, family = stats::poisson(link = "identity"), data = day7,
                       start = c(mean(y), 0, 0))
    full <- stats::glm(y ~ cell, family = stats::poisson(), data = day7)
    row <- out$taxa$time == "7" & out$taxa$taxon == taxon
    # The statistic: the drop in deviance divided by the dispersion of the
    # null model estimated from its deviance.
    d0 <- stats::deviance(null)
    expect_equal(out$taxa$statistic[row],
                 stats::df.residual(null) * (d0 - stats::deviance(full)) / d0, tolerance = 1e-5)
    m <- tapply(y, day7$cell, mean)
    expect_equal(out$taxa$estimate[row],
                 unname(m[["1 1"]] - max(0, m[["1 0"]] + m[["0 1"]] - m[["0 0"]])))
  }
  expect_equal(unique(out$taxa$scale), "response")
  expect_null(out$diagnostics$size)
})

test_that("AD_taxa expects a taxon that the stressors eliminate to be absent", {
  d <- community_counts()
  noise <- c(-1, 1, 0, 0)[d$rep]
  both <- d$a == "1" & d$b == "1"
  one <- xor(d$a == "1", d$b == "1")
  # `kept` falls from 20 to 15 with each stressor alone and is absent with
  # both: AD expects 15 + 15 - 20 = 10, a deviation. `gone` falls to 2 with
  # each stressor alone: AD expects 2 + 2 - 20 < 0, truncated at zero, and
  # its absence with both follows AD_taxa.
  d$kept <- ifelse(both, 0, ifelse(one, 15, 20) + noise)
  d$gone <- ifelse(both, 0, ifelse(one, 2, 20 + noise))
  out <- test_community(d, c("kept", "gone"), community_design(d), family = "poisson",
                        null_models = "AD_taxa", overall = FALSE, n_boot = 99, seed = 1,
                        quiet = TRUE)
  st <- out$taxa
  expect_equal(st$statistic[st$taxon == "gone"], c(0, 0))
  expect_equal(st$estimate[st$taxon == "gone"], c(0, 0))
  # The null model bends the means of the reference and the stressors alone
  # towards a zero combination; almost all its deviance is the interaction.
  expect_true(all(st$statistic[st$taxon == "kept"] > 0.9 * 13))
  expect_equal(st$estimate[st$taxon == "kept"], c(-10, -10))
})

test_that("negative binomial AD_taxa fixes the size of each taxon at its estimate", {
  skip_if_not_installed("MASS")
  d <- community_counts()
  out <- test_community(d, c("sp1", "sp2"), community_design(d), null_models = "AD_taxa",
                        n_boot = 99, seed = 1, quiet = TRUE)
  size <- out$diagnostics$size
  d$cell <- paste(d$a, d$b, d$day)
  for (taxon in c("sp1", "sp2")) {
    fit <- MASS::glm.nb(stats::as.formula(paste(taxon, "~ cell")), data = d)
    expect_equal(size$size[size$taxon == taxon], fit$theta, tolerance = 1e-3)
  }
  expect_false("mvabund" %in% names(out$settings$versions))
  expect_true(all(out$tests$status == "ok"))
})

test_that("MU and AD_taxa can be tested together", {
  skip_if_not_installed("mvabund")
  d <- community_counts()
  both <- test_community(d, c("sp1", "sp2"), community_design(d),
                         null_models = c("MU", "AD_taxa"), overall = FALSE, n_boot = 99,
                         seed = 1, quiet = TRUE)
  expect_equal(both$tests$null_model, c("MU", "AD_taxa", "MU", "AD_taxa"))
  expect_equal(both$tests$p_adjusted[both$tests$null_model == "MU"],
               stats::p.adjust(both$tests$p_value[both$tests$null_model == "MU"], "BH"))
  expect_output(print(both), "AD_taxa")
  expect_setequal(names(both$settings$statistic), c("MU", "AD_taxa"))
})

test_that("step-down p-values follow the successive maxima", {
  observed <- c(a = 5, b = 1, c = 3)
  boot <- rbind(c(6, 0, 1), c(1, 2, 2), c(2, 0.5, 4), c(0, 3, 0))
  # Order a, c, b: successive maxima from the smallest (b) up.
  q_b <- boot[, 2]
  q_c <- pmax(boot[, 3], q_b)
  q_a <- pmax(boot[, 1], q_c)
  p <- c(a = (1 + sum(q_a >= 5)) / 5, c = (1 + sum(q_c >= 3)) / 5, b = (1 + sum(q_b >= 1)) / 5)
  p <- cummax(p)
  expect_equal(step_down_p(observed, boot), p[c("a", "b", "c")])
})

test_that("test_community() checks its input", {
  d <- community_counts()
  des <- community_design(d)
  expect_error(test_community(d, "sp1", des), "two or more taxa", class = "must_error_input")
  expect_error(test_community(d, c("sp1", "other"), des), class = "must_error_input")
  expect_error(test_community(d, c("sp1", "mass"), des), "'mass' must hold counts",
               class = "must_error_input")
  expect_error(test_community(d, c("sp1", "mass"), des), "family = \"tweedie\"",
               class = "must_error_input")
  negative <- d
  negative$sp2[[3L]] <- -1
  expect_error(test_community(negative, c("sp1", "sp2"), des), "row 3",
               class = "must_error_input")
  expect_error(test_community(negative, c("sp1", "sp2"), des, family = "tweedie"), "row 3",
               class = "must_error_input")
  expect_error(test_community(d, c("sp1", "a"), des), "design column",
               class = "must_error_input")
  expect_error(test_community(d, c("sp1", "sp2"), des, family = "gaussian"),
               class = "must_error_input")
  expect_error(test_community(d, c("sp1", "sp2"), des, resampling = "case"),
               class = "must_error_input")
  expect_error(test_community(d, c("sp1", "sp2"), des, family = "tweedie",
                              resampling = "perm.resid"), "pit.trap",
               class = "must_error_input")
  expect_error(test_community(d, c("sp1", "sp2"), des, null_models = "AD_taxa",
                              resampling = "perm.resid"), "pit.trap",
               class = "must_error_input")
  expect_error(test_community(d, c("sp1", "sp2"), des, null_models = "AD"), "AD_taxa",
               class = "must_error_input")
  expect_error(test_community(d, c("sp1", "sp2"), des, null_models = c("MU", "MU")),
               class = "must_error_input")
  expect_error(test_community(d, c("sp1", "sp2"), des, n_boot = 10),
               class = "must_error_input")
  expect_error(test_community(d, c("sp1", "sp2"), list()), "describe_design",
               class = "must_error_input")
  expect_error(test_community(d, c("sp1", "sp2"), des, effort = "volume"),
               class = "must_error_input")
  expect_error(test_community(d, c("sp1", "sp2"), des, effort = c("mass", "mass", "mass")),
               "`effort`", class = "must_error_input")
  expect_error(test_community(d, c("sp1", "sp2"), des, effort = "sp2"), "taxon or a design",
               class = "must_error_input")
  zero <- d
  zero$volume <- 1
  zero$volume[[2L]] <- 0
  expect_error(test_community(zero, c("sp1", "sp2"), des, effort = "volume"), "volume",
               class = "must_error_input")
  local_mocked_bindings(is_installed = function(pkg) FALSE)
  expect_error(test_community(d, c("sp1", "sp2"), des), "mvabund")
})
