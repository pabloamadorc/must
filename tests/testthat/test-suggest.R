known <- function(...) {
  # Every criterion of the guidelines known, with the values given in `...`.
  args <- list(similar_action = FALSE, specific_action = FALSE, mortality = FALSE,
               stressor_effect = FALSE, sensitivity_distribution = FALSE,
               correlation = "weak", effects = c(a = 0.2, b = 0.3))
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(suggest_null_models, args)
}

step_of <- function(x) {
  reason <- x$models$reason[x$models$status != "excluded"]
  as.integer(sub("^Step ([0-9]).*", "\\1", reason))
}

test_that("each step of Fig. 3 of Schaefer and Piggott (2018) gives its null model", {
  # Step 1: similar mode of action and known stressor-effect relationship.
  s1 <- known(similar_action = TRUE, stressor_effect = TRUE, correlation = "positive",
              mortality = TRUE)
  expect_identical(s1$null_models, "CA")
  expect_identical(step_of(s1), 1L)
  # Step 2: one stressor acting specifically, mortality, and a known
  # stressor-effect relationship or sensitivity distribution.
  s2 <- known(specific_action = TRUE, mortality = TRUE, sensitivity_distribution = TRUE)
  expect_identical(s2$null_models, "SAM")
  expect_identical(step_of(s2), 2L)
  expect_identical(known(specific_action = TRUE, mortality = TRUE,
                         stressor_effect = TRUE)$null_models, "SAM")
  expect_identical(known(specific_action = TRUE, mortality = FALSE,
                         stressor_effect = TRUE)$null_models, "AD")
  # Step 3: highly positively correlated sensitivities, or an effect close to
  # the effect limit.
  s3 <- known(correlation = "positive")
  expect_identical(s3$null_models, "DO")
  expect_identical(step_of(s3), 3L)
  expect_identical(known(effects = c(a = 0.95, b = 0.02))$null_models, "DO")
  # Step 4: highly negatively correlated and the sum below the limit, even
  # for mortality.
  s4 <- known(correlation = "negative", mortality = TRUE)
  expect_identical(s4$null_models, "AD")
  expect_identical(step_of(s4), 4L)
  # Step 5: mortality.
  s5 <- known(mortality = TRUE)
  expect_identical(s5$null_models, "MU")
  expect_identical(step_of(s5), 5L)
  expect_identical(known(mortality = TRUE, correlation = "negative",
                         effects = c(a = 0.6, b = 0.7))$null_models, "MU")
  # Step 6: the sum above the effect limit.
  s6 <- known(effects = c(a = 0.6, b = 0.7))
  expect_identical(s6$null_models, "MU")
  expect_identical(step_of(s6), 6L)
  expect_identical(known(correlation = "negative", effects = c(a = 0.6, b = 0.7))$null_models,
                   "MU")
  # Step 7: otherwise.
  s7 <- known()
  expect_identical(s7$null_models, "AD")
  expect_identical(step_of(s7), 7L)
  expect_identical(s7$models$status[s7$models$null_model == "AD"], "suggested")
  expect_true(all(s7$models$status[s7$models$null_model != "AD"] == "excluded"))
  # A sum equal to the effect limit fails steps 4 and 6.
  expect_identical(known(effects = c(a = 0.5, b = 0.5))$null_models, "AD")
  expect_identical(known(effects = c(a = 0.5, b = 0.5), correlation = "negative",
                         mortality = TRUE)$null_models, "MU")
})

test_that("unknown criteria keep every null model reached with one of their values", {
  everything <- suggest_null_models()
  expect_identical(everything$null_models, c("AD", "MU", "DO", "SAM", "CA"))
  expect_true(all(everything$models$status == "possible"))
  values <- list(similar_action = c(TRUE, FALSE), specific_action = c(TRUE, FALSE),
                 mortality = c(TRUE, FALSE), stressor_effect = c(TRUE, FALSE),
                 sensitivity_distribution = c(TRUE, FALSE),
                 correlation = c("positive", "negative", "weak"))
  # With one criterion unknown, the models kept are the union of the models
  # given by each of its values.
  base <- list(similar_action = FALSE, specific_action = TRUE, mortality = TRUE,
               stressor_effect = FALSE, sensitivity_distribution = TRUE,
               correlation = "negative", effects = c(a = 0.3, b = 0.4))
  for (cn in names(values)) {
    unknown <- base
    unknown[[cn]] <- NA
    union <- unique(unlist(lapply(values[[cn]], function(v) {
      args <- base
      args[[cn]] <- v
      do.call(suggest_null_models, args)$null_models
    })))
    expect_setequal(do.call(suggest_null_models, unknown)$null_models, union)
  }
  # The example of section 4: mode of action and correlation unknown.
  s <- suggest_null_models(similar_action = NA, correlation = NA, mortality = FALSE,
                           stressor_effect = FALSE, effects = c(a = 0.3, b = 0.5))
  expect_identical(s$null_models, c("AD", "DO"))
  expect_identical(s$decisive$criterion, rep("correlation", 3))
  expect_identical(s$decisive$null_models, c("DO", "AD", "AD"))
  expect_identical(unique(s$models$status[s$models$null_model %in% c("AD", "DO")]),
                   "possible")
  expect_match(s$models$reason[s$models$null_model == "MU"], "not mortality")
})

test_that("effects given by the user decide the criteria that depend on them", {
  near <- function(...) known(...)$criteria$value[[7L]]
  total <- function(...) known(...)$criteria$value[[8L]]
  expect_identical(near(effects = c(a = 0.9, b = 0.1)), "yes")
  expect_identical(near(effects = c(a = 0.9, b = 0.1), near_limit = 0.95), "no")
  expect_identical(near(effects = c(a = NA, b = 0.95)), "yes")
  expect_identical(near(effects = c(a = NA, b = 0.2)), "unknown")
  expect_identical(total(effects = c(a = NA, b = 0.2)), "unknown")
  expect_identical(total(effects = c(a = 0.3, b = 0.4)), "below")
  expect_identical(total(effects = c(a = 0.6, b = 0.7)), "above")
  expect_identical(total(effects = c(a = 0.25, b = 0.75)), "equal")
  # An effect that increases the response is kept, with its direction.
  up <- known(effects = c(a = -0.5, b = 0.3))
  expect_identical(up$effects$direction, c("increase", "decrease"))
  expect_identical(up$criteria$value[9:10], c("yes", "yes"))
  expect_identical(up$null_models, "AD")
  # Without effects, the criteria that depend on them are unknown.
  none <- known(effects = NULL)
  expect_identical(none$criteria$value[7:8], c("unknown", "unknown"))
  expect_identical(none$criteria$source[7:8], c("unknown", "unknown"))
  expect_null(none$effects)
  expect_identical(none$null_models, c("AD", "MU", "DO"))
})

test_that("the directions of the effects give the types of Piggott et al. (2015)", {
  expect_identical(effect_type(c("decrease", "decrease")), "double negative")
  expect_identical(effect_type(c("increase", "increase")), "double positive")
  expect_identical(effect_type(c("decrease", "increase")), "opposing")
  expect_identical(effect_type(c("decrease", "none")), "negative neutral")
  expect_identical(effect_type(c("none", "increase")), "positive neutral")
  expect_identical(effect_type(c("none", "none")), "double neutral")
  expect_identical(effect_type(c("decrease", "decrease", "none")), "negative neutral")
  expect_identical(effect_type(c("decrease", "decrease", "decrease")), "all negative")
  expect_identical(effect_type(c("increase", "unknown", "decrease")), "opposing")
  expect_true(is.na(effect_type(c("increase", "unknown"))))
  expect_identical(effect_directions(c(0.3, -0.2, 0, NA), NULL, NULL, use_points = TRUE),
                   c("decrease", "increase", "none", "unknown"))
  expect_identical(effect_directions(c(0.3, -0.2, 0.1, 0.5), c(0.1, -0.5, -0.1, NA),
                                     c(0.5, -0.1, 0.3, NA), use_points = FALSE),
                   c("decrease", "increase", "none", "unknown"))
  # Whether some stressor increases the response or two have opposite effects.
  dc <- function(...) unlist(direction_criteria(c(...))[c("increase", "opposing")])
  expect_identical(dc("decrease", "none"), c(increase = FALSE, opposing = FALSE))
  expect_identical(dc("increase", "none"), c(increase = TRUE, opposing = FALSE))
  expect_identical(dc("increase", "decrease"), c(increase = TRUE, opposing = TRUE))
  expect_identical(dc("increase", "unknown"), c(increase = TRUE, opposing = NA))
  expect_identical(dc("none", "unknown"), c(increase = NA, opposing = FALSE))
  expect_identical(dc("unknown", "unknown"), c(increase = NA, opposing = NA))
})

test_that("concentration and stress addition need effects in the same direction", {
  # Step 1 with opposite effects, and step 2 with an increase, do not apply.
  ca <- known(similar_action = TRUE, stressor_effect = TRUE, effects = c(a = 0.3, b = -0.2))
  expect_false("CA" %in% ca$null_models)
  expect_match(ca$models$reason[ca$models$null_model == "CA"], "opposite effects")
  expect_identical(known(similar_action = TRUE, stressor_effect = TRUE,
                         effects = c(a = -0.3, b = -0.2))$null_models, "CA")
  sam <- known(specific_action = TRUE, mortality = TRUE, sensitivity_distribution = TRUE,
               effects = c(a = 0.3, b = -0.2))
  expect_identical(sam$null_models, "MU")
  expect_match(sam$models$reason[sam$models$null_model == "SAM"], "increases the response")
  expect_identical(known(specific_action = TRUE, mortality = TRUE,
                         sensitivity_distribution = TRUE,
                         effects = c(a = 0.3, b = 0))$null_models, "SAM")
  # With the effects unknown, the stress addition model needs no increase,
  # and then the effects are not opposite either, so step 1 applies first.
  open <- suggest_null_models(similar_action = TRUE, stressor_effect = TRUE,
                              specific_action = TRUE, mortality = TRUE)
  expect_identical(open$null_models, c("AD", "MU", "DO", "CA"))
  expect_identical(open$decisive$null_models[open$decisive$criterion == "opposing"],
                   c("AD, MU, DO", "CA"))
  expect_false("SAM" %in% open$null_models)
})

test_that("effects estimated from a Gaussian fit match hand calculations", {
  d <- two_by_two_gaussian()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  s <- suggest_null_models(fit, similar_action = FALSE, mortality = FALSE,
                           stressor_effect = FALSE, correlation = "weak")
  # Cell means 10, 8, 7 and 3; residual variance RSS / (N - 4) = 0.625.
  v <- 0.625 / 5
  q <- stats::qt(0.975, 16)
  r <- c(8, 7) / 10
  se <- sqrt(v * (1 / c(8, 7)^2 + 1 / 10^2))
  e <- s$effects
  expect_equal(e$stressor, c("a", "b"))
  expect_equal(e$estimate, 1 - r, tolerance = 1e-6)
  expect_equal(e$lower, 1 - r * exp(q * se), tolerance = 1e-5)
  expect_equal(e$upper, 1 - r * exp(-q * se), tolerance = 1e-5)
  expect_equal(e$df, c(16, 16))
  # The sum of the effects is 2 - (8 + 7) / 10, with its interval from the log
  # scale of the sum of the ratios.
  total <- (8 + 7) / 10
  se_log <- sqrt(v * ((8 + 7)^2 / 10^4 + 2 / 10^2)) / total
  expect_equal(s$sets$sum, 0.5, tolerance = 1e-6)
  expect_equal(c(s$sets$sum_lower, s$sets$sum_upper),
               2 - total * exp(c(1, -1) * q * se_log), tolerance = 1e-5)
  expect_false(s$sets$near_limit)
  expect_identical(s$sets$sum_vs_limit, "below")
  expect_identical(s$null_models, "AD")
  expect_identical(s$criteria$source[7:8], c("estimated", "estimated"))
  # The point estimates alone, without intervals.
  p <- suggest_null_models(fit, level = NULL)
  expect_true(all(is.na(p$effects$lower)))
  expect_identical(p$sets$sum_vs_limit, "below")
})

test_that("an interval that spans a threshold leaves its criterion unknown", {
  d <- two_by_two_gaussian(means = c(10, 4.6, 4.8, 3), dev = c(-1.5, -0.75, 0, 0.75, 1.5))
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  s <- suggest_null_models(fit, similar_action = FALSE, mortality = FALSE,
                           stressor_effect = FALSE, correlation = "weak")
  # The sum of the effects is 1.06 with an interval that spans 1.
  expect_lt(s$sets$sum_lower, 1)
  expect_gt(s$sets$sum_upper, 1)
  expect_true(is.na(s$sets$sum_vs_limit))
  expect_identical(s$criteria$value[[8L]], "uncertain")
  expect_identical(s$null_models, c("AD", "MU"))
  # With a single set and time, whether the sum is above the limit decides.
  expect_identical(s$decisive$criterion, c("sum", "sum"))
  expect_identical(s$decisive$null_models, c("AD", "MU"))
  # With the point estimates the sum is above the limit.
  expect_identical(suggest_null_models(fit, similar_action = FALSE, mortality = FALSE,
                                       stressor_effect = FALSE, correlation = "weak",
                                       level = NULL)$null_models, "MU")
})

test_that("a cell with only zeros gives an effect without interval", {
  d <- two_by_two_counts()
  d$count[d$a == "0" & d$b == "1"] <- 0
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  s <- suggest_null_models(fit, similar_action = FALSE, mortality = FALSE,
                           stressor_effect = FALSE, correlation = "weak")
  e <- s$effects
  expect_equal(e$estimate[e$stressor == "b"], 1, tolerance = 1e-6)
  expect_true(is.na(e$lower[e$stressor == "b"]))
  expect_identical(e$note, c("", "boundary"))
  expect_false(is.na(e$lower[e$stressor == "a"]))
  # Whether the effect of b is close to the limit is not known.
  expect_true(is.na(s$sets$near_limit))
  expect_true(is.na(s$sets$sum_lower))
  expect_identical(s$null_models, c("AD", "MU", "DO"))
})

test_that("effects are not defined without an effect limit", {
  d <- two_by_two_gaussian(means = c(1, -2, 0.5, -1))
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  s <- suggest_null_models(fit, mortality = FALSE)
  expect_identical(s$sets$status, "not_defined")
  expect_match(s$sets$reason, "negative values")
  expect_identical(s$criteria$value[7:8], c("unknown", "unknown"))
  expect_identical(s$settings$effects, "unknown")
  g <- two_by_two_gaussian()
  gfit <- fit_responses(g, "y", describe_design(g, c("a", "b"), control = c(a = "0", b = "0"),
                                                unit = "id"),
                        families = list(gaussian = stats::gaussian()),
                        residual_checks = FALSE, quiet = TRUE)
  expect_match(suggest_null_models(gfit, sensitivity = "y")$sets$reason, "no upper")
})

test_that("mortality is complemented to survival with sensitivity", {
  d <- expand.grid(a = c("0", "1"), b = c("0", "1"), rep = 1:5, stringsAsFactors = FALSE)
  d$id <- seq_len(nrow(d))
  d$n <- 20
  p <- c(0.9, 0.6, 0.5, 0.3)[match(paste(d$a, d$b), c("0 0", "1 0", "0 1", "1 1"))]
  d$alive <- round(d$n * p + c(-1, 0, 1, 0, 0)[d$rep])
  d$dead <- d$n - d$alive
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, c("alive", "dead"), des, trials = c(alive = "n", dead = "n"),
                       residual_checks = FALSE, quiet = TRUE)
  alive <- suggest_null_models(fit, response = "alive", mortality = TRUE,
                               similar_action = FALSE, specific_action = FALSE,
                               correlation = "weak")
  dead <- suggest_null_models(fit, response = "dead", sensitivity = "dead", mortality = TRUE,
                              similar_action = FALSE, specific_action = FALSE,
                              correlation = "weak")
  expect_equal(alive$effects$estimate, 1 - c(0.6, 0.5) / 0.9, tolerance = 1e-4)
  expect_equal(dead$effects[, c("estimate", "lower", "upper")],
               alive$effects[, c("estimate", "lower", "upper")], tolerance = 1e-6)
  expect_identical(dead$null_models, "MU")
  # The interactions must then be estimated on the same scale.
  expect_error(estimate_interactions(fit, null_models = dead), class = "must_error_input")
  est <- estimate_interactions(fit, null_models = dead, sensitivity = "dead",
                               inference = "delta")
  expect_identical(est$settings$null_models, "MU")
})

test_that("a proportion that the stressors increase is judged on its complement", {
  d <- expand.grid(a = c("0", "1"), b = c("0", "1"), rep = 1:5, stringsAsFactors = FALSE)
  d$id <- seq_len(nrow(d))
  d$n <- 40
  p <- c(0.5, 0.8, 0.85, 0.95)[match(paste(d$a, d$b), c("0 0", "1 0", "0 1", "1 1"))]
  d$alive <- round(d$n * p + c(-1, 0, 1, 0, 0)[d$rep])
  d$dead <- d$n - d$alive
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, c("alive", "dead"), des, trials = c(alive = "n", dead = "n"),
                       residual_checks = FALSE, quiet = TRUE)
  s <- suggest_null_models(fit, response = "alive", mortality = FALSE, similar_action = FALSE,
                           specific_action = FALSE, correlation = "weak")
  e <- s$effects
  expect_identical(e$scale, c("response", "response", "complement", "complement"))
  expect_identical(e$direction, c("increase", "increase", "decrease", "decrease"))
  # Survival rises from 0.5 to 0.8 and 0.85: mortality falls by 60% and 70%.
  expect_equal(e$estimate, c(1 - c(0.8, 0.85) / 0.5, 1 - c(0.2, 0.15) / 0.5), tolerance = 1e-4)
  expect_identical(s$sets$scale, "complement")
  expect_identical(s$sets$type, "double positive")
  expect_identical(s$sets$sum_vs_limit, "above")
  expect_identical(s$null_models, "MU")
  expect_identical(s$settings$complement, "alive")
  expect_identical(s$settings$suggested_sensitivity, "alive")
  expect_match(s$criteria$detail[[8L]], "complement")
  expect_output(print(s), "sensitivity = \"alive\"")
  # The interactions must then be estimated on the complement.
  expect_error(estimate_interactions(fit, null_models = s, inference = "delta"), "complement")
  est <- estimate_interactions(fit, null_models = s, sensitivity = "alive", inference = "delta")
  expect_identical(est$settings$null_model_choice$complement, "alive")
  # Analysed as survival through `sensitivity`, mortality is then left as it is.
  dead <- suggest_null_models(fit, response = "dead", sensitivity = "dead", mortality = FALSE,
                              similar_action = FALSE, specific_action = FALSE,
                              correlation = "weak")
  expect_identical(dead$settings$complement, "dead")
  expect_identical(dead$settings$suggested_sensitivity, character())
  expect_equal(dead$sets$sum, s$sets$sum, tolerance = 1e-6)
  # Opposite effects keep the scale as analysed.
  d$alive2 <- round(d$n * c(0.5, 0.8, 0.3, 0.5)[match(paste(d$a, d$b),
                                                       c("0 0", "1 0", "0 1", "1 1"))])
  fit2 <- fit_responses(d, "alive2", des, trials = c(alive2 = "n"), residual_checks = FALSE,
                        quiet = TRUE)
  mixed <- suggest_null_models(fit2, mortality = FALSE)
  expect_identical(mixed$sets$scale, "response")
  expect_identical(mixed$sets$type, "opposing")
  expect_identical(mixed$settings$complement, character())
})

test_that("the complement also applies to a binomial model fitted outside must", {
  skip_if_not_installed("emmeans")
  d <- expand.grid(a = c("0", "1"), b = c("0", "1"), rep = 1:5, stringsAsFactors = FALSE)
  d$id <- seq_len(nrow(d))
  d$n <- 40
  p <- c(0.5, 0.8, 0.85, 0.95)[match(paste(d$a, d$b), c("0 0", "1 0", "0 1", "1 1"))]
  d$alive <- round(d$n * p + c(-1, 0, 1, 0, 0)[d$rep])
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  model <- stats::glm(cbind(alive, n - alive) ~ a * b, family = stats::binomial(), data = d)
  s <- suggest_null_models(model, design = des, mortality = FALSE, similar_action = FALSE,
                           specific_action = FALSE, correlation = "weak")
  expect_identical(s$sets$scale, "complement")
  expect_true(s$settings$suggested_sensitivity)
  expect_equal(s$effects$estimate[s$effects$scale == "complement"], 1 - c(0.2, 0.15) / 0.5,
               tolerance = 1e-4)
  expect_error(estimate_interactions(model, des, null_models = s, inference = "delta"),
               class = "must_error_input")
  est <- estimate_interactions(model, des, null_models = s, sensitivity = TRUE,
                               inference = "delta")
  expect_identical(est$settings$null_models, "MU")
})

test_that("estimate_interactions() takes the suggestion and keeps its reasons", {
  d <- two_by_two_gaussian()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  s <- suggest_null_models(fit, similar_action = FALSE, mortality = FALSE)
  expect_identical(s$null_models, c("AD", "DO"))
  est <- estimate_interactions(fit, null_models = s, n_sim = 500, seed = 1)
  expect_identical(est$settings$null_models, c("AD", "DO"))
  choice <- est$settings$null_model_choice
  expect_identical(choice$null_models, s$null_models)
  expect_identical(choice$criteria, s$criteria)
  expect_identical(choice$responses, "y")
  expect_output(print(est), "from suggest_null_models")
  expect_null(estimate_interactions(fit, n_sim = 500, seed = 1)$settings$null_model_choice)
  # The same estimates as with the codes.
  codes <- estimate_interactions(fit, null_models = c("AD", "DO"), n_sim = 500, seed = 1)
  expect_identical(est$estimates, codes$estimates)
})

test_that("a model fitted outside must gives the same effects", {
  skip_if_not_installed("emmeans")
  d <- two_by_two_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  model <- stats::glm(count ~ a * b, family = stats::poisson(), data = d)
  inside <- suggest_null_models(fit, mortality = FALSE)
  outside <- suggest_null_models(model, design = des, mortality = FALSE)
  cols <- c("stressor", "estimate", "lower", "upper")
  expect_equal(outside$effects[, cols], inside$effects[, cols], tolerance = 1e-5)
  expect_identical(outside$null_models, inside$null_models)
  est <- estimate_interactions(model, des, null_models = outside, n_sim = 500, seed = 1)
  expect_identical(est$settings$null_models, outside$null_models)
  expect_error(estimate_interactions(model, des, null_models = outside, inference = "delta"),
               "setdiff")
})

test_that("the effects of each `by` group use the control shared with it", {
  d <- shared_regimes()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id",
                         by = "regime", shared = c(regime = "both"))
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  s <- suggest_null_models(fit, mortality = FALSE)
  e <- s$effects
  expect_identical(e$regime, c("press", "press", "pulse", "pulse"))
  expect_equal(e$estimate, 1 - c(9, 6, 8, 7) / 10, tolerance = 1e-6)
  expect_identical(nrow(s$sets), 2L)
})

test_that("three stressors use every effect of the set", {
  d <- three_way_gaussian()
  des <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                         unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  s <- suggest_null_models(fit, mortality = FALSE, correlation = "weak")
  three <- s$sets[s$sets$order == 3L, ]
  expect_identical(nrow(three), 1L)
  m <- tapply(d$y, paste(d$n, d$s, d$f), mean)
  single <- 1 - m[c("1 0 0", "0 1 0", "0 0 1")] / m[["0 0 0"]]
  expect_equal(three$sum, sum(single), tolerance = 1e-6)
  expect_identical(s$models$design[s$models$null_model == "DO"], "some sets")
  expect_match(s$models$design_note[s$models$null_model == "DO"], "two stressors only")
})

test_that("the design says whether concentration addition can be computed", {
  g <- expand.grid(a = c("0", "1", "2", "4"), b = c("0", "1", "2"), rep = 1:3,
                   stringsAsFactors = FALSE)
  g$id <- seq_len(nrow(g))
  g$y <- 50 * exp(-0.2 * as.numeric(g$a)) * exp(-0.3 * as.numeric(g$b)) + c(-1, 0, 1)[g$rep]
  dd <- describe_design(g, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  s <- suggest_null_models(dd, similar_action = TRUE, stressor_effect = TRUE,
                           effects = c(a = 0.3, b = 0.4))
  expect_identical(s$null_models, "CA")
  # Without the effects, their directions are unknown: opposite effects would
  # rule concentration addition out.
  unknown <- suggest_null_models(dd, similar_action = TRUE, stressor_effect = TRUE)
  expect_identical(unknown$null_models, c("AD", "MU", "DO", "CA"))
  ca <- s$models[s$models$null_model == "CA", ]
  expect_identical(ca$design, "yes")
  expect_match(ca$design_note, "4 levels of 'a' and 3 levels of 'b'")
  d <- two_by_two_gaussian()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  ca2 <- suggest_null_models(des)$models
  expect_identical(ca2$design[ca2$null_model == "CA"], "no")
  expect_match(ca2$design_note[ca2$null_model == "CA"], "at least two doses")
  expect_true(all(is.na(suggest_null_models()$models$design)))
})

test_that("invalid input gives informative errors", {
  d <- two_by_two_gaussian()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  expect_error(suggest_null_models(mortality = "yes"), class = "must_error_input")
  expect_error(suggest_null_models(similar_action = c(TRUE, FALSE)),
               class = "must_error_input")
  expect_error(suggest_null_models(correlation = "high"), class = "must_error_input")
  expect_error(suggest_null_models(effects = c(0.2, 0.3)), class = "must_error_input")
  expect_error(suggest_null_models(effects = c(a = 0.2)), class = "must_error_input")
  expect_error(suggest_null_models(effects = c(a = 1.2, b = 0.3)), "cannot exceed 1")
  expect_error(suggest_null_models(des, effects = c(a = 0.2, z = 0.3)), "'z'")
  expect_error(suggest_null_models(fit, effects = c(a = 0.2, b = 0.3)), "not both")
  expect_error(suggest_null_models(near_limit = 0), class = "must_error_input")
  expect_error(suggest_null_models(level = 1), class = "must_error_input")
  expect_error(suggest_null_models(stats::lm(y ~ a * b, data = d)), "`design`")
  expect_error(suggest_null_models(des, response = "y"), class = "must_error_input")
  expect_error(suggest_null_models(fit, response = "z"), class = "must_error_input")
  expect_error(suggest_null_models(des, sensitivity = "y"), class = "must_error_input")
  expect_error(suggest_null_models(fit, design = des), class = "must_error_input")
  expect_error(suggest_null_models(fit, other = 1), class = "must_error_input")
  est <- estimate_interactions(fit, inference = "delta")
  expect_error(suggest_null_models(est), "estimate_interactions")
})

test_that("the print method shows the criteria, the models and the reasons", {
  s <- suggest_null_models(similar_action = FALSE, mortality = FALSE, stressor_effect = FALSE,
                           effects = c(salinity = 0.3, warming = 0.5))
  expect_snapshot(print(s))
})
