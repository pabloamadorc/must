test_that("the reduced rules reproduce the 45 rows of Piggott et al. (2015, Table 1)", {
  # Typed from Table 1 (p. 1541): for each interaction type, the rows for an
  # interaction coded ++, +, 0, - and --.
  types <- list(
    c("+", "+"), c("-", "-"), c("+", "-"), c("-", "+"), c("0", "-"),
    c("-", "0"), c("+", "0"), c("0", "+"), c("0", "0")
  )
  classes <- list(
    c("+S", "+S", "AD", "+A", "-S"),  # double positive
    c("+S", "-A", "AD", "-S", "-S"),  # double negative
    c("+S", "-A", "AD", "+A", "-S"),  # opposing
    c("+S", "-A", "AD", "+A", "-S"),  # opposing
    c("+S", "-A", "AD", "-S", "-S"),  # negative neutral
    c("+S", "-A", "AD", "-S", "-S"),  # negative neutral
    c("+S", "+S", "AD", "+A", "-S"),  # positive neutral
    c("+S", "+S", "AD", "+A", "-S"),  # positive neutral
    c("+S", "+S", "AD", "-S", "-S")   # double neutral
  )
  codes <- c("++", "+", "0", "-", "--")
  for (k in seq_along(types)) {
    got <- vapply(codes, function(i) piggott_class(types[[k]], i), character(1))
    expect_equal(unname(got), classes[[k]], label = paste(types[[k]], collapse = ""))
  }
})

test_that("strong interactions follow the legend of Table 1", {
  # Mixture positive although both stressors are negative
  expect_equal(strong_interaction_code("+", singles = c(-1, -2), combined = 0.5,
                                       interaction = 3.5), "++")
  # Interaction larger than both single effects in absolute value
  expect_equal(strong_interaction_code("+", singles = c(-1, -2), combined = -0.5,
                                       interaction = 2.5), "++")
  expect_equal(strong_interaction_code("+", singles = c(-1, -2), combined = -2.5,
                                       interaction = 0.5), "+")
  expect_equal(strong_interaction_code("-", singles = c(1, 2), combined = -0.5,
                                       interaction = -3.5), "--")
  expect_equal(strong_interaction_code("0", singles = c(1, 2), combined = 9,
                                       interaction = 6), "0")
  # An infinite interaction (MU with a combination of zeros) is strong
  expect_equal(strong_interaction_code("-", singles = c(-0.1, 0.2), combined = -Inf,
                                       interaction = -Inf), "--")
})

test_that("the rules extend to any number of single effects", {
  expect_equal(piggott_class(c("+", "0", "-"), "+"), "-A")
  expect_equal(piggott_class(c("+", "0", "0"), "+"), "+S")
  expect_equal(piggott_class(c("+", "0", "-"), "-"), "+A")
  expect_equal(piggott_class(c("-", "-", "0"), "-"), "-S")
  expect_equal(piggott_class(c("-", "-", "0"), "0"), "AD")
  expect_equal(strong_interaction_code("+", singles = c(1, 2, 3), combined = 3.5,
                                       interaction = 1), "++")
  expect_equal(strong_interaction_code("+", singles = c(1, 2, 3), combined = 2.5,
                                       interaction = 1), "+")
  expect_true(is.na(piggott_class(c("+", NA), "+")))
})

test_that("a strong negative deviation is classified as negative synergism", {
  d <- two_by_two_gaussian(means = c(20, 15, 15, 2), dev = c(-0.2, -0.1, 0, 0.1, 0.2))
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("AD", "MU", "DO"), n_sim = 4000, seed = 3)
  cls <- classify_interactions(est)
  ad <- cls[cls$null_model == "AD", ]
  expect_equal(ad$code_a, "-")
  expect_equal(ad$code_b, "-")
  expect_equal(ad$code_interaction, "--")
  expect_equal(ad$class, "-S")
  expect_match(ad$meaning, "more negative")
  expect_equal(cls$status[cls$null_model == "DO"], "not_applicable")
})

test_that("sets of three stressors are classified from their net deviation", {
  d <- three_way_gaussian()
  des <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                         unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = "AD", inference = "delta")
  cls <- classify_interactions(est)
  expect_equal(nrow(cls), 4L)
  expect_true(all(cls$status == "ok"))
  three <- cls[cls$order == 3L, ]
  expect_equal(three$estimand, "net")
  expect_true(is.na(three$code_a))
  expect_equal(lengths(strsplit(three$single_codes, "; ", fixed = TRUE)), 3L)
  expect_match(three$pattern, "^n . \\| s . \\| f . \\| i ")
  e <- est$estimates
  piece <- e[e$order == 3L & e$component == "interaction" & e$estimand == "net", ]
  single <- e[e$order == 3L & e$component == "single", ]
  codes <- vapply(seq_len(nrow(single)), function(i) {
    effect_code(single$estimate[[i]], single$p_adjusted[[i]], 0.05)
  }, character(1))
  code_i <- strong_interaction_code(
    effect_code(piece$estimate, piece$p_adjusted, 0.05), single$estimate,
    e$estimate[e$order == 3L & e$component == "combined"], piece$estimate
  )
  expect_equal(three$class, piggott_class(codes, code_i))
})

test_that("borderline classes are flagged from the Monte Carlo error", {
  near <- data.frame(estimate = 1, p_value = 0.049, p_adjusted = 0.049, mc_se = 0.002)
  far <- data.frame(estimate = 1, p_value = 0.01, p_adjusted = 0.01, mc_se = 0.001)
  expect_true(is_borderline(near, alpha = 0.05))
  expect_false(is_borderline(far, alpha = 0.05))
  expect_false(is_borderline(transform(near, mc_se = NA_real_), alpha = 0.05))
})

test_that("combinations with only zeros are classified, and not_applicable rows are kept", {
  d <- zero_combination_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "y", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("AD", "MU", "DO"), n_sim = 2000, seed = 1,
                               sensitivity = "y")
  cls <- classify_interactions(est)
  ad <- cls[cls$null_model == "AD", ]
  expect_equal(ad$status, "ok")
  expect_equal(ad$class, "-S")
  expect_equal(cls$status[cls$null_model == "MU"], "not_applicable")
  est2 <- estimate_interactions(fit, null_models = c("AD", "MU", "DO"), n_sim = 2000, seed = 1)
  mu <- classify_interactions(est2)
  mu <- mu[mu$null_model == "MU", ]
  expect_equal(mu$code_interaction, "--")
  expect_equal(mu$class, "-S")
})

test_that("Bliss interactions are classified against the product of relative effects", {
  d <- two_by_two_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("MU", "BL"), n_sim = 2000, seed = 1)
  cls <- classify_interactions(est)
  bl <- cls[cls$null_model == "BL", ]
  mu <- cls[cls$null_model == "MU", ]
  expect_equal(bl$status, "ok")
  expect_equal(bl$scale, "response / control")
  # Same tests as MU, so the same codes and class
  expect_equal(bl[c("code_a", "code_b", "code_interaction", "class")],
               mu[c("code_a", "code_b", "code_interaction", "class")], ignore_attr = TRUE)
  expect_match(bl$meaning, "Bliss independence expectation")
  sam <- classify_interactions(estimate_interactions(fit, null_models = "SAM", n_sim = 2000,
                                                     seed = 1))
  expect_equal(sam$status, "ok")
  expect_match(sam$meaning, "stress addition expectation")
})

test_that("classify_interactions() checks its input", {
  expect_error(classify_interactions(list()), class = "must_error_input")
})
