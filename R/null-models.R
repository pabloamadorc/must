# Numerical core of the null models (internal) --------------------------------
#
# A contrast set with m stressors uses the cells indexed by the subsets T of
# its stressors, in the order of all_subsets(m): control first, then the m
# single-stressor cells, then pairs, ..., and the combination of all m last.
# `M` is a matrix of means on the response scale with one column per subset
# and one row per evaluation (one row for the point estimate, one per
# simulation draw).

# Weight vectors of the linear contrasts. Applied to f(M) with f = identity
# they give the simple-addition (AD) quantities; with f = log, the
# multiplicative (MU) ones.
contrast_weights <- function(m) {
  subsets <- all_subsets(m)
  size <- lengths(subsets)
  n <- length(subsets)
  single <- lapply(seq_len(m), function(i) {
    w <- numeric(n)
    w[size == 0L] <- -1
    w[vapply(subsets, function(t) identical(t, i), logical(1))] <- 1
    w
  })
  combined <- numeric(n)
  combined[size == 0L] <- -1
  combined[size == m] <- 1
  # Net deviation: y_K - sum(y_i) + (m - 1) y_0, i.e. the combined effect
  # minus the sum of the single effects.
  net <- numeric(n)
  net[size == m] <- 1
  net[size == 1L] <- -1
  net[size == 0L] <- m - 1
  # Pure m-th order interaction (inclusion-exclusion).
  pure <- (-1)^(m - size)
  list(single = single, combined = combined, net = net, pure = pure, size = size)
}

null_model_scale <- function(null_model) {
  switch(null_model, AD = "response", MU = "log(response)", DO = "response")
}

# Transformation of the means for each null model.
transform_means <- function(M, null_model) {
  if (identical(null_model, "MU")) log(M) else M
}

# Dominance interaction for two stressors (Madge Pimentel et al. 2025, eq. 9
# with factorial stressors): min(y_AB, y_0) - min(y_A, y_B).
dominance_interaction <- function(M) {
  pmin(M[, 4L], M[, 1L]) - pmin(M[, 2L], M[, 3L])
}

# Expected combined response (response scale) under each null model, from the
# point estimates of the means `mu` (vector in subset order).
expected_combined <- function(mu, m, null_model) {
  y0 <- mu[[1L]]
  singles <- mu[1L + seq_len(m)]
  switch(
    null_model,
    AD = y0 + sum(singles - y0),
    MU = if (all(mu > 0)) y0 * prod(singles / y0) else NA_real_,
    DO = if (min(singles) < y0) min(singles) else NA_real_
  )
}

# Summaries of a simulated distribution: standard error, percentile interval
# and two-sided p-value 2 * min(P(x <= 0), P(x >= 0)), bounded below by
# 1 / n_sim.
summarise_draws <- function(draws, level) {
  alpha <- (1 - level) / 2
  q <- stats::quantile(draws, probs = c(alpha, 1 - alpha), names = FALSE, type = 7)
  n <- length(draws)
  p <- 2 * min(mean(draws <= 0), mean(draws >= 0))
  list(
    se = stats::sd(draws),
    lower = q[[1L]],
    upper = q[[2L]],
    p_value = min(1, max(p, 1 / n))
  )
}
