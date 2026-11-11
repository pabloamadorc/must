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
  switch(null_model, AD = "response", MU = "log(response)", DO = "response",
         BL = , SAM = , CA = "response / control")
}

# Null models whose effects are responses relative to the reference cell
# (w = y / y0): single and combined effects w - 1 and interactions on that
# scale.
relative_null_models <- function() {
  c("BL", "SAM", "CA")
}

# Stress addition model (Liess et al. 2016) -----------------------------------
#
# Each individual tolerates a general stress drawn from a beta(p, p)
# distribution. A stressor that leaves a relative response (survival) w
# exerts the general stress whose beta quantile is the mortality 1 - w; the
# general stresses of independent stressors add up, and the predicted
# survival is the share of individuals that tolerate their sum (zero when
# it reaches 1). Responses at or above the reference (w >= 1) exert no
# stress.
sam_stress <- function(w, shape) {
  stats::qbeta(pmin(pmax(1 - w, 0), 1), shape, shape)
}

sam_survival <- function(stress, shape) {
  1 - stats::pbeta(pmin(stress, 1), shape, shape)
}

# Expected relative response of the combination from the relative responses
# of the single stressors (matrix: one row per evaluation, one column per
# stressor).
sam_expected <- function(W, shape) {
  W <- as.matrix(W)
  stress <- matrix(sam_stress(W, shape), nrow = nrow(W))
  sam_survival(rowSums(stress), shape)
}

# Concentration addition (Loewe additivity) -------------------------------------
#
# Each stressor has a dose-response curve relative to the reference cell,
# here a two-parameter log-logistic w(x) = 1 / (1 + (x / EC50)^b), that is
# plogis(alpha + beta * log(x)) with beta = -b < 0 and alpha = b * log(EC50),
# fitted to the relative responses at the doses tested alone. Under
# concentration addition, the combination of doses x_i has the relative
# response w at which sum_i x_i / EC_i(w) = 1, with EC_i(w) the dose of
# stressor i alone that gives w (Loewe 1953; Tang et al. 2015, eq. 2).

# Least-squares fit of plogis(alpha + beta * u) to the rows of W (one row
# per evaluation, one column per dose; u = log(dose)): Gauss-Newton with
# Levenberg damping, vectorised over rows. Returns the matrix (alpha, beta).
fit_loglogistic <- function(W, u) {
  W <- as.matrix(W)
  n <- nrow(W)
  U <- matrix(u, nrow = n, ncol = length(u), byrow = TRUE)
  ok <- is.finite(W)
  # Start: least squares of the logits of the bounded responses on log dose.
  z <- stats::qlogis(pmin(pmax(W, 0.02), 0.98))
  z[!ok] <- 0
  wt <- ok * 1
  sw <- rowSums(wt)
  ub <- rowSums(wt * U) / sw
  zb <- rowSums(wt * z) / sw
  sxx <- rowSums(wt * (U - ub)^2)
  beta <- rowSums(wt * (U - ub) * (z - zb)) / sxx
  beta[!is.finite(beta)] <- 0
  alpha <- zb - beta * ub
  W0 <- W
  W0[!ok] <- 0
  sse <- function(a, b) {
    p <- stats::plogis(a + b * U)
    rowSums(wt * (W0 - p)^2)
  }
  current <- sse(alpha, beta)
  lambda <- rep(1e-3, n)
  for (iter in seq_len(100L)) {
    p <- stats::plogis(alpha + beta * U)
    g <- p * (1 - p) * wt
    r <- W0 - p
    a11 <- rowSums(g^2)
    a12 <- rowSums(g^2 * U)
    a22 <- rowSums(g^2 * U^2)
    b1 <- rowSums(g * r)
    b2 <- rowSums(g * r * U)
    d11 <- a11 * (1 + lambda) + 1e-12
    d22 <- a22 * (1 + lambda) + 1e-12
    det <- d11 * d22 - a12^2
    step_a <- (d22 * b1 - a12 * b2) / det
    step_b <- (d11 * b2 - a12 * b1) / det
    step_a[!is.finite(step_a)] <- 0
    step_b[!is.finite(step_b)] <- 0
    trial <- sse(alpha + step_a, beta + step_b)
    better <- is.finite(trial) & trial <= current
    alpha[better] <- alpha[better] + step_a[better]
    beta[better] <- beta[better] + step_b[better]
    current[better] <- trial[better]
    lambda <- ifelse(better, lambda / 10, lambda * 10)
    if (all(abs(step_a) + abs(step_b) < 1e-10 | lambda > 1e10)) {
      break
    }
  }
  cbind(alpha = alpha, beta = beta)
}

# Relative response of a combination under concentration addition. `curves`
# is a list with one (alpha, beta) matrix per stressor (rows: evaluations)
# and `x` the doses of the combination. Curves that do not decrease with dose
# (beta >= 0) contribute no dose equivalent. The root of
# sum_i x_i exp((t - alpha_i) / -beta_i) = 1 in t = logit(w) is found by
# bisection, vectorised over rows.
ca_expected <- function(curves, x) {
  n <- nrow(curves[[1L]])
  h <- function(t) {
    total <- numeric(n)
    for (i in seq_along(curves)) {
      a <- curves[[i]][, 1L]
      b <- curves[[i]][, 2L]
      active <- is.finite(a) & is.finite(b) & b < 0
      term <- numeric(n)
      term[active] <- x[[i]] * exp((t[active] - a[active]) / -b[active])
      total <- total + term
    }
    total - 1
  }
  lower <- rep(-60, n)
  upper <- rep(60, n)
  for (iter in seq_len(60L)) {
    mid <- (lower + upper) / 2
    above <- h(mid) > 0
    upper[above] <- mid[above]
    lower[!above] <- mid[!above]
  }
  w <- stats::plogis((lower + upper) / 2)
  decreasing <- lapply(curves, function(cv) is.finite(cv[, 2L]) & cv[, 2L] < 0)
  w[!Reduce(`|`, decreasing)] <- 1
  w
}

# Bliss independence (Bliss 1939), with the measures of Beppler et al. (2016)
# and Tekin et al. (2018). Each mean is taken relative to the reference cell,
# w_T = y_T / y_0. Single and combined effects are w - 1; the net interaction
# (deviation from additivity, DA) is w_K - prod(w_i); the emergent interaction
# of three stressors (E3, Beppler et al. 2016, eq. 2.4) is
#   w_XYZ - w_X w_YZ - w_Y w_XZ - w_Z w_XY + 2 w_X w_Y w_Z,
# the part of the three-way interaction that the pairs do not explain. `M` has
# one column per subset in all_subsets(m) order (columns not needed may be
# NA) and `sp` is a component specification of component_specs().
bliss_contrast <- function(M, m, sp) {
  W <- M / M[, 1L]
  size <- lengths(all_subsets(m))
  full <- length(size)
  if (!identical(sp$component, "interaction")) {
    return(W[, which(sp$weights == 1)] - 1)
  }
  if (identical(sp$estimand, "pure")) {
    if (m != 3L) {
      abort_internal("the emergent Bliss interaction is defined for three stressors only.")
    }
    # Columns: 1 control; 2-4 X, Y, Z; 5 XY; 6 XZ; 7 YZ; 8 XYZ.
    return(W[, 8L] - W[, 2L] * W[, 7L] - W[, 3L] * W[, 6L] - W[, 4L] * W[, 5L] +
             2 * W[, 2L] * W[, 3L] * W[, 4L])
  }
  W[, full] - apply(W[, size == 1L, drop = FALSE], 1L, prod)
}

# Stress addition contrasts: relative single and combined effects as under
# Bliss independence; the net interaction is w_K minus the survival that the
# stress addition model predicts from the single stressors.
sam_contrast <- function(M, m, sp) {
  if (!identical(sp$component, "interaction")) {
    return(bliss_contrast(M, m, sp))
  }
  W <- M / M[, 1L]
  size <- lengths(all_subsets(m))
  W[, length(size)] - sam_expected(W[, size == 1L, drop = FALSE], sp$shape)
}

# Log-logistic curves of the stressors of a concentration-addition contrast,
# from the relative means `W` of the cells of sp$cols (one row per
# evaluation).
ca_curves <- function(W, sp) {
  lapply(seq_along(sp$ca$curve_cols), function(i) {
    fit_loglogistic(W[, sp$ca$curve_cols[[i]], drop = FALSE], sp$ca$u[[i]])
  })
}

# Net concentration-addition interaction from the means of the cells of
# sp$cols: the reference first, the combination last.
ca_contrast <- function(D, sp) {
  D <- as.matrix(D)
  W <- D / D[, 1L]
  W[, ncol(W)] - ca_expected(ca_curves(W, sp), sp$ca$x)
}

# Rescaled Bliss interaction (Segre et al. 2005; Beppler et al. 2016): the
# interaction divided by the absolute value of the same expression with the
# combined response replaced by 0 when it is negative (lethal synergy) and by
# the smallest single-stressor response when it is positive (complete
# buffering), so that -1 and 1 mark those two extremes.
bliss_rescaled <- function(value, M, m) {
  W <- M / M[, 1L]
  size <- lengths(all_subsets(m))
  wk <- W[, length(size)]
  smallest <- apply(W[, size == 1L, drop = FALSE], 1L, min)
  out <- ifelse(value < 0, value / abs(value - wk),
                ifelse(value > 0, value / abs(value - wk + smallest), 0))
  out[!is.finite(out)] <- NA_real_
  out
}

# Descriptive classes of Beppler et al. (2016, section 2.2.4) from the
# rescaled interaction, with their conservative cut-offs of -0.5 and 0.5. They
# ignore the uncertainty of the estimate and are reported only as a description.
bliss_class <- function(rescaled) {
  ifelse(is.na(rescaled), NA_character_,
         ifelse(rescaled < -0.5, "synergy", ifelse(rescaled > 0.5, "antagonism", "additive")))
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
    DO = if (min(singles) < y0) min(singles) else NA_real_,
    BL = if (y0 > 0) y0 * prod(singles / y0) else NA_real_,
    abort_internal(sprintf("no expected combination for null model '%s' here.", null_model))
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
