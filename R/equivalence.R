# Tests of equivalence (internal) ----------------------------------------------
#
# An interaction is equivalent to its null model when its relative deviation
# d = (observed - expected) / reference lies within (-margin, margin). Two
# one-sided tests (Schuirmann 1987): "d <= -margin" and "d >= margin". Each
# uses the signed root of the distance from the estimated means of the cells
# to the closest means at which d equals the bound, in the metric of their
# covariance, on a scale on which the variance of a mean does not depend on
# the mean (variance-stabilizing; Anscombe 1948). On that scale the
# log-likelihood of the means is nearly quadratic, so the squared distance is
# close to the likelihood-ratio statistic, and a standard error does not
# shrink as an estimate moves away from the bound, as it does for counts on
# the scale of the model. When the bound is linear in the means (AD with a
# constant variance) the distance is the t statistic of observed - expected
# -/+ margin * reference, the t test of a ratio to a control (Hauschke,
# Kieser & Hothorn 1999). The test is the same with either inference.

# Variance-stabilizing scale of the cell means for the error distribution of
# the model: h() maps means to it, h_inv() back and dh() is the derivative of
# h(). NULL when the family is not one of these, or when the model was fitted
# to a transformed response (the variance is then stable on the scale of the
# model already): the tests use the scale of the linear predictor.
stabilizing_scale <- function(family, theta = NA_real_, transformed = FALSE) {
  if (isTRUE(transformed) || !inherits(family, "family") || !is.character(family$family)) {
    return(NULL)
  }
  name <- tolower(family$family[[1L]])
  if (identical(name, "gaussian")) {
    return(list(h = function(mu) mu, h_inv = function(z) z,
                dh = function(mu) rep(1, length(mu))))
  }
  if (name %in% c("poisson", "quasipoisson", "nbinom1", "genpois")) {
    return(list(h = function(mu) sqrt(mu), h_inv = function(z) pmax(z, 0)^2,
                dh = function(mu) 1 / (2 * sqrt(mu))))
  }
  if (name %in% c("binomial", "quasibinomial", "betabinomial", "beta")) {
    return(list(h = function(p) asin(sqrt(p)),
                h_inv = function(z) sin(pmin(pmax(z, 0), pi / 2))^2,
                dh = function(p) 1 / (2 * sqrt(p * (1 - p)))))
  }
  if (name %in% c("gamma", "lognormal")) {
    return(list(h = function(mu) log(mu), h_inv = function(z) exp(z),
                dh = function(mu) 1 / mu))
  }
  if (identical(name, "nbinom2") || startsWith(name, "negative binomial")) {
    if (!startsWith(name, "nbinom2")) {
      theta <- suppressWarnings(as.numeric(sub("^negative binomial\\((.*)\\)$", "\\1", name)))
    }
    if (length(theta) == 1L && is.finite(theta) && theta > 0) {
      return(list(h = function(mu) asinh(sqrt(mu / theta)),
                  h_inv = function(z) theta * sinh(pmax(z, 0))^2,
                  dh = function(mu) 1 / (2 * sqrt(theta * mu * (1 + mu / theta)))))
    }
  }
  NULL
}

# Size parameter of a negative binomial (variance mu + mu^2 / theta) fitted
# with glmmTMB with a single dispersion; NA otherwise.
nbinom2_theta <- function(model) {
  if (!inherits(model, "glmmTMB")) {
    return(NA_real_)
  }
  disp <- tryCatch(glmmTMB::fixef(model)$disp, error = function(e) NULL)
  if (length(disp) != 1L) {
    return(NA_real_)
  }
  theta <- tryCatch(suppressWarnings(stats::sigma(model)), error = function(e) NA_real_)
  if (length(theta) == 1L && is.finite(theta)) theta else NA_real_
}

# p-value of the equivalence test of the relative deviation `rel_fun` (a
# function of a matrix of means of the cells `cols`, one row per evaluation):
# the larger of the one-sided p-values against -margin and against margin,
# with a t distribution with `df` degrees of freedom (normal if infinite). NA
# if a closest point cannot be found.
equivalence_p_signed <- function(rel_fun, pred, cols, df, margin) {
  mu <- pred$mu[cols]
  V_eta <- pred$V_eta[cols, cols, drop = FALSE]
  st <- pred$stable
  dh <- if (is.null(st)) NA_real_ else st$dh(mu)
  if (!is.null(st) && all(is.finite(dh)) && all(is.finite(pred$dmu[cols]))) {
    z0 <- st$h(mu)
    J <- dh * pred$dmu[cols]
    V <- V_eta * outer(J, J)
    to_rel <- function(Z) rel_fun(st$h_inv(Z))
  } else {
    z0 <- pred$eta[cols]
    V <- V_eta
    to_rel <- function(Z) rel_fun(pred$to_mean(Z))
  }
  lower <- signed_distance(z0, V, function(Z) to_rel(Z) + margin)
  upper <- signed_distance(z0, V, function(Z) to_rel(Z) - margin)
  if (!is.finite(lower) || !is.finite(upper)) {
    return(NA_real_)
  }
  p_lower <- if (is.finite(df)) stats::pt(lower, df, lower.tail = FALSE) else stats::pnorm(lower, lower.tail = FALSE)
  p_upper <- if (is.finite(df)) stats::pt(upper, df) else stats::pnorm(upper)
  max(p_lower, p_upper)
}

# Signed root of the distance from `z0` to the surface g(z) = 0 in the metric
# of the covariance `V` (the Mahalanobis distance to its closest point),
# positive when g(z0) > 0. `g` takes a one-row matrix. The closest point is
# found by projecting z0 on successive tangent planes of the surface
# (Gauss-Newton); when the surface curves too much for those to converge, by
# Newton's method on the conditions of the closest point, z - z0 + lambda V
# grad g(z) = 0 and g(z) = 0, which use the curvature of g. NA if neither
# converges. At the closest point z - z0 = -lambda V grad g, so the distance
# is |lambda| sqrt(grad g' V grad g), without inverting V.
signed_distance <- function(z0, V, g, max_iter = 50L) {
  at <- function(z) g(matrix(z, nrow = 1L))
  g0 <- at(z0)
  if (!is.finite(g0)) {
    return(NA_real_)
  }
  if (g0 == 0) {
    return(0)
  }
  se <- sqrt(pmax(diag(V), 0))
  tol <- 1e-8 * max(se, 1e-12)
  # Central differences with steps on the scale of the standard errors (the
  # values may be close to zero on the scale of the linear predictor).
  steps <- 1e-6 * pmax(abs(z0), se, 1e-6)
  gradient <- function(z, h = steps) {
    vapply(seq_along(z), function(j) {
      up <- z
      down <- z
      up[[j]] <- z[[j]] + h[[j]]
      down[[j]] <- z[[j]] - h[[j]]
      (at(up) - at(down)) / (2 * h[[j]])
    }, numeric(1))
  }
  multiplier <- function(z, gr) {
    (at(z) + sum(gr * (z0 - z))) / sum(gr * drop(V %*% gr))
  }
  distance <- function(z) {
    gr <- gradient(z)
    s <- sum(gr * drop(V %*% gr))
    if (!is.finite(s) || s <= 0) {
      return(NA_real_)
    }
    sign(g0) * abs(multiplier(z, gr)) * sqrt(s)
  }
  converged <- function(change, z) {
    change <= tol && isTRUE(abs(at(z)) <= 1e-8 * max(1, abs(g0)))
  }
  z <- z0
  for (k in seq_len(max_iter)) {
    gr <- gradient(z)
    s <- sum(gr * drop(V %*% gr))
    if (!all(is.finite(gr)) || !is.finite(s) || s <= 0) {
      break
    }
    z_new <- z0 - multiplier(z, gr) * drop(V %*% gr)
    if (!all(is.finite(z_new)) || !is.finite(at(z_new))) {
      break
    }
    change <- max(abs(z_new - z))
    z <- z_new
    if (converged(change, z)) {
      return(distance(z))
    }
  }
  n <- length(z0)
  residual <- function(z, lambda) {
    c(z - z0 + lambda * drop(V %*% gradient(z)), at(z))
  }
  h_outer <- 1e-4 * pmax(abs(z0), se, 1e-6)
  # Residuals small enough for the closest point (`loose` when the steps can
  # no longer reduce them, at the precision of the numerical derivatives).
  solved <- function(r, loose = 1) {
    max(abs(r[seq_len(n)])) <= loose * 1e-7 * max(se, 1e-12) &&
      abs(r[[n + 1L]]) <= loose * 1e-8 * max(1, abs(g0))
  }
  newton <- function(z) {
    gr <- gradient(z)
    lambda <- multiplier(z, gr)
    r <- residual(z, lambda)
    for (k in seq_len(max_iter)) {
      if (!all(is.finite(r)) || !is.finite(lambda)) {
        return(NA_real_)
      }
      if (solved(r)) {
        return(distance(z))
      }
      gr <- gradient(z)
      H <- vapply(seq_len(n), function(j) {
        up <- z
        down <- z
        up[[j]] <- z[[j]] + h_outer[[j]]
        down[[j]] <- z[[j]] - h_outer[[j]]
        (gradient(up) - gradient(down)) / (2 * h_outer[[j]])
      }, numeric(n))
      H <- (H + t(H)) / 2
      jacobian <- rbind(cbind(diag(n) + lambda * V %*% H, drop(V %*% gr)), c(gr, 0))
      delta <- tryCatch(solve(jacobian, -r), error = function(e) NULL)
      if (is.null(delta) || !all(is.finite(delta))) {
        return(NA_real_)
      }
      # Halve the step until the residuals decrease.
      alpha <- 1
      repeat {
        z_try <- z + alpha * delta[seq_len(n)]
        lambda_try <- lambda + alpha * delta[[n + 1L]]
        r_try <- residual(z_try, lambda_try)
        if (all(is.finite(r_try)) && sum(r_try^2) < sum(r^2)) {
          break
        }
        alpha <- alpha / 2
        if (alpha < 1 / 1024) {
          return(if (solved(r, loose = 10)) distance(z) else NA_real_)
        }
      }
      z <- z_try
      lambda <- lambda_try
      r <- r_try
    }
    if (solved(r)) distance(z) else NA_real_
  }
  out <- if (all(is.finite(z)) && is.finite(at(z))) newton(z) else NA_real_
  if (is.na(out)) {
    out <- newton(z0)
  }
  if (is.na(out)) {
    # Far from the estimate (several standard errors on the linearised
    # scale), where the surface may reach the edge of the parameter space,
    # the linearised distance is used; nearer, there is no result.
    gr <- gradient(z0)
    s <- sum(gr * drop(V %*% gr))
    linear <- if (is.finite(s) && s > 0) g0 / sqrt(s) else NA_real_
    out <- if (is.finite(linear) && abs(linear) >= 4) linear else NA_real_
  }
  out
}
