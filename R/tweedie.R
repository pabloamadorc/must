# Tweedie distributions with a power between 1 and 2 (internal) ---------------
#
# For 1 < p < 2 a Tweedie variable with mean mu and dispersion phi is a
# compound Poisson-gamma sum: N ~ Poisson(lambda) amounts, each gamma with
# shape alpha and scale gamma, where lambda = mu^(2 - p) / (phi (2 - p)),
# alpha = (2 - p) / (p - 1) and gamma = phi (p - 1) mu^(p - 1) (Dunn & Smyth
# 2005). Its mean is mu, its variance phi mu^p and P(Y = 0) = exp(-lambda).
# The distribution function is the Poisson mixture of gamma distribution
# functions, summed here over the numbers of amounts that hold all but 1e-12
# of the Poisson probability. With more than `tweedie_lambda_max()` expected
# amounts the series grows long and the distribution is close to the gamma
# with its mean and variance (shape mu^(2 - p) / phi, scale phi mu^(p - 1)),
# which replaces it: they differ by less than 0.003 in probability, and
# P(Y = 0) is below exp(-1000).

tweedie_lambda_max <- function() 1000

tweedie_parameters <- function(mu, phi, power) {
  list(
    lambda = mu^(2 - power) / (phi * (2 - power)),
    alpha = (2 - power) / (power - 1),
    scale = phi * (power - 1) * mu^(power - 1),
    gamma_shape = mu^(2 - power) / phi,
    gamma_scale = phi * mu^(power - 1)
  )
}

# Terms of the series for each element: numbers of amounts n from `lo` to
# `hi` (at least one), with their log Poisson weights.
tweedie_terms <- function(lambda) {
  lo <- pmax(1, stats::qpois(1e-12, lambda))
  hi <- pmax(lo, stats::qpois(1e-12, lambda, lower.tail = FALSE))
  len <- hi - lo + 1
  element <- rep.int(seq_along(lambda), len)
  start <- cumsum(c(0, len[-length(len)]))
  n <- lo[element] + seq_along(element) - 1 - start[element]
  list(element = element, n = n, log_weight = stats::dpois(n, lambda[element], log = TRUE),
       index = split(seq_along(element), element))
}

recycle_tweedie <- function(x, mu, phi, power) {
  n <- length(x)
  list(x = x, mu = rep_len(mu, n), phi = rep_len(phi, n), power = rep_len(power, n))
}

# Distribution function at q >= 0.
ptweedie_series <- function(q, mu, phi, power) {
  a <- recycle_tweedie(q, mu, phi, power)
  prm <- tweedie_parameters(a$mu, a$phi, a$power)
  out <- exp(-prm$lambda)
  long <- prm$lambda > tweedie_lambda_max()
  out[long] <- stats::pgamma(a$x[long], shape = prm$gamma_shape[long],
                             scale = prm$gamma_scale[long])
  pos <- which(a$x > 0 & !long)
  if (length(pos) > 0L) {
    out[pos] <- out[pos] + tweedie_sums(a$x[pos], prm, pos, tweedie_terms(prm$lambda[pos]),
                                        seq_along(pos))$cdf
  }
  pmin(out, 1)
}

# Sums of the series for the elements `active` (positions in `pos`) at the
# values x (one per active element): the distribution function without the
# probability of zero, and the density.
tweedie_sums <- function(x, prm, pos, terms, active) {
  idx <- unlist(terms$index[active], use.names = FALSE)
  e <- terms$element[idx]
  at <- match(e, active)
  shape <- terms$n[idx] * prm$alpha[pos][e]
  scale <- prm$scale[pos][e]
  w <- terms$log_weight[idx]
  cdf <- rowsum(exp(w + stats::pgamma(x[at], shape = shape, scale = scale, log.p = TRUE)), at,
                reorder = TRUE)
  dens <- rowsum(exp(w + stats::dgamma(x[at], shape = shape, scale = scale, log = TRUE)), at,
                 reorder = TRUE)
  list(cdf = as.vector(cdf), density = as.vector(dens))
}

# Quantile function: zero up to P(Y = 0), above it the root of the
# distribution function, by Newton steps kept inside a bracket (bisection
# when a step leaves it).
qtweedie_series <- function(p, mu, phi, power, tol = 1e-10, max_iter = 200L) {
  a <- recycle_tweedie(p, mu, phi, power)
  prm <- tweedie_parameters(a$mu, a$phi, a$power)
  p0 <- exp(-prm$lambda)
  out <- numeric(length(p))
  long <- prm$lambda > tweedie_lambda_max()
  out[long] <- stats::qgamma(a$x[long], shape = prm$gamma_shape[long],
                             scale = prm$gamma_scale[long])
  pos <- which(a$x > p0 & !long)
  if (length(pos) == 0L) {
    return(out)
  }
  terms <- tweedie_terms(prm$lambda[pos])
  target <- a$x[pos] - p0[pos]
  k <- length(pos)
  sd <- sqrt(a$phi[pos] * a$mu[pos]^a$power[pos])
  lo <- numeric(k)
  hi <- a$mu[pos] + 10 * sd
  for (i in seq_len(200L)) {
    low <- which(tweedie_sums(hi, prm, pos, terms, seq_len(k))$cdf < target)
    if (length(low) == 0L) {
      break
    }
    hi[low] <- 2 * hi[low]
  }
  y <- a$mu[pos] + sd * stats::qnorm(pmin(a$x[pos], 1 - 1e-12))
  outside <- !is.finite(y) | y <= lo | y >= hi
  y[outside] <- (lo[outside] + hi[outside]) / 2
  active <- seq_len(k)
  for (iter in seq_len(max_iter)) {
    s <- tweedie_sums(y[active], prm, pos, terms, active)
    gap <- s$cdf - target[active]
    below <- gap < 0
    lo[active[below]] <- y[active[below]]
    hi[active[!below]] <- y[active[!below]]
    done <- abs(gap) < tol | (hi[active] - lo[active]) <= 1e-13 * y[active]
    step <- y[active] - gap / s$density
    bad <- !is.finite(step) | step <= lo[active] | step >= hi[active]
    step[bad] <- (lo[active[bad]] + hi[active[bad]]) / 2
    y[active[!done]] <- step[!done]
    active <- active[!done]
    if (length(active) == 0L) {
      break
    }
  }
  out[pos] <- y
  out
}

# Unit deviance of a Tweedie distribution with power 1 < p < 2.
tweedie_unit_deviance <- function(y, mu, power) {
  mu <- rep_len(mu, length(y))
  d <- 2 * mu^(2 - power) / (2 - power)
  pos <- y > 0
  if (any(pos)) {
    yp <- y[pos]
    mp <- mu[pos]
    d[pos] <- 2 * (yp^(2 - power) / ((1 - power) * (2 - power)) -
                     yp * mp^(1 - power) / (1 - power) + mp^(2 - power) / (2 - power))
  }
  pmax(d, 0)
}

# Family object for stats::glm.fit(): Tweedie with power `power` and a log
# link, as statmod::tweedie(var.power = power, link.power = 0).
tweedie_glm_family <- function(power) {
  force(power)
  structure(list(
    family = "Tweedie",
    link = "log",
    linkfun = function(mu) log(mu),
    linkinv = function(eta) pmax(exp(eta), .Machine$double.eps),
    mu.eta = function(eta) pmax(exp(eta), .Machine$double.eps),
    valideta = function(eta) TRUE,
    variance = function(mu) mu^power,
    validmu = function(mu) all(mu >= 0),
    dev.resids = function(y, mu, wt) wt * tweedie_unit_deviance(y, mu, power),
    aic = function(y, n, mu, wt, dev) NA_real_,
    initialize = expression({
      n <- rep(1, nobs)
      mustart <- y + 0.1 * (y == 0)
    }),
    power = power
  ), class = "family")
}

# Probability integral transform residuals: F(y) for positive values, a
# uniform draw below P(Y = 0) for zeros.
tweedie_pit <- function(y, mu, phi, power) {
  u <- ptweedie_series(y, mu, phi, power)
  zero <- y == 0
  u[zero] <- stats::runif(sum(zero)) * u[zero]
  pmin(pmax(u, 0), 1 - 1e-10)
}
