# Small deterministic data sets used by several test files.

# Balanced 2 x 2 factorial with Gaussian-like values: cell means are known
# exactly because the deviations within each cell sum to zero.
two_by_two_gaussian <- function(means = c(10, 8, 7, 3), dev = c(-1, -0.5, 0, 0.5, 1)) {
  cells <- data.frame(a = c("0", "1", "0", "1"), b = c("0", "0", "1", "1"),
                      mean = means, stringsAsFactors = FALSE)
  d <- cells[rep(seq_len(4), each = length(dev)), ]
  d$y <- d$mean + rep(dev, 4)
  d$id <- seq_len(nrow(d))
  rownames(d) <- NULL
  d
}

# Balanced 2 x 2 factorial with counts; cell totals 100, 50, 70 and 25.
two_by_two_counts <- function() {
  d <- data.frame(
    a = rep(c("0", "1", "0", "1"), each = 5),
    b = rep(c("0", "0", "1", "1"), each = 5),
    stringsAsFactors = FALSE
  )
  d$count <- c(18, 20, 22, 19, 21, 9, 11, 10, 12, 8, 15, 13, 14, 16, 12, 4, 6, 5, 7, 3)
  d$id <- seq_len(nrow(d))
  d
}

# Full 2 x 2 x 2 factorial, three replicates per cell, Gaussian-like values.
three_way_gaussian <- function() {
  grid <- expand.grid(n = c("0", "1"), s = c("0", "1"), f = c("0", "1"),
                      stringsAsFactors = FALSE)
  na <- grid$n == "1"
  sa <- grid$s == "1"
  fa <- grid$f == "1"
  grid$mean <- 20 - 3 * na - 4 * sa - 2 * fa + 1.5 * na * sa - 1 * na * fa +
    0.5 * sa * fa - 2.5 * na * sa * fa
  d <- grid[rep(seq_len(8), each = 3), ]
  d$y <- d$mean + rep(c(-0.6, 0.1, 0.5), 8)
  d$id <- seq_len(nrow(d))
  rownames(d) <- NULL
  d
}

# 2 x 2 factorial with four tanks per cell sampled at three times.
repeated_counts <- function() {
  tanks <- expand.grid(a = c("0", "1"), b = c("0", "1"), rep = 1:4,
                       stringsAsFactors = FALSE)
  tanks$tank <- paste0("T", seq_len(nrow(tanks)))
  d <- tanks[rep(seq_len(nrow(tanks)), times = 3), ]
  d$day <- rep(c(7, 14, 21), each = nrow(tanks))
  base <- 30 * ifelse(d$a == "1", 0.6, 1) * ifelse(d$b == "1", 0.5, 1)
  tank_effect <- rep(c(0.9, 1.1, 1, 1.05), each = 4)[match(d$tank, tanks$tank)]
  d$count <- round(base * tank_effect * rep(c(1, 0.9, 1.1), each = nrow(tanks)))
  rownames(d) <- NULL
  d
}

# Twelve stressors in the layout of a number-of-factors experiment: 20
# controls, 8 units of each stressor alone and 139 combinations of two or
# three stressors with a single unit each. Gaussian-like values with a small
# spread: with 153 cell means the optimiser needs more than glmmTMB's default
# 300 iterations.
many_single_cells <- function() {
  k <- 12L
  label <- function(on) {
    v <- rep("0", k)
    v[on] <- "1"
    v
  }
  rows <- c(
    rep(list(label(integer(0))), 20L),
    rep(lapply(seq_len(k), label), each = 8L),
    lapply(utils::combn(k, 2L, simplify = FALSE), label),
    lapply(utils::combn(k, 3L, simplify = FALSE)[seq(1L, 220L, by = 3L)], label)
  )
  d <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(d) <- paste0("s", seq_len(k))
  d$id <- seq_len(nrow(d))
  z <- stats::qnorm(stats::ppoints(nrow(d)))[order(sin(seq_len(nrow(d)) * 7.3))]
  d$y <- 0.65 + 0.025 * z
  d
}

# 2 x 2 design with jars inside rooms. `rooms_per_level` rooms hold each
# temperature (temperature is applied to whole rooms); nutrients vary within
# rooms. Two jars per nutrient level and room.
jars_in_rooms <- function(rooms_per_level = 1L) {
  d <- expand.grid(temperature = c("15", "20"), room_rep = seq_len(rooms_per_level),
                   nutrients = c("no", "yes"), jar = 1:2, stringsAsFactors = FALSE)
  d$room <- paste0("R", d$temperature, "-", d$room_rep)
  d$id <- seq_len(nrow(d))
  room_effect <- c(0.4, -0.3, 0.2, -0.1)[match(d$room, sort(unique(d$room)))]
  d$y <- 10 - 2 * (d$temperature == "20") - 1 * (d$nutrients == "yes") +
    room_effect + c(-0.2, 0.2)[d$jar]
  d
}

# 2 x 2 factorial, four tanks per cell sampled on six consecutive days, with
# first-order autoregressive noise (rho = 0.7) within each tank.
autocorrelated_gaussian <- function() {
  withr::with_seed(11, {
    tanks <- expand.grid(a = c("0", "1"), b = c("0", "1"), rep = 1:4,
                         stringsAsFactors = FALSE)
    tanks$tank <- paste0("T", seq_len(nrow(tanks)))
    noise <- vapply(seq_len(nrow(tanks)), function(i) {
      as.numeric(stats::arima.sim(list(ar = 0.7), n = 6))
    }, numeric(6))
  })
  d <- tanks[rep(seq_len(nrow(tanks)), times = 6), ]
  d$day <- rep(1:6, each = nrow(tanks))
  d$y <- 10 - 2 * (d$a == "1") - (d$b == "1") +
    0.5 * noise[cbind(d$day, match(d$tank, tanks$tank))]
  rownames(d) <- NULL
  d
}

# 2 x 2 counts in which the combination has only zeros (the response is
# eliminated); the other cells average 20, 15 and 12.
zero_combination_counts <- function() {
  d <- data.frame(a = rep(c("0", "1", "0", "1"), each = 5),
                  b = rep(c("0", "0", "1", "1"), each = 5), stringsAsFactors = FALSE)
  d$y <- c(18, 22, 20, 19, 21, 14, 16, 15, 13, 17, 11, 13, 12, 10, 14, 0, 0, 0, 0, 0)
  d$id <- seq_len(nrow(d))
  d
}

# 2 x 2 survival experiment: `alive` of `n` animals per jar.
survival_jars <- function(n = 10L) {
  d <- data.frame(a = rep(c("0", "1", "0", "1"), each = 5),
                  b = rep(c("0", "0", "1", "1"), each = 5), stringsAsFactors = FALSE)
  d$alive <- c(9, 10, 9, 8, 9, 7, 8, 6, 7, 7, 8, 7, 8, 9, 8, 1, 2, 0, 1, 1)
  d$n <- n
  d$id <- seq_len(nrow(d))
  d
}

# Dose-response data whose cell means follow concentration addition exactly:
# log-logistic curves with EC50 = 2 and slope 1.5 (a) and EC50 = 5 and slope 2
# (b).
ca_factorial <- function(a = c(0, 1, 2, 4), b = c(0, 2, 5, 10)) {
  curve_a <- function(x) 1 / (1 + (x / 2)^1.5)
  curve_b <- function(x) 1 / (1 + (x / 5)^2)
  ca <- function(xa, xb) {
    if (xa == 0) return(curve_b(xb))
    if (xb == 0) return(curve_a(xa))
    eq <- function(w) xa / (2 * ((1 - w) / w)^(1 / 1.5)) + xb / (5 * ((1 - w) / w)^(1 / 2)) - 1
    stats::uniroot(eq, c(1e-9, 1 - 1e-9), tol = 1e-13)$root
  }
  g <- expand.grid(a = a, b = b)
  g$w <- mapply(ca, g$a, g$b)
  d <- g[rep(seq_len(nrow(g)), each = 4), ]
  d$y <- 100 * d$w * rep(c(0.97, 0.99, 1.01, 1.03), nrow(g))
  d$id <- seq_len(nrow(d))
  rownames(d) <- NULL
  d
}

# The same two stressors under two exposure regimes (pulse and press) with a
# control shared by both (regime "both"); four units per cell.
shared_regimes <- function() {
  cells <- data.frame(
    regime = c("both", "pulse", "pulse", "pulse", "press", "press", "press"),
    a = c("0", "1", "0", "1", "1", "0", "1"),
    b = c("0", "0", "1", "1", "0", "1", "1"),
    mean = c(10, 8, 7, 3, 9, 6, 5),
    stringsAsFactors = FALSE
  )
  d <- cells[rep(seq_len(7), each = 4), ]
  d$y <- d$mean + rep(c(-0.6, -0.2, 0.3, 0.5), 7)
  d$id <- seq_len(nrow(d))
  rownames(d) <- NULL
  d
}

# Full 2 x 2 x 2 factorial, three tanks per cell sampled on days 0, 10 and 20.
three_way_repeated <- function() {
  tanks <- expand.grid(n = c("0", "1"), s = c("0", "1"), f = c("0", "1"), rep = 1:3,
                       stringsAsFactors = FALSE)
  tanks$tank <- paste0("T", seq_len(nrow(tanks)))
  d <- tanks[rep(seq_len(nrow(tanks)), times = 3), ]
  d$day <- rep(c(0, 10, 20), each = nrow(tanks))
  na <- d$n == "1"
  sa <- d$s == "1"
  fa <- d$f == "1" & d$day > 0
  tank_effect <- round(sin(seq_len(nrow(tanks)) * 1.7), 2)[match(d$tank, tanks$tank)]
  d$y <- 20 - 3 * na - 4 * sa - 2 * fa + 1.5 * na * sa + tank_effect +
    round(sin(seq_len(nrow(d)) * 2.3) * 0.8, 2)
  rownames(d) <- NULL
  d
}
