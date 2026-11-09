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
