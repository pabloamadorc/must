# plot() of must_interactions: what each type draws (returned invisibly),
# the labels and classes, the restored settings and the errors. Figures go
# to a null device.

local_null_device <- function(env = parent.frame()) {
  grDevices::pdf(NULL)
  device <- grDevices::dev.cur()
  withr::defer(grDevices::dev.off(device), envir = env)
}

counts_estimates <- function(null_models = c("AD", "MU", "DO")) {
  d <- two_by_two_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  estimate_interactions(fit, null_models = null_models, n_sim = 500, seed = 1)
}

repeated_estimates <- function() {
  d <- repeated_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank", time = "day")
  fit <- fit_responses(d, "count", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  estimate_interactions(fit, null_models = c("AD", "MU"), inference = "delta")
}

three_way_estimates <- function() {
  d <- three_way_gaussian()
  des <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                         unit = "id")
  fit <- fit_responses(d, "y", des, families = list(gaussian = stats::gaussian()),
                       residual_checks = FALSE, quiet = TRUE)
  estimate_interactions(fit, null_models = c("AD", "DO"), n_sim = 500, seed = 1)
}

test_that("the default plot shows how each interaction arises from the single effects", {
  est <- counts_estimates()
  local_null_device()
  out <- expect_invisible(plot(est))
  expect_identical(out$component, c("single", "single", "expected", "combined", "interaction"))
  expect_identical(out$stressor, c("a", "b", NA, NA, NA))
  # The first null model estimated that gives single effects
  expect_true(all(out$null_model == "AD"))
  e <- est$estimates[est$estimates$null_model == "AD", ]
  expect_equal(out$estimate[out$component == "single"], e$estimate[e$component == "single"])
  component <- function(x) out$estimate[out$component == x]
  expected <- component("combined") - component("interaction")
  expect_equal(out$estimate[out$component == "expected"], expected)
  expect_true(is.na(out$lower[out$component == "expected"]))
  cls <- classify_interactions(est)
  code <- unname(plot_codes[cls$class[cls$null_model == "AD"]])
  expect_identical(out$class_code, c("", "", "", "", code))
  expect_true(all(out$page == 1L))
  # MU on its own scale
  mu <- plot(est, null_model = "MU")
  expect_equal(mu$estimate[mu$component == "interaction"],
               est$estimates$estimate[est$estimates$null_model == "MU" &
                                         est$estimates$component == "interaction"])
})

test_that("pure interactions are drawn without an expected combined effect and without class", {
  est <- three_way_estimates()
  local_null_device()
  out <- plot(est, estimand = "pure")
  expect_identical(out$component, c("single", "single", "single", "combined", "interaction"))
  expect_identical(out$class_code[out$component == "interaction"], "n/a")
  # Net interactions of the four sets, three per page
  net <- plot(est, per_page = 3)
  expect_equal(sum(net$component == "interaction"), 4L)
  expect_identical(sort(unique(net$page)), 1:2)
  # Without `per_page`, as many as fit in the height of the device
  grDevices::pdf(NULL, width = 7, height = 4)
  small <- plot(est)
  grDevices::dev.off()
  expect_gt(max(small$page), 1L)
  expect_error(plot(counts_estimates("AD"), estimand = "pure"), "three or more stressors",
               class = "must_error_input")
})

test_that("type = \"interactions\" draws one row per interaction, on either scale", {
  est <- repeated_estimates()
  local_null_device()
  out <- plot(est, type = "interactions")
  e <- est$estimates
  int <- e[e$null_model == "AD" & e$component == "interaction", ]
  expect_equal(nrow(out), nrow(int))
  # Sampling times in numeric order, then their average
  expect_identical(out$time, c("7", "14", "21", "overall"))
  rel <- plot(est, type = "interactions", null_model = "MU", scale = "relative")
  mu <- e[e$null_model == "MU" & e$component == "interaction", ]
  expect_equal(rel$estimate, mu$relative[match(rel$time, mu$time)])
  expect_equal(rel$lower, mu$relative_lower[match(rel$time, mu$time)])
  grouped <- plot(est, type = "interactions", group_by = "set", time = c("21", "7"))
  expect_identical(grouped$time, c("7", "21"))
  # The dominance model draws interactions, which have no class
  do <- plot(counts_estimates(), type = "interactions", null_model = "DO")
  expect_identical(do$class_code, "n/a")
})

test_that("type = \"models\" draws the relative deviation under each null model", {
  est <- counts_estimates()
  local_null_device()
  out <- plot(est, type = "models")
  expect_identical(out$null_model, c("AD", "MU", "DO"))
  int <- est$estimates[est$estimates$component == "interaction", ]
  expect_equal(out$estimate, int$relative[match(out$null_model, int$null_model)])
  expect_identical(out$class_code[out$null_model == "DO"], "n/a")
  cls <- classify_interactions(est)
  expect_identical(out$class_code[out$null_model == "MU"],
                   unname(plot_codes[cls$class[cls$null_model == "MU"]]))
  two <- plot(est, type = "models", null_model = c("MU", "AD"), classes = "none")
  expect_identical(two$null_model, c("MU", "AD"))
  # One block per set and time; four per page by default
  rep_est <- repeated_estimates()
  expect_true(all(plot(rep_est, type = "models")$page == 1L))
  pages <- plot(rep_est, type = "models", per_page = 3)
  expect_identical(pages$page, rep(1:2, c(6L, 2L)))
})

test_that("type = \"time\" draws each set over the sampling times, without the average", {
  est <- repeated_estimates()
  local_null_device()
  out <- plot(est, type = "time", null_model = "MU")
  expect_false("overall" %in% out$time)
  expect_identical(unique(out$component), c("single", "expected", "combined", "interaction"))
  expect_equal(sum(out$component == "interaction"), 3L)
  int <- out[out$component == "interaction", ]
  expect_true(all(int$class_code %in% c(plot_codes, "NC")))
  expect_true(all(out$page == 1L))
  expect_error(plot(est, type = "time", time = "overall"), "other than \"overall\"",
               class = "must_error_input")
  expect_error(plot(counts_estimates(), type = "time"), "sampling times",
               class = "must_error_input")
})

test_that("combinations with only zeros and infinite estimates are drawn", {
  d <- zero_combination_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "y", des, families = list(poisson = stats::poisson()),
                       residual_checks = FALSE, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("AD", "MU"), n_sim = 500, seed = 1)
  local_null_device()
  out <- plot(est, null_model = "MU")
  expect_true(any(out$status == "boundary"))
  expect_true(any(is.infinite(out$estimate)))
  expect_silent(plot(est, type = "interactions", null_model = "MU"))
  expect_silent(plot(est, type = "models"))
})

test_that("sets are labelled by the levels of their stressors, or as given", {
  est <- three_way_estimates()
  d <- est$estimates[est$estimates$component == "interaction" &
                       est$estimates$estimand == "net" & est$estimates$null_model == "AD", ]
  opts <- list(label_style = "full", stressor_labels = NULL, labels = NULL)
  full <- plot_set_labels(d, opts)
  expect_identical(full[d$order == 3L], "n: 1 + s: 1 + f: 1")
  levels_only <- plot_set_labels(d, modifyList(opts, list(label_style = "levels")))
  expect_identical(levels_only[d$order == 3L], "1 + 1 + 1")
  named <- modifyList(opts, list(stressor_labels = c(n = "Nitrate"),
                                 labels = stats::setNames("All three", d$set_id[d$order == 3L])))
  out <- plot_set_labels(d, named)
  expect_identical(out[d$order == 3L], "All three")
  expect_match(out[d$order == 2L & grepl("n=1", d$levels) & grepl("s=1", d$levels)],
               "^Nitrate: 1 \\+ s: 1$")
  short <- list(label_style = "short", stressor_labels = c(n = "nitrate"), labels = NULL)
  expect_match(plot_set_labels(d, short)[d$order == 3L], "^nitr: 1 \\+ s: 1 \\+ f: 1$")
  # The background, in brackets
  bg <- d[1L, ]
  bg$stressors <- "s x f"
  bg$background <- "n=1"
  expect_identical(plot_set_labels(bg, opts), "s: 1 + f: 1 [n: 1]")
})

test_that("the rows, wrapped labels and headers of a page do not overlap", {
  line <- 0.16
  n_text <- c(1, 2, 2, 1, 1)
  n_head <- c(2, 0, 0, 1, 0)
  pos <- forest_layout(n_text, n_head, height = 3, line = line)
  half <- function(n) n * line / 2 / (3 / pos$top)
  top <- pos$y + half(n_text)
  bottom <- pos$y - half(n_text)
  expect_true(all(bottom[-5] >= top[-1]))
  heads <- which(n_head > 0)
  expect_true(all(pos$y_head[heads] - half(n_head[heads]) >= top[heads]))
  expect_true(pos$y_head[[4]] + half(n_head[[4]]) <= bottom[[3]])
  expect_true(pos$y_head[[1]] + half(n_head[[1]]) <= pos$top)
  # Too many rows for the device: the regular spacing
  crowded <- forest_layout(rep(1, 40), c(1, rep(0, 39)), height = 2, line = line)
  expect_equal(-diff(crowded$y), rep(1, 39))
})

test_that("labels are wrapped between stressors rather than inside them", {
  expect_identical(
    plot_wrap("insecticide: chlorpyrifos + herbicide: diuron + nutrients: added", 42),
    "insecticide: chlorpyrifos +\nherbicide: diuron + nutrients: added"
  )
  # A separator that ends a line may overhang it
  expect_identical(plot_wrap("species: T. gravida \u00b7 temperature: 9", 20),
                   "species: T. gravida \u00b7\ntemperature: 9")
  # A part longer than a line is wrapped at its spaces
  expect_identical(plot_wrap("a very long level name: indeed", 12),
                   "a very long\nlevel name:\nindeed")
  expect_identical(plot_wrap(c("", NA), 10), c("", ""))
})

test_that("plot() rejects arguments it does not use, naming plot() in the error", {
  est <- counts_estimates()
  local_null_device()
  err <- expect_error(plot(est, col = "red"), "Unknown argument\\(s\\): 'col'",
                      class = "must_error_input")
  expect_identical(conditionCall(err)[[1L]], as.name("plot"))
  expect_error(plot(est, type = "forest"), class = "must_error_input")
  expect_error(plot(est, null_model = "DO"), "dominance model", class = "must_error_input")
  expect_error(plot(est, type = "time", null_model = "DO"), "dominance model",
               class = "must_error_input")
  expect_error(plot(est, null_model = "BL"), class = "must_error_input")
  expect_error(plot(est, type = "models", null_model = c("AD", "BL")), "estimated",
               class = "must_error_input")
  expect_error(plot(est, response = "other"), class = "must_error_input")
  expect_error(plot(est, scale = "relative"), "applies to type = \"interactions\"",
               class = "must_error_input")
  expect_error(plot(est, group_by = "set"), "applies to type = \"interactions\"",
               class = "must_error_input")
  expect_error(plot(est, type = "models", classes = "names"), "as codes",
               class = "must_error_input")
  expect_error(plot(est, ylim = c(0, 1)), "applies to type = \"time\"",
               class = "must_error_input")
  expect_error(plot(est, xlim = c(1, 0)), "increasing order", class = "must_error_input")
  expect_error(plot(est, per_page = 1.5), "whole number", class = "must_error_input")
  expect_error(plot(est, alpha = 1), class = "must_error_input")
  expect_error(plot(est, set_id = "S09"), "Available: 'S01'", class = "must_error_input")
  expect_error(plot(est, context = "site=A"), "have no `context`", class = "must_error_input")
  expect_error(plot(est, labels = c(S09 = "x")), "set identifiers", class = "must_error_input")
  expect_error(plot(est, labels = "x"), "unique names", class = "must_error_input")
  expect_error(plot(est, stressor_labels = c(z = "Z")), "must be stressors",
               class = "must_error_input")
  rep_est <- repeated_estimates()
  expect_error(plot(rep_est, type = "time", estimand = "pure"), "net interactions",
               class = "must_error_input")
  expect_error(plot(rep_est, type = "time", per_page = 2), "one set per page",
               class = "must_error_input")
  expect_error(plot(rep_est, type = "time", xlim = c(0, 1)), "ylim",
               class = "must_error_input")
  expect_error(plot(rep_est, time = "99"), "Available: '7', '14', '21', 'overall'",
               class = "must_error_input")
})

test_that("plot() restores the graphical parameters and the prompt between pages", {
  est <- repeated_estimates()
  local_null_device()
  old <- graphics::par(c("mar", "mai", "mfrow", "cex", "las", "oma", "mgp"))
  ask <- grDevices::devAskNewPage()
  for (type in c("effects", "interactions", "models", "time")) {
    plot(est, type = type, main = type)
    expect_equal(graphics::par(names(old)), old)
    expect_identical(grDevices::devAskNewPage(), ask)
  }
})
