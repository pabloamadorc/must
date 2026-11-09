test_that("a full 2 x 2 factorial has one identifiable set and no missing cells", {
  d <- two_by_two_gaussian()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  expect_s3_class(des, "must_design")
  expect_equal(nrow(des$cells), 4L)
  expect_equal(des$cells$n_units, rep(5L, 4))
  expect_equal(nrow(des$missing_cells), 0L)
  expect_equal(nrow(des$sets), 1L)
  expect_true(des$sets$net_identifiable)
  expect_true(is.na(des$sets$pure_identifiable))
  expect_false(des$sets$shared)
  expect_equal(des$sets$min_units, 5L)
  expect_false(any(des$checks$severity == "problem"))
})

test_that("a ray design identifies only the diagonal sets, which share the control", {
  # Single compounds at 2 and 20 and the mixture only at matching doses
  ray <- data.frame(a = c(0, 2, 20, 0, 0, 2, 20), b = c(0, 0, 0, 2, 20, 2, 20))
  ray <- ray[rep(seq_len(nrow(ray)), each = 3), ]
  des <- describe_design(ray, c("a", "b"), control = c(a = "0", b = "0"))
  expect_equal(des$sets$levels, c("a=2 | b=2", "a=20 | b=20"))
  expect_true(all(des$sets$net_identifiable))
  expect_setequal(des$missing_cells$cell, c("a=2 | b=20", "a=20 | b=2"))
  control_use <- des$cell_use$n_sets[des$cell_use$cell == "a=0 | b=0"]
  expect_equal(control_use, 2L)
  expect_true(all(des$sets$shared))
  expect_true("shared_cells" %in% des$checks$check)
})

test_that("three stressors without pairwise cells: net identifiable, pure not", {
  d <- data.frame(n = c(0, 1, 0, 0, 1), s = c(0, 0, 1, 0, 1), f = c(0, 0, 0, 1, 1))
  d <- d[rep(seq_len(5), each = 2), ]
  des <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"))
  expect_equal(nrow(des$sets), 1L)
  expect_true(des$sets$net_identifiable)
  expect_false(des$sets$pure_identifiable)
  expect_match(des$sets$missing_for_pure, "n=1 | s=1 | f=0", fixed = TRUE)
  expect_true("pure_not_identifiable" %in% des$checks$check)
})

test_that("a full 2 x 2 x 2 factorial gives three pairs and one pure triple", {
  des <- describe_design(three_way_gaussian(), c("n", "s", "f"),
                         control = c(n = "0", s = "0", f = "0"), unit = "id")
  expect_equal(sort(des$sets$order), c(2L, 2L, 2L, 3L))
  expect_true(all(des$sets$net_identifiable))
  expect_true(des$sets$pure_identifiable[des$sets$order == 3L])
  # The control is used by the four sets
  expect_equal(des$cell_use$n_sets[des$cell_use$cell == "n=0 | s=0 | f=0"], 4L)
})

test_that("treatments that change within a unit are a problem and block fitting", {
  tv <- data.frame(u = rep(1:4, each = 2), t = rep(1:2, 4),
                   a = c(0, 0, 1, 1, 0, 0, 1, 1), b = c(0, 0, 0, 0, 1, 1, 1, 0),
                   y = 1:8)
  des <- describe_design(tv, c("a", "b"), control = c(a = "0", b = "0"), unit = "u", time = "t")
  expect_true("problem" %in% des$checks$severity)
  expect_error(fit_responses(tv, "y", des, residual_checks = FALSE, quiet = TRUE),
               class = "must_error_input")
})

test_that("repeated measures are detected and reported", {
  d <- repeated_counts()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank", time = "day")
  expect_equal(unique(des$cells$n_times), 3L)
  expect_equal(unique(des$cells$n_units), 4L)
  msg <- des$checks$message[des$checks$check == "repeated_measures"]
  expect_match(msg, "measured repeatedly")
})

test_that("unbalanced sampling over time is flagged", {
  d <- repeated_counts()
  d <- d[!(d$tank == "T1" & d$day == 21), ]
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "tank", time = "day")
  expect_true("unbalanced_time" %in% des$checks$check)
})

test_that("sets are formed within `by` groups", {
  d <- two_by_two_gaussian()
  d2 <- d
  d2$y <- d2$y + 1
  both <- rbind(cbind(d, pop = "urban"), cbind(d2, pop = "rural"))
  both$id <- seq_len(nrow(both))
  des <- describe_design(both, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "id", by = "pop")
  expect_equal(nrow(des$sets), 2L)
  expect_setequal(des$sets$pop, c("urban", "rural"))
})

test_that("results do not depend on the order of the rows", {
  d <- three_way_gaussian()
  shuffled <- d[rev(seq_len(nrow(d))), ]
  a <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"), unit = "id")
  b <- describe_design(shuffled, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                       unit = "id")
  expect_equal(a$cells, b$cells)
  expect_equal(a$sets, b$sets)
})

test_that("factors and text give the same design", {
  d <- two_by_two_gaussian()
  f <- d
  f$a <- factor(f$a)
  f$b <- factor(f$b)
  a <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"))
  b <- describe_design(f, c("a", "b"), control = c(a = "0", b = "0"))
  expect_equal(a$sets, b$sets)
})

test_that("invalid specifications give classed, informative errors", {
  d <- two_by_two_gaussian()
  expect_error(describe_design(d, c("a", "b")), class = "must_error_input")
  expect_error(describe_design(d, c("a", "b"), control = c(a = "0")), class = "must_error_input")
  expect_error(describe_design(d, c("a", "b"), control = c(a = "none", b = "0")),
               "not in the data", class = "must_error_input")
  expect_error(describe_design(d, "a", control = c(a = "0")), class = "must_error_input")
  expect_error(describe_design(d, c("a", "zz"), control = c(a = "0", zz = "0")),
               class = "must_error_input")
  expect_error(describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), time = "y"),
               "unit", class = "must_error_input")
  d$a[3] <- NA
  expect_error(describe_design(d, c("a", "b"), control = c(a = "0", b = "0")),
               "missing values", class = "must_error_input")
  e <- two_by_two_gaussian()
  names(e)[names(e) == "b"] <- "cell"
  expect_error(describe_design(e, c("a", "cell"), control = c(a = "0", cell = "0")),
               "Rename", class = "must_error_input")
})

test_that("a design without the needed cells has no identifiable set", {
  d <- data.frame(a = c(0, 1, 1), b = c(0, 0, 1))
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"))
  expect_false(des$sets$net_identifiable)
  expect_equal(des$sets$missing_for_net, "a=0 | b=1")
  expect_true("problem" %in% des$checks$severity)
})

test_that("print() summarises the design", {
  d <- two_by_two_gaussian()
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  expect_snapshot(print(des))
})

test_that("a stressor applied to one group per level is confounded with the groups", {
  d <- jars_in_rooms(rooms_per_level = 1L)
  des <- describe_design(d, c("temperature", "nutrients"),
                         control = c(temperature = "15", nutrients = "no"),
                         unit = "id", groups = "room")
  expect_equal(des$groups$structure, "confounded")
  expect_equal(des$groups$between, "temperature")
  expect_equal(des$groups$groups_per_combination, "1")
  expect_false(des$groups$random)
  expect_true(des$sets$confounded)
  check <- des$checks[des$checks$check == "groups", ]
  expect_equal(check$severity, "warning")
  expect_match(check$message, "cannot be separated from differences among groups of 'room'")
  expect_match(check$message, "No random intercept")
})

test_that("replicated whole plots are reported and not flagged as confounded", {
  d <- jars_in_rooms(rooms_per_level = 2L)
  des <- describe_design(d, c("temperature", "nutrients"),
                         control = c(temperature = "15", nutrients = "no"),
                         unit = "id", groups = "room")
  expect_equal(des$groups$structure, "whole_units")
  expect_equal(des$groups$groups_per_combination, "2")
  expect_true(des$groups$random)
  expect_false(des$sets$confounded)
  msg <- des$checks$message[des$checks$check == "groups"]
  expect_match(msg, "the replication is the number of groups")
})

test_that("blocks holding the treatment combinations are complete or incomplete", {
  d <- two_by_two_gaussian()
  d$block <- rep(paste0("B", 1:5), 4)
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"),
                         unit = "id", groups = "block")
  expect_equal(des$groups$structure, "blocks")
  expect_equal(des$groups$between, "")
  expect_true(des$groups$random)
  expect_match(des$checks$message[des$checks$check == "groups"], "complete blocks")
  incomplete <- d[!(d$block == "B1" & d$a == "1" & d$b == "1"), ]
  des2 <- describe_design(incomplete, c("a", "b"), control = c(a = "0", b = "0"),
                          unit = "id", groups = "block")
  expect_match(des2$checks$message[des2$checks$check == "groups"], "incomplete blocks")
  expect_equal(des2$groups$cells_per_group, "3-4")
})

test_that("a combination given to a single group is confounded even if each level is not", {
  # One chamber per treatment combination: each level of each stressor is in
  # two chambers, but every combination is in one
  d <- expand.grid(a = c("0", "1"), b = c("0", "1"), jar = 1:3, stringsAsFactors = FALSE)
  d$chamber <- paste0("C", d$a, d$b)
  d$id <- seq_len(nrow(d))
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id",
                         groups = "chamber")
  expect_equal(des$groups$between, "a, b")
  expect_equal(des$groups$structure, "confounded")
  expect_true(des$sets$confounded)
})

test_that("only the sets that involve the stressor applied to whole groups are confounded", {
  # Nitrate applied to whole blocks (one block per level); flow within blocks
  d <- expand.grid(nitrate = c("ambient", "high"), flow = c("ambient", "low"),
                   sediment = c("ambient", "added"), rep = 1:2, stringsAsFactors = FALSE)
  d$block <- paste0("B-", d$nitrate)
  d$id <- seq_len(nrow(d))
  des <- describe_design(d, c("nitrate", "flow", "sediment"),
                         control = c(nitrate = "ambient", flow = "ambient", sediment = "ambient"),
                         unit = "id", groups = "block")
  flagged <- stats::setNames(des$sets$confounded, des$sets$stressors)
  expect_false(flagged[["flow x sediment"]])
  expect_true(all(flagged[names(flagged) != "flow x sediment"]))
  expect_match(des$checks$message[des$checks$check == "groups"], "3 contrast set")
})

test_that("groups that repeat the experimental unit are not used", {
  d <- two_by_two_gaussian()
  d$tank <- paste0("T", d$id)
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id",
                         groups = "tank")
  expect_equal(des$groups$structure, "units")
  expect_false(des$groups$random)
  expect_equal(des$checks$severity[des$checks$check == "groups"], "warning")
})

test_that("invalid `groups` give classed errors", {
  d <- two_by_two_gaussian()
  d$block <- rep(1:5, 4)
  expect_error(describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), groups = "zz"),
               class = "must_error_input")
  expect_error(describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), groups = "a"),
               "two roles", class = "must_error_input")
  expect_error(describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), groups = 1),
               class = "must_error_input")
  d$block[2] <- NA
  expect_error(describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), groups = "block"),
               "missing values", class = "must_error_input")
})

test_that("messages show only the stressors away from control, with their `by` group", {
  d <- two_by_two_gaussian()
  d2 <- d
  both <- rbind(cbind(d, pop = "urban"), cbind(d2, pop = "rural"))
  both$id <- seq_len(nrow(both))
  both <- both[!(both$pop == "rural" & both$a == "1" & both$b == "1" & both$id > 36), ]
  des <- describe_design(both, c("a", "b"), control = c(a = "0", b = "0"), unit = "id",
                         by = "pop")
  msg <- des$checks$message[des$checks$check == "replication"]
  expect_match(msg, "'pop=rural: a=1 + b=1'", fixed = TRUE)
})

test_that("background = \"all\" adds sets with other stressors away from control", {
  d <- three_way_gaussian()
  base <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                          unit = "id")
  all_bg <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                            unit = "id", background = "all")
  expect_true(all(base$sets$background == ""))
  expect_equal(nrow(all_bg$sets), 7L)
  bg <- all_bg$sets[nzchar(all_bg$sets$background), ]
  expect_setequal(bg$background, c("n=1", "s=1", "f=1"))
  expect_true(all(bg$order == 2L))
  # The reference cell of n x s at f = 1 is f alone
  sc <- all_bg$set_cells[[bg$set_id[bg$background == "f=1"]]]
  expect_equal(sc$subset[[1L]], "reference")
  expect_equal(sc$cell[[1L]], "n=0 | s=0 | f=1")
  # The cell with f alone is now used by the sets that take it as reference
  expect_gt(all_bg$cell_use$n_sets[all_bg$cell_use$cell == "n=0 | s=0 | f=1"],
            base$cell_use$n_sets[base$cell_use$cell == "n=0 | s=0 | f=1"])
  expect_output(print(all_bg), "other stressors at control or not")
  limited <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                             unit = "id", background = "all", max_order = 2)
  expect_equal(nrow(limited$sets), 6L)
  expect_error(describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                               background = "some"), class = "must_error_input")
})

test_that("groups labelled within other groups are detected", {
  # Split plot: nitrate applied to whole plots, plots numbered 1-2 within blocks
  d <- expand.grid(block = c("A", "B", "C"), plot = c("1", "2"), sub = 1:2,
                   stringsAsFactors = FALSE)
  d$nitrate <- ifelse((d$plot == "1") == (d$block %in% c("A", "C")), "high", "low")
  d$warming <- ifelse(d$sub == 1, "no", "yes")
  d$id <- seq_len(nrow(d))
  des <- describe_design(d, c("nitrate", "warming"),
                         control = c(nitrate = "low", warming = "no"), unit = "id",
                         groups = c("block", "plot"))
  nested <- des$checks[des$checks$check == "nested_labels", ]
  expect_equal(nrow(nested), 1L)
  expect_equal(nested$severity, "warning")
  expect_match(nested$message, "'nitrate'")
  d$plot_id <- paste(d$block, d$plot)
  ok <- describe_design(d, c("nitrate", "warming"),
                        control = c(nitrate = "low", warming = "no"), unit = "id",
                        groups = c("block", "plot_id"))
  expect_false("nested_labels" %in% ok$checks$check)
  expect_equal(ok$groups$between[ok$groups$group == "plot_id"], "nitrate")
})

test_that("reused unit labels get a hint in the problem message", {
  d <- two_by_two_gaussian()
  d$tank <- rep(1:5, 4)
  des <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "tank")
  msg <- des$checks$message[des$checks$severity == "problem"]
  expect_match(msg, "paste(treatment, tank)", fixed = TRUE)
})
