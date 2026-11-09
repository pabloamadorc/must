test_that("stressor levels put the control first and sort numbers numerically", {
  expect_equal(stressor_levels(c(20, 2, 200, 0), "0"), c("0", "2", "20", "200"))
  expect_equal(stressor_levels(c("high", "ambient", "low"), "ambient"),
               c("ambient", "high", "low"))
  expect_equal(stressor_levels(factor(c("b", "a"), levels = c("b", "a")), "a"), c("a", "b"))
})

test_that("cell and context labels name every stressor", {
  lv <- data.frame(t = c("20", "25"), n = c("0", "1"), stringsAsFactors = FALSE)
  expect_equal(make_cell_label(lv, c("t", "n")), c("t=20 | n=0", "t=25 | n=1"))
  expect_equal(make_cell_label(list(t = "20", n = "1"), c("t", "n")), "t=20 | n=1")
  expect_equal(make_context_label(data.frame(pop = "urban"), "pop"), "pop=urban")
  expect_equal(make_context_label(data.frame(x = 1:2), NULL), c("", ""))
})

test_that("subsets are ordered by size, control first and full set last", {
  s <- all_subsets(3)
  expect_equal(length(s), 8L)
  expect_equal(s[[1]], integer(0))
  expect_equal(s[2:4], list(1L, 2L, 3L))
  expect_equal(s[[8]], 1:3)
})

test_that("small helpers format and bind as expected", {
  expect_equal(format_range(c(3, 3)), "3")
  expect_equal(format_range(c(2, 5, 3)), "2-5")
  expect_equal(backtick("a b"), "`a b`")
  template <- data.frame(x = numeric())
  expect_identical(bind_rows_base(list(NULL), template), template)
  expect_equal(nrow(bind_rows_base(list(data.frame(x = 1), data.frame(x = 2)))), 2L)
})

test_that("helpers that paste return nothing for empty input", {
  expect_identical(backtick(character()), character())
  expect_identical(group_column(character()), character())
  expect_identical(short_cells(character(), c(a = "0")), character())
})
