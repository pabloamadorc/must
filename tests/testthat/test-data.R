test_that("polazzo_zooplankton matches its documentation", {
  d <- polazzo_zooplankton
  expect_equal(dim(d), c(120L, 8L))
  expect_named(d, c("mesocosm", "day", "insecticide", "herbicide", "nutrients",
                    "total_density", "daphnia", "cyclopoida"))
  expect_false(anyNA(d))
  expect_setequal(unique(d$day), c(-5L, 7L, 14L, 30L, 50L))
  # Full factorial, three mesocosms per combination, one treatment per mesocosm
  cells <- unique(d[c("mesocosm", "insecticide", "herbicide", "nutrients")])
  expect_equal(nrow(cells), 24L)
  expect_true(all(table(paste(cells$insecticide, cells$herbicide, cells$nutrients)) == 3L))
  expect_true(all(d[c("total_density", "daphnia", "cyclopoida")] >= 0))
})
