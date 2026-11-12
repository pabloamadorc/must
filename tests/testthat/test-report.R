# Interactions estimated from a 2 x 2 factorial of counts, five jars per cell.
report_fixture <- function(null_models = c("AD", "MU"), ...) {
  d <- two_by_two_counts()
  design <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "count", design, residual_checks = FALSE, seed = 1, quiet = TRUE)
  estimate_interactions(fit, null_models = null_models, n_sim = 500, seed = 2, ...)
}

methods_text <- function(r) paste(r$methods, collapse = " ")

decision <- function(r, argument, column = "value", step = NULL) {
  rows <- r$decisions$argument == argument
  if (!is.null(step)) rows <- rows & r$decisions$step == step
  r$decisions[[column]][rows]
}

test_that("the report describes the design, the model and the inference that were run", {
  est <- report_fixture()
  r <- report_interactions(est)
  expect_s3_class(r, "must_report")
  expect_named(r, c("methods", "decisions", "references"))
  expect_named(r$decisions, c("step", "decision", "argument", "value", "default", "source"))
  text <- methods_text(r)
  expect_match(text, "2 stressors, a (0 and 1, with 0 as control)", fixed = TRUE)
  expect_match(text, "full factorial design (4 treatment combinations)", fixed = TRUE)
  expect_match(text, "5 experimental units (id) per combination", fixed = TRUE)
  expect_match(text, "one contrast set, a pairwise interaction", fixed = TRUE)
  # The model selected by fit_responses().
  family <- est$settings$fit$responses[[1L]]$family_name
  expect_true(family %in% c("poisson", "nbinom2"))
  expect_match(text, sprintf("The model selected was %s.", family_label(family, "log")), fixed = TRUE)
  expect_match(text, "glmmTMB", fixed = TRUE)
  # The settings of estimate_interactions().
  expect_match(text, "simple addition (AD) and the multiplicative model (MU)", fixed = TRUE)
  expect_match(text, "from 500 draws of the fixed effects", fixed = TRUE)
  expect_match(text, "at the 95% level", fixed = TRUE)
  expect_match(text, "false discovery rate (Benjamini & Hochberg 1995)", fixed = TRUE)
  expect_match(text, "seed 2.", fixed = TRUE)
  expect_match(text, "Madge Pimentel et al. (2025)", fixed = TRUE)
  # Every work cited is in the references, and only those.
  expect_true(any(grepl("^Madge Pimentel, I\\.", r$references)))
  expect_true(any(grepl("^Brooks, M\\. E\\.", r$references)))
  expect_true(any(grepl("^Krinsky, I\\.", r$references)))
  expect_false(any(grepl("^Liess", r$references)))
  expect_identical(r$references, r$references[order(tolower(strip_markup(r$references)),
                                                    method = "radix")])
})

test_that("the table of decisions records the values used and where they came from", {
  est <- report_fixture()
  r <- report_interactions(est)
  expect_identical(decision(r, "n_sim"), "500")
  expect_identical(decision(r, "n_sim", "default"), "10000")
  expect_identical(decision(r, "n_sim", "source"), "set")
  expect_identical(decision(r, "level", "source"), "default")
  expect_identical(decision(r, "null_models"), "AD, MU")
  expect_identical(decision(r, "seed", "source"), "set")
  expect_identical(decision(r, "unit"), "id")
  expect_identical(decision(r, "residual_checks"), "FALSE")
  selected <- r$decisions[r$decisions$decision == "Distribution selected (count)", ]
  expect_identical(selected$source, "chosen")
  expect_true(all(c("must", "R", "glmmTMB") %in%
                    sub("^Version of ", "", r$decisions$decision[r$decisions$step == "software"])))
  # A seed drawn at random is chosen, not set.
  d <- two_by_two_counts()
  design <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"))
  fit <- fit_responses(d, "count", design, residual_checks = FALSE, quiet = TRUE)
  expect_false(fit$settings$seed_given)
  drawn <- report_interactions(estimate_interactions(fit, n_sim = 200))
  expect_identical(decision(drawn, "seed", "source", step = "estimate_interactions()"), "chosen")
  # Arguments changed from their defaults are "set" in the table and in the text.
  other <- report_interactions(report_fixture(level = 0.9, p_adjust = "holm",
                                              correction = "none", inference = "delta"))
  expect_identical(decision(other, "p_adjust"), "holm")
  expect_identical(decision(other, "correction", "source"), "set")
  text <- methods_text(other)
  expect_match(text, "at the 90% level", fixed = TRUE)
  expect_match(text, "method of Holm (1979)", fixed = TRUE)
  expect_match(text, "No correction for small samples was applied", fixed = TRUE)
  expect_match(text, "delta method", fixed = TRUE)
  # The delta method draws no random numbers: no seed in the text or the table.
  expect_false(grepl("seed", text, fixed = TRUE))
  expect_length(decision(other, "seed", step = "estimate_interactions()"), 0L)
})

test_that("the reasons for the null models come from suggest_null_models() or the user", {
  d <- two_by_two_counts()
  design <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(d, "count", design, residual_checks = FALSE, seed = 1, quiet = TRUE)
  choice <- suggest_null_models(fit, similar_action = FALSE, specific_action = TRUE,
                                mortality = FALSE, stressor_effect = FALSE)
  est <- estimate_interactions(fit, null_models = choice, n_sim = 300, seed = 1)
  r <- report_interactions(est, justification = "Both stressors reduce growth.")
  text <- methods_text(r)
  expect_match(text, "guidelines of Schäfer & Piggott (2018)", fixed = TRUE)
  expect_match(text, "the stressors do not act similarly", fixed = TRUE)
  expect_match(text, "They did not state whether the distribution of sensitivities is known",
               fixed = TRUE)
  expect_match(text, "evaluated from the effects estimated by the models", fixed = TRUE)
  expect_match(text, "Both stressors reduce growth.", fixed = TRUE)
  expect_identical(decision(r, "null_models", "source"), "suggest_null_models()")
  expect_identical(decision(r, "similar_action"), "no")
  expect_identical(decision(r, "similar_action", "source"), "set")
  expect_identical(decision(r, "correlation", "source"), "default")
  expect_identical(decision(r, "justification"), "given")
  expect_true("estimated" %in% r$decisions$source[r$decisions$step == "suggest_null_models()"])
  expect_true(any(grepl("^Schäfer, R\\. B\\., & Piggott", r$references)))
  # The default null model, without a justification, is explained as the default.
  default <- report_interactions(report_fixture(null_models = "AD"))
  expect_match(methods_text(default), "default null model of must", fixed = TRUE)
  # Null models chosen by the user without a reason get none.
  plain <- report_interactions(report_fixture())
  expect_false(grepl("default null model", methods_text(plain), fixed = TRUE))
  expect_identical(decision(plain, "justification"), "none")
})

test_that("classes and agreement are described with the significance level used", {
  est <- report_fixture(null_models = c("AD", "MU", "DO"))
  classes <- classify_interactions(est, alpha = 0.1)
  agreement <- compare_null_models(est, alpha = 0.01)
  expect_identical(attr(classes, "alpha"), 0.1)
  expect_identical(attr(agreement, "alpha"), 0.01)
  expect_identical(attr(agreement, "null_models"), c("AD", "MU", "DO"))
  r <- report_interactions(est, classes = classes, agreement = agreement)
  text <- methods_text(r)
  expect_match(text, "directional scheme of Piggott et al. (2015)", fixed = TRUE)
  expect_match(text, "(alpha = 0.1)", fixed = TRUE)
  expect_match(text, "Dominance, which does not predict the combination", fixed = TRUE)
  expect_match(text, "null models compared (AD, MU and DO)", fixed = TRUE)
  expect_match(text, "adjusted p < 0.01", fixed = TRUE)
  expect_identical(decision(r, "alpha", step = "classify_interactions()"), "0.1")
  expect_identical(decision(r, "alpha", "source", step = "compare_null_models()"), "set")
  expect_true(any(grepl("^Piggott, J\\. J\\.", r$references)))
  # Outputs that lost their attributes are still described.
  bare <- classes
  attr(bare, "alpha") <- NULL
  r2 <- report_interactions(est, classes = bare)
  expect_identical(decision(r2, "alpha", "source", step = "classify_interactions()"), "recorded")
})

test_that("format() gives plain text or markdown", {
  r <- report_interactions(report_fixture(), classes = classify_interactions(report_fixture()))
  text <- format(r, width = 60L)
  expect_identical(text[[1L]], "Methods")
  expect_true(all(c("Decisions", "References", "describe_design()", "software") %in% text))
  long <- nchar(text) > 60L & !grepl("<http", text, fixed = TRUE)
  expect_false(any(long))
  expect_false(any(grepl("*", text, fixed = TRUE)))
  md <- format(r, format = "markdown")
  expect_true("| Step | Decision | Argument | Value | Default | Source |" %in% md)
  expect_true(any(startsWith(md, "| `fit_responses()` | Distribution selected (count) |")))
  expect_true("**References**" %in% md)
  expect_true(any(grepl("*Ecology and Evolution*", md, fixed = TRUE)))
  expect_identical(md[[1L]], r$methods[[1L]])
  expect_output(print(r), "<must_report>", fixed = TRUE)
  expect_error(format(r, format = "html"), class = "must_error_input")
})

test_that("repeated measures, baselines, groups and margins are described", {
  d <- three_way_repeated()
  design <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                            unit = "tank", time = "day", baseline = list(f = 0))
  fit <- fit_responses(d, "y", design, residual_checks = FALSE, seed = 3, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("AD", "MU", "BL"), n_sim = 300, seed = 1,
                               margin = 0.25)
  text <- methods_text(report_interactions(est))
  expect_match(text, "sampled at 3 times (day 0, 10 and 20)", fixed = TRUE)
  expect_match(text, "not tested at the samplings before it was applied (f: day 0)", fixed = TRUE)
  expect_match(text, "3 pairwise interactions and 1 three-way interaction", fixed = TRUE)
  expect_match(text, "one mean per treatment combination and sampling time", fixed = TRUE)
  expect_match(text, "a random intercept for each tank", fixed = TRUE)
  expect_match(text, "random effects at zero", fixed = TRUE)
  expect_match(text, "emergent interaction of Beppler et al. (2016)", fixed = TRUE)
  expect_match(text, "across contrast sets and sampling times", fixed = TRUE)
  expect_match(text, "lay between -0.25 and 0.25", fixed = TRUE)
  expect_match(text, sprintf("experimental units (%s)", format(unique(est$estimates$df[
    est$estimates$component == "interaction" & is.finite(est$estimates$df)]))), fixed = TRUE)
  expect_identical(decision(report_interactions(est), "baseline"), "f: 0")
  # A baseline common to every stressor is described once.
  common <- describe_design(d, c("n", "s", "f"), control = c(n = "0", s = "0", f = "0"),
                            unit = "tank", time = "day", baseline = 0)
  common_fit <- fit_responses(d, "y", common, residual_checks = FALSE, seed = 3, quiet = TRUE)
  common_report <- report_interactions(estimate_interactions(common_fit, n_sim = 300, seed = 1))
  expect_match(common_report$methods[[1L]],
               "The sampling before the stressors were applied (day 0) was not tested.",
               fixed = TRUE)
  expect_identical(decision(common_report, "baseline"), "0")
  # Stressors applied to whole groups of units.
  rooms <- jars_in_rooms(rooms_per_level = 2L)
  room_design <- describe_design(rooms, c("temperature", "nutrients"),
                                 control = c(temperature = "15", nutrients = "no"),
                                 unit = "id", groups = "room")
  room_fit <- fit_responses(rooms, "y", room_design, residual_checks = FALSE, seed = 1,
                            quiet = TRUE)
  room_text <- methods_text(report_interactions(estimate_interactions(room_fit, n_sim = 300,
                                                                      seed = 1)))
  expect_match(room_text, "grouped by room (4 groups), and temperature was applied to whole groups",
               fixed = TRUE)
  expect_match(room_text, "a random intercept for each room", fixed = TRUE)
})

test_that("proportions, sensitivity and combinations with only zeros are described", {
  jars <- survival_jars()
  design <- describe_design(jars, c("a", "b"), control = c(a = "0", b = "0"), unit = "id")
  fit <- fit_responses(jars, "alive", design, trials = c(alive = "n"), residual_checks = FALSE,
                       seed = 1, quiet = TRUE)
  est <- estimate_interactions(fit, null_models = c("AD", "BL", "SAM"), n_sim = 300, seed = 1,
                               sensitivity = "alive")
  r <- report_interactions(est)
  text <- methods_text(r)
  expect_match(text, "successes out of a number of trials (alive out of n)", fixed = TRUE)
  expect_match(text, "binomial (logit link) and beta-binomial (logit link)", fixed = TRUE)
  expect_match(text, "applied to its complement (1 - p)", fixed = TRUE)
  expect_match(text, "both shape parameters equal to 3.2", fixed = TRUE)
  expect_match(text, "relative to the reference combination (BL and SAM)", fixed = TRUE)
  expect_identical(decision(r, "sensitivity"), "alive")
  expect_identical(decision(r, "trials"), "alive: n")
  zeros <- zero_combination_counts()
  zero_design <- describe_design(zeros, c("a", "b"), control = c(a = "0", b = "0"))
  zero_fit <- fit_responses(zeros, "y", zero_design, residual_checks = FALSE, seed = 1,
                            quiet = TRUE)
  zero_text <- methods_text(report_interactions(estimate_interactions(zero_fit, n_sim = 300,
                                                                      seed = 1)))
  expect_match(zero_text, "The mean of the treatment combination with only zeros was fixed",
               fixed = TRUE)
})

test_that("a model fitted outside must is described by its class and emmeans", {
  skip_if_not_installed("emmeans")
  d <- two_by_two_counts()
  design <- describe_design(d, c("a", "b"), control = c(a = "0", b = "0"))
  model <- stats::glm(count ~ a * b, family = stats::quasipoisson, data = d)
  est <- estimate_interactions(model, design, null_models = c("AD", "MU"), n_sim = 300, seed = 1)
  r <- report_interactions(est)
  text <- methods_text(r)
  expect_match(text, "modelled outside must, with a model of class glm", fixed = TRUE)
  expect_match(text, "estimated marginal means from emmeans", fixed = TRUE)
  expect_match(text, "residual degrees of freedom of the model (16)", fixed = TRUE)
  expect_false(grepl("glmmTMB", text, fixed = TRUE))
  expect_identical(decision(r, "fit"), "glm: quasipoisson (log)")
  expect_identical(decision(r, "df"), "residual degrees of freedom of the model")
  expect_true("emmeans" %in% names(est$settings$versions))
  expect_true(any(grepl("^Lenth", r$references)))
})

test_that("interactions estimated by earlier versions of must are still reported", {
  est <- report_fixture()
  est$settings$fit <- NULL
  est$settings$versions <- NULL
  est$settings$seed_given <- NULL
  r <- report_interactions(est)
  text <- methods_text(r)
  expect_match(text, "The analysis was run with the R package must (", fixed = TRUE)
  expect_match(text, "Poisson (log link) for count", fixed = TRUE)
  expect_identical(decision(r, "seed", "source"), "recorded")
  expect_false("fit_responses()" %in% r$decisions$step)
})

test_that("report_interactions() checks its inputs", {
  est <- report_fixture()
  expect_error(report_interactions(est$estimates), class = "must_error_input")
  expect_error(report_interactions(est, classes = est$estimates), class = "must_error_input")
  other <- report_fixture(null_models = c("AD", "BL"))
  expect_error(report_interactions(est, classes = classify_interactions(other)),
               "the same responses and null models", class = "must_error_input")
  expect_error(report_interactions(est, agreement = compare_null_models(other)),
               class = "must_error_input")
  expect_error(report_interactions(est, justification = 1), class = "must_error_input")
  expect_error(report_interactions(est, justification = NA_character_),
               class = "must_error_input")
  expect_error(report_interactions(est, justification = "  "), class = "must_error_input")
  # Several strings are joined into one text.
  r <- report_interactions(est, justification = c("First.", "Second."))
  expect_match(methods_text(r), "First. Second.", fixed = TRUE)
})

test_that("the pieces of the report are written correctly", {
  expect_identical(family_label(c("nbinom2", "gamma_log", "beta"), c("log", "log", "logit")),
                   c("negative binomial (NB2, log link)", "Gamma (log link)", "beta (logit link)"))
  expect_identical(settings_family_label("poisson (log)"), "Poisson (log link)")
  expect_identical(default_text(estimate_interactions.must_fit, "n_sim"), "10000")
  expect_identical(default_text(estimate_interactions.must_fit, "seed"), "none")
  expect_identical(default_text(describe_design, "control"), "required")
  expect_identical(show_value(c(0.5, 1)), "0.5, 1")
  expect_identical(or_list(c("a", "b", "c")), "a, b or c")
  expect_identical(sort_times(c("20", "5", "10")), c("5", "10", "20"))
  expect_identical(strip_markup("*Ecology* and _Evolution_ of dry_mass"),
                   "Ecology and Evolution of dry_mass")
  refs <- new_citations()
  expect_identical(refs$cite("holm1979"), "Holm 1979")
  expect_identical(refs$citet("holm1979"), "Holm (1979)")
  expect_length(refs$formatted(), 1L)
  expect_error(refs$cite("nobody2000"), class = "must_error_internal")
})

test_that("manuals are cited with their organization as text and their DOI", {
  title <- "R: A Language and Environment for Statistical Computing"
  expect_identical(
    bib_text(utils::person("R Foundation for Statistical Computing",
                           comment = c(ROR = "05qewa988"))),
    "R Foundation for Statistical Computing"
  )
  expect_identical(bib_text("Vienna, Austria"), "Vienna, Austria")
  expect_null(bib_text(NULL))
  old <- utils::bibentry("Manual", title = title, author = utils::person("R Core Team"),
                         organization = "R Foundation for Statistical Computing",
                         address = "Vienna, Austria", year = "2024",
                         url = "https://www.R-project.org/")
  expect_identical(bib_reference(old), list(
    who = "R Core Team", year = "2024",
    full = paste("R Core Team (2024). *R: A Language and Environment for Statistical",
                 "Computing*. R Foundation for Statistical Computing, Vienna, Austria.",
                 "<https://www.R-project.org/>")
  ))
  with_doi <- utils::bibentry("Manual", title = title, author = utils::person("R Core Team"),
                              year = "2026", doi = "10.32614/R.manuals",
                              url = "https://www.R-project.org/")
  expect_match(bib_reference(with_doi)$full, "Computing*. <https://doi.org/10.32614/R.manuals>",
               fixed = TRUE)
  # The citation of the R that runs the tests (from R 4.6, with a ROR
  # identifier for its organization and a DOI)
  r <- report_reference("R")
  expect_match(r$full, "R Foundation for Statistical Computing, Vienna, Austria.",
               fixed = TRUE)
  expect_false(grepl("ROR|05qewa988", r$full))
})

test_that("must is cited with both of its authors", {
  ref <- report_reference("must")
  expect_identical(ref$who, "Amador Crespo & Rico")
  expect_match(ref$full, "^Amador Crespo, P\\., & Rico, A\\. \\(2026\\)\\. \\*must: ")
})

test_that("the organization of the citation of R from R 4.6 is cited by its name", {
  new <- utils::bibentry(
    "Manual", title = "R: A Language and Environment for Statistical Computing",
    author = utils::person("R Core Team", comment = c(ROR = "02zz1nj61")),
    organization = utils::person("R Foundation for Statistical Computing",
                                 comment = c(ROR = "05qewa988")),
    address = "Vienna, Austria", year = "2026", doi = "10.32614/R.manuals",
    url = "https://www.R-project.org/"
  )
  if (!inherits(new$organization, "person")) {
    skip("before R 4.6, bibentry() turns a person in other fields into text")
  }
  expect_identical(bib_reference(new)$full, paste(
    "R Core Team (2026). *R: A Language and Environment for Statistical Computing*.",
    "R Foundation for Statistical Computing, Vienna, Austria.",
    "<https://doi.org/10.32614/R.manuals>"
  ))
})
