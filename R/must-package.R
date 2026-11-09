#' must: multiple-stressor interactions under alternative null models
#'
#' Whether a combined effect is called additive, synergistic or antagonistic
#' depends on the null model that defines "no interaction" and on the scale
#' on which effects are combined. must makes these choices explicit and
#' separates three questions:
#'
#' 1. **What the design can identify.** [describe_design()] reports missing
#'    treatment combinations, cells shared by several contrasts, the
#'    replication unit and repeated measures, and which interactions are
#'    identifiable, before any model is fitted.
#' 2. **Which model fits the data.** [fit_responses()] fits a generalised
#'    linear mixed model to each response and chooses its error distribution
#'    by convergence, residual checks and AIC, independently of the null
#'    models.
#' 3. **What "non-additive" means.** [estimate_interactions()] tests the
#'    simple addition, multiplicative and dominance null models from the
#'    predictions of that model, and distinguishes the net deviation of a
#'    combination of three or more stressors from its pure highest-order
#'    interaction. By default the tests are corrected for small samples;
#'    `correction = "none"` gives the asymptotic procedure of Madge Pimentel
#'    et al. (2025). [classify_interactions()] assigns the directional classes
#'    of Piggott et al. (2015), and `summary()` and `plot()` give an overview.
#'
#' @references
#' Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S.
#' J., Matthaei, C. D., Orr, J. A., Piggott, J. J., & Schäfer, R. B. (2025).
#' Hypothesis-driven research on multiple stressors: An analytical framework
#' for stressor interactions. *Ecology and Evolution*, 15, e71959.
#' \doi{10.1002/ece3.71959}
#'
#' Orr, J. A., et al. (2020). Towards a unified study of multiple stressors:
#' Divisions and common goals across research disciplines. *Proceedings of
#' the Royal Society B*, 287, 20200421. \doi{10.1098/rspb.2020.0421}
#'
#' Piggott, J. J., Townsend, C. R., & Matthaei, C. D. (2015).
#' Reconceptualizing synergism and antagonism among multiple stressors.
#' *Ecology and Evolution*, 5(7), 1538–1547. \doi{10.1002/ece3.1465}
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
