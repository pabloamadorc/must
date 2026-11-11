# Report the methods and the decisions of an analysis

Writes, from the output of
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md),
the methods of the statistical analysis as paragraphs of text and a
table with every decision behind it: the design, the model of each
response, the null models and why they were chosen, the inference and,
when given, the classification of the interactions and the agreement
among null models. The text cites the methods used and lists their
references. Every value comes from the settings of the objects, so the
description matches the analysis that was run.

## Usage

``` r
report_interactions(x, classes = NULL, agreement = NULL, justification = NULL)
```

## Arguments

- x:

  The output of
  [`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md).

- classes:

  Optional: the output of
  [`classify_interactions()`](https://pabloamadorc.github.io/must/reference/classify_interactions.md)
  for `x`, to describe the classification of the interactions.

- agreement:

  Optional: the output of
  [`compare_null_models()`](https://pabloamadorc.github.io/must/reference/compare_null_models.md)
  for `x`, to describe the agreement among the null models.

- justification:

  Optional text, one or more sentences, that justifies the choice of
  null models; it is added to the methods as given. When the null models
  came from
  [`suggest_null_models()`](https://pabloamadorc.github.io/must/reference/suggest_null_models.md),
  the criteria behind them are described without it.

## Value

An object of class `must_report`, a list with `methods` (the paragraphs,
a character vector), `decisions` (a data frame with columns `step`,
`decision`, `argument`, `value`, `default` and `source`) and
`references` (the works cited, formatted, a character vector).
[`format()`](https://rdrr.io/r/base/format.html) gives the whole report
as lines of plain text or, with `format = "markdown"`, of markdown for
an R Markdown or Quarto document (in a chunk with `results = "asis"`);
[`print()`](https://rdrr.io/r/base/print.html) shows the plain text.

## Details

The methods have up to five paragraphs: the design (stressors, units,
times and the interactions the design identifies); the software and the
model of each response (for a model fitted with
[`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md),
the candidates for its error distribution, how one was chosen and the
random effects; for a model fitted outside must, its class and family);
the null models and why they were chosen (the criteria of
[`suggest_null_models()`](https://pabloamadorc.github.io/must/reference/suggest_null_models.md),
the `justification` given or, for the default, why it is the default);
how the interactions were estimated and tested; and, when `classes` or
`agreement` are given, how the interactions were classified and
compared. They are written in English, in the past tense, for the
methods section of a manuscript, and are meant to be edited: the
description of the experiment itself (organisms, exposure, measurements)
is left to the authors.

The table of decisions has one row per decision, with the function
(`step`) and argument that control it, the value used, the default and
its `source`: `"default"`, `"set"` (given by the user), `"chosen"` (by
must: the distribution selected by
[`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md),
a seed drawn at random), `"estimated"` (criteria of
[`suggest_null_models()`](https://pabloamadorc.github.io/must/reference/suggest_null_models.md)
evaluated from the effects estimated), `"suggest_null_models()"` (the
null models suggested) or `"recorded"` (versions of the software). The
versions are those recorded when the models were fitted and the
interactions estimated; objects created by earlier versions of must lack
them and the summary of the models, and the report then says less.

## See also

[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md),
[`suggest_null_models()`](https://pabloamadorc.github.io/must/reference/suggest_null_models.md),
[`classify_interactions()`](https://pabloamadorc.github.io/must/reference/classify_interactions.md),
[`compare_null_models()`](https://pabloamadorc.github.io/must/reference/compare_null_models.md).

## Examples

``` r
jars <- expand.grid(
  salinity = c("ambient", "high"),
  warming = c("no", "yes"),
  replicate = 1:5,
  stringsAsFactors = FALSE
)
jars$jar <- seq_len(nrow(jars))
mu <- 40 * ifelse(jars$salinity == "high", 0.5, 1) *
  ifelse(jars$warming == "yes", 0.5, 1)
jars$algae <- stats::qpois(c(0.2, 0.35, 0.5, 0.65, 0.8)[jars$replicate], mu)
design <- describe_design(jars, c("salinity", "warming"),
  control = c(salinity = "ambient", warming = "no"), unit = "jar")
fits <- fit_responses(jars, "algae", design, residual_checks = FALSE,
  seed = 1, quiet = TRUE)
est <- estimate_interactions(fits, null_models = c("AD", "MU"), n_sim = 2000,
  seed = 1)
report <- report_interactions(est, classes = classify_interactions(est))
report
#> <must_report>
#> Methods
#> 
#> The experiment combined 2 stressors, salinity (ambient and high, with ambient
#> as control) and warming (no and yes, with no as control), in a full factorial
#> design (4 treatment combinations), with 5 experimental units (jar) per
#> combination. The design identified one contrast set, a pairwise interaction.
#> 
#> The analysis was run with the R package must 0.2.0 (Amador Crespo 2026) in R
#> 4.6.2 (R Core Team 2026). The response (algae) was modelled with a generalized
#> linear model with one mean per treatment combination (a cell-means model),
#> fitted by maximum likelihood with glmmTMB 1.1.15.2 (Brooks et al. 2017). The
#> error distribution was chosen among candidates suited to the values of the
#> response (Poisson (log link) and negative binomial (NB2, log link)): the
#> candidate that converged with a positive-definite Hessian and had the lowest
#> AIC was selected, except that the Poisson was used only if the negative
#> binomial did not converge. The model selected was Poisson (log link). For the
#> inference, the selected model was refitted by restricted maximum likelihood
#> (REML).
#> 
#> Each interaction was tested against simple addition (AD) and the multiplicative
#> model (MU). Following Madge Pimentel et al. (2025), each interaction was
#> computed from the means predicted by the model (post-estimation inference), as
#> the deviation of the combination of stressors from the expectation of the null
#> model, as a difference on the scale of the response (AD) and as a difference on
#> the log scale (MU). Uncertainty was obtained from 2,000 draws of the fixed
#> effects from their estimated distribution (Krinsky & Robb 1986), which gave
#> percentile intervals and two-sided p-values. To correct for small samples, the
#> covariance came from the REML refit; the distributions stayed normal, because
#> the model has no dispersion parameter or random effects. Intervals were
#> computed at the 95% level, and p-values were adjusted for the false discovery
#> rate (Benjamini & Hochberg 1995) within each response, null model and estimate,
#> across contrast sets. Random numbers were generated with seed 1.
#> 
#> Interactions were classified with the directional scheme of Piggott et al.
#> (2015): each single, combined and interaction effect was coded as positive,
#> negative or not significant from its adjusted p-value (alpha = 0.05), and each
#> net interaction was classed as additive (no deviation detected), as a positive
#> or negative synergism or as a positive or negative antagonism (their Table 1).
#> 
#> Decisions
#> 
#> describe_design()
#>   Stressors and their control levels: salinity (control ambient), warming
#>     (control no) [stressors, control, set]
#>   Experimental unit: jar [unit, set]
#>   Sampling time: none [time, default]
#>   Groups of units (blocks, rooms, whole plots): none [groups, default]
#>   Factors within which interactions are estimated: none [by, default]
#>   Level of the stressors outside each set: control [background, default]
#>   Highest number of stressors in a set: 2 (all stressors) [max_order, default]
#> 
#> fit_responses()
#>   Candidate distributions (algae): Poisson (log link) and negative binomial
#>     (NB2, log link) (count, detected) [families, default]
#>   Observations that share a dispersion parameter: common [dispersion, default]
#>   Random effects: automatic: none [random, default]
#>   Limits of the optimiser: scaled iteration limits [control, default]
#>   Residual checks with DHARMa: FALSE [residual_checks, set]
#>   Distribution selected (algae): Poisson (log link) [chosen]
#> 
#> estimate_interactions()
#>   Null models: AD, MU [null_models, set]
#>   Inference: simulation [inference, default]
#>   Number of simulation draws: 2000 [n_sim, set]
#>   Small-sample correction: small_sample [correction, default]
#>   Confidence level: 0.95 [level, default]
#>   Adjustment of p-values: BH [p_adjust, default]
#>   Responses in which an increase means sensitivity: none [sensitivity, default]
#>   Equivalence margin: none [margin, default]
#>   Seed of the simulations: 1 [seed, set]
#> 
#> report_interactions()
#>   Justification of the null models: none [justification, default]
#> 
#> classify_interactions()
#>   Significance level of the classes: 0.05 [alpha, default]
#> 
#> software
#>   Version of must: 0.2.0 [recorded]
#>   Version of R: 4.6.2 [recorded]
#>   Version of glmmTMB: 1.1.15.2 [recorded]
#> 
#> References
#> 
#> Amador Crespo, P. (2026). must: Multiple-Stressor Interactions Under
#>   Alternative Null Models. R package version 0.2.0.
#>   <https://github.com/pabloamadorc/must>
#> 
#> Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A
#>   practical and powerful approach to multiple testing. Journal of the Royal
#>   Statistical Society: Series B, 57(1), 289–300.
#>   <https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>
#> 
#> Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg, C. W.,
#>   Nielsen, A., Skaug, H. J., Mächler, M., & Bolker, B. M. (2017). glmmTMB
#>   balances speed and flexibility among packages for zero-inflated generalized
#>   linear mixed modeling. The R Journal, 9(2), 378–400.
#>   <https://doi.org/10.32614/RJ-2017-066>
#> 
#> Krinsky, I., & Robb, A. L. (1986). On approximating the statistical properties
#>   of elasticities. The Review of Economics and Statistics, 68(4), 715–719.
#>   <https://doi.org/10.2307/1924536>
#> 
#> Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S. J.,
#>   Matthaei, C. D., Orr, J. A., Piggott, J. J., & Schäfer, R. B. (2025).
#>   Hypothesis-driven research on multiple stressors: An analytical framework for
#>   stressor interactions. Ecology and Evolution, 15, e71959.
#>   <https://doi.org/10.1002/ece3.71959>
#> 
#> Piggott, J. J., Townsend, C. R., & Matthaei, C. D. (2015). Reconceptualizing
#>   synergism and antagonism among multiple stressors. Ecology and Evolution,
#>   5(7), 1538–1547. <https://doi.org/10.1002/ece3.1465>
#> 
#> R Core Team (2026). R: A Language and Environment for Statistical Computing. R
#>   Foundation for Statistical Computing, Vienna, Austria.
#>   <https://doi.org/10.32614/R.manuals>
#> 
report$decisions[, c("decision", "value", "source")]
#>                                            decision
#> 1                Stressors and their control levels
#> 2                                 Experimental unit
#> 3                                     Sampling time
#> 4      Groups of units (blocks, rooms, whole plots)
#> 5   Factors within which interactions are estimated
#> 6           Level of the stressors outside each set
#> 7              Highest number of stressors in a set
#> 8                   Candidate distributions (algae)
#> 9    Observations that share a dispersion parameter
#> 10                                   Random effects
#> 11                          Limits of the optimiser
#> 12                      Residual checks with DHARMa
#> 13                    Distribution selected (algae)
#> 14                                      Null models
#> 15                                        Inference
#> 16                       Number of simulation draws
#> 17                          Small-sample correction
#> 18                                 Confidence level
#> 19                           Adjustment of p-values
#> 20 Responses in which an increase means sensitivity
#> 21                               Equivalence margin
#> 22                          Seed of the simulations
#> 23                 Justification of the null models
#> 24                Significance level of the classes
#> 25                                  Version of must
#> 26                                     Version of R
#> 27                               Version of glmmTMB
#>                                                                         value
#> 1                            salinity (control ambient), warming (control no)
#> 2                                                                         jar
#> 3                                                                        none
#> 4                                                                        none
#> 5                                                                        none
#> 6                                                                     control
#> 7                                                           2 (all stressors)
#> 8  Poisson (log link) and negative binomial (NB2, log link) (count, detected)
#> 9                                                                      common
#> 10                                                            automatic: none
#> 11                                                    scaled iteration limits
#> 12                                                                      FALSE
#> 13                                                         Poisson (log link)
#> 14                                                                     AD, MU
#> 15                                                                 simulation
#> 16                                                                       2000
#> 17                                                               small_sample
#> 18                                                                       0.95
#> 19                                                                         BH
#> 20                                                                       none
#> 21                                                                       none
#> 22                                                                          1
#> 23                                                                       none
#> 24                                                                       0.05
#> 25                                                                      0.2.0
#> 26                                                                      4.6.2
#> 27                                                                   1.1.15.2
#>      source
#> 1       set
#> 2       set
#> 3   default
#> 4   default
#> 5   default
#> 6   default
#> 7   default
#> 8   default
#> 9   default
#> 10  default
#> 11  default
#> 12      set
#> 13   chosen
#> 14      set
#> 15  default
#> 16      set
#> 17  default
#> 18  default
#> 19  default
#> 20  default
#> 21  default
#> 22      set
#> 23  default
#> 24  default
#> 25 recorded
#> 26 recorded
#> 27 recorded
# Markdown, for an R Markdown document (in a chunk with results = "asis")
cat(format(report, format = "markdown"), sep = "\n")
#> The experiment combined 2 stressors, salinity (ambient and high, with ambient as control) and warming (no and yes, with no as control), in a full factorial design (4 treatment combinations), with 5 experimental units (jar) per combination. The design identified one contrast set, a pairwise interaction.
#> 
#> The analysis was run with the R package must 0.2.0 (Amador Crespo 2026) in R 4.6.2 (R Core Team 2026). The response (algae) was modelled with a generalized linear model with one mean per treatment combination (a cell-means model), fitted by maximum likelihood with glmmTMB 1.1.15.2 (Brooks et al. 2017). The error distribution was chosen among candidates suited to the values of the response (Poisson (log link) and negative binomial (NB2, log link)): the candidate that converged with a positive-definite Hessian and had the lowest AIC was selected, except that the Poisson was used only if the negative binomial did not converge. The model selected was Poisson (log link). For the inference, the selected model was refitted by restricted maximum likelihood (REML).
#> 
#> Each interaction was tested against simple addition (AD) and the multiplicative model (MU). Following Madge Pimentel et al. (2025), each interaction was computed from the means predicted by the model (post-estimation inference), as the deviation of the combination of stressors from the expectation of the null model, as a difference on the scale of the response (AD) and as a difference on the log scale (MU). Uncertainty was obtained from 2,000 draws of the fixed effects from their estimated distribution (Krinsky & Robb 1986), which gave percentile intervals and two-sided p-values. To correct for small samples, the covariance came from the REML refit; the distributions stayed normal, because the model has no dispersion parameter or random effects. Intervals were computed at the 95% level, and p-values were adjusted for the false discovery rate (Benjamini & Hochberg 1995) within each response, null model and estimate, across contrast sets. Random numbers were generated with seed 1.
#> 
#> Interactions were classified with the directional scheme of Piggott et al. (2015): each single, combined and interaction effect was coded as positive, negative or not significant from its adjusted p-value (alpha = 0.05), and each net interaction was classed as additive (no deviation detected), as a positive or negative synergism or as a positive or negative antagonism (their Table 1).
#> 
#> | Step | Decision | Argument | Value | Default | Source |
#> |---|---|---|---|---|---|
#> | `describe_design()` | Stressors and their control levels | `stressors, control` | salinity (control ambient), warming (control no) | required | set |
#> | `describe_design()` | Experimental unit | `unit` | jar | none | set |
#> | `describe_design()` | Sampling time | `time` | none | none | default |
#> | `describe_design()` | Groups of units (blocks, rooms, whole plots) | `groups` | none | none | default |
#> | `describe_design()` | Factors within which interactions are estimated | `by` | none | none | default |
#> | `describe_design()` | Level of the stressors outside each set | `background` | control | control | default |
#> | `describe_design()` | Highest number of stressors in a set | `max_order` | 2 (all stressors) | all stressors | default |
#> | `fit_responses()` | Candidate distributions (algae) | `families` | Poisson (log link) and negative binomial (NB2, log link) (count, detected) | chosen from the values | default |
#> | `fit_responses()` | Observations that share a dispersion parameter | `dispersion` | common | common | default |
#> | `fit_responses()` | Random effects | `random` | automatic: none | automatic | default |
#> | `fit_responses()` | Limits of the optimiser | `control` | scaled iteration limits | scaled iteration limits | default |
#> | `fit_responses()` | Residual checks with DHARMa | `residual_checks` | FALSE | TRUE | set |
#> | `fit_responses()` | Distribution selected (algae) |  | Poisson (log link) |  | chosen |
#> | `estimate_interactions()` | Null models | `null_models` | AD, MU | AD | set |
#> | `estimate_interactions()` | Inference | `inference` | simulation | simulation | default |
#> | `estimate_interactions()` | Number of simulation draws | `n_sim` | 2000 | 10000 | set |
#> | `estimate_interactions()` | Small-sample correction | `correction` | small_sample | small_sample | default |
#> | `estimate_interactions()` | Confidence level | `level` | 0.95 | 0.95 | default |
#> | `estimate_interactions()` | Adjustment of p-values | `p_adjust` | BH | BH | default |
#> | `estimate_interactions()` | Responses in which an increase means sensitivity | `sensitivity` | none | none | default |
#> | `estimate_interactions()` | Equivalence margin | `margin` | none | none | default |
#> | `estimate_interactions()` | Seed of the simulations | `seed` | 1 | drawn at random | set |
#> | `report_interactions()` | Justification of the null models | `justification` | none | none | default |
#> | `classify_interactions()` | Significance level of the classes | `alpha` | 0.05 | 0.05 | default |
#> | software | Version of must |  | 0.2.0 |  | recorded |
#> | software | Version of R |  | 4.6.2 |  | recorded |
#> | software | Version of glmmTMB |  | 1.1.15.2 |  | recorded |
#> 
#> **References**
#> 
#> Amador Crespo, P. (2026). *must: Multiple-Stressor Interactions Under Alternative Null Models*. R package version 0.2.0. <https://github.com/pabloamadorc/must>
#> 
#> Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300. <https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>
#> 
#> Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg, C. W., Nielsen, A., Skaug, H. J., Mächler, M., & Bolker, B. M. (2017). glmmTMB balances speed and flexibility among packages for zero-inflated generalized linear mixed modeling. *The R Journal*, 9(2), 378–400. <https://doi.org/10.32614/RJ-2017-066>
#> 
#> Krinsky, I., & Robb, A. L. (1986). On approximating the statistical properties of elasticities. *The Review of Economics and Statistics*, 68(4), 715–719. <https://doi.org/10.2307/1924536>
#> 
#> Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S. J., Matthaei, C. D., Orr, J. A., Piggott, J. J., & Schäfer, R. B. (2025). Hypothesis-driven research on multiple stressors: An analytical framework for stressor interactions. *Ecology and Evolution*, 15, e71959. <https://doi.org/10.1002/ece3.71959>
#> 
#> Piggott, J. J., Townsend, C. R., & Matthaei, C. D. (2015). Reconceptualizing synergism and antagonism among multiple stressors. *Ecology and Evolution*, 5(7), 1538–1547. <https://doi.org/10.1002/ece3.1465>
#> 
#> R Core Team (2026). *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing, Vienna, Austria. <https://doi.org/10.32614/R.manuals>
#> 
```
