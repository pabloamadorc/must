# Classify interactions (Piggott et al. 2015)

Codes the single, combined and interaction effects of each identifiable
contrast set as positive, negative or not significant, and assigns the
directional class of Piggott et al. (2015, Table 1) under the
simple-addition (AD), multiplicative (MU), Bliss independence (BL),
stress addition (SAM) and concentration addition (CA) null models.

## Usage

``` r
classify_interactions(x, alpha = 0.05)
```

## Arguments

- x:

  A `must_interactions` object from
  [`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md).

- alpha:

  Significance level applied to the adjusted p-values (`p_adjusted`).
  Default 0.05.

## Value

A data frame with one row per response, contrast set, sampling time and
null model: `response`, `context` and the `by` columns, `time`,
`set_id`, `order`, `stressors`, `levels`, `background`, `null_model`,
`scale`, `estimand` (`"net"`), `stressor_a`, `code_a`, `stressor_b`,
`code_b` (two-stressor sets only), `single_codes` (the code of every
single effect, `"salinity +; warming 0"`), `code_combined`,
`code_interaction` (`"+"`, `"-"`, `"0"`, or `"++"` / `"--"` for the
strong cases defined by Piggott et al.), `pattern`, `class` (`"AD"`,
`"+S"`, `"-S"`, `"+A"`, `"-A"`), `meaning`, `verdict` (with an
equivalence `margin` in
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md):
`"deviation"`, `"negligible"`, `"equivalent"` or `"inconclusive"`; `NA`
otherwise), `borderline` (some adjusted p-value used is within two Monte
Carlo standard errors of `alpha`, so that another seed could change the
class), `shared` (the set shares cells with other sets), `confounded`
(the set is confounded with differences among groups of units, see
[`describe_design()`](https://pabloamadorc.github.io/must/reference/describe_design.md)),
`min_units`, `status` (`"ok"`, `"not_applicable"` or `"not_estimable"`)
and `reason`. The attribute `alpha` keeps the significance level, for
[`report_interactions()`](https://pabloamadorc.github.io/must/reference/report_interactions.md).

## Details

Each effect is coded `"0"` if its adjusted p-value is at least `alpha`,
and `"+"` or `"-"` after the sign of its estimate otherwise. A
significant interaction is coded `"++"` (`"--"`) when the combined
effect is larger (smaller) than every single effect, or when the
absolute interaction is larger than every absolute single effect (legend
of Table 1 in Piggott et al. 2015), comparing point estimates on the
scale of the null model. Table 1 then reduces to:

|             |                       |       |
|-------------|-----------------------|-------|
| Interaction | Single effects        | Class |
| `0`         | any                   | AD    |
| `++`        | any                   | +S    |
| `--`        | any                   | -S    |
| `+`         | none negative         | +S    |
| `+`         | at least one negative | -A    |
| `-`         | none positive         | -S    |
| `-`         | at least one positive | +A    |

Under MU, BL, SAM and CA, +S and -S mean "more positive" or "more
negative" than the expectation of that null model, not than the sum of
the single effects. The BL, SAM and CA classes are computed from the
relative effects (`w - 1`) and the interaction of each model; the
descriptive `rescaled_class` of
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)
is not used.

**Three or more stressors.** The net deviation of a set of k stressors
is the combined effect minus the sum of the k single effects (AD, or the
same on the log scale for MU), the comparison that Table 1 classifies;
the same rules are applied with k single effects. This generalization is
a design decision of this package: Piggott et al. (2015) define the
classes for two stressors. The pure k-way interaction is not classified.

**Verdicts.** With an equivalence margin, the net interaction gets a
verdict from its two adjusted p-values, both compared with `alpha`:
`"deviation"` (deviation detected, not equivalent), `"negligible"` (a
deviation detected but within the margin), `"equivalent"` (no deviation
detected and within the margin) or `"inconclusive"` (neither). The class
`"AD"` covers both `"equivalent"` and `"inconclusive"`: only the first
supports additivity (see "Equivalence" in
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)).
Combinations with only zeros have no test of equivalence, so their
verdict is `"deviation"` or `"inconclusive"`. DO interactions get a
verdict too, although they have no class.

**Combinations with only zeros** (`status = "boundary"` in
[`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md))
are classified from the sign of their estimate (which may be infinite
under MU) and their p-value.

## Assumptions

- AD means that no deviation from the null model was detected. It is not
  evidence of additivity: with few replicates most interactions will be
  classed AD. With a `margin`, `verdict` tells equivalent from
  inconclusive.

- The class of a set with `confounded = TRUE` may reflect differences
  among groups of units (rooms, blocks) as much as the stressors.

- The classification is defined for null models that predict the
  combination from the single effects (AD, MU, BL, SAM and CA). The
  dominance model gets `status = "not_applicable"`, as do rows of
  `"not_applicable"` estimates (see `sensitivity` in
  [`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)).

## References

Piggott, J. J., Townsend, C. R., & Matthaei, C. D. (2015).
Reconceptualizing synergism and antagonism among multiple stressors.
*Ecology and Evolution*, 5(7), 1538–1547.
[doi:10.1002/ece3.1465](https://doi.org/10.1002/ece3.1465)

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
est <- estimate_interactions(fits, n_sim = 2000, seed = 1)
classify_interactions(est)
#>   response context time set_id order          stressors
#> 1    algae         <NA>    S01     2 salinity x warming
#>                        levels background null_model    scale estimand
#> 1 salinity=high | warming=yes                    AD response      net
#>   stressor_a code_a stressor_b code_b          single_codes code_combined
#> 1   salinity      -    warming      - salinity -; warming -             -
#>   code_interaction         pattern class
#> 1                + a - | b - | i +    -A
#>                                                            meaning verdict
#> 1 Combined effect less negative than the sum of the single effects    <NA>
#>   borderline shared confounded min_units status reason
#> 1      FALSE  FALSE      FALSE         5     ok       
```
