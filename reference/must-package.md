# must: multiple-stressor interactions under alternative null models

Whether a combined effect is called additive, synergistic or
antagonistic depends on the null model that defines "no interaction" and
on the scale on which effects are combined. must makes these choices
explicit and separates three questions:

## Details

1.  **What the design can identify.**
    [`describe_design()`](https://pabloamadorc.github.io/must/reference/describe_design.md)
    reports missing treatment combinations, cells shared by several
    contrasts, the replication unit and repeated measures, and which
    interactions are identifiable, before any model is fitted.

2.  **Which model fits the data.**
    [`fit_responses()`](https://pabloamadorc.github.io/must/reference/fit_responses.md)
    fits a generalised linear mixed model to each response and chooses
    its error distribution by convergence, residual checks and AIC,
    independently of the null models.

3.  **What "non-additive" means.**
    [`estimate_interactions()`](https://pabloamadorc.github.io/must/reference/estimate_interactions.md)
    tests each combination of stressors against one or more null models
    with the same criteria, from the predictions of that model or of one
    fitted with another package: simple addition (the default),
    multiplicative, dominance, Bliss independence, stress addition and
    concentration addition. It distinguishes the net deviation of a
    combination of three or more stressors from its pure highest-order
    interaction. By default the tests are corrected for small samples;
    `correction = "none"` gives the asymptotic procedure of Madge
    Pimentel et al. (2025).
    [`classify_interactions()`](https://pabloamadorc.github.io/must/reference/classify_interactions.md)
    assigns the directional classes of Piggott et al. (2015), and
    [`summary()`](https://rdrr.io/r/base/summary.html) and
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) give an
    overview.

[`suggest_null_models()`](https://pabloamadorc.github.io/must/reference/suggest_null_models.md)
applies the guidelines of Schäfer & Piggott (2018) to choose the null
models, and
[`report_interactions()`](https://pabloamadorc.github.io/must/reference/report_interactions.md)
writes the methods of the analysis and a table of every decision behind
it.

## References

Madge Pimentel, I., Albini, D., Beermann, A. J., Leese, F., Macaulay, S.
J., Matthaei, C. D., Orr, J. A., Piggott, J. J., & Schäfer, R. B.
(2025). Hypothesis-driven research on multiple stressors: An analytical
framework for stressor interactions. *Ecology and Evolution*, 15,
e71959. [doi:10.1002/ece3.71959](https://doi.org/10.1002/ece3.71959)

Orr, J. A., et al. (2020). Towards a unified study of multiple
stressors: Divisions and common goals across research disciplines.
*Proceedings of the Royal Society B*, 287, 20200421.
[doi:10.1098/rspb.2020.0421](https://doi.org/10.1098/rspb.2020.0421)

Piggott, J. J., Townsend, C. R., & Matthaei, C. D. (2015).
Reconceptualizing synergism and antagonism among multiple stressors.
*Ecology and Evolution*, 5(7), 1538–1547.
[doi:10.1002/ece3.1465](https://doi.org/10.1002/ece3.1465)

Schäfer, R. B., & Piggott, J. J. (2018). Advancing understanding and
prediction in multiple stressor research through a mechanistic basis for
null models. *Global Change Biology*, 24(5), 1817–1826.
[doi:10.1111/gcb.14073](https://doi.org/10.1111/gcb.14073)

## See also

Useful links:

- <https://pabloamadorc.github.io/must/>

- <https://github.com/pabloamadorc/must>

- Report bugs at <https://github.com/pabloamadorc/must/issues>

## Author

**Maintainer**: Pablo Amador Crespo <pamadorc@gmail.com>
([ORCID](https://orcid.org/0000-0002-3793-8794))
