# Robustness of the null models across mechanisms

Simulated with null_model_robustness.R on 2026-11-09, with must 0.1.0.9000 and R version 4.3.3 (2024-02-29): 100 replicates per scenario and 1000 Krinsky-Robb draws per test (57 minutes).

Rates in %, over the combinations tested: 2 x 2 factorial and 3 x 3 dose designs, survival (binomial) and counts (quasi-Poisson), moderate and weak effects; p < 0.05 without adjustment in each combination. CA is tested in the dose design only. Consensus: an interaction only if every null model flags it with the same sign.

## False interactions (no interaction, by the mechanism that generated the data)

|  | AD | MU | DO | BL | SAM | CA | consensus |
|---|---|---|---|---|---|---|---|
| AD | 4.8 | 10.6 | 20.2 | 10.6 | 28.2 | 8.1 | 0.3 |
| MU | 6.7 | 4.9 | 12.1 | 4.9 | 32.5 | 9.7 | 0.2 |
| DO | 19.7 | 14.8 | 3.2 | 14.8 | 44.5 | 25.0 | 2.3 |
| SAM | 60.2 | 76.4 | 85.8 | 76.4 | 4.0 | 34.7 | 0.2 |
| CA | 9.7 | 13.9 | 37.9 | 13.9 | 23.4 | 2.6 | 0.2 |

## Power against interactions that no mechanism explains

|  | AD | MU | DO | BL | SAM | CA | consensus |
|---|---|---|---|---|---|---|---|
| synergy | 79.0 | 97.4 | 99.4 | 97.4 | 10.8 | 51.2 | 8.3 |
| antagonism | 37.6 | 36.9 | 21.6 | 36.9 | 64.8 | 51.8 | 21.5 |

## Direction of the false interactions

Share of the false interactions with a negative estimate (synergy, for these responses), in %.

|  | AD | MU | DO | BL | SAM | CA | consensus |
|---|---|---|---|---|---|---|---|
| AD | 55.2 | 90.8 | 99.6 | 90.8 | 0.0 | 24.4 | 0.0 |
| MU | 20.9 | 49.1 | 97.4 | 49.1 | 0.0 | 1.4 | 0.0 |
| DO | 3.1 | 5.1 | 37.0 | 5.1 | 0.0 | 0.0 | 0.0 |
| SAM | 99.7 | 100.0 | 100.0 | 100.0 | 14.8 | 99.8 | 100.0 |
| CA | 76.7 | 98.8 | 99.9 | 98.8 | 0.0 | 31.6 | 0.0 |

## Rule

Among the null models defined in two-level factorials and with an average power of at least 80% of the highest, the one with the lowest average rate of false interactions (the five mechanisms with the same weight); within one percentage point, the more powerful one.

Result: **AD**.

## Summary by null model

| Null model | False interactions, mean | False interactions, maximum | Power, mean |
|---|---|---|---|
| AD | 20.2 | 60.2 | 58.3 |
| MU | 24.1 | 76.4 | 67.1 |
| DO | 31.9 | 85.8 | 60.5 |
| BL | 24.1 | 76.4 | 67.1 |
| SAM | 26.5 | 44.5 | 37.8 |
| CA | 16.0 | 34.7 | 51.5 |
| consensus | 0.6 | 2.3 | 14.9 |

## By design, response and size of the effects

Mean rate of false interactions (five mechanisms) / mean power (synergy and antagonism), in %.

| Scenario | AD | MU | DO | BL | SAM | CA | consensus |
|---|---|---|---|---|---|---|---|
| dose_3x3 counts moderate | 19.6 / 50.8 | 30.7 / 71.5 | 40.3 / 63.6 | 30.7 / 71.5 | 28.8 / 39.2 | 15.0 / 51.2 | 0.4 / 12.7 |
| dose_3x3 counts weak | 12.8 / 36.2 | 17.3 / 54.5 | 21.6 / 53.4 | 17.3 / 54.5 | 15.7 / 26.5 | 4.4 / 22.8 | 0.4 / 3.6 |
| dose_3x3 survival moderate | 34.2 / 93.1 | 38.1 / 95.1 | 53.1 / 85.4 | 38.1 / 95.1 | 57.4 / 54.8 | 32.4 / 88.2 | 0.6 / 39.1 |
| dose_3x3 survival weak | 26.2 / 71.6 | 28.5 / 71.4 | 37.9 / 63.3 | 28.5 / 71.4 | 27.5 / 49.5 | 12.3 / 43.7 | 0.6 / 18.7 |
| factorial_2x2 counts moderate | 14.4 / 38.5 | 20.8 / 57.0 | 25.4 / 52.5 | 20.8 / 57.0 | 13.6 / 22.6 | - | 1.0 / 3.4 |
| factorial_2x2 counts weak | 6.4 / 29.5 | 7.4 / 43.5 | 6.0 / 48.5 | 7.4 / 43.5 | 17.1 / 17.5 | - | 1.2 / 0.0 |
| factorial_2x2 survival moderate | 30.4 / 87.0 | 30.4 / 85.0 | 45.4 / 64.5 | 30.4 / 85.0 | 38.9 / 51.5 | - | 0.6 / 22.2 |
| factorial_2x2 survival weak | 17.8 / 59.5 | 19.4 / 59.0 | 25.2 / 53.0 | 19.4 / 59.0 | 13.3 / 40.6 | - | 0.3 / 19.4 |

