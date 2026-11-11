# the print method shows the criteria, the models and the reasons

    Code
      print(s)
    Output
      <must_null_models> after the guidelines of Schäfer and Piggott (2018), Fig. 3
      Criteria:
        The stressors act similarly (same mode of action)                       no
        One stressor acts specifically on physiology                            unknown
        The effect is mortality                                                 no
        The stressor-effect relationship is known                               no
        The sensitivity distribution is known                                   unknown
        Correlation of the sensitivities                                        unknown
        The effect of one stressor is close to the effect limit (at least 90%)  no
        The sum of the effects compared with the effect limit                   below
        A stressor alone increases the response                                 no
        Stressors alone with opposite effects                                   no
      Effects of the stressors alone (given): salinity = 0.3 (decrease),
        warming = 0.5 (decrease); type: double negative.
      Null models:
        AD (simple addition): possible. Step 4: the sensitivities are highly
          negatively correlated and the sum of the effects is below the
          effect limit; step 7: none of the earlier steps applies.
        MU (multiplicative): excluded. The multiplicative model needs
          mortality as the effect (step 5) or the sum of the effects above
          the effect limit (step 6), but the effect is not mortality and the
          sum of the effects does not exceed the effect limit.
        DO (dominance): possible. Step 3: the sensitivities are highly
          positively correlated, or the effect of one stressor is close to
          the effect limit.
        SAM (stress addition): excluded. The stress addition model needs one
          stressor that acts specifically on physiology, mortality as the
          effect, stressors that do not increase the response and a known
          stressor-effect relationship or sensitivity distribution (step 2),
          but the effect is not mortality.
        CA (concentration addition): excluded. Concentration addition needs
          stressors that act similarly, with effects in the same direction,
          and a known stressor-effect relationship (step 1), but the
          stressors do not act similarly and their stressor-effect
          relationship is not known.
      Unknown criteria that would change the choice:
        Correlation of the sensitivities (`correlation`): positive -> DO;
          negative -> AD; weak -> AD
      For estimate_interactions(): null_models = c("AD", "DO")

