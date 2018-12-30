# Sensitivity-analysis-on-heterogeneity-effect-using-data-simulation

this project is a research on treatment effect estimates of highway safety countermeasure.

The purpose of this project is to:
1. Explore the sensitivity of matching methods to unobserved heterogeneity.
2. Using data simulation based on observational data to ocntrol the effect of unobserved covaraites.
3. Comparing the performance between Bayesian propensity score matching and typical propensity score matching methods.
 
The code inventory includes: 
1. sensitivity_analysis_data_simulation.R: the main file of simulation analysis.
2. jags_simulation.R: function of data simulation, and Bayesian modeling using JAGS.
3. matching_asmd: compute the absolute standardized mean difference between treatment and control group before and after matching.
4. visualization.R: visualize the distribution of simulation anlaysis results.
