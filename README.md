# Sensitivity-analysis-on-heterogeneity-effect-using-data-simulation

this project is a research on treatment effect estimates of highway safety countermeasure.

The purpose of this project is to:
1. Explore the sensitivity of matching methods to unobserved heterogeneity.
2. Using data simulation based on observational data to ocntrol the effect of unobserved covaraites.
3. Comparing the performance between Bayesian propensity score matching and typical propensity score matching methods.

The case study:
1. In this case, the safety effect (crash reduction effect) of roadway rumble strips on PA two-lane highways is to be estimated.
2. The data includes 151 miles (326 segments) roads in the treatment group, and 6225.9 miles (13,091 segments) road in the control group.
3. 10-year (2003 to 2012) crash records on each roads are included.
4. Outcome variable: total crash frequency. Treatment: rumble strips (binary). Confounders: before-period crash history, AADT, paved width, speed limit, access density, roadside hazard rating, degree of curvature, etc.

The reserach study:
1. A no-treatment analysis to validate the effectiveness of propensity score matching methods. (not included, refer to my dissertation)
2. A simulation-based analysis: 
  1) Use the real treatment and confounder data to simulate outcome, 
  2) Add simulated unobserved variable
  3) Compare the matching estimate to simulated effect and test the sensitivity of the methods to unobserved effect.
  
The packages:
1. Propensity score matching: MatchIT in R
2. Bayesian modeling and simulation: JAGS and Rjags in R

The code inventory includes: 
1. sensitivity_analysis_data_simulation.R: the main file of simulation analysis.
2. jags_simulation.R: function of data simulation, and Bayesian modeling using JAGS.
3. matching_asmd: compute the absolute standardized mean difference between treatment and control group before and after matching.
4. visualization.R: visualize the distribution of simulation anlaysis results.
