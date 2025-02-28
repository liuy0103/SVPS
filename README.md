# [AISTATS 2025] Selecting the Number of Communities for Weighted Degree-Corrected Stochastic Block Models

This repository contains the R code for numerical experiments in the paper **Selecting the Number of Communities for Weighted Degree-Corrected Stochastic Block Models**, as published in **AISTATS 2025**.

### Synthetic Networks

For Simulation 1 (Poisson), run `sim_poisson.R`; for Simulation 2 (Binomial), run `sim_binom.R`; for Simulation 3 (Negative Binomial), run `sim_nbinom.R`; for additional experiments in Appendix D, run `additional_poisson.R`. Please adjust the DCSBM parameters $\rho$ and $r$ accordingly to reproduce each setup in the paper.

### Les Miserables Network

The data is originally compiled in Knuth (1993) and downloaded from https://websites.umich.edu/~mejn/netdata/. For the estimated number of communities in the network by different methods, run `lesmis_plot.R`; for the plots of node clustering by SCORE and RSC, run `lesmis_plot.R`.

## References
D. E. Knuth, *The Stanford GraphBase: A Platform for Combinatorial Computing*. Addison-Wesley, Reading, MA (1993).