# kSZ 4pt analysis using SPT x SPIRE
## Paper:
* Raghunathan S., SPT-3G and SPTpol Collaboration 2024; arXiv:[2403.02337](https://arxiv.org/abs/2403.02337).

## Overview:
* This work uses a quadratic estimator, similar to CMB-lensing, to reconstruct the velocity-induced correlations on large-scales of the small-scale kSZ signal.
* The method is based on [Smith & Ferroro 2016](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.021301)(arXiv:[1607.01769](https://arxiv.org/abs/1607.01769)).
 * Given that the reionisation process is expected to be patchy, it should give rise to non-Gaussian kSZ signals.
 * The trispectrum measurement will help us to differentiate between the contributions from the reionisation and the post-reionisation (due to haloes in the local Universe) kSZ signals which are hard to de distinguished in the power spectrum space.

## Results:
* [read_files_and_make_plots.ipynb](https://github.com/sriniraghunathan/kSZ_4pt_SPT_SPIRE/blob/main/read_files_and_make_plots.ipynb): Reproduce plots in the paper.
* [K_map_data_sims_meanfield.npy](https://github.com/sriniraghunathan/kSZ_4pt_SPT_SPIRE/blob/main/results/K_map_data_sims_meanfield.npy): Reconstructed $\hat{K}(\hat{n})$ map from simulations.
  * Note that the data map is currently not released.
  * Used for Figure 1 of the paper.
* [CL_kk_data_sims_amberksz_nzero.npy](https://github.com/sriniraghunathan/kSZ_4pt_SPT_SPIRE/blob/main/results/CL_kk_data_sims_amberksz_nzero.npy): $C_{L}^{KK}$ measurements from simulations and data.
  * Used for Figure 2 of the paper.
* ...: Posterior distributions and priors.
  * Used for Figure 3 of the paper.
* [sims_data_for_likelihood_with_and_without_systematic_for_plotting.npy](https://github.com/sriniraghunathan/kSZ_4pt_SPT_SPIRE/blob/main/results/sims_data_for_likelihood_with_and_without_systematic_for_plotting.npy): $C_{L}^{KK}$ measurements from simulations for diffenent analysis choices.
  * Used for (Appendix) Figure 4 of the paper.
* ...: Posterior distributions and priors.
  * Used for (Appendix) Figure 5 of the paper.
