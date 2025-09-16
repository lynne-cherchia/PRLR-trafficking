# PRLR-trafficking

Simulation code for ODE model of PRLR trafficking. This repository will accompany a manuscript describing the use of sFCS and simulations of the trafficking model to characterize PRLR organization.

## JAK/STAT Model

The folder "baseModelOutputs" contains .mat files holding predictions of the base model created by [Mortlock et al.](https://link.springer.com/article/10.1007/s12195-020-00647-8) 

The folder "baseModelOutputs_t24h" contains a .mat file holding predictions of the base model created by Mortlock et al. for BCL over a 24 h time course.

The folder "error_ranked_parameters" contains .mat files holding the error-ranked parameters and initial values found by [Simoni et al.](https://academic.oup.com/ib/article/14/2/37/6562848) for the base JAK/STAT model. In the revised trafficking model, we use the best fitted (lowest error) parameter and initial values sets:

    lowest_error_params = error_ranked_params(1,:);
    lowest_error_initvals = error_ranked_initvals(1,:);

The folder "traffickingParams_initValues" contains .mat files holding parameters and initial values for simulating the trafficking model under the varied initial conditions (all species RJ internal, all RJ at cell surface, etc.)

"downstreamRateTuning" simulates and plots results for varying signal transduction parameters (k6, k10, and k21) in the trafficking model.

"meshgridSims_250rates_plotHeatmaps" simulates and plots heatmaps for varying PRLR-complex trafficking parameters (krecfree, kintfree, krecbound, and kintbound) in the trafficking model; "meshgridSims_250rates_plotTimeCourses" also simulates varying the trafficking parameters and plots 0-6 h time courses of the receptor localization ratio. These files are written for use on a high-performance computing (HPC) cluster as the large meshgrid requires large amounts of RAM.


## Technical details and requirements

All MATLAB scripts were written in MATLAB R2024a on macOS Sonoma (version 14.4) with an Apple M3 Pro chip and 18 GB RAM. "downstreamRateTuning" was executed locally; meshgrid simulations were executed on USC's CARC HPC cluster.


Simulating the trafficking model can take hours-days depending on the number of parameter values tested as well as computer specification and number of CPU cores utilized.
