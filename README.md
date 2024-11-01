# Master's Thesis

This repo contains the implented code used for my master's thesis.

The repo is comprised of three folders, WGS, WGS_particle and WGS_particle_reaction.

WGS contains the Manifest and Project files used in Julia, which can be activated by using the activate("WGS") function.

WGS_particle is the diffusion model described in the thesis. The model is implemented in WGSParticle_lab.jl, and all functions used for the implementation are in functions_WGSParticle.jl. The plots folder contains generated plots, used in the thesis.

The folder WGS_particle_reaction, contains the reaction model described in the thesis, and the trained surrogate models. The reaction model is in WGSParticleReaction.jl, and all functions are in functions_WGSParticleReaction.jl. The surrogate models can be found in the SurrogatesWGS.jl, but are used in WGSParticleReaction.jl. The folders parity_plots, appx_figures and figures_result contain generated plots, which are used in thesis. In the folder csv_files, all csv files can be found. This includes sampled parameter, training and testin data, and calulated quality metrics. For the sample csv files, init depicts the first approach, beyond is the second, and random is the third. The file all_data_m2.csv is the training data used, and test_data_m2.csv is the test data. The file all_errors.csv, is the indices of all failed simulations, and cci_t0_beyond.csv are the intital concentrations used for the second approach. 