# PFTC6-AT-cleaning

This repo houses scripts for downloading and cleaning Dataset II: Assimilation-Temperature response data for the following publication:

Halbritter et al., 2025. Plant trait, carbon flux, reflectance and climate data from global change experiments and gradients in Norway. Scientific Data.

The script "code/main.R" is the main script for downloading and cleaning assimilation-temperature response data. Data is automatically downloaded from the OSF repo associated with this paper ([link](https://osf.io/fcbw4/)). Other scripts included are helper functions which are called by the main script as needed. The script outputs cleaned data in a single .csv file, "data/AT_clean.csv".

Last modified, 9 April 2025, Josef Garen.
