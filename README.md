# MVGC2

This is Version 2 (early development) of the MVGC Multivariate Granger Causality MATLAB software suite.

The MVGC2 project will enhance and extend the existing open-source MVGC Multivariate
Granger Causality MATLAB toolbox (see github.com/SacklerCentre/MVGC1). This project
will update the current MVGC: (i) to take advantage of new state-of-the-art state-space and
spectral methods developed in-house; (ii) to implement the MVGC functionality in Python
(MVGC2-P); and (iii) to integrate the MVGC functionality with the popular neuroimaging
software suites EEGLAB (MVGC2-EEGLAB) and MNE (MVGC2-MNE). The enhanced algorithms, as well
as improving on efficiency and accuracy, also address some well-known problems associated
with standard Granger causality inference from neuroimaging data. The Python implementation
and EEGLAB/MNE integration facilitate users to deploy the MVGC functionality within their
standard toolchains.

Currently, this repository contains the MVGC2 standalone MATLAB toolbox. To use the toolbox,
first inspect the config.m script in the root directory, and edit to taste.

Make sure to run the startup.m script (e.g., by starting MATLAB in th MVGC2 root directory).

WARNING: Do NOT add the entire MVGC2 directory tree to your MATLAB path - this will cause
things to break! MVGC2 paths are set appropriately by the startup script.
