# 
This repository contains the code to 1) generate synthetic data in the correct format 2) perform the variable using either the Broken Adaptive Ridge (BAR) method or the Minimum Information Criterion (MIC) method

This repository contains three main files. The first file named *genRecurandTermData.R* contains the code to generate data required for the simulation study of either Scenario 1 or Scenario 2. It will save the synthetic data in .csv format. After pre-specifying how many datasets to generate and running the code, the second file named *BAR.R* contains the code to perform BAR method on the newly created datasets. This file requires C++ code to successfully run the R program, hence the user needs to install *Rcpp* and *RcppArmadillo*. 

Similarly, the third file named *MIC.R* contains the code to perform the MIC method. 
