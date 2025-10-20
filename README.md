# Parametric Modal Regression with Error-Contaminated Covariates

This repository contains all code and data to reproduce the results from the manuscript "Parametric Modal Regression with Error Contaminated Covariates"
by Yanfei He, Jianhong Shi, and Weixing Song.

##  Overview

This study develops parametric modal regression methods for datasets with error-contaminated covariates. The repository includes:
- Simulation studies validating the proposed methods
- Two real-data case studies (AD data and Dietary data)
- Complete reproducible research files

##  Environment

- **R version**: 4.4.0
- **Processor**: Intel(R) Core(TM) i5-8500 CPU @ 3.00GHz
- **RAM**: 8.00 GB
- **Matrix products**: Default

##  Required R Packages
```r
MASS, pracma, matrixcalc, MultiRNG, gofgamma, dcov, VGAM, latex2exp

##  Repository Structure

├── ReadME.txt: ReadME file 
├── simulation_time.xlsx               # Simulation time recorded for all the simulation studies
├── case_study/ 
│   ├── AD_data/
│   │   ├── estimate.R                 # Main analysis (Tables 6-8, Figures 7-9)
│   │   ├── test_Cvm_KS.R              # P-value calculations for CvM and KS tests
│   │   ├── test_MCCM1.R               # P-value calculations for MCCM1 test
│   │   ├── gammaTest.R                # Subroutine on Gamma test from package gofgamma
│   │   ├── resampling.R               # Subroutine for bootstrap critical values 
│   │   ├── p-value.xlsx               # Estimated p-values for gamma test
│   │   └── BJadni.xls                 # AD dataset
│   └── Dietary_data/
│       ├── Diet1_estimate.Hatx.R      # Right panels of Figures 9-11
│       ├── Diet1_estimate.R           # Main analysis (Tables 4-6, left panels of Figures 9-11)
│       ├── Diet1_test_Cvm_KS.R        # P-value calculations for CvM and KS tests
│       ├── Diet1_test_MCCM1.R         # P-value calculations for MCCM1 test
│       ├── Diet2_estimate.R           # Table 7 and Figure 12
│       ├── Diet2_test_Cvm_KS.R        # P-value calculations for CvM and KS tests
│       ├── Diet2_test_MCCM1.R         # P-value calculations for MCCM1 test
│       ├── gammaTest.R                # Subroutine on Gamma test from package gofgamma
│       ├── resampling.R               # Subroutine on bootstrap critical values
│       └── wishreg.xls                # Dietary dataset
└── simulation/
    ├── simulation_1/
    │   ├── Table_1.R                  # Table 1 results
    │   ├── Table_2.R                  # Table 2 results
    │   ├── Table_3.R                  # Table 3 results
    │   ├── Figure2.dependent.R        # Dependent case figures
    │   └── Figure2.independent.R      # Independent case figures
    ├── simulation_2/
    │   ├── Figure_3.R                 # Figure 3
    │   ├── Figure_4.R                 # Figure 4
    │   ├── Figure_5.R                 # Figure 5
    │   ├── Figure_6.R                 # Figure 6
    │   ├── Figure_7_Model_2.1.R       # Figure 7 Model 2.1
    │   ├── Figure_7_Model_2.2.R       # Figure 7 Model 2.2
    │   ├── Figure_7_Model_2.3.R       # Figure 7 Model 2.3
    │   ├── Figure_7_8_print.R         # Print Figures 7-8
    │   ├── Figure_8_CvM_KS.R          # Figure 8 CvM & KS tests
    │   ├── Figure_8_DC.R              # Figure 8 Distance Correlation
    │   ├── Figure_8_MCCM1.R           # Figure 8 MCCM1 tests
    │   └── Figure_8_MCCM2.R           # Figure 8 MCCM2 tests
    ├── simulation_B1/
    │   ├── Table_B.13.R               # Table B.13 results
    │   ├── Table_B.14.R               # Table B.14 results
    │   ├── Table_B.15.R               # Table B.15 results
    │   ├── Figure_B.13.dependent.R    # Figure B.11 dependent
    │   └── Figure_B.13.independent.R  # Figure B.11 independent
    └── simulation_B2/
        ├── Figure_B.14.R              # Figure B.14
        ├── Figure_B.15.R              # Figure B.15
        ├── Figure_B.16.R              # Figure B.16
        ├── Figure_B.17.R              # Figure B.17
        ├── Figure_B.18_B.1.R          # Figure B.18 Model B.1
        ├── Figure_B.18_B.2.R          # Figure B.18 Model B.2
        ├── Figure_B.18_B.3.R          # Figure B.18 Model B.3
        ├── Figure_B.18_B.4.R          # Figure B.18 Model B.4
        ├── Figure_18_print.R          # Print Figure 18
        ├── Figure_18_CvM_KS.R         # Figure 18 CvM & KS tests
        ├── Figure_18_DC.R             # Figure 18 Distance Correlation
        ├── Figure_18_MCCM2.R          # Figure 18 MCCM2 tests
        ├── Model_B.3_CvM_KS.R         # Model B.3 CvM & KS tests
        ├── Model_B.3_DC.R             # Model B.3 Distance Correlation
        └── Model_B.3_MCCM2.R          # Model B.3 MCCM2 tests


##  Output

Each subfolder "results" within the "simulation" folder contains all the intermediate results ready for check.

##  Support

For questions or comments about this code, please contact:

Weixing Song - weixing@ksu.edu
