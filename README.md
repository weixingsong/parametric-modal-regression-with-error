# Parametric-Modal-Regression-with-Error-Contaminated-Covariates
#  by Yanfei He, Jianhong Shi and Weixing Song  
  
This repository contains all the R-codes and data for the simulation studies and real data applicaitons 
in our manuscript submitted to the Biometrical Journal. 


--1. File (diagnostics1) Description
       File M1.1: Corresponds to Simulation 1 in Section 6 of the paper, generating results for Tables 1–3 and Figures 2–6.
       File misM2.1: Corresponds to Model 2.1 in Table 4.
       File misM2.2: Corresponds to Model 2.2 in Table 4.
       File misM2.3: Corresponds to Model 2.3 in Table 4 and Table 5.

--2. File (diagnostics2) Description
      File B1: Corresponds to Simulation B1 in Appendix B, generating results for Tables B.13–B.15 and Figures B.11–B.15.
      File B2.misM.B.1: Corresponds to Model B.1 in Table B.16.
      File B2.misM.B.2: Corresponds to Model B.2 in Table B.16.
      File B2.misM.B.3: Corresponds to Model B.3 in Table B.16 and Table B.17.
      File B2.misM.B.4: Corresponds to Model B.4 in Table B.16 and Table B.18.


--3.File (Dietary Data) Description: Dietary Data Analysis:

   Data File
   wishreg.xls: The Dietary Data dataset (Excel format).

   Analysis Scripts
   1.Diet.R, 1.Diet.Hatx.R, :
   This script tests is used to generate the results for Tables 6-8 and figure 7-9.
   1.Diet.test.ks.R & 1.Diet.test.mccl.R:
   These scripts  whether the data follows a Gamma distribution assumption,  including the CvM, KS in 1.Diet.test.ks.R and the proposed MCCM1 tests in 1.Diet.test.mccl.R.

   2.Diet.R, 2.Diet.Hatx.R, :
   This script tests is used to generate the results for Tables 9 and figure 10.
   2.Diet.test.ks.R & 2.Diet.test.mccl.R:
  These scripts  whether the data follows a Gamma distribution assumption,  including the CvM, KS in 2.Diet.test.ks.R and the proposed MCCM1 tests in 2.Diet.test.mccl.R.

--4.File (AD data) Description: Alzheimer’s Disease (AD) Data Analysis:

   Data File
   BJadni.xls: The Alzheimer’s Disease dataset (Excel format).

   Analysis Scripts
   AD.estimate.R:
   This script tests is used to generate the results for Tables 10–12, involving parameter estimation and hypothesis testing.
   AD.test.R & AD.test.mccl.R:
   These scripts  whether the data follows a Gamma distribution assumption,  including the CvM, KS in AD.test.R and the proposed MCCM1 tests in AD.test.mccl.R.


