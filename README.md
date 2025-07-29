# Parametric-Modal-Regression-with-Error-Contaminated-Covariates
###  by Yanfei He, Jianhong Shi and Weixing Song  
  
This repository contains all the R-codes and data for the simulation studies and real data applicaitons in our manuscript submitted to the Biometrical Journal. 

The folder "Simulation Codes" contains two subfolders

  -- Folder "Diagnostic 1"

     M1.1: Tables 1–3 and Figures 2–6 in Simulation 1, Section 6.
     misM2.1: Table 4 in Model 2.1.
     misM2.2: Table 4 in Model 2.2.
     misM2.3: Table 4 and 5 in Model 2.3.

  -- Folder "Diagnostics2"
  
     B1: Tables B.13–B.15 and Figures B.11–B.15 in Simulation B1 in Appendix B.
     B2.misM.B.1: Table B.16 in Model B.1.
     B2.misM.B.2: Table B.16 in Model B.2.
     B2.misM.B.3: Table B.16 and Table B.17 in Model B.3.
     B2.misM.B.4: Table B.16 and Table B.18 in Model B.4.


  -- Folder "Dietary Data"
  
     Data File: 
     
     wishreg.xls: The Dietary Data dataset (Excel format).
       
     Analysis Scripts
        Diet1.R, Diet1.Hatx.R:  Tables 6-8 and figure 7-9.
        
        Diet.test1.ks.R & Diet1.test.mccl.R:   
        Test whether the data follows a Gamma distribution assumption, including the 
        CvM, KS test from Diet1.test.ks.R and the proposed MCCM1 tests in 
        Diet1.test.mccl.R.

        Diet2.R, Diet2.Hatx.R: Tables 9 and Figure 10
   
        Diet2.test.ks.R & Diet2.test.mccl.R:
        
         Test whether the data follows a Gamma distribution assumption. The 
         CvM, KS tests are included in Diet2.test.ks.R and the proposed MCCM1 test 
         is included in Diet2.test.mccl.R.

  -- Folder "AD data": 

     Data File:
     
     BJadni.xls: The Alzheimer’s Disease dataset (Excel format).

     Analysis Scripts
     
     AD.estimate.R:  Tables 10–12. 
   
     AD.test.R & AD.test.mccl.R:
     
     Test whether the data follows a Gamma distribution assumption. The CvM, KS are 
     included in AD.test.R and the proposed MCCM1 test is included in AD.test.mccl.R.


