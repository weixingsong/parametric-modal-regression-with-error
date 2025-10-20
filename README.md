Supplementary information / reproducible research files for the manuscript 
Title: "Parametric Modal Regression with Error Contaminated Covariates"

Authors: Yanfei He, Jianhong Shi, Weixing Song
Code was written by Yanfei He, Jianhong Shi and Weixing Song
In case of questions or comments please contact  weixing@ksu.edu.

The code was written/evaluated in R with the following software versions:
R version 4.4.0
Processor: Inter(R) Core(TM) i5-8500 CPU @ 3.00GHz 3.00GHZ 
Computer RAM: 8.00GB

Matrix products: default

attached  packages:
MASS  pracma  matrixcalc  MultiRNG  gofgamma  dcov  VGAM   latex2exp           

This folder contains the following  files that can be used to reproduce all table and figures of the manuscript.
It contains two subfolders containing the following files:

./case_study/:
This folder contains the following  two files：
     ./AD data/
     This subfolder folder contains the following code(.R) and AD data(.xls)
         estimate.R
         An R script that contains the code of the analysis reported in the paper (section 6).  One can obtain the result 
         of Table6-8 and figure 7-9

         test_Cvm_KS.R and test_MCCM1.R
         One can obtain the pvalue.
             
         gammaTest.R and resampling.R
         Subroutine required by test_Cvm_KS.R.R

         BJadni.xls
         An XLS sheet containing AD data.
        

     ./Dietary data/
     This subfolder folder contains the following code(.R) and Dietary data(.xls)
         Diet1_estimate.Hatx.R
         One can obtain the right-panel of Figure 7-9 

         Diet1_estimate.R
         An R script that contains the code of the analysis reported in the paper (section 6).  One can obtain the result 
         of Table 6-8 and the left-panel of Figure 7-9.

         Diet1_test_Cvm_KS.R and Diet1_test_MCCM1.R
          One can obtain the pvalue.

         Diet2_estimate.R
         One can obtain the result of Table 9 and Figure 10.

         Diet2_test_Cvm_KS.R and Diet2_test_MCCM1.R
         One can obtain the pvalue.
             
         gammaTest.R and resampling.R
         Subroutine required by Diet2_test_Cvm_KS.R

         wishreg.xls
         An XLS sheet containing Dietary data.
    
    
./simulation/
This folder contains the following four files, each corresponding to a different simulation in the manuscript.
   ./simulation 1/
    This subfolder folder contains the following code(.R) and a file.
        Table 1.R
        All results of Table 1 can be generated using the following parameter combinations: n = {100, 200}, B = {100, 200}, 
        su2 = {0.25, 1}.  In the code, Btnumsum and NBtnumsum correspond to the proposed method and the Naive method 
        in the manuscript, respectively.

       Table 2.R,Table 3.R, figure2.dependent.R and figure2.independent.R
       Running the above will generate the corresponding tables and figures.  In the code, Btnumsum(Sdnumsum),  TBtnumsum(TSdnumsum) 
       and NBtnumsum(NSdnumsum) correspond to the proposed method(MCCL_1) , MCCL_2 and the Naive method in the manuscript, respectively.

   
   ./simulation 2/
   Running the code in the above file will generate the corresponding tables by varying the value of \beta_4 and the sample size n.
    
   ./simulation B1/
    This subfolder folder contains the following code(.R) and a file. 
       Table B.13.R
       All results of Table B.13 can be generated using the following parameter combinations: n = {100, 200}, B = {100, 200}, 
       su2 = {0.25, 1}.  In the code, Btnumsum and NBtnumsum correspond to the proposed method and the Naive method 
       in the manuscript, respectively.

      Table B.14.R,Table B.15.R, figure B.11.dependent.R and figure B.11.independent.R
      Running the above will generate the corresponding tables.
   
   ./simulation B2/
    Running the in the above file will generate the corresponding tables by varying the sample size n.
    
      results.R
      Each folder's result file contains the results generated from within that folder. Specifically, these results correspond to 
      the tables and figures presented in the manuscript's online supplementary information.
