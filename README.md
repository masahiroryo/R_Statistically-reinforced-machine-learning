# R_Statistically reinforced machine learning
  
## **Ryo M. and M.C. Rillig (2017) Statistically-reinforced machine learning for nonlinear patterns and variable interactions. Ecosphere.**


**Corresponding Author: Masahiro Ryo (masahiroryo@gmail.com)**  
  Do not hesitate to contact me if you get in any troubles in using these scripts :)

## **How to use?**
1. Download all the following scripts and csv file into a folder
2. Run SML_application.R, which is the main body of the script
3. If you want to try with your own data, adjust "read.csv(...)"
     
     
    
### *data_example.csv*   
 This is an example dataset generated with "Data_generator.R"   
 The sample size was set to 300.   
   
### *Data_generator.R*   
 This is used to generate data_example.csv   
 It can change sample size, predictors, associations, the ratio of missing values and so on.   
 (This is not called in the main script, "SML_application.R")   
    
### *SML_application.R*   
 This is the main script. To run this script, "data_example.csv" and "Functions.R" must be contained in the same folder.   
   
### *Functions.R*  
 This is called in SML_application.R to define two functions for permutation-based random forest modeling.  
 The script is originally available from Hapfelmeier and Ulm (2013), and this version allows the user to perform parallel computing on Windows OS. 
