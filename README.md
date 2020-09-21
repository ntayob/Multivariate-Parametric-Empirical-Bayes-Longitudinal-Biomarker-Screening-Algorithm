# Multivariate-Parametric-Empirical-Bayes-Longitudinal-Biomarker-Screening-Algorithm

We have provided an example set of code to implement the multivariate parametric empirical Bayes (mPEB) proposed in "A multivariate parametric empirical Bayes screening approach for early detection of hepatocellular carcinoma using multiple longitudinal biomarkers". We propose a parametric empirical Bayes approach for hepatocellular carcinoma (HCC) screening using multiple longitudinal biomarkers. The goal of HCC screening is to detect the cancer at an early disease stage when there are multiple treatment options available since HCC cases that are diagnosed at a later stage, once patients start exhibiting symptoms, have poor survival. We have developed a general screening algorithm for multiple biomarkers that incorporates longitudinal screening history by defining patient specific thresholds. We propose a minimal model for the joint biomarkers to maintain a robust and flexible algorithm.

The example datasets are generated based on Simulation Scenario A in the manuscript. We provide a training dataset (data_training.csv) to estimate parameters of the mPEB algorithm and a validation dataset (data_validation.csv) to evaluate the performance of the algorithm.

To implement the analysis, the following steps are required

1. Run code in "Analyze Training Data.R". This will produce the following .csv files that will be used in the Matlab code:
a. mu.csv
b. sigma_theta.csv
c. sigma.csv
d. mu_star.csv
e. sigma_theta_star.csv
f. sigma_star.csv

2. Run code in "three_biomarker_sims_cluster.m" to implement the mPEB algorithm. The functions "probability_cases.m" and "probability_controls.m" are the objective function and contrainst used by fmincon to find local optimal solutions. The output from Matlab (data_validation_output.txt) is the number of positive biomarkers at each screening occasion for varying population level false-positive rates (f0).

3. Run code in "Analyze Validation Data.R" to import results from Matlab and estimate the patient-level sensitivity (for multiple time-periods of clinical interest) when screening-level specificity is fixed.
