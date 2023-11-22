# devMSMs
<br>
Scientists who study humans are fundamentally interested in questions of causation, yet conceptual, methodological, and practical barriers have historically prevented their use of methods for causal inference developed in other fields. More specifically, scientists, clinicians, educators, and policymakers alike are often interested in *causal processes*, involving questions about when (timing) and to what extent (dose) different factors influence human functioning and development, in order to inform our scientific understanding and improve people's lives.  
<br>   
   
Marginal structural models (MSMs; Robins et al., 2000), orginating in epidemiology and public health, represent one under-utilized tool for improving causal inference with longitudinal observational data. In brief, MSMs leverage inverse-probability-of-treatment-weights (IPTW) and the potential outcomes framework to attenuate associations between confounders and an exposure (e.g., experience, characteristic, evevent --from biology to the broader environment) over time, and uncover causal relations between the time-varying exposure and future outcome, given certain assumptions. 
 <br>   
   
*devMSMs* is an R package accompanying our tutorial paper, *Investigating Causal Questions in Human Development using Marginal Structural Models: A Tutorial Introduction to the devMSMs Package in R* (*insert preprint link here*), for implementing MSMs with longitudinal data to answer causal questions about the dose and timing effects of a given exposure on a future outcome. 

Core features of this package include: 
  
- flexible functions with built-in guidance, drawing on established expertise and best practices for implementing longitudinal IPTW weighting and outcome modeling to answer substantive causal questions about dose and timing
  
- accommodation of data in the form of either a complete dataframe or multiple imputation 
   
- a recommended workflow for using the *devMSMs* functions with longitudinal data 
  
- step-by-step user guidance to *deveMSMs* worflow in the form of vignettes and a <a href="https://github.com/istallworthy/devMSMs/blob/main/examplePipelineRevised.Rmd">R markdown template file</a> for users new to the MSM technique and R programming

- an accompanying suite of <a href="https://github.com/istallworthy/devMSMsHelpers">helper functions</a> to assist users in preparing and inspecting their data prior to use of *devMSMs*

- conceptual introduction and example empirical application in accompanying tutorial paper

<br>

## Overview
The package contains 6 core functions for conducting longitudinal confounder adjustment and outcome modeling of longitudinal data with time-varying exposures.  
  
<br>  

<img width="742" alt="Screen Shot 2023-11-08 at 4 00 24 PM" src="https://github.com/istallworthy/devMSMs/assets/31548151/cbab3b78-ffa8-4ffc-9b97-d082a9f145b5">  
<br> 
<br> 
<br>   
  

<img width="658" alt="Screen Shot 2023-11-09 at 4 54 22 PM" src="https://github.com/istallworthy/devMSMs/assets/31548151/4b62730e-3975-4002-a06b-205c93d04224">    

<br> 
 
## Installation
*devMSMs* can be installed in R Studio from Github using the *devtools* package:   
<br> 
`library(devtools)`  
`install_github("istallworthy/devMSMs")`  
`library(devMSMs)`  

The helper functions can be installed from the accompanying *devMSMsHelpers* repo:   
`install_github("istallworthy/devMSMsHelpers`  
`library(devMSMsHelpers)`
  
<br> 
  
## Recommended Workflow

<img width="651" alt="Screen Shot 2023-11-08 at 4 00 44 PM" src="https://github.com/istallworthy/devMSMs/assets/31548151/91216ef3-e4c3-4c8a-a522-f683105bfc66">

<br> 

## Additional Resources
Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity score for a continuous treatment: Application to the efficacy of political advertisements. The Annals of Applied Statistics, 12(1), 156–177. https://doi.org/10.1214/17-AOAS1101  

Hirano, K., & Imbens, G. W. (2004). The Propensity Score with Continuous Treatments. In Applied Bayesian Modeling and Causal Inference from Incomplete-Data Perspectives (pp. 73–84). John Wiley & Sons, Ltd. https://doi.org/10.1002/0470090456.ch7  

Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal Structural Models and Causal Inference in Epidemiology. Epidemiology, 11(5), 550–560.  

Thoemmes, F., & Ong, A. D. (2016). A Primer on Inverse Probability of Treatment Weighting and Marginal Structural Models. https://doi.org/10.1177/2167696815621645  






