# devMSMs

*devMSMs* is an R package for implementing marginal structural models (MSMs) with longitudinal data to answer causal questions in developmental science. This package has an accompanying tutorial paper, *Investigating Causal Questions in Human Development using Marginal Structural Models: A Tutorial Introduction to the devMSMs Package in R*, that is forthcoming. We suggest reading this paper for a brief conceptual introdcution to MSMs prior to implementing this package. We have provided additional resources on this page for more in-depth material on this method and its application.  

Core features of this package include: 
- step-by-step user guidance for users new to the MSM technique and R programming
- flexible functions for implementing IPTW weighting and outcome modeling to answer substantive causal questions about dose and timing
  
## Overview
The package contains 6 core functions for conducting confounder adjustment and outcome modeling of longitudinal data with time-varying exposures.  

  
<img width="742" alt="Screen Shot 2023-11-08 at 4 00 24 PM" src="https://github.com/istallworthy/devMSMs/assets/31548151/cbab3b78-ffa8-4ffc-9b97-d082a9f145b5">

    
  [Table 1.pdf](https://github.com/istallworthy/devMSMs/files/13313050/Table.1.pdf)


## Installation
devMSMs can be installed in R Studio from Github using the *devtools* package:  
`library(devtools)`  
`install_github("istallworthy/devMSMs")`  
`library(devMSMs)`  

The helper functions can be installed from the accompanying *devMSMsHelpers* repo:   
`install_github("istallworthy/devMSMsHelpers`  
`library(devMSMsHelpers)`
  

## Recommended Workflow

<img width="651" alt="Screen Shot 2023-11-08 at 4 00 44 PM" src="https://github.com/istallworthy/devMSMs/assets/31548151/91216ef3-e4c3-4c8a-a522-f683105bfc66">

  
## Additional Resources
Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity score for a continuous treatment: Application to the efficacy of political advertisements. The Annals of Applied Statistics, 12(1), 156–177. https://doi.org/10.1214/17-AOAS1101  

Hirano, K., & Imbens, G. W. (2004). The Propensity Score with Continuous Treatments. In Applied Bayesian Modeling and Causal Inference from Incomplete-Data Perspectives (pp. 73–84). John Wiley & Sons, Ltd. https://doi.org/10.1002/0470090456.ch7  

Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal Structural Models and Causal Inference in Epidemiology. Epidemiology, 11(5), 550–560.  

Thoemmes, F., & Ong, A. D. (2016). A Primer on Inverse Probability of Treatment Weighting and Marginal Structural Models. https://doi.org/10.1177/2167696815621645  






