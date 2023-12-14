
<!-- README.md is generated from README.Rmd. Please edit that file -->

# devMSMs: Implementing Marginal <img src="https://github.com/istallworthy/devMSMs/assets/31548151/0a1c69d7-1984-4ebe-8d16-0948d670a8fb" align="right" width="100"/> Structural Models with Longitudinal Data

<!-- badges: start -->
<!-- badges: end -->

<br>


Those who study and work with humans are fundamentally interested in questions of causation. More specifically, scientists, clinicians, educators, and policymakers alike are often interested in *causal processes* involving questions about when (timing) and to what extent (dose) different factors influence human functioning and development, in order to inform our scientific understanding and improve people's lives. However, for many, conceptual, methodological, and practical barriers have prevented the use of methods for causal inference developed in other fields.
<br>

The goal of this *devMSMs* package and accompanying tutorial paper, *Investigating Causal Questions in Human Development Using Marginal Structural Models: A Tutorial Introduction to the devMSMs Package in R* (*insert preprint link here*), is to provide a set of tools for implementing marginal structural models (**MSMs**; Robins et al., 2000).   
  
MSMs orginated in epidemiology and public health and represent one under-utilized tool for improving causal inference with longitudinal observational data, given certain assumptions. In brief, MSMs leverage inverse-probability-of-treatment-weights (IPTW) and the potential outcomes framework. MSMs first focus on the problem of confounding, using IPTW to attenuate associations between measured confounders and an exposure (e.g., experience, characteristic, event --from biology to the broader environment) over time. A weighted model can then be fitted relating a time-varying exposure and a future outcome. Finally, the model-predicted effects of different exposure histories that vary in dose and timing can be evaluated and compared as counterfactuals, to reveal putative causal effects.
<br>

We employ the term *exposure* (sometimes referred to as "treatment" in other literatures) to encompass a variety of environmental factors, individual characteristics, or experiences that constitute the putative causal events within a causal model. Exposures may be distal or proximal, reflecting a developing child’s experience within different environments at many levels (Bronfenbrenner & Ceci, 1994), ranging from the family (e.g., parenting), home (e.g., economic strain), school (e.g., teacher quality), neighborhood (e.g., diversity), to the greater politico-cultural-economic context (e.g., inequality). Exposures could also reflect factors internal to the child, including neurodevelopmental (e.g., risk markers), physiological (e.g., stress), and behavioral (e.g., anxiety) patterns to which the child’s development is exposed.  
<br>


## Core Features 
Core features of *devMSMs* include:  
  
- flexible functions with built-in user guidance, drawing on established expertise and best practices for implementing longitudinal IPTW weighting and outcome modeling, to answer substantive causal questions about dose and timing
  
- functions that accept complete or imputed data to accommodate missingness often found in human studies  
   
- a novel recommended workflow, based on expertise from several disciplines, for using the *devMSMs* functions with longitudinal data (see *Workflows* vignettes)

- an accompanying simulated longitudinal dataset, based on the real-world, Family Life Project (FLP) study of human development, for getting to know the package functions

- an accompanying suite of <a href="https://github.com/istallworthy/devMSMsHelpers">helper functions</a> to assist users in preparing and inspecting their data prior to the implementation of *devMSMs* (see the <a href="https://istallworthy.github.io/devMSMs/articles/Preliminary_Steps.html"> Preliminary Steps vignette</a>)  
   
- executable, step-by-step user guidance for implementing the *deveMSMs* worflow and preliminary steps in the form of vignettes geared toward users of all levels of R programming experience, along with a <a href="https://github.com/istallworthy/devMSMs/blob/main/examplePipelineRevised.Rmd">R markdown template file</a> 

- a brief conceptual introduction, example empirical application, and additional resources in the accompanying tutorial paper

<br>

## Overview
The package contains 6 core functions for implementing the two phases of the MSM process: longitudinal confounder adjustment and outcome modeling of longitudinal data with time-varying exposures.
<br>  
<img width="820" alt="overview" src="https://github.com/istallworthy/devMSMs/assets/31548151/db0efeb1-c1e5-4c71-8062-01a5f34d56ab">
<br>  
<br>
    
Below is a summary of the terms used in the *devMSMs* vignettes and functions. More details and examples can be found in the accompanying manuscript. 
<br>
  
<img width="513" alt="term summary table" src="https://github.com/istallworthy/devMSMs/assets/31548151/e8839e88-20dc-4e55-9f88-b936876be75b">

<br> 
 
## Installation
*devMSMs* can be installed in R Studio from Github using the *devtools* package: 
<br>  
`library(devtools)`  
`install_github("istallworthy/devMSMs")`  
`library(devMSMs)`  

<br> 

The helper functions can be installed from the accompanying *devMSMsHelpers* repo: 
<br>   
`install_github("istallworthy/devMSMsHelpers")`  
`library(devMSMsHelpers)`
  
<br> 
  
## Recommended Workflow
We propose a recommended workflow for using *devMSMs* to answer causal questions with longituinal data. We suggest using the vignettes in the order they appear in the Articles tab. After reading the accompanying manuscript, We recommend first reviewing the <a href="https://istallworthy.github.io/devMSMs/articles/Terminology.html">Terminology</a> and <a href="https://istallworthy.github.io/devMSMs/articles/Data_Requirements.html">Data Requirements</a> vignettes as you begin preparing your data. We then recommend downloading the <a href="https://github.com/istallworthy/devMSMs/blob/main/examplePipelineRevised.Rmd">R markdown template file</a> which contains all the code desribed in the <a href="https://istallworthy.github.io/devMSMs/articles/Specify_Core_Inputs.html">Specify Core Inputs</a>, <a href="https://istallworthy.github.io/devMSMs/articles/Preliminary_Steps.html">Preliminary Steps</a>, and *Workflows* vignettes (for binary (TBA) or <a href="https://istallworthy.github.io/devMSMs/articles/Workflow_Continuous_Exposure.html">continuous</a> exposures) for implementing the steps below.
  
<img width="800" alt="Screen Shot 2023-11-08 at 4 00 44 PM" src="https://github.com/istallworthy/devMSMs/assets/31548151/91216ef3-e4c3-4c8a-a522-f683105bfc66">

<br> 

## Citation & Bug Reports
Please cite your use *devMSMs* using the following citation: 
<br> 
Stallworthy I, Greifer N, DeJoseph M, Padrutt E, Berry D (2023). *devMSMs*: Implementing Marginal Structural Models with Longitudinal Data. R package version 0.0.0.9000, https://istallworthy.github.io/devMSMs/.  

<br>   

Please report any bugs at the following link: https://github.com/istallworthy/devMSMs/issues 

<br> 

## Additional Resources
Austin, P. C. (2011). An Introduction to Propensity Score Methods for Reducing the Effects of Confounding in Observational Studies. Multivariate Behavioral Research, 46(3), 399–424. https://doi.org/10.1080/00273171.2011.568786  

Blackwell, M. (2013). A Framework for Dynamic Causal Inference in Political Science. American Journal of Political Science, 57(2), 504–520. https://doi.org/10.1111/j.1540-5907.2012.00626.x  
  
Cole, S. R., & Hernán, M. A. (2008). Constructing Inverse Probability Weights for Marginal Structural Models. American Journal of Epidemiology, 168(6), 656–664. https://doi.org/10.1093/aje/kwn164  
  
Eronen, M. I. (2020). Causal discovery and the problem of psychological interventions. New Ideas in Psychology, 59, 100785. https://doi.org/10.1016/j.newideapsych.2020.100785  
   
Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity score for a continuous treatment: Application to the efficacy of political advertisements. The Annals of Applied Statistics, 12(1), 156–177. https://doi.org/10.1214/17-AOAS1101  
  
Greifer N (2023). cobalt: Covariate Balance Tables and Plots. R package version 4.5.2, https://github.com/ngreifer/cobalt, https://ngreifer.github.io/cobalt/  

Greifer N (2023). WeightIt: Weighting for Covariate Balance in Observational Studies. https://ngreifer.github.io/WeightIt/, https://github.com/ngreifer/WeightIt  
   
Haber, N. A., Wood, M. E., Wieten, S., & Breskin, A. (2022). DAG With Omitted Objects Displayed (DAGWOOD): A framework for revealing causal assumptions in DAGs. Annals of Epidemiology, 68, 64–71. https://doi.org/10.1016/j.annepidem.2022.01.001  
  
Hirano, K., & Imbens, G. W. (2004). The Propensity Score with Continuous Treatments. In Applied Bayesian Modeling and Causal Inference from Incomplete-Data Perspectives (pp. 73–84). John Wiley & Sons, Ltd. https://doi.org/10.1002/0470090456.ch7  

Kainz, K., Greifer, N., Givens, A., Swietek, K., Lombardi, B. M., Zietz, S., & Kohn, J. L. (2017). Improving Causal Inference: Recommendations for Covariate Selection and Balance in Propensity Score Methods. Journal of the Society for Social Work and Research, 8(2), 279–303. https://doi.org/10.1086/691464  

Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal Structural Models and Causal Inference in Epidemiology. Epidemiology, 11(5), 550–560.  

Thoemmes, F., & Ong, A. D. (2016). A Primer on Inverse Probability of Treatment Weighting and Marginal Structural Models. https://doi.org/10.1177/2167696815621645  


