# devMSMs
<br>
Scientists who study humans are fundamentally interested in questions of causation, yet conceptual, methodological, and practical barriers have historically prevented their use of methods for causal inference developed in other fields. More specifically, scientists, clinicians, educators, and policymakers alike are often interested in *causal processes* involving questions about when (timing) and to what extent (dose) different factors influence human functioning and development, in order to inform our scientific understanding and improve people's lives.  
<br>
Marginal structural models (MSMs; Robins et al., 2000), orginating in epidemiology and public health, represent one under-utilized tool for improving causal inference with longitudinal observational data, given certain assumptions. In brief, MSMs leverage inverse-probability-of-treatment-weights (IPTW) and the potential outcomes framework. MSMs first focus on the problem of confounding, using IPTW to attenuate associations between measured confounders and an exposure (e.g., experience, characteristic, event --from biology to the broader environment) over time. A weighted model can then be fitted relating a time-varying exposure and a future outcome. Finally, the model-predicted effects of different exposure histories that vary in dose and timing can be evaluated and compared as counterfactuals, to reveal putative causal effects.    
 <br>
*devMSMs* is an R package accompanying our tutorial paper, *Investigating Causal Questions in Human Development using Marginal Structural Models: A Tutorial Introduction to the devMSMs Package in R* (*insert preprint link here*). Together, they offer conceptual and practice guidance for implementing MSMs with longitudinal data to answer causal questions about the dose and timing effects of a given exposure on a future outcome. 

Core features of *devMSMs* include:  
  
- flexible functions with built-in user guidance, drawing on established expertise and best practices for implementing longitudinal IPTW weighting and outcome modeling, to answer substantive causal questions about dose and timing
  
- functions that accept complete or imputed data to accommodate missingness often found in human studies  
   
- a novel recommended workflow, based on expertise from several disciplines, for using the *devMSMs* functions with longitudinal data (see *Workflows* vignettes)
  
- step-by-step user guidance for implementing the *deveMSMs* worflow in the form of vignettes applied to simulated data, geared toward users of all levels of R programming experience, along with a <a href="https://github.com/istallworthy/devMSMs/blob/main/examplePipelineRevised.Rmd">R markdown template file</a> 
  
- an accompanying suite of <a href="https://github.com/istallworthy/devMSMsHelpers">helper functions</a> to assist users in preparing and inspecting their data prior to the implementaiton of *devMSMs* (see the <a href="https://istallworthy.github.io/devMSMs/articles/Preliminary_Steps.html"> Preliminary steps vignette</a>)

- a conceptual introduction, an example empirical application, and additional resources in the accompanying tutorial paper

<br>

## Overview
The package contains 6 core functions for conducting longitudinal confounder adjustment and outcome modeling of longitudinal data with time-varying exposures.
<br>  

<img width="742" alt="Screen Shot 2023-11-08 at 4 00 24 PM" src="https://github.com/istallworthy/devMSMs/assets/31548151/cbab3b78-ffa8-4ffc-9b97-d082a9f145b5">  
<br>  
    
Below is a summary of the terms used in the *devMSMs* vignettes and functions. More details and examples can be found in the accompanying manuscript. 
<br>
  
<img width="742" alt="Screen Shot 2023-11-08 at 4 00 24 PM" src="https://github.com/istallworthy/devMSMs/blob/main/man/figures/term%20summary%20table.png">  

<br> 
 
## Installation
*devMSMs* can be installed in R Studio from Github using the *devtools* package:   
<br>  
`library(devtools)`  
`install_github("istallworthy/devMSMs")`  
`library(devMSMs)`  

The helper functions can be installed from the accompanying *devMSMsHelpers* repo:   
<br>   
`install_github("istallworthy/devMSMsHelpers")`  
`library(devMSMsHelpers)`
  
<br> 
  
## Recommended Workflow
We propose a recommended workflow for using *devMSMs* to answer causal questions with longituinal data. Please see the *Workflows* vignettes and accompanying manuscript for more details.  
  
<img width="651" alt="Screen Shot 2023-11-08 at 4 00 44 PM" src="https://github.com/istallworthy/devMSMs/assets/31548151/91216ef3-e4c3-4c8a-a522-f683105bfc66">

<br> 

## Additional Resources
Austin, P. C. (2011). An Introduction to Propensity Score Methods for Reducing the Effects of Confounding in Observational Studies. Multivariate Behavioral Research, 46(3), 399–424. https://doi.org/10.1080/00273171.2011.568786  

Blackwell, M. (2013). A Framework for Dynamic Causal Inference in Political Science. American Journal of Political Science, 57(2), 504–520. https://doi.org/10.1111/j.1540-5907.2012.00626.x  
  
Cole, S. R., & Hernán, M. A. (2008). Constructing Inverse Probability Weights for Marginal Structural Models. American Journal of Epidemiology, 168(6), 656–664. https://doi.org/10.1093/aje/kwn164  
  
Eronen, M. I. (2020). Causal discovery and the problem of psychological interventions. New Ideas in Psychology, 59, 100785. https://doi.org/10.1016/j.newideapsych.2020.100785  
   
Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity score for a continuous treatment: Application to the efficacy of political advertisements. The Annals of Applied Statistics, 12(1), 156–177. https://doi.org/10.1214/17-AOAS1101  
  
Haber, N. A., Wood, M. E., Wieten, S., & Breskin, A. (2022). DAG With Omitted Objects Displayed (DAGWOOD): A framework for revealing causal assumptions in DAGs. Annals of Epidemiology, 68, 64–71. https://doi.org/10.1016/j.annepidem.2022.01.001  
  
Hirano, K., & Imbens, G. W. (2004). The Propensity Score with Continuous Treatments. In Applied Bayesian Modeling and Causal Inference from Incomplete-Data Perspectives (pp. 73–84). John Wiley & Sons, Ltd. https://doi.org/10.1002/0470090456.ch7  

Kainz, K., Greifer, N., Givens, A., Swietek, K., Lombardi, B. M., Zietz, S., & Kohn, J. L. (2017). Improving Causal Inference: Recommendations for Covariate Selection and Balance in Propensity Score Methods. Journal of the Society for Social Work and Research, 8(2), 279–303. https://doi.org/10.1086/691464  

Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal Structural Models and Causal Inference in Epidemiology. Epidemiology, 11(5), 550–560.  

Thoemmes, F., & Ong, A. D. (2016). A Primer on Inverse Probability of Treatment Weighting and Marginal Structural Models. https://doi.org/10.1177/2167696815621645  






