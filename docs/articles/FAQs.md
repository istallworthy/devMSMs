# devMSMs FAQs

## Getting Started

### 1. Who was *devMSMs* designed for?

We designed the package to make MSMs accessible to researchers at all
stages of careers in developmental science, psychology, and/or
neuroscience. However, the package could be useful to professionals from
a wide variety of domains who are interested in causal processes.
*devMSMs* is particularly useful for those interested in causal
questions about the accumulation and timing of causal effects on an
outcome.

### 2. What background knowledge do I need to use *devMSMs?*

The package was designed to make implementing MSMs accessible to
researchers who do not have training in causal inference methods.
However, we strongly recommend users read the
[manuscript](https://academic.oup.com/chidev/article/97/2/331/8516421?utm_source=authortollfreelink&utm_campaign=chidev&utm_medium=email&guestAccessKey=)
accompanying this package (which also provides recommendations for
further reading) to gain a working understanding of the underlying
concepts and assumptions of MSMs prior to using the package. It is also
helpful for users to have working knowledge of R programming (basic
coding, data cleaning, familiarity with the
[Rmarkdown](https://rmarkdown.rstudio.com/articles_intro.html) format).

### 3. What main terms do I need to know to use *devMSMs?*

We strongly recommend users read the
[manuscript](https://academic.oup.com/chidev/article/97/2/331/8516421?utm_source=authortollfreelink&utm_campaign=chidev&utm_medium=email&guestAccessKey=)
accompanying this package to gain a working understanding of the
underlying concepts and terms (see Table 1 in manuscript) relevant to
MSMs prior to using the package.

- *Exposure* (“treatment” in other literatures) = hypothesized cause of
  your outcome (can encompass various environmental factors, individual
  characteristics, or experiences that constitute putative causal
  events; may be distal or proximal, reflecting a developing child’s
  experience within different environments at many levels, from the
  family, home, school, and neighborhood to the greater
  politico-cultural-economic context; could reflect factors internal to
  the child, including neurodevelopmental, physiological, and
  psychological patterns; should be time-varying)

- *Outcome* = hypothesized effect of your cause (measured at a single
  timepoint at the last exposure timepoint or ideally following it)

- *Confounder* = common cause of exposure and outcome (can be
  time-invariant or time-varying)

Please also see the [Terminology
vignette](https://istallworthy.github.io/devMSMs/articles/Terminology.html)
for more information.

### 4. What do I need to get started with *devMSMs?*

To get started, you need a specific causal question/set of questions
that you want to test and longitudinal data measuring your time-varying
exposure (hypothesized cause), outcome (hypothesized effect), and
time-invariant and time-varying confounders (causes of both exposure and
outcome).

### 5. How do I formulate my causal `question(s)` for *devMSMs*?

As detailed in the
[manuscript](https://academic.oup.com/chidev/article/97/2/331/8516421?utm_source=authortollfreelink&utm_campaign=chidev&utm_medium=email&guestAccessKey=)
accompanying this package, MSMs operate within the potential outcomes
framework (Imbens & Rubin, 2015) which asks, *what would have happened
(i.e., outcome) under circumstances (i.e., exposure) other than what
actually occurred?* We probe causal effects of exposure by comparing
predicted outcomes under different exposure circumstances. In this vein,
using MSMs requires that you express your causal question about the
effect of exposure on outcome in terms of *comparisons of outcomes that
follow from different exposure histories (i.e., sequences).*

It helps to make each causal question as specific as possible, in the
form of a specific hypothesis, in terms of the effects of different
*histories of exposure* on an outcome. For example, I could consider the
question: *Are there causal effects of economic strain experien*ce*d
over infancy, toddlerhood, and early childhood on children’s behavior
problems?* From this, for example, I could formulate a more specific
hypothesis that *chronic exposure to economic strain* over this period
would result in *worse behavior problems* compared to *no exposure to
economic strain.* I can now think about the exposure history comparisons
that could test this hypothesis. Here I could compare outcomes following
*high exposure at all timepoints (“h”-“h”-“h”)* to those following *low
exposure at all timepoints (“l”-“l”-“l”)*. The package requires you to
identify your causal question ahead of time and identify specific
exposure history `comparison(s)` to test that question. Please refer to
the
[manuscript](https://academic.oup.com/chidev/article/97/2/331/8516421?utm_source=authortollfreelink&utm_campaign=chidev&utm_medium=email&guestAccessKey=)
accompanying this package for more extensive consideration of how to
implement this research question for *devMSMs,* and our customizable
[project planning
template](https://docs.google.com/document/d/1Lj-uWObhoorDg87nYKyqaG4M3NoV0uLYeueR1_yJeKE/edit?tab=t.0)
for more information.

### 6. What are all the things I should consider while planning a project using *devMSMs?*

Please refer to our customizable [project planning
template](https://docs.google.com/document/d/1Lj-uWObhoorDg87nYKyqaG4M3NoV0uLYeueR1_yJeKE/edit?tab=t.0)
for starting a project using existing data or planning new data
collection efforts.

## Data

### 7. What kind of data do I need to use *devMSMs*?

*devMSMs* was designed for use with longitudinal datasets that have at
least two (though ideally more) timepoints of data per person and
include a measure of exposure (measured at multiple timepoints), at
least one confounder (ideally multiple time invariant and time-varying),
and one outcome measure (measured at one timepoint). We note that it can
sometimes be challenging to create IPTWs when there is very little
variation in exposure (especially with binary exposures). Please see the
[Data Requirements & Preparation
vignette](https://istallworthy.github.io/devMSMs/articles/Data_Requirements.html)
for more information. We strongly recommend using large longitudinal
datasets such as the ABCD, HBCD, ECLS-K, FLP, Future of Families, and
RAPID datasets.

### 8. Can I use *devMSMs* with accelerated longitudinal data (i.e., data are not collected at the exact same timepoints across individuals)?

The package requires that you identify a focal set of timepoints that
apply to everyone in your sample. These timepoints (or a subset of them)
need to be the timepoints at which your exposure, confounders, and
outcome are measured for everyone. For accelerated longitudinal designs,
you may need to collapse across timepoints to accommodate data from all
individuals. For example, if some people have three timepoints of data
collected at 1, 3, and 5 months and others at 2, 4, and 6 months, you
might need to identify 1-2, 3-4, and 5-6 months as your three timepoints
that accommodate everyone.

### 9. How does *devMSMs* handle missing data?

Missing data are ubiquitous in developmental science. The package
requires complete data and is designed to accept multiple imputed
datasets (either output from the
[*mice*](https://cran.r-project.org/web/packages/mice/index.html)
package in “*mids”* format) or in list form. *devMSMs* conducts all
steps on each imputed dataset separately in accordance with best
practices (Granger et al., 2019; Leyrat et al., 2021), pooling estimates
only at the final stage to fit just one outcome model used to compare
histories. Please see the [helper
functions](https://github.com/istallworthy/devMSMsHelpers) and [Data
Requirements & Preparation
vignette](https://istallworthy.github.io/devMSMs/articles/Data_Requirements.html)
for more information.

### 10. How is the *devMSMs* workflow implemented with multiple imputed datasets?

Given existing work demonstrating its superiority, *devMSMs* implements
the ‘within’ approach for imputed data (Granger et al., 2019; Leyrat et
al., 2021), conducting all steps on each imputed dataset before pooling
estimates using Rubin’s rules to create one final set of average
predictions and history comparisons.

### 11. How should I format my longitudinal data for *devMSMs?*

The package requires complete data in wide format (i.e., one row per ID,
with exposures, outcome, and confounders as columns) with no missing
values, in the form of multiple imputed datasets or a single complete
dataset. Please see the [Data Requirements & Preparation
vignette](https://istallworthy.github.io/devMSMs/articles/Data_Requirements.html)
for more information.

### 12. Is there a minimum sample size requirement for using *devMSMs*?

Power analysis techniques for MSMs are still under development (e.g.,
Shook-Sa & Hudgens, 2022). Beyond our recommendation to use large
longitudinal datasets (such as ABCD, HBCD, ECLS-K, FLP, Future of
Families, and RAPID), we cannot offer specific sample size requirements.
The package also does provide ways for researchers to examine how well
represented each of their user-specified exposure histories are in their
data (see Workflow vignettes) prior to conducting MSMs.

## Inputs

Exposure

### 13. What kind of exposure can I use with *devMSMs*?

We use the term “exposure” broadly to refer to any measurable construct
with a binary or continuous distribution that is hypothesized to
causally impact some measurable outcome. Exposure can be internal or
external to the person. These should ideally be measured at multiple
timepoints. We note that it can sometimes be challenging to create IPTWs
when there is very little variation in exposure (especially with binary
exposures). Please see the [Specify Core Inputs
vignette](https://istallworthy.github.io/devMSMs/articles/Specify_Core_Inputs.html)
for more information.

### 14. What if I have multiple different exposures?

The package requires you to pick a single exposure, measured at multiple
timepoints, with respect to a single outcome, measured at a final
timepoint. If you have causal questions about multiple exposures and
outcomes, we suggest fitting separate MSMs for each exposure-outcome
pair.

### 15. How do I specify exposure histories of high and low exposure when my exposure is continuously distributed?

To specify exposure histories (i.e., sequences of high vs low levels of
exposure) with continuously distributed exposures, the user is required
to specify a priori exemplar values of what could be considered high
versus low levels of the continuous variable (e.g., 75th and 25th
quantile values). The package also does provide ways for researchers to
examine how well represented each of their user-specified exposure
histories are in their data (see Workflow vignettes). In the
accompanying
[manuscript](https://academic.oup.com/chidev/article/97/2/331/8516421?utm_source=authortollfreelink&utm_campaign=chidev&utm_medium=email&guestAccessKey=)
and Workflow vignettes, we provide guidance on how to decide these
values with an emphasis on a priori decision-making and transparency.
This practice is similar to probing simple slopes of an interaction. We
also note that the full variation of continuous exposures is preserved
throughout the process of creating and assessing IPTW weights and
fitting the weighted outcome model and exemplar high and low values are
only used in the final step of comparing exposure histories.

Outcome

### 16. What kind of outcome can I use with *devMSMs?*

We use the term “outcome” broadly to refer to any measurable construct
that is hypothesized to be caused by exposure. This should ideally be
measured after the exposures. Given the flexibility of the outcome model
(i.e., can accept any family/link functions accepted by the
[`glm()`](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/glm)
function from the *stats* package), there are no specific distributional
requirements for the outcome variable. Please see the [Specify Core
Inputs
vignette](https://istallworthy.github.io/devMSMs/articles/Specify_Core_Inputs.html)
for more information.

### 17. What if my outcome is measured at multiple timepoints?

Right now, *devMSMs* can only support an outcome model and history
comparisons with

respect to outcome at a single timepoint (ideally at or following the
final exposure timepoint). However, outcome measured at timepoints prior
to that final timepoint can be included by the user as a time-varying
confounder.

### 18. What if I have multiple different outcomes?

The package requires you to pick a single exposure, measured at multiple
timepoints, with respect to a single outcome, measured at a final
timepoint. If you have causal questions about multiple different
exposures and outcomes, we suggest fitting separate MSMs for each
exposure-outcome pair.

Confounders

### 19. How do I identify confounders of my exposure and outcome?

The process of identifying confounders should draw on existing theory
and

domain-specific literature. We also strongly suggest constructing a
directed acyclic graph

(DAG), or a diagram of directed arrows representing the possible causal
relations in your project. Please refer to these [additional DAG
resources](https://docs.google.com/presentation/d/1-43rAAZu2nFzZ1xuVadT1tPhbtn2vlYR/edit?slide=id.g34560c88db0_0_192#slide=id.g34560c88db0_0_192).
It is important to include all possible confounders of exposure and
outcome. Possible mediators and colliders for a given exposure timepoint
should *not* be balanced. However, it is okay if a confounder at a later
exposure timepoint is also a mediator of a prior exposure timepoint.
*devMSMs* has mechanisms to ensure confounders are not included in the
IPTW process in ways that they could act as mediators or colliders (see
subsequent Q&A).

### 20. What kinds of confounders can I include with *devMSMs*?

One of the biggest assumptions of this approach is that all possible
confounders (i.e., common causes) of exposure and outcome are measured
and included in the IPTW process. We suggest including as many time
invariant and time-varying confounders (ideally measured at multiple
timepoints), which can be any kind of measurable construct. When
possible, we suggest including multiple different measures of the same
construct. It is important to include all possible confounders of
exposure and outcome, even if they could also function as colliders or
mediators at some timepoints. *devMSMs* has mechanisms to ensure any
confounders are not included in the IPTW process in ways that they could
act as mediators or colliders (see subsequent Q&A). Please see the
[Specify Core Inputs
vignette](https://istallworthy.github.io/devMSMs/articles/Specify_Core_Inputs.html)
for more information.

### 21. How can I be sure my confounders (common *causes* of exposure and outcome) are not colliders (common *effects* of exposure and outcome) when they’re used for IPTW?

It can be challenging to disentangle confounders from colliders solely
based on the literature. Sometimes a given variable can function as a
confounder for exposure at some timepoints and a collider of exposure at
other timepoints, making it all the more confusing (see the accompanying
[manuscript](https://academic.oup.com/chidev/article/97/2/331/8516421?utm_source=authortollfreelink&utm_campaign=chidev&utm_medium=email&guestAccessKey=)).
Conditioning or balancing on a collider is problematic as it opens up a
backdoor path between exposure and outcome leading to bias. To avoid
this, the package has default settings to ensure that when creating IPTW
weights at each exposure timepoint, no user-specified confounders are
included in such a way that they could function as a collider. As a
default, the
[`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)
function (used to create and assess IPTWs) only includes time-varying
confounders at lagged timepoints relative to each exposure timepoint,
ensuring that colliders (which would have to be contemporaneous with a
given exposure timepoint or in the future) are not included. This
setting can be overridden by the user in cases when they are certain a
variable measured contemporaneously with the exposure functions as a
confounder and not a collider. Please see the [Customize Weights
Formulas
vignette](https://istallworthy.github.io/devMSMs/articles/Customize_Balancing_Formulas.html)
for more information.

### 22. How can I be sure my confounders (common cause of exposure and outcome) are not mediators (an intermediary variable along the causal path through which exposure affects outcome) when they’re used for IPTW?

It can be challenging to disentangle confounders from mediators solely
based on the literature. Often, a given variable functions as a
confounder for the exposure at some timepoints and a mediator of the
exposure at other timepoints (see the accompanying
[manuscript](https://academic.oup.com/chidev/article/97/2/331/8516421?utm_source=authortollfreelink&utm_campaign=chidev&utm_medium=email&guestAccessKey=)).
Conditioning or balancing on a mediator is problematic as it partially
blocks the effects of exposure at that timepoint on the outcome. To
avoid this, the package has default settings to ensure that when
creating IPTW weights at each exposure timepoint, no user-specified
confounders are included in such a way that they could function as a
mediator. As a default, the
[`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)
function (used to create and assess IPTWs) only includes time-varying
confounders at lagged timepoints relative to each exposure timepoint,
ensuring that confounders that function as mediators (which would have
to be contemporaneous with a given exposure timepoint or in the future)
are not included. This setting can be overridden by the user in cases
when they are certain a variable measured contemporaneously to the
exposure functions as a confounder rather than a mediator.

## Inverse-Probability-of-Treatment-Weighting (IPTW)

Fitting IPTWs

### 23. How does each person end up with a single IPTW weight based on their confounding at each exposure timepoint?

For each person, IPTW weights are created at each exposure timepoint
which are then multiplied together to create a single IPTW weight per
person that reflects the inverse of their probability of receiving their
sequence of exposures given all confounders (Cole & Hernan, 2008).
Please see the `weightitMSM()`
[documentation](https://search.r-project.org/CRAN/refmans/WeightIt/html/weightitMSM.html)
for more information and the
[manuscript](https://academic.oup.com/chidev/article/97/2/331/8516421?utm_source=authortollfreelink&utm_campaign=chidev&utm_medium=email&guestAccessKey=)
accompanying this package.

### 24. How are IPTW weights stabilized and how do I implement additional stabilization?

For continuous exposures, weights are automatically stabilized within
exposure timepoint drawing on the `weightitMSM()` function. With binary
exposures, the user needs to specify they wish to stabilize (see Binary
Exposure Workflow). Otherwise, simply specifying *stabilize = TRUE* to
the
[`createWeights()`](https://istallworthy.github.io/devMSMs/reference/createWeights.md)
function stabilizes the weights across exposure timepoints.

### 25. What other customized user specifications for fitting IPTWs are available using *devMSMs*?

Please see the *WeightIt* package
[documentation](https://ngreifer.github.io/WeightIt/) for more
information about the additional specifications available when fitting
IPTW weights via *devMSMs.* The
[`createWeights()`](https://istallworthy.github.io/devMSMs/reference/createWeights.md)
function in *devMSMs* can take any input accepted by the
[`weightitmsm()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.html)
function in *WeightIt*.

Assessing IPTWs

### 26. How do I know which weighting method to use when fitting IPTWs with *devMSMs*?

The best IPTW weighting method is the one that most successfully
attenuates all confounding at all exposure timepoints in your data. The
package draws on the [*WeightIt*](https://ngreifer.github.io/WeightIt/)
package which accommodates a whole (and growing) variety of weighting
methods (e.g., covariate balancing propensity score, Bayesian additive
regression trees, GLM, generalized boosted modeling) see their
documentation for more detail). *devMSMs* allows you to fit weights
using all of those methods, assess how well they attenuate confounding,
and select the one that works best for your data. Please see the
Workflow vignettes for more detail.

### 27. How will I know if IPTWs fitted through *devMSMs* successfully attenuate all measured confounding for all exposure timepoints in my data?

The package comes with functions for empirically assessing how well
IPTWs balance or attenuate each confounder at each timepoint, and allow
you to select the best-performing weights out of several weighting
methods. Please see the Workflow vignettes for more information.

### 28. How is confounder balance assessed for the IPTW weights in *devMSMs*?

Please see [Jackson, 2016](https://pubmed.ncbi.nlm.nih.gov/27479649/)
for more information about assessing balance for time-varying exposures.

### 29. What if there are confounders that cannot be successfully balanced using IPTW with *devMSMs*?

With real-world data, it is possible that some confounders with respect
to exposure at some timepoints cannot be fully balanced using the
best-performing IPTW weights. You could also consider squaring or
interacting the confounder with others. If the imbalanced confounder is
*time invariant*, you could consider adding it as a covariate to the
final outcome model (see Workflow vignettes). If the imbalanced
confounder is time-varying, we recommend including it as a point in your
discussion section where you could talk about the magnitude of its
remaining relation to the exposure (from the balance assessment output,
see Workflow vignettes), at which timepoint, and the extent to which you
believe it could still be a confounder.

## General

### 30. What assumptions do I make when using *devMSMs*?

The assumptions of using MSMs with IPTWs are as follows:

- No unmeasured confounders: assumes that the only confounders that
  could bias our estimates are those that we *(a)* measured in our
  study, and *(b)* included in our correctly specified balancing model.
  This assures that individuals and their potential outcomes are
  exchangeable, and potential outcomes are equally likely across
  exposure levels (i.e., like random assignment).

- Positivity: assumption holds that each individual has a non-zero
  probability of receiving each level of the exposure, given their
  confounders and prior exposure history.

- Consistency: assumption states that the observed outcome accurately
  reflects the potential outcome for the exposure received.
  Conceptually, this means that the observed outcome is what was caused
  by the exposure condition. For example, there was no subgroup of
  individuals who actually experienced high levels of economic strain
  yet have invalid observed scores because they understood the construct
  by a different name or maintained a different scaling. Additionally,
  it means that the potential outcomes for a given child do not depend
  on the exposure status of other units in the sample.

### 31. How can I be sure that I have measured all possible confounders of my exposure and outcome?

Unfortunately, you can never be absolutely certain you have measured all
possible confounders as this assumption cannot be empirically
verifiable. As with all scientific endeavours, we can only do our best
and justify our efforts. In this case, we must do our due diligence
using theory and existing literature to think of and measure as many
confounders as possible.

### 32. Can I use *devMSMs* to do causal mediation?

We are in the process of formulating a vignette for one approach to
causal mediation using the package. Stay tuned!

### 33. What other packages do *devMSMs* functions draw on?

*devMSMs* would not be possible without several other amazing existing
packages, namely [*WeightIt*](https://ngreifer.github.io/WeightIt/) (for
fitting IPTW weights) and [*cobalt*](https://ngreifer.github.io/cobalt/)
(for assessing IPTW weights) by our collaborator Noah Greifer, and
[*marginaleffects*](https://marginaleffects.com/) *(*for comparing
exposure histories) by Vincent Arel-Bundock. We highly suggest you check
out their work!

### 34. Where can I find more background reading on MSMs?

In addition to the
[manuscript](https://academic.oup.com/chidev/article/97/2/331/8516421?utm_source=authortollfreelink&utm_campaign=chidev&utm_medium=email&guestAccessKey=)
accompanying this package, this repository of [additional
resources](https://drive.google.com/drive/folders/1naueDwrm5LUa4QAw_e0IThrXFKRruN_3)
is a good place to start!

### 35. Who should I contact if I have additional questions about implementing *devMSMs* for my project?

Contact Isa at istall@seas.upenn.edu

## References

Granger, E., Sergeant, J. C., & Lunt, M. (2019). Avoiding pitfalls when
combining multiple imputation and propensity scores. *Statistics in
Medicine*, *38*(26), 5120–5132. <https://doi.org/10.1002/sim.8355>

Imbens, G. W., & Rubin, D. B. (2015). *Causal Inference in Statistics,
Social, and Biomedical Sciences*. Cambridge University Press.

Jackson, J. W. (2016). Diagnostics for Confounding of Time-varying and
Other Joint Exposures. *Epidemiology*, *27*(6), 859.
<https://doi.org/10.1097/EDE.0000000000000547>

Leyrat, C., Carpenter, J. R., Bailly, S., & Williamson, E. J. (2021).
Common Methods for Handling Missing Data in Marginal Structural Models:
What Works and Why. *American Journal of Epidemiology*, *190*(4),
663–672. <https://doi.org/10.1093/aje/kwaa225>

Shook-Sa, B. E., & Hudgens, M. G. (2020). *Power and Sample Size for
Marginal Structural Models* (arXiv:2003.05979). arXiv.
<https://doi.org/10.48550/arXiv.2003.05979>

Stallworthy, I. C., DeJoseph, M. L., Padrutt, E. R., Greifer, N., &
Berry, D. (2026). Investigating causal questions about temporal and
cumulative developmental effects: An introduction to the devMSMs package
in R. *Child Development*, *97*(2), 331–351.
<https://doi.org/10.1093/chidev/aacaf024>
