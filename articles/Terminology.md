# Terminology

This vignette contains some key definitions for using the [Data
Requirements](https://istallworthy.github.io/devMSMs/articles/Data_Requirements.html),
[Specifying Core
Inputs](https://istallworthy.github.io/devMSMs/articles/Specify_Core_Inputs.html),
[Recommended Preliminary
Steps](https://istallworthy.github.io/devMSMs/articles/Preliminary_Steps.html),
and *Workflows* vignettes. A more thorough introduction to these
concepts can be found in the forthcoming manuscript: *\[insert
link\]*.  
  

Key Terminology

| Term                          | Definition                                                                                                                                                                                                                                                     |
|-------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Exposure**                  | Exposure or experience that constitutes the causal event of interest and is measured at at least two time points, with at least one time point occurring prior to the outcome.                                                                                 |
| **Outcome**                   | Any developmental construct measured at least once at a final outcome time point upon which the exposure is theorized to have causal effects.                                                                                                                  |
| **Exposure Time Points**      | Time points in development when the exposure was measured, at which weights formulas will be created.                                                                                                                                                          |
| **Exposure Epochs**           | *(optional)* Further delineation of exposure time points into meaningful units of developmental time, each of which could encompass multiple exposure time points, that together constitute exposure main effects in the outcome model and exposure histories. |
| **Exposure Histories**        | Sequences of relatively high (`'h'`) or low (`'l'`) levels of exposure at each exposure time point or exposure epoch.                                                                                                                                          |
| **Exposure Dosage**           | Total cumulative exposure epochs/time points during which an individual experienced high (or low) levels of exposure, across an entire exposure history.                                                                                                       |
| **Confounder**                | Pre-exposure variable that represents a common cause of exposure at a given time point and outcome; adjusting for all of which successfully blocks all backdoor paths.                                                                                         |
| **Time-varying confounder**   | A confounder that often changes over time (even if it is not measured at every time point), and is affected by prior exposure, either directly or indirectly.                                                                                                  |
| **Time invariant confounder** | A confounder that occurs only at a single time point, prior to the exposure and remains stable and/or is not possibly affected by exposure.                                                                                                                    |
| **Collider**                  | A variable that represents a common effect of exposure at a given time point and outcome; adjusting for which introduces bias.                                                                                                                                 |

Summary of terminology used in *devMSMs*

  

## Weights Formula

Weights formula are formula, or mathematical equations, created at each
exposure time point, regressing exposure on confounders for the purpose
of creating and/or assessing balancing weights. The general form is
given by:  
  

    ESETA1.t ~ time-invariant confounders + time-varying confounders + lagged exposures + lagged outcomes

  

Below is a table that summarizes the three formula `type`s that will be
used.

Table 2.

[TABLE]

Summary of the formula terminology used in the *devMSMs* package for
creating balancing weights and assessing their performance in devMSMs
(*Workflows* vignettes Steps 4-6). The different kinds of formulas vary
with respect to their inclusion of time-varying confounders and the
`type` specification to create them using the
[`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)
function.

  

## 7 core functions of *devMSMs*

There are 7 key functions that are used in the *devMSMs* workflow.

| Function                                                                          | Purpose                                                                                              | Required.Input                                                        | Output                                                                                                     |
|-----------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------|
| [`initMSM()`](https://istallworthy.github.io/reference/initMSM)                   | Creates MSM object with core variables for use in all other functions.                               | `exposure`, `time invariant confounders`, `time-varying confounders`) | msmObject                                                                                                  |
| [`createFormulas()`](https://istallworthy.github.io/reference/createFormulas.md)  | Creates balancing formulas at each exposure time point relating exposure to confounders.             | `msmObject`, `type` (`'full'`, `'short'`, `'update'`)                 | List of balancing formulas for each exposure time point.                                                   |
| [`assessBalance()`](https://istallworthy.github.io/reference/assessBalance)       | Assesses confounder balance at each exposure time point according to best practices (Jackson, 2016). | `data`, `msmObject`                                                   | List of balance statistics. Tables & love plots displaying balance statistics at each exposure time point. |
| [`createWeights()`](https://istallworthy.github.io/reference/createWeights)       | Creates balancing weights for each person (wrapper for `weightit::weightitMSM`).                     | `data`, `msmObject`                                                   | Weights object. Histogram & summary of balancing weights.                                                  |
| [`trimWeights()`](https://istallworthy.github.io/reference/trimWeights)           | Trims balancing weights to account for heavy tails.                                                  | `data`, `msmObject`                                                   | Trimmed weights object.Histogram & summary of trimmed weight.                                              |
| [`fitModel()`](https://istallworthy.github.io/reference/fitModel)                 | Fits user-selected, weighted marginal model relating exposure to outcome.                            | `data`, `msmObject`, `weights`, `outcome`, `model` (`m0`-`m3`)        | Fitted model object. Omnibus test & table of model evidence.                                               |
| [`compareHistories()`](https://istallworthy.github.io/reference/compareHistories) | Estimate, compare, & visualize model-predicted outcomes as a function of exposure history.           | `data`, `msmObject`, `fitted model`                                   | Tables of estimated values and comparisons. Boxplot.                                                       |

Summary of the 7 core functions of the *devMSMs* package.
