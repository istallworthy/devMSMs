# Customize Weights Formulas

This vignette provides guidance on how to customize weights formulas for
use in *devMSMs*. Please first review the [Specify Core
Inputs](https://istallworthy.github.io/devMSMs/articles/Specify_Core_Inputs.html)
and any of the *Workflows* vignettes.

The code contained in this vignette is also available, integrated code
from the other vignettes, in the ExampleWorkflow.rmd.

## Core Inputs

We specify core inputs to the *devMSMs* functions.

``` r
exposure <- c("ESETA1.6", "ESETA1.15", "ESETA1.24", "ESETA1.35", "ESETA1.58")

tv_conf <- c(
  "SAAmylase.6", "SAAmylase.15", "SAAmylase.24",
  "MDI.6", "MDI.15",
  "RHasSO.6", "RHasSO.15", "RHasSO.24", "RHasSO.35",
  "WndNbrhood.6", "WndNbrhood.24", "WndNbrhood.35",
  "IBRAttn.6", "IBRAttn.15", "IBRAttn.24",
  "B18Raw.6", "B18Raw.15", "B18Raw.24",
  "HOMEETA1.6", "HOMEETA1.15", "HOMEETA1.24", "HOMEETA1.35", 
  "InRatioCor.6", "InRatioCor.15", "InRatioCor.24", "InRatioCor.35",
  "CORTB.6", "CORTB.15", "CORTB.24",
  "EARS_TJo.24", "EARS_TJo.35",
  "LESMnPos.24", "LESMnPos.35",
  "LESMnNeg.24", "LESMnNeg.35",
  "StrDif_Tot.35", 
  "fscore.35"
)

ti_conf <- c(
  "state", "BioDadInHH2", "PmAge2", "PmBlac2", "TcBlac2", "PmMrSt2", "PmEd2", "KFASTScr",
  "RMomAgeU", "RHealth", "HomeOwnd", "SWghtLB", "SurpPreg", "SmokTotl", "DrnkFreq",
  "peri_health", "caregiv_health", "gov_assist"
)

data(sim_data_wide, package = "devMSMs")

obj <- devMSMs::initMSM(
  data = sim_data_wide, 
  exposure = exposure, 
  tv_conf = tv_conf, 
  ti_conf = ti_conf
)
```

  

The
[`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)
function has 3 different modes: “full” (Steps 1 & 4), “short” (Step 2),
and “update” (Step 3), corresponding to different steps in the
*Workflows* vignettes, and can be customized in each mode. We suggest
staying consistent with any customization throughout use of package
workflow. We recommend being sure to include all measured confounders
(and no mediators or colliders) for each exposure time point in the
weights formulas. Please refer to the accompanying manuscript for more
details.

Below, we illustrate the customization options for each of the three
modes of the `createFormuls()` function using information from the
package dataset as an example. These data are simulated based on data
from the Family Life Project (FLP), a longitudinal study following 1,292
families representative of two geographic areas (three counties in North
Carolina and three counties in Pennsylvania) with high rural child
poverty (Vernon-Feagans et al., 2013; Burchinal et al., 2008). We take
the example exposure of economic strain (ESETA1) measured at 6, 15, 24,
35, and 58 months in relation to the outcome of behavior problems
(StrDif_Tot) measured at 58 months.

  

## Customize Full Weights Formulas (Steps 1 & 4)

The default use of
[`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)
with `type` = “full” creates weights formulas at each exposure time
point. Each formula contains all time invariant confounders as well as
all time-varying confounders measured at time points prior to that
exposure time point. For example, for the 24-month weights formula, we
would have time-varying confounders only at 6 and 15 months.

### Time-Varying Confounders Measured Contemporaneously With Exposure Time Point

However, if users wish to also include any time-varying confounders
measured contemporaneously with the exposure time point (e.g., also
include 24-month time-varying confounders in the 24-month weights
formula), they can specify those time-varying confounders in the
optional argument to
[`initMSM()`](https://istallworthy.github.io/devMSMs/reference/initMSM.md),
`concur_conf`.

For example, to include 24-month maternal depression in the 24-month
weights formula, we do the following:

We first specify confounders to retain contemporaneously when creating
our MSM object.

``` r
obj <- devMSMs::initMSM(
  data = sim_data_wide, 
  exposure = exposure, 
  tv_conf = tv_conf, 
  ti_conf = ti_conf,
  concur_conf = "B18Raw.24"
)
```

``` r
full_formulas <- createFormulas(obj  = obj, 
                                type = "full")

print(full_formulas)
#> USER ALERT: Please manually inspect the full balancing formula below:
#> At time point 6, the full formula for ESETA1.6 is: 
#> ESETA1.6 ~ state + BioDadInHH2 + PmAge2 + PmBlac2 + TcBlac2 + PmMrSt2 + PmEd2 + KFASTScr + RMomAgeU + RHealth + HomeOwnd + SWghtLB + SurpPreg + SmokTotl + DrnkFreq + peri_health + caregiv_health + gov_assist
#> USER ALERT: Please manually inspect the full balancing formula below:
#> At time point 15, the full formula for ESETA1.15 is: 
#> ESETA1.15 ~ state + BioDadInHH2 + PmAge2 + PmBlac2 + TcBlac2 + PmMrSt2 + PmEd2 + KFASTScr + RMomAgeU + RHealth + HomeOwnd + SWghtLB + SurpPreg + SmokTotl + DrnkFreq + peri_health + caregiv_health + gov_assist + SAAmylase.6 + MDI.6 + RHasSO.6 + WndNbrhood.6 + IBRAttn.6 + B18Raw.6 + HOMEETA1.6 + InRatioCor.6 + CORTB.6 + ESETA1.6
#> USER ALERT: Please manually inspect the full balancing formula below:
#> At time point 24, the full formula for ESETA1.24 is: 
#> ESETA1.24 ~ state + BioDadInHH2 + PmAge2 + PmBlac2 + TcBlac2 + PmMrSt2 + PmEd2 + KFASTScr + RMomAgeU + RHealth + HomeOwnd + SWghtLB + SurpPreg + SmokTotl + DrnkFreq + peri_health + caregiv_health + gov_assist + SAAmylase.6 + SAAmylase.15 + MDI.6 + MDI.15 + RHasSO.6 + RHasSO.15 + WndNbrhood.6 + IBRAttn.6 + IBRAttn.15 + B18Raw.6 + B18Raw.15 + B18Raw.24 + HOMEETA1.6 + HOMEETA1.15 + InRatioCor.6 + InRatioCor.15 + CORTB.6 + CORTB.15 + ESETA1.6 + ESETA1.15
#> USER ALERT: Please manually inspect the full balancing formula below:
#> At time point 35, the full formula for ESETA1.35 is: 
#> ESETA1.35 ~ state + BioDadInHH2 + PmAge2 + PmBlac2 + TcBlac2 + PmMrSt2 + PmEd2 + KFASTScr + RMomAgeU + RHealth + HomeOwnd + SWghtLB + SurpPreg + SmokTotl + DrnkFreq + peri_health + caregiv_health + gov_assist + SAAmylase.6 + SAAmylase.15 + SAAmylase.24 + MDI.6 + MDI.15 + RHasSO.6 + RHasSO.15 + RHasSO.24 + WndNbrhood.6 + WndNbrhood.24 + IBRAttn.6 + IBRAttn.15 + IBRAttn.24 + B18Raw.6 + B18Raw.15 + B18Raw.24 + HOMEETA1.6 + HOMEETA1.15 + HOMEETA1.24 + InRatioCor.6 + InRatioCor.15 + InRatioCor.24 + CORTB.6 +      CORTB.15 + CORTB.24 + EARS_TJo.24 + LESMnPos.24 + LESMnNeg.24 + ESETA1.6 + ESETA1.15 + ESETA1.24
#> USER ALERT: Please manually inspect the full balancing formula below:
#> At time point 58, the full formula for ESETA1.58 is: 
#> ESETA1.58 ~ state + BioDadInHH2 + PmAge2 + PmBlac2 + TcBlac2 + PmMrSt2 + PmEd2 + KFASTScr + RMomAgeU + RHealth + HomeOwnd + SWghtLB + SurpPreg + SmokTotl + DrnkFreq + peri_health + caregiv_health + gov_assist + SAAmylase.6 + SAAmylase.15 + SAAmylase.24 + MDI.6 + MDI.15 + RHasSO.6 + RHasSO.15 + RHasSO.24 + RHasSO.35 + WndNbrhood.6 + WndNbrhood.24 + WndNbrhood.35 + IBRAttn.6 + IBRAttn.15 + IBRAttn.24 + B18Raw.6 + B18Raw.15 + B18Raw.24 + HOMEETA1.6 + HOMEETA1.15 + HOMEETA1.24 + HOMEETA1.35 + InRatioCor.6 +      InRatioCor.15 + InRatioCor.24 + InRatioCor.35 + CORTB.6 + CORTB.15 + CORTB.24 + EARS_TJo.24 + EARS_TJo.35 + LESMnPos.24 + LESMnPos.35 + LESMnNeg.24 + LESMnNeg.35 + StrDif_Tot.35 + fscore.35 + ESETA1.6 + ESETA1.15 + ESETA1.24 + ESETA1.35
```

Now `B18Raw.24` is included in the 24-month weights formula in addition
to the 35-and 58-month formulas.

## Fully Custom Formulas

Alternatively, a user could completely customize the full formulas by
manually creating formulas for each exposure time point using the
optional `custom` argument.

To create custom formulas, the user creates a list with an entry for
each exposure time point, each of which contains the weights formula for
exposure at that time point.
[`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)
will conduct checks to ensure there is a formula for each time point and
provides cautionary warnings for unusual variables (e.g., those measured
after the exposure time point). However, the user is responsible for
ensuring the validity of these formulas.

Of note, if the user creates a custom “full” formula, they will also
have to create custom formulas when using
[`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)
to make “short” and “updated” formulas.

Below, we illustrate an abridged example for specifying the `custom`
arguments for the exposure, ESETA1, at 6-, 15-, 24-, 35-, and 58-month
time points.

``` r
custom <- list(
  ESETA1.6 ~ BioDadInHH2 + DrnkFreq + gov_assist,
  ESETA1.15 ~ BioDadInHH2 + DrnkFreq + gov_assist,
  ESETA1.24 ~ BioDadInHH2 + DrnkFreq + gov_assist,
  ESETA1.35 ~ BioDadInHH2 + DrnkFreq + gov_assist,
  ESETA1.58 ~ BioDadInHH2 + DrnkFreq + gov_assist
)
```

[`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)
then checks and retains these formulas for further use with the package.

``` r
full_formulas <- createFormulas(obj = obj, 
                                custom = custom)

print(full_formulas)
#> USER ALERT: Please manually inspect the slightly cleaned custom formula below:
#> At time point 6, the custom formula for ESETA1.6 is: 
#> ESETA1.6 ~ BioDadInHH2 + DrnkFreq + gov_assist
#> USER ALERT: Please manually inspect the slightly cleaned custom formula below:
#> At time point 15, the custom formula for ESETA1.15 is: 
#> ESETA1.15 ~ BioDadInHH2 + DrnkFreq + gov_assist
#> USER ALERT: Please manually inspect the slightly cleaned custom formula below:
#> At time point 24, the custom formula for ESETA1.24 is: 
#> ESETA1.24 ~ BioDadInHH2 + DrnkFreq + gov_assist
#> USER ALERT: Please manually inspect the slightly cleaned custom formula below:
#> At time point 35, the custom formula for ESETA1.35 is: 
#> ESETA1.35 ~ BioDadInHH2 + DrnkFreq + gov_assist
#> USER ALERT: Please manually inspect the slightly cleaned custom formula below:
#> At time point 58, the custom formula for ESETA1.58 is: 
#> ESETA1.58 ~ BioDadInHH2 + DrnkFreq + gov_assist
```

  

## Customize Short Weights Formulas (Step 2)

The default use of
[`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)
with `type` = “short” creates shortened weights formulas at each
exposure time point. Each formula contains all time invariant
confounders as well as only those time-varying confounders measured at
the time point directly prior (*t*-1 lag) to that exposure time point.
For example, for the 24-month weights formula, we would only include
time-varying confounders at 15 months. Please see the *Workflows*
vignettes and accompanying manuscript for more detail.

When creating short formulas, the user can specify the optional
`concur_conf` and `custom` arguments detailed above. Of note,
[`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)
does not check that formulas in the `custom` field meet criteria for the
shortened formula, and thus the user is responsible for ensuring the
validity of these formulas.

When creating short formulas, the user can also specify an optional
`keep_conf` argument. Users would specify this second argument with any
time-varying confounders to always retain in the formulas, in lagged
form, overriding the *t*-1 lag default. Users may wish to use this
argument if they have time-varying confounders that are not highly
consistent over time and have strong reasons to include them in the
initial phase of selecting the optimal weighting method (see *Workflows*
vignettes for more details). Note that unless these variables are also
specified in `concuf_conf`, they will only be retained in lagged form.

For example, to always retain income in all formulas in lagged form,
overriding the *t*-1 lag default via `keep_conf`.

As shown below, the
[`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)
function now retains income at 6 months in all 15-, 24-, 35-, and
58-month formulas, overriding the *t*-1 lag rule for that variable.
Without this specification, income at 6 months would have been omitted
in formulas for exposure at 24, 35, and 58 months.

``` r
short_formulas <- createFormulas(obj = obj, 
                                 type = "short", 
                                 keep_conf = c("InRatioCor.6"))
print(short_formulas)
#> USER ALERT: Please manually inspect the short balancing formula below that includes time-varying confounders at t-1 only:
#> At time point 6, the short formula for ESETA1.6 is: 
#> ESETA1.6 ~ state + BioDadInHH2 + PmAge2 + PmBlac2 + TcBlac2 + PmMrSt2 + PmEd2 + KFASTScr + RMomAgeU + RHealth + HomeOwnd + SWghtLB + SurpPreg + SmokTotl + DrnkFreq + peri_health + caregiv_health + gov_assist
#> USER ALERT: Please manually inspect the short balancing formula below that includes time-varying confounders at t-1 only:
#> At time point 15, the short formula for ESETA1.15 is: 
#> ESETA1.15 ~ state + BioDadInHH2 + PmAge2 + PmBlac2 + TcBlac2 + PmMrSt2 + PmEd2 + KFASTScr + RMomAgeU + RHealth + HomeOwnd + SWghtLB + SurpPreg + SmokTotl + DrnkFreq + peri_health + caregiv_health + gov_assist + SAAmylase.6 + MDI.6 + RHasSO.6 + WndNbrhood.6 + IBRAttn.6 + B18Raw.6 + HOMEETA1.6 + InRatioCor.6 + CORTB.6 + ESETA1.6
#> USER ALERT: Please manually inspect the short balancing formula below that includes time-varying confounders at t-1 only:
#> At time point 24, the short formula for ESETA1.24 is: 
#> ESETA1.24 ~ state + BioDadInHH2 + PmAge2 + PmBlac2 + TcBlac2 + PmMrSt2 + PmEd2 + KFASTScr + RMomAgeU + RHealth + HomeOwnd + SWghtLB + SurpPreg + SmokTotl + DrnkFreq + peri_health + caregiv_health + gov_assist + SAAmylase.15 + MDI.15 + RHasSO.15 + IBRAttn.15 + B18Raw.15 + B18Raw.24 + HOMEETA1.15 + InRatioCor.15 + CORTB.15 + ESETA1.15 + InRatioCor.6
#> USER ALERT: Please manually inspect the short balancing formula below that includes time-varying confounders at t-1 only:
#> At time point 35, the short formula for ESETA1.35 is: 
#> ESETA1.35 ~ state + BioDadInHH2 + PmAge2 + PmBlac2 + TcBlac2 + PmMrSt2 + PmEd2 + KFASTScr + RMomAgeU + RHealth + HomeOwnd + SWghtLB + SurpPreg + SmokTotl + DrnkFreq + peri_health + caregiv_health + gov_assist + SAAmylase.24 + RHasSO.24 + WndNbrhood.24 + IBRAttn.24 + HOMEETA1.24 + InRatioCor.24 + CORTB.24 + EARS_TJo.24 + LESMnPos.24 + LESMnNeg.24 + ESETA1.24 + InRatioCor.6
#> USER ALERT: Please manually inspect the short balancing formula below that includes time-varying confounders at t-1 only:
#> At time point 58, the short formula for ESETA1.58 is: 
#> ESETA1.58 ~ state + BioDadInHH2 + PmAge2 + PmBlac2 + TcBlac2 + PmMrSt2 + PmEd2 + KFASTScr + RMomAgeU + RHealth + HomeOwnd + SWghtLB + SurpPreg + SmokTotl + DrnkFreq + peri_health + caregiv_health + gov_assist + RHasSO.35 + WndNbrhood.35 + HOMEETA1.35 + InRatioCor.35 + EARS_TJo.35 + LESMnPos.35 + LESMnNeg.35 + StrDif_Tot.35 + fscore.35 + ESETA1.35 + InRatioCor.6
```

  

## Customize Updated Weights Formulas (Step 3)

The default use of
[`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)
with `type` = “update” creates updated weights formulas at each exposure
time point. Each formula contains all time invariant confounders;
time-varying confounders measured at the time point directly prior
(*t*-1 lag) to that exposure time point; as well as any time-varying
confounders at greater lags that were not successfully balanced using
the short weights formulas.

When creating updated formulas, the user can specify the optional
`concur_conf`, `keep_conf`, and `custom` arguments detailed above. Of
note,
[`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)
does not check that formulas in the `custom` field meet criteria for the
updated formula or identify any imbalanced confounders, and the user is
responsible for ensuring the validity of these formulas.  

## References

Burchinal, M., Howes, C., Pianta, R., Bryant, D., Early, D., Clifford,
R., & Barbarin, O. (2008). Predicting Child Outcomes at the End of
Kindergarten from the Quality of Pre-Kindergarten Teacher–Child
Interactions and Instruction. Applied Developmental Science, 12(3),
140–153. https://doi.org/10.1080/10888690802199418

Vernon-Feagans, L., Cox, M., Willoughby, M., Burchinal, M.,
Garrett-Peters, P., Mills-Koonce, R., Garrett-Peiers, P., Conger, R. D.,
& Bauer, P. J. (2013). The Family Life Project: An Epidemiological and
Developmental Study of Young Children Living in Poor Rural Communities.
Monographs of the Society for Research in Child Development, 78(5),
i–150.
