# Fit outcome model

Fits weighted marginal outcome model as a generalized linear model of
the user's choosing, relating exposure main effects to outcome using
IPTW weights.

## Usage

``` r
fitModel(
  data,
  obj,
  weights = NULL,
  outcome,
  model = "m0",
  int_order = NA,
  covariates = NULL,
  family = gaussian(),
  link = "identity",
  verbose = FALSE,
  save.out = FALSE
)

# S3 method for class 'devMSM_models'
print(x, i = NA, save.out = FALSE, ...)
```

## Arguments

- data:

  data in wide format as: a data frame, list of imputed data frames, or
  `mids` object from the `mice` package

- obj:

  initialized MSM object from
  [`initMSM()`](https://istallworthy.github.io/devMSMs/reference/initMSM.md)

- weights:

  list of IPTW weights output from
  [`createWeights()`](https://istallworthy.github.io/devMSMs/reference/createWeights.md)

- outcome:

  name of outcome variable with ".timepoint" suffix. See
  [`initMSM()`](https://istallworthy.github.io/devMSMs/reference/initMSM.md)
  for details on suffix

- model:

  character indicating one of the following outcome models:

  - "m0" (exposure main effects)

  - "m1" (exposure main effects & covariates)

  - "m2" (exposure main effects & their interactions)

  - "m3" (exposure main effects, their interactions, & covariates)

- int_order:

  integer specification of highest order exposure main effects
  interaction, required for interaction models ("m2", "m3")

- covariates:

  list of characters reflecting variable names of covariates, required
  for covariate models ("m1", "m3")

- family:

  (optional) family function specification for
  [`WeightIt::glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.html)
  model. Note that this should be specified as as a function, not a
  character string unless you have a multinomial outcome, in which case
  set this to "multinomial", or if you have an ordinal outcome with more
  than 2 levels, in which case set this to "ordinal". (default is
  gaussian)

- link:

  (optional) link function specification for
  [`WeightIt::glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.html)
  model.

- verbose:

  (optional) TRUE or FALSE indicator for printing output to console.
  default is FALSE.

- save.out:

  (optional) Either logical or a character string. If `TRUE`, it will
  output the result to a default file name within `home_dir` set in
  [`initMSM()`](https://istallworthy.github.io/devMSMs/reference/initMSM.md).
  You can load the data with `x <- readRDS(file)`. To use a non-default
  file name, specify a character string with the file name. It will save
  relative to `home_dir`. There might be naming conflicts where two
  objects get saved to the same file. In these cases, users should
  specify a custom name. default is FALSE.

- x:

  devMSM_models object from `fitModel`

- i:

  For multiply imputed datasets, `i` selects which imputation to print
  results for. Default is `i = 1`. With `i = TRUE`, all imputed datasets
  will be looped over. With `i = NULL`, will average over all imputed
  datasets and summarize that. Ignored for non-imputed data.

- ...:

  ignored

## Value

list containing
[`WeightIt::glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.html)
model output. It is the length of the number of datasets (1 for a
data.frame or the number of imputed datasets)

## See also

[`WeightIt::glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.html)
for more on family/link specifications.

## Examples

``` r
library(devMSMs)
data <- data.frame(
  ID = 1:50,
  A.1 = rnorm(n = 50),
  A.2 = rnorm(n = 50),
  A.3 = rnorm(n = 50),
  B.1 = rnorm(n = 50),
  B.2 = rnorm(n = 50),
  B.3 = rnorm(n = 50),
  C = rnorm(n = 50),
  D.3 = rnorm(n = 50)
)
obj <- initMSM(
  data,
  exposure = c("A.1", "A.2", "A.3"),
  ti_conf = c("C"),
  tv_conf = c("B.1", "B.2", "B.3", "D.3")
)
f <- createFormulas(obj, type = "short")
w <- createWeights(data = data, formulas = f)

fit_m0 <- fitModel(
  data = data, weights = w, 
  outcome = "D.3", model = "m0"
)
print(fit_m0)
#> Please inspect the Wald test to determine if the exposures collectively predict significant variation in the outcome compared to a model without exposure terms.
#> We strongly suggest only conducting history comparisons if the test is significant.
#> 
#> Wald test
#> Variance: HC0 robust (adjusted for estimation of weights)
#> 
#> Model 1: D.3 ~ A.1 + A.2 + A.3
#> Model 2: D.3 ~ 1
#> 
#>   Res.Df Df  Chisq Pr(>Chisq)
#> 1     46                     
#> 2     49  3 1.7185     0.6328
#> 
#> 
#> The marginal model, m0, is summarized below:
#> +-------------+--------+-----------------+-------+
#> |             | (1)                              |
#> +-------------+--------+-----------------+-------+
#> |             | Est.   | CI              | p     |
#> +=============+========+=================+=======+
#> | (Intercept) | -0.265 | [-0.546, 0.015] | 0.063 |
#> +-------------+--------+-----------------+-------+
#> | A.1         | 0.086  | [-0.142, 0.314] | 0.461 |
#> +-------------+--------+-----------------+-------+
#> | A.2         | -0.053 | [-0.325, 0.220] | 0.705 |
#> +-------------+--------+-----------------+-------+
#> | A.3         | 0.149  | [-0.140, 0.437] | 0.312 |
#> +-------------+--------+-----------------+-------+
#> | Num.Obs.    | 50     |                 |       |
#> +-------------+--------+-----------------+-------+ 

fit_m1 <- fitModel(
  data = data, weights = w, 
  outcome = "D.3", model = "m1", 
  covariates = c("C")
)
print(fit_m1)
#> Please inspect the Wald test to determine if the exposures collectively predict significant variation in the outcome compared to a model without exposure terms.
#> We strongly suggest only conducting history comparisons if the test is significant.
#> 
#> Wald test
#> Variance: HC0 robust (adjusted for estimation of weights)
#> 
#> Model 1: D.3 ~ A.1 + A.2 + A.3 + C
#> Model 2: D.3 ~ C
#> 
#>   Res.Df Df  Chisq Pr(>Chisq)
#> 1     45                     
#> 2     48  3 1.6466     0.6489
#> 
#> 
#> The marginal model, m1, is summarized below:
#> +-------------+--------+------------------+-------+
#> |             | (1)                               |
#> +-------------+--------+------------------+-------+
#> |             | Est.   | CI               | p     |
#> +=============+========+==================+=======+
#> | (Intercept) | -0.273 | [-0.533, -0.012] | 0.040 |
#> +-------------+--------+------------------+-------+
#> | A.1         | 0.065  | [-0.159, 0.289]  | 0.568 |
#> +-------------+--------+------------------+-------+
#> | A.2         | -0.064 | [-0.353, 0.225]  | 0.665 |
#> +-------------+--------+------------------+-------+
#> | A.3         | 0.149  | [-0.129, 0.427]  | 0.294 |
#> +-------------+--------+------------------+-------+
#> | C           | 0.334  | [0.074, 0.594]   | 0.012 |
#> +-------------+--------+------------------+-------+
#> | Num.Obs.    | 50     |                  |       |
#> +-------------+--------+------------------+-------+ 

fit_m2 <- fitModel(
  data = data, weights = w, 
  outcome = "D.3", model = "m2", 
  int_order = 2
)
print(fit_m2)
#> Please inspect the Wald test to determine if the exposures collectively predict significant variation in the outcome compared to a model without exposure terms.
#> We strongly suggest only conducting history comparisons if the test is significant.
#> 
#> Wald test
#> Variance: HC0 robust (adjusted for estimation of weights)
#> 
#> Model 1: D.3 ~ A.1 + A.2 + A.3 + A.1 + A.2 + A.3 + A.1:A.2 + A.1:A.3 + A.2:A.3
#> Model 2: D.3 ~ 1
#> 
#>   Res.Df Df  Chisq Pr(>Chisq)
#> 1     43                     
#> 2     49  6 2.8312     0.8297
#> 
#> 
#> The marginal model, m2, is summarized below:
#> +-------------+--------+-----------------+-------+
#> |             | (1)                              |
#> +-------------+--------+-----------------+-------+
#> |             | Est.   | CI              | p     |
#> +=============+========+=================+=======+
#> | (Intercept) | -0.268 | [-0.556, 0.020] | 0.068 |
#> +-------------+--------+-----------------+-------+
#> | A.1         | 0.072  | [-0.173, 0.316] | 0.566 |
#> +-------------+--------+-----------------+-------+
#> | A.2         | -0.028 | [-0.307, 0.250] | 0.843 |
#> +-------------+--------+-----------------+-------+
#> | A.3         | 0.142  | [-0.143, 0.427] | 0.329 |
#> +-------------+--------+-----------------+-------+
#> | A.1 × A.2   | -0.088 | [-0.395, 0.219] | 0.574 |
#> +-------------+--------+-----------------+-------+
#> | A.1 × A.3   | 0.010  | [-0.376, 0.397] | 0.958 |
#> +-------------+--------+-----------------+-------+
#> | A.2 × A.3   | -0.095 | [-0.423, 0.233] | 0.571 |
#> +-------------+--------+-----------------+-------+
#> | Num.Obs.    | 50     |                 |       |
#> +-------------+--------+-----------------+-------+ 

fit_m3 <- fitModel(
  data = data, weights = w, 
  outcome = "D.3", model = "m3",
  int_order = 2, covariates = c("C")
)
print(fit_m3)
#> Please inspect the Wald test to determine if the exposures collectively predict significant variation in the outcome compared to a model without exposure terms.
#> We strongly suggest only conducting history comparisons if the test is significant.
#> 
#> Wald test
#> Variance: HC0 robust (adjusted for estimation of weights)
#> 
#> Model 1: D.3 ~ A.1 + A.2 + A.3 + A.1 + A.2 + A.3 + C + A.1:A.2 + A.1:A.3 + A.1:C + A.2:A.3 + A.2:C + A.3:C
#> Model 2: D.3 ~ C
#> 
#>   Res.Df Df  Chisq Pr(>Chisq)
#> 1     39                     
#> 2     48  9 10.134     0.3397
#> 
#> 
#> The marginal model, m3, is summarized below:
#> +-------------+--------+------------------+-------+
#> |             | (1)                               |
#> +-------------+--------+------------------+-------+
#> |             | Est.   | CI               | p     |
#> +=============+========+==================+=======+
#> | (Intercept) | -0.232 | [-0.486, 0.023]  | 0.075 |
#> +-------------+--------+------------------+-------+
#> | A.1         | 0.049  | [-0.223, 0.322]  | 0.723 |
#> +-------------+--------+------------------+-------+
#> | A.2         | -0.046 | [-0.333, 0.240]  | 0.751 |
#> +-------------+--------+------------------+-------+
#> | A.3         | 0.126  | [-0.141, 0.392]  | 0.354 |
#> +-------------+--------+------------------+-------+
#> | C           | 0.358  | [0.066, 0.649]   | 0.016 |
#> +-------------+--------+------------------+-------+
#> | A.1 × A.2   | 0.067  | [-0.245, 0.380]  | 0.673 |
#> +-------------+--------+------------------+-------+
#> | A.1 × A.3   | 0.195  | [-0.156, 0.546]  | 0.276 |
#> +-------------+--------+------------------+-------+
#> | A.1 × C     | -0.076 | [-0.318, 0.167]  | 0.540 |
#> +-------------+--------+------------------+-------+
#> | A.2 × A.3   | -0.164 | [-0.480, 0.151]  | 0.308 |
#> +-------------+--------+------------------+-------+
#> | A.2 × C     | -0.343 | [-0.684, -0.001] | 0.049 |
#> +-------------+--------+------------------+-------+
#> | A.3 × C     | -0.108 | [-0.333, 0.117]  | 0.347 |
#> +-------------+--------+------------------+-------+
#> | Num.Obs.    | 50     |                  |       |
#> +-------------+--------+------------------+-------+ 

data <- data.frame(
  ID = 1:50,
  A.1 = rnorm(n = 50),
  A.2 = rnorm(n = 50),
  A.3 = rnorm(n = 50),
  B.1 = rnorm(n = 50),
  B.2 = rnorm(n = 50),
  B.3 = rnorm(n = 50),
  C = rnorm(n = 50),
  D.3 = c(rep(c("A", "B", "C"), 16), "A", "B")
)
obj <- initMSM(
  data,
  exposure = c("A.1", "A.2", "A.3"),
  ti_conf = c("C"),
  tv_conf = c("B.1", "B.2", "B.3", "D.3")
)
f <- createFormulas(obj, type = "short")
w <- createWeights(data = data, formulas = f)

fit_m0 <- fitModel(
  data = data, weights = w, 
  outcome = "D.3", model = "m0", family = "multinomial"
)
print(fit_m0)
#> Please inspect the Wald test to determine if the exposures collectively predict significant variation in the outcome compared to a model without exposure terms.
#> We strongly suggest only conducting history comparisons if the test is significant.
#> 
#> Wald test
#> Variance: HC0 robust (adjusted for estimation of weights)
#> 
#> Model 1: D.3 ~ A.1 + A.2 + A.3
#> Model 2: D.3 ~ 1
#> 
#>   Res.Df Df Chisq Pr(>Chisq)
#> 1     42                    
#> 2     48  6 7.964     0.2408
#> 
#> 
#> The marginal model, m0, is summarized below:
#> +-------------+----------+--------+-----------------+-------+
#> |                        | (1)                              |
#> +-------------+----------+--------+-----------------+-------+
#> |             | response | Est.   | CI              | p     |
#> +=============+==========+========+=================+=======+
#> | (Intercept) | B        | -0.245 | [-0.989, 0.499] | 0.519 |
#> +-------------+----------+--------+-----------------+-------+
#> |             | C        | -0.253 | [-0.971, 0.464] | 0.489 |
#> +-------------+----------+--------+-----------------+-------+
#> | A.1         | B        | 0.315  | [-0.582, 1.212] | 0.491 |
#> +-------------+----------+--------+-----------------+-------+
#> |             | C        | -0.763 | [-1.814, 0.289] | 0.155 |
#> +-------------+----------+--------+-----------------+-------+
#> | A.2         | B        | 0.147  | [-0.594, 0.889] | 0.697 |
#> +-------------+----------+--------+-----------------+-------+
#> |             | C        | 0.295  | [-0.392, 0.982] | 0.400 |
#> +-------------+----------+--------+-----------------+-------+
#> | A.3         | B        | 0.092  | [-0.753, 0.937] | 0.831 |
#> +-------------+----------+--------+-----------------+-------+
#> |             | C        | -1.030 | [-2.093, 0.034] | 0.058 |
#> +-------------+----------+--------+-----------------+-------+
#> | Num.Obs.    |          | 50     |                 |       |
#> +-------------+----------+--------+-----------------+-------+ 

```
