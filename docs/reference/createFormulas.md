# Create balancing formulas

Creates balancing formulas relating exposure to all relevant
time-varying and time invariant confounders at each exposure time point
to be used to create IPTW weights.

## Usage

``` r
createFormulas(
  obj,
  type = c("full", "short", "update"),
  custom = NULL,
  keep_conf = NULL,
  bal_stats = NULL,
  verbose = FALSE,
  save.out = FALSE
)

# S3 method for class 'devMSM_formulas'
print(x, ...)
```

## Arguments

- obj:

  initialized MSM object from
  [`initMSM()`](https://istallworthy.github.io/devMSMs/reference/initMSM.md)

- type:

  type of formula to create from 'full' (includes all lagged
  time-varying confounders), 'short' (includes time-varying confounders
  at t-1 lag only), or 'update' (adds to 'short' formulas any imbalanced
  time-varying confounders at lags great than t-1)

- custom:

  (optional) custom list of formulas at each exposure time point
  (default is to create automatically according to type)

- keep_conf:

  (optional) For 'short' formulas only, list of variable names
  reflecting confounders that should be included always.

- bal_stats:

  list of balance statistics from
  [`assessBalance()`](https://istallworthy.github.io/devMSMs/reference/assessBalance.md)

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

  devMSM_formulas object from `createFormulas()`

- ...:

  ignored

## Value

a list containing balancing formulas. It is the length of the number of
exposure variables.

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

# Full Formulas
f <- createFormulas(obj, type = "full")
print(f)
#> USER ALERT: Please manually inspect the full balancing formula below:
#> At time point 1, the full formula for A.1 is: 
#> A.1 ~ C
#> 
#> USER ALERT: Please manually inspect the full balancing formula below:
#> At time point 2, the full formula for A.2 is: 
#> A.2 ~ C + B.1 + A.1
#> 
#> USER ALERT: Please manually inspect the full balancing formula below:
#> At time point 3, the full formula for A.3 is: 
#> A.3 ~ C + B.1 + B.2 + A.1 + A.2

# Short Formulas
f <- createFormulas(obj, type = "short")
print(f)
#> USER ALERT: Please manually inspect the short balancing formula below that includes time-varying confounders at t-1 only:
#> At time point 1, the short formula for A.1 is: 
#> A.1 ~ C
#> 
#> USER ALERT: Please manually inspect the short balancing formula below that includes time-varying confounders at t-1 only:
#> At time point 2, the short formula for A.2 is: 
#> A.2 ~ C + B.1 + A.1
#> 
#> USER ALERT: Please manually inspect the short balancing formula below that includes time-varying confounders at t-1 only:
#> At time point 3, the short formula for A.3 is: 
#> A.3 ~ C + B.2 + A.2

# Update Formulas
w <- createWeights(data = data, formulas = f)
b <- assessBalance(data = data, weights = w)
f <- createFormulas(obj, type = "update", bal_stats = b)
print(f)
#> USER ALERT: Please manually inspect the updated balancing formula below that includes time-varying confounders at t-1 and those greater at further lags that remained imbalanced:
#> At time point 1, the update formula for A.1 is: 
#> A.1 ~ C
#> 
#> USER ALERT: Please manually inspect the updated balancing formula below that includes time-varying confounders at t-1 and those greater at further lags that remained imbalanced:
#> At time point 2, the update formula for A.2 is: 
#> A.2 ~ C + B.1 + A.1
#> 
#> USER ALERT: Please manually inspect the updated balancing formula below that includes time-varying confounders at t-1 and those greater at further lags that remained imbalanced:
#> At time point 3, the update formula for A.3 is: 
#> A.3 ~ C + B.2 + A.2
```
