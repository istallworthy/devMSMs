# Creates IPTW balancing weights

Creates IPTW balancing weights at each user-specified exposure time
point using balancing formulas that relate exposure at each time point
to all relevant confounders.

## Usage

``` r
createWeights(
  data,
  formulas,
  method = "glm",
  verbose = FALSE,
  save.out = FALSE,
  ...
)

# S3 method for class 'devMSM_weights'
print(x, i = 1, ...)

# S3 method for class 'devMSM_weights'
plot(x, i = 1, save.out = FALSE, ...)

# S3 method for class 'devMSM_weights'
summary(object, i = 1, ...)
```

## Arguments

- data:

  data in wide format as: a data frame, list of imputed data frames, or
  `mids` object from the `mice` package

- formulas:

  list of balancing formulas at each time point output from
  [`createFormulas()`](https://istallworthy.github.io/devMSMs/reference/createFormulas.md)

- method:

  character string of `weightitMSM()` balancing method abbreviation
  (default is generalized linear models propensity score weighting
  `"glm"`)

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

- ...:

  arguments passed to
  [`WeightIt::weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.html)
  or `summary.weightitMSM()`

- x, object:

  `devMSM_weights` object from `createWeights()`

- i:

  For multiply imputed datasets, `i` selects which imputation to print
  results for. Default is `i = 1`. With `i = TRUE`, all imputed datasets
  will be looped over. With `i = NULL`, will average over all imputed
  datasets and summarize that. Ignored for non-imputed data.

## Value

a list containing
[`WeightIt::weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.html)
output. It is the length of the number of datasets (1 for a data.frame
or the number of imputed datasets).

## See also

[`WeightIt::weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.html)

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
print(w)
#> 
#> For the `glm` weighting method, the median weight value is 0.97 (SD = 0.47; range = 0.38-3).
plot(w)


# Methods from `WeightIt::weightitMSM`
w <- createWeights(data = data, formulas = f,
                   method = "glm")

w <- createWeights(data = data, formulas = f,
                   method = "cbps")
w <- createWeights(data = data, formulas = f,
                   method = "gbm")
w <- createWeights(data = data, formulas = f,
                   method = "bart")
w <- createWeights(data = data, formulas = f,
                   method = "super")
#> Loading required package: nnls
#> Warning: All algorithms have zero weight
#> Warning: All metalearner coefficients are zero, predictions will all be equal to 0
#> Warning: All algorithms have zero weight
#> Warning: All metalearner coefficients are zero, predictions will all be equal to 0
```
