# Initial step in `devMSMs` workflow

Initial step in `devMSMs` workflow

## Usage

``` r
initMSM(
  data,
  exposure,
  epoch = NULL,
  tv_conf,
  ti_conf = NULL,
  concur_conf = NULL,
  home_dir = NULL,
  sep = "[\\._]"
)
```

## Arguments

- data:

  data in wide format as: a data frame, list of imputed data frames, or
  `mids` object from the `mice` package

- exposure:

  names of exposure variables with ".timepoint" suffix

- epoch:

  (optional) group set of exposure variables into categories. Provide a
  character vector corresponding to the category for each exposure

- tv_conf:

  list of time-varying confounders with ".timepoint" suffix, should
  include exposure and outcome variables (at least time-varying exposure
  variables required here)

- ti_conf:

  list of time invariant confounders. Can be left as NULL for none.

- concur_conf:

  (optional) list of variable names reflecting time-varying confounders
  to retain in formulas contemporaneously (default is none)

- home_dir:

  (optional) directory for saving output. Either an absolute path or a
  relative path with respect to
  [`getwd()`](https://rdrr.io/r/base/getwd.html)

- sep:

  (optional) The seperator between the variable and the time period. The
  variable names will be split by the last occurance of `sep` with the
  second string containing the time. This uses regex notation, so `.`
  must be `\\.`

## Value

object of class `devMSM` that contains the initialized information.

## Details

By `.timepoint` suffix, we mean that time-varying and exposure variable
names must end with either `.#` or a `_#`. This allows us to extract the
time-period when the variable is measured and will allow us to properly
create the formulae (omitting future mediators)

## Examples

``` r
data <- data.frame(
  A.1 = rnorm(n = 50),
  A.2 = rnorm(n = 50),
  A.3 = rnorm(n = 50),
  B.1 = rnorm(n = 50),
  B.2 = rnorm(n = 50),
  B.3 = rnorm(n = 50),
  D.3 = rnorm(n = 50),
  L.1 = sample(c(0, 1), size = 50, replace = TRUE),
  C   = rnorm(n = 50)
)
obj <- initMSM(
  data = data,
  exposure = c("A.1", "A.2", "A.3"),
  tv_conf = c("B.1", "B.2", "B.3", "D.3"),
  ti_conf = "C"
)

obj
#> Exposure (continuous): A.1, A.2, A.3
#> Variable and their encodings:
#>  var     type time
#>  A.1 exposure    1
#>  A.2 exposure    2
#>  A.3 exposure    3
#>  B.1  tv_conf    1
#>  B.2  tv_conf    2
#>  B.3  tv_conf    3
#>  D.3  tv_conf    3
#>    C  ti_conf   -1
```
