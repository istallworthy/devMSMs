#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import stats
#' @import utils
#' @importFrom ggplot2 .data
## usethis namespace: end
NULL

#' @name devMSM_common_docs
#' 
#' @param data data in wide format as: a data frame, list of imputed data
#'  frames, or `mids` object from the `mice` package
#' @param obj initialized MSM object from [initMSM()]
#' @param formulas list of balancing formulas at each time point output from [createFormulas()]
#' @param weights list of IPTW weights output from [createWeights()]
#' @param fit list of model outputs from [fitModel()]
#' @param bal_stats list of balance statistics from [assessBalance()]
#' 
#' @param verbose (optional) TRUE or FALSE indicator for printing output to console.
#' default is FALSE.
#' @param save.out (optional) Either logical or a character string. If `TRUE`, 
#' it will output the result to a default file name within `home_dir` set in `initMSM()`. You can load the data with `x <- readRDS(file)`. 
#' To use a non-default file name, specify a character string with the file name. It will save relative to `home_dir`. 
#' There might be naming conflicts where two objects get saved to the same file. In these cases, users should specify a custom name.
#' default is FALSE. 
#' 
#' @param i For multiply imputed datasets, `i` selects which imputation to print results for.
#'  Default is `i = 1`. With `i = TRUE`, all imputed datasets will be looped over. With `i = NULL`, will average over all imputed datasets and summarize that. Ignored for non-imputed data.
#' @param t Which exposure variable to use. Can either be an index from 1 to the number of exposures or 
#'  a string containing the exposure variable name (e.g. `"A.3"`). With `t = TRUE`, all exposure 
#'  variables will be looped over.
NULL
#Note: see https://github.com/r-lib/roxygen2/issues/1159#issuecomment-708231897 for how this is used

