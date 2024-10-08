% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assessBalance.R
\name{assessBalance}
\alias{assessBalance}
\alias{print.devMSM_bal_stats}
\alias{summary.devMSM_bal_stats}
\alias{plot.devMSM_bal_stats}
\title{Assesses confounder balancing}
\usage{
assessBalance(
  data,
  obj,
  weights = NULL,
  balance_thresh = NULL,
  imp_conf = NULL,
  verbose = FALSE,
  save.out = FALSE
)

\method{print}{devMSM_bal_stats}(x, i = NA, t = TRUE, save.out = FALSE, ...)

\method{summary}{devMSM_bal_stats}(object, i = NA, save.out = FALSE, ...)

\method{plot}{devMSM_bal_stats}(x, i = NA, t = TRUE, save.out = FALSE, ...)
}
\arguments{
\item{data}{data in wide format as: a data frame, list of imputed data
frames, or \code{mids} object from the \code{mice} package}

\item{obj}{initialized MSM object from \code{\link[=initMSM]{initMSM()}}}

\item{weights}{(optional) list of IPTW weights output from \code{\link[=createWeights]{createWeights()}}}

\item{balance_thresh}{(optional) one or two numbers between 0 and 1
indicating a single balancing threshold or thresholds for more and less
important confounders, respectively (default = 0.1)}

\item{imp_conf}{(optional) list of variable names reflecting important
confounders, required if two balance thresholds are supplied}

\item{verbose}{(optional) TRUE or FALSE indicator for printing output to console.
default is FALSE.}

\item{save.out}{(optional) Either logical or a character string. If \code{TRUE},
it will output the result to a default file name within \code{home_dir} set in \code{initMSM()}. You can load the data with \code{x <- readRDS(file)}.
To use a non-default file name, specify a character string with the file name. It will save relative to \code{home_dir}.
There might be naming conflicts where two objects get saved to the same file. In these cases, users should specify a custom name.
default is FALSE.}

\item{x}{devMSM_bal_stats object from \code{assessBalance}}

\item{i}{For multiply imputed datasets, \code{i} selects which imputation to print results for.
Default is \code{i = 1}. With \code{i = TRUE}, all imputed datasets will be looped over. With \code{i = NULL}, will average over all imputed datasets and summarize that. Ignored for non-imputed data.}

\item{t}{Which exposure variable to use. Can either be an index from 1 to the number of exposures or
a string containing the exposure variable name (e.g. \code{"A.3"}). With \code{t = TRUE}, all exposure
variables will be looped over.}

\item{...}{ignored}

\item{object}{devMSM_bal_stats object from \code{\link[=assessBalance]{assessBalance()}}}
}
\value{
a list containing balance statistics as a dataframe. It is the length
of the number of datasets (1 for a data.frame or the number of imputed datasets)
}
\description{
Draws on functions from the cobalt package to quantify the relations between
exposure and confounders at each exposure time point according to the
guidelines from Jackson, 2016 on how to assess balance for time-varying
exposures.
}
\examples{
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

# Prebalance
b <- assessBalance(data = data, obj = obj)
print(b)

# returns ggplot of balance stats for all exposure variables
plots <- plot(b, t = TRUE)
# can plot only specific exposure time periods
plot(b, t = "A.3")
plot(b, t = 3)

# Weighted
f <- createFormulas(obj, type = "short")
w <- createWeights(data = data, formulas = f)
bw <- assessBalance(data = data, weights = w)
print(bw)
plot(bw)

}
\seealso{
\link{cobalt} package, \url{https://cran.r-project.org/web/packages/cobalt/index.html};
Jackson, 2016 for more on assessing balance for time-varying exposures,
\url{https://pubmed.ncbi.nlm.nih.gov/27479649/}
}
