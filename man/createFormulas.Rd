% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createFormulas.R
\name{createFormulas}
\alias{createFormulas}
\title{Create balancing formulas}
\usage{
createFormulas(
  exposure,
  exposure_time_pts,
  outcome,
  type,
  ti_confounders,
  tv_confounders,
  bal_stats = NULL,
  concur_conf = NULL,
  keep_conf = NULL,
  home_dir = NULL,
  custom = NULL,
  verbose = TRUE,
  save.out = TRUE
)
}
\arguments{
\item{exposure}{name of exposure variable}

\item{exposure_time_pts}{list of integers at which weights will be
created/assessed that correspond to time points when exposure wass measured}

\item{outcome}{name of outcome variable with ".timepoint" suffix}

\item{type}{type of formula to create from 'full' (includes all lagged
time-varying confounders), 'short' (includes time-varying confounders at
t-1 lag only), or 'update' (adds to 'short' formulas any imbalanced
time-varying confounders at lags great than t-1)}

\item{ti_confounders}{list of time invariant confounders (at least one
required)}

\item{tv_confounders}{list of time-varying confounders with ".timepoint"
suffix, should include exposure and outcome variables (at least
time-varying exposure variables required here)}

\item{bal_stats}{list of balance statistics from assessBalance(), required
for 'update' type}

\item{concur_conf}{(optional) list of variable names reflecting time-varying
confounders to retain in formulas contemporaneously (default is none)}

\item{keep_conf}{(optional) list of variable names reflecting confounders to
always retain in formulas (default depends on type)}

\item{home_dir}{path to home directory (required if 'save.out' = TRUE)}

\item{custom}{(optional) custom list of formulas at each exposure time point
(default is to create automatically according to type)}

\item{verbose}{(optional) TRUE or FALSE indicator for user output (default is
TRUE)}

\item{save.out}{(optional) TRUE or FALSE indicator to save output and
intermediary output locally (default is TRUE)}
}
\value{
list of balancing formulas at each exposure time point
}
\description{
Creates balancing formulas relating exposure to all relevant time-varying and
time invariant confounders at each exposure time point to be used to create
IPTW weights.
}
\examples{
#Full Formulas
f <- createFormulas(exposure = "A",
                    exposure_time_pts = c(1, 2, 3),
                    outcome = "D.3",
                    tv_confounders = c("A.1", "A.2", "A.3"),
                    ti_confounders = "C",
                    type = "full",
                    save.out = FALSE)

f <- createFormulas(exposure = "A",
                    exposure_time_pts = c(1, 2, 3),
                    outcome = "D.3",
                    tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
                    ti_confounders = "C",
                    type = "full",
                    save.out = FALSE)

#Short Formulas
f <- createFormulas(exposure = "A",
                    exposure_time_pts = c(1, 2, 3),
                    outcome = "D.3",
                    tv_confounders = c("A.1", "A.2", "A.3"),
                    ti_confounders = "C",
                    type = "short",
                    save.out = FALSE)
f <- createFormulas(exposure = "A",
                    exposure_time_pts = c(1, 2, 3),
                    outcome = "D.3",
                    tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
                    ti_confounders = "C",
                    type = "short",
                    save.out = FALSE)

c <- list("short_form-1" = as.formula(A.1 ~ C),
          "short_form-2" = as.formula(A.2 ~ A.1 + B.1 + C),
          "short_form-3" = as.formula(A.3 ~ A.2 + B.2 + C))

f <- createFormulas(exposure = "A",
                    exposure_time_pts = c(1, 2, 3),
                    outcome = "D.3",
                    tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
                    ti_confounders = "C",
                    type = "short",
                    custom = c,
                    save.out = FALSE)

#Update Formulas
test <- data.frame(ID = 1:50,
                   A.1 = rnorm(n = 50),
                   A.2 = rnorm(n = 50),
                   A.3 = rnorm(n = 50),
                   B.1 = rnorm(n = 50),
                   B.2 = rnorm(n = 50),
                   B.3 = rnorm(n = 50),
                   C = rnorm(n = 50),
                   D.3 = rnorm(n = 50))
test[, c("A.1", "A.2", "A.3")] <- lapply(test[, c("A.1", "A.2", "A.3")], as.numeric)

w <- createWeights(data = test,
                   exposure = "A",
                   outcome = "D.3",
                   formulas = f,
                   save.out = FALSE)

b <- assessBalance(data = test,
                   exposure = "A",
                   exposure_time_pts = c(1, 2, 3),
                   outcome = "D.3",
                   type = "weighted",
                   weights = w,
                   formulas = f,
                   save.out = FALSE)
f <- createFormulas(exposure = "A",
                    exposure_time_pts = c(1, 2, 3),
                    outcome = "D.3",
                    tv_confounders = c("A.1", "A.2", "A.3"),
                    ti_confounders = "C",
                    type = "update",
                    bal_stats = b,
                    save.out = FALSE)
f <- createFormulas(exposure = "A",
                    exposure_time_pts = c(1, 2, 3),
                    outcome = "D.3",
                    tv_confounders = c("A.1", "A.2", "A.3", "B.1", "B.2", "B.3"),
                    ti_confounders = "C",
                    type = "update",
                    bal_stats = b,
                    save.out = FALSE)

}
