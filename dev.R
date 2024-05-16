library(devtools)
devtools::load_all()

# TODO: Discuss how things should be saved?
# TODO: make sample_data.png a markdown table 

data <- data.frame(
  A.1 = rnorm(n = 50),
  A.2 = rnorm(n = 50),
  A.3 = rnorm(n = 50),
  A.14 = rnorm(n = 50),
  B.1 = rnorm(n = 50),
  B.2 = rnorm(n = 50),
  B.3 = rnorm(n = 50),
  L.1 = sample(c(0,1), size = 50, replace = TRUE),
  C   = rnorm(n = 50),
  F   = sample(c("a", "b"), size = 50, replace = TRUE)
)
data$Y.14 <- data$A.1 * 0.2 + data$A.2 * 0.1 + data$A.3 * 0.05 + rnorm(nrow(data))

# Possible TODO: sweights = ?, (passed to balance functions, WeightIt, and outcome model) (low priority for now)
# Possible TODO: Give basename "A.*"
# Without epochs
# obj <- initMSM(
#   data, 
#   exposure = c("A.1", "A.2", "A.3"),
#   tv_conf = c("B.1", "B.2", "B.3"),
#   ti_conf = c("C", "F"),
#   concur_conf = "B.1"
# )

# With epochs
obj <- initMSM(
  data, 
  exposure = c("A.1", "A.2", "A.3", "A.14"),
  epoch = c("Infant", "Infant", "Toddlerhood", "Teenager"),
  tv_conf = c("B.1", "B.2", "B.3"),
  ti_conf = c("C", "F"),
  concur_conf = "B.1"
)

# Temp act like mids
data <- lapply(1:2, \(i) data)

# Formula ----
f <- createFormulas(obj = obj, type = "full")
print(f)
f <- createFormulas(obj = obj, type = "short")
print(f)

# Weights ----
w <- createWeights(data = data, obj = obj, formulas = f, method = "cbps")
print(w)
w <- createWeights(data = data, obj = obj, formulas = f, method = "glm")
print(w, i = 1)
print(w, i = TRUE)

# overwrite default options of `WeightIt::weightitMSM`
w <- createWeights(
  data = data, obj = obj, formulas = f, 
  method = "bart", stabilize = FALSE
)
plot(w, i = 1)
plot(w, i = 2)

# Trim weights ----
t <- trimWeights(w, at = 0.975, lower = FALSE)
print(t, i = 1)
print(t, i = TRUE)


# Assess balance ----
# Current confounders argument
b <- assessBalance(data = data, obj = obj)

print(b)
print(b, i = 1)
print(b, i = 1, t = 3)
print(b, i = 1, t = TRUE)
print(b, i = 1, t = c("A.1", "A.2", "A.3"))
summary(b)
plot(b, t = TRUE)
plot(b, t = 3)
plot(b, t = c(1, 2, 3))
plot(b, t = c("A.1", "A.2"))


bw <- assessBalance(data = data, obj = obj, weights = w)
bw
print(bw)
plot(bw)
plot(bw, t = 1)


# Update formula ----
f = createFormulas(obj = obj, bal_stats = b, type = "update")
print(f)

# Custom formula ----
custom = list(
  A.1 ~ C + B.1 + B.1:C, 
  A.2 ~ C + A.1, 
  A.3 ~ C + A.2,
  A.14 ~ C + A.2 + A.1 + B.1 + B.1:C
)
try({ fc = createFormulas(obj = obj, custom = custom) })


# Fit Model ----
fit_m0 <- fitModel(
  data = data, obj = obj, weights = w, 
  outcome = "Y.14", model = "m0"
)
print(fit_m0)
(fit_m1 <- fitModel(
  data = data, obj = obj, weights = w, 
  outcome = "Y.14", model = "m1", 
  covariates = c("C")
))
(fit_m2 <- fitModel(
  data = data, obj = obj, weights = w, 
  outcome = "Y.14", model = "m2", 
  int_order = 2
))
(fit_m3 <- fitModel(
  data = data, obj = obj, weights = w, 
  outcome = "Y.14", model = "m3",
  int_order = 2, covariates = c("C")
))


# Compare histories ----
comp = compareHistories(
  obj, fit = fit_m0,
  hi_lo_cut = c(0.3, 0.6)
)
print(comp)
plot(comp)
summary(comp, "preds")
summary(comp, "comps")

comp2 = compareHistories(
  obj, fit = fit_m0, 
  reference = "l-l-l",
  comparison = c("h-h-h", "h-h-l") 
) 
print(comp2)
plot(comp2)
summary(comp2, "preds")
summary(comp2, "comps")

