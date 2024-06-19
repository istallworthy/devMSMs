library(devtools)
devtools::load_all()

# TODO: Is it better to return `ggplot` objects from `plot()` and `tinytable` objects from `print()` and tell users how to save?

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
  F   = sample(c("a", "b", "c"), size = 50, replace = TRUE)
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

# If want to simulate `mids` object
data <- lapply(1:2, \(i) data)

# With epochs
obj <- initMSM(
  data, 
  exposure = c("A.1", "A.2", "A.3", "A.14"),
  epoch = c("Infant", "Infant", "Toddlerhood", "Teenager"),
  tv_conf = c("B.1", "B.2", "B.3"),
  ti_conf = c("C", "F"),
  concur_conf = "B.1",
  home_dir = "../temp-out/"
)


# Formula ----
f <- createFormulas(obj = obj, type = "full", save.out = TRUE)
print(f)
f <- createFormulas(obj = obj, type = "short", save.out = TRUE)
print(f)

f <- createFormulas(obj = obj, type = "short", save.out = "test_custom_name.rds")

# Weights ----
w <- createWeights(data = data, obj = obj, formulas = f, method = "cbps", save.out = TRUE)
print(w)
w <- createWeights(data = data, obj = obj, formulas = f, method = "glm")
print(w, i = 1)
print(w, i = TRUE)

# overwrite default options of `WeightIt::weightitMSM`
w <- createWeights(
  data = data, obj = obj, formulas = f, 
  method = "bart", stabilize = FALSE, 
  save.out = TRUE
)
plot(w, i = 1, save.out = TRUE)
plot(w, i = 2, save.out = TRUE)
p <- plot(w, i = TRUE, save.out = TRUE)
plot(p[[1]])

# TODO: Should we document that you can use WeightIt internal functions? Can think of: 
# - plots with `i = TRUE`, (p[[i]] extracts `ggplot2` object)
# - weights (w[[i]] extracts `weightitMSM`)
# - fits (f[[i]] extracts `glm_weightit`)
# - balance statistics (b[[i]] extracts `list` of `data.frame`)
#
# here's an example of WeightIt functions being called:
summary(w[[1]])
plot(summary(w[[1]]))
plot(summary(w[[1]]), time = 1)


# Trim weights ----
t <- trimWeights(obj = obj, weights = w, at = 0.975, lower = FALSE, save.out = TRUE)
print(t, i = 1)
print(t, i = TRUE)
plot(t, save.out = TRUE)


# Assess balance ----
# Current confounders argument
b <- assessBalance(data = data, obj = obj, save.out = TRUE)

print(b)
print(b, i = 1, save.out = TRUE)
print(b, i = 1, t = 3)
print(b, i = 1, t = TRUE)
print(b, i = 1, t = c("A.1", "A.2", "A.3"))
summary(b, save.out = TRUE)
plot(b, t = TRUE, save.out = TRUE)
plot(b, t = 3)
plot(b, t = c(1, 2, 3))
plot(b, t = c("A.1", "A.2"))

# b/c of tinytable, can use `.tex`/`.html`/`.png`/`.pdf`
summary(b, save.out = "summ_b.html") 
summary(b, save.out = "summ_b.tex") 
print(b, save.out = "print_b.html") 
print(b, save.out = "print_b.tex") 
print(b, save.out = "print_b.pdf") 


bw <- assessBalance(data = data, obj = obj, weights = w)
bw
print(bw, save.out = TRUE)
plot(bw, save.out = TRUE)
plot(bw, t = 1, save.out = TRUE)


# Update formula ----
f = createFormulas(obj = obj, bal_stats = b, type = "update", save.out = TRUE)
print(f)

# Custom formula ----
custom = list(
  A.1 ~ C + B.1 + B.1:C, 
  A.2 ~ C + A.1, 
  A.3 ~ C + A.2,
  A.14 ~ C + A.2 + A.1 + B.1 + B.1:C
)
fc = createFormulas(obj = obj, custom = custom)


# Fit Model ----
fit_m0 <- fitModel(
  data = data, obj = obj, weights = w, 
  outcome = "Y.14", model = "m0",
  save.out = TRUE
)
print(fit_m0, i = 2, save.out = TRUE)
print(fit_m0, i = 2, save.out = "fit_m0.html")
print(fit_m0, i = 2, save.out = "fit_m0.tex")

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
  hi_lo_cut = c(0.3, 0.6), 
  save.out = TRUE
)
print(comp, save.out = TRUE)
print(comp, save.out = "comp.md")
plot(comp, save.out = TRUE)
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
