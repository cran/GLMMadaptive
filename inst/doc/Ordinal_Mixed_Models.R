## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library("GLMMadaptive")
library("lattice")

## ---- sim_data, eval = TRUE---------------------------------------------------
set.seed(1234)
n <- 300 # number of subjects
K <- 8 # number of measurements per subject
t_max <- 15 # maximum follow-up time

# we construct a data frame with the design: 
# everyone has a baseline measurement, and then measurements at random follow-up times
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects
# we exclude the intercept from the design matrix of the fixed effects because in the
# CR model we have K intercepts (the alpha_k coefficients in the formulation above)
X <- model.matrix(~ sex * time, data = DF)[, -1]
Z <- model.matrix(~ time, data = DF)

thrs <- c(-1.5, 0, 0.9) # thresholds for the different ordinal categories
betas <- c(-0.25, 0.24, -0.05) # fixed effects coefficients
D11 <- 0.48 # variance of random intercepts
D22 <- 0.1 # variance of random slopes

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
# linear predictor
eta_y <- drop(X %*% betas + rowSums(Z * b[DF$id, , drop = FALSE]))
# linear predictor for each category under forward CR formulation
# for the backward formulation, check the note below
eta_y <- outer(eta_y, thrs, "+")
# marginal probabilities per category
mprobs <- cr_marg_probs(eta_y)
# we simulate ordinal longitudinal data
DF$y <- unname(apply(mprobs, 1, sample, x = ncol(mprobs), size = 1, replace = TRUE))
DF$y <- factor(DF$y, levels = 1:4, labels = c("none", "mild", "moderate", "severe"))

## ---- set_up_data, eval = TRUE------------------------------------------------
cr_vals <- cr_setup(DF$y)
cr_data <- DF[cr_vals$subs, ]
cr_data$y_new <- cr_vals$y
cr_data$cohort <- cr_vals$cohort

## ---- random_intercepts, eval = TRUE------------------------------------------
fm <- mixed_model(y_new ~ cohort + sex + time, random = ~ 1 | id, 
                  data = cr_data, family = binomial())

fm

## ---- relax_CR_assumption, eval = TRUE----------------------------------------
gm <- mixed_model(y_new ~ cohort * sex + time, random = ~ 1 | id, 
                  data = cr_data, family = binomial())

gm

## ---- LRT_CR_assumption, eval = TRUE------------------------------------------
anova(fm, gm)

## ---- effect_plot_data, eval = TRUE-------------------------------------------
nDF <- with(cr_data, expand.grid(cohort = levels(cohort), sex = levels(sex), 
                                 time = seq(0, 10, length.out = 55)))

plot_data <- effectPlotData(fm, nDF)

## ---- CR_probs_plot, eval = TRUE, fig.align = "center", fig.width = 8.5, fig.height = 7.5----
expit <- function (x) exp(x) / (1 + exp(x))
my_panel_bands <- function(x, y, upper, lower, fill, col, subscripts, ..., font, 
                           fontface) {
    upper <- upper[subscripts]
    lower <- lower[subscripts]
    panel.polygon(c(x, rev(x)), c(upper, rev(lower)), col = fill, border = FALSE, ...)
}

xyplot(expit(pred) ~ time | sex, group = cohort, data = plot_data, 
       upper = expit(plot_data$upp), low = expit(plot_data$low), type = "l",
       panel = function (x, y, ...) {
           panel.superpose(x, y, panel.groups = my_panel_bands, ...)
           panel.xyplot(x, y, lwd = 2,  ...)
       }, xlab = "Follow-up time", ylab = "Continuation Ratio Probabilities")

## ---- effect_plot_data_marg, eval = TRUE--------------------------------------
plot_data_m <- effectPlotData(fm, nDF, CR_cohort_varname = "cohort", 
                              direction = "forward")

## ---- CR_probs_plot_marg, eval = TRUE, fig.align = "center", fig.width = 8.5, fig.height = 7.5----
key <- list(space = "top", rep = FALSE,
            text = list(levels(DF$y)[1:2]),
            lines = list(lty = c(1, 1), lwd = c(2, 2), col = c("#0080ff", "#ff00ff")),
            text = list(levels(DF$y)[3:4]),
            lines = list(lty = c(1, 1), lwd = c(2, 2), col = c("darkgreen", "#ff0000")))

xyplot(expit(pred) ~ time | sex, group = ordinal_response, data = plot_data_m, 
       upper = expit(plot_data_m$upp), low = expit(plot_data_m$low), type = "l",
       panel = function (x, y, ...) {
           panel.superpose(x, y, panel.groups = my_panel_bands, ...)
           panel.xyplot(x, y, lwd = 2, ...)
       }, xlab = "Follow-up time", ylab = "Marginal Probabilities", key = key)

## ---- effect_plot_data_marg2, eval = TRUE-------------------------------------
plot_data_m2 <- effectPlotData(fm, nDF, CR_cohort_varname = "cohort", 
                               direction = "forward", marginal = TRUE, cores = 2)

## ---- CR_probs_plot_marg2, eval = TRUE, fig.align = "center", fig.width = 8.5, fig.height = 7.5----
xyplot(expit(pred) ~ time | sex, group = ordinal_response, data = plot_data_m2, 
       upper = expit(plot_data_m2$upp), low = expit(plot_data_m2$low), type = "l",
       panel = function (x, y, ...) {
           panel.superpose(x, y, panel.groups = my_panel_bands, ...)
           panel.xyplot(x, y, lwd = 2,  ...)
       }, xlab = "Follow-up time", 
       ylab = "Marginal Probabilities\nalso w.r.t Random Effects", 
       key = key)

