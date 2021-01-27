## ---- setup, include = FALSE, message = FALSE, warning = FALSE----------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library("GLMMadaptive")
library("DHARMa")

## ---- simulate_data-----------------------------------------------------------
set.seed(123)
n <- 300 # number of subjects
K <- 8 # number of measurements per subject
t_max <- 5 # maximum follow-up time

# we construct a data frame with the design: 
# everyone has a baseline measurement, and then measurements at random follow-up times
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects non-zero part
X <- model.matrix(~ sex * time, data = DF)
Z <- model.matrix(~ 1, data = DF)
# design matrices for the fixed and random effects zero part
X_zi <- model.matrix(~ sex, data = DF)
Z_zi <- model.matrix(~ 1, data = DF)

betas <- c(1.5, 0.05, 0.05, -0.03) # fixed effects coefficients non-zero part
shape <- 2 # shape/size parameter of the negative binomial distribution
gammas <- c(-1.5, 0.5) # fixed effects coefficients zero part
D11 <- 0.5 # variance of random intercepts non-zero part
D22 <- 0.4 # variance of random intercepts zero part

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
# linear predictor non-zero part
eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, 1, drop = FALSE]))
# linear predictor zero part
eta_zi <- as.vector(X_zi %*% gammas + rowSums(Z_zi * b[DF$id, 2, drop = FALSE]))
# we simulate negative binomial longitudinal data
DF$y <- rnbinom(n * K, size = shape, mu = exp(eta_y))
# we set the extra zeros
DF$y[as.logical(rbinom(n * K, size = 1, prob = plogis(eta_zi)))] <- 0

## ---- resids_plot_fun---------------------------------------------------------
resids_plot <- function (object, y, nsim = 1000,
                         type = c("subject_specific", "mean_subject"),
                         integerResponse = NULL) {
    if (!inherits(object, "MixMod"))
        stop("this function works for 'MixMod' objects.\n")
    type <- match.arg(type)
    if (is.null(integerResponse)) {
        integer_families <- c("binomial", "poisson", "negative binomial",
                              "zero-inflated poisson", "zero-inflated negative binomial", 
                              "hurdle poisson", "hurdle negative binomial")
        numeric_families <- c("hurdle log-normal", "beta", "hurdle beta", "Gamma")
        if (object$family$family %in% integer_families) {
            integerResponse <- TRUE
        } else if (object$family$family %in% numeric_families) {
            integerResponse <- FALSE
        } else {
            stop("non build-in family object; you need to specify the 'integerResponse',\n",
                 "\targument indicating whether the outcome variable is integer or not.\n")
        }
    }
    sims <- simulate(object, nsim = nsim, type = type)
    fits <- fitted(object, type = type)
    dharmaRes <- DHARMa::createDHARMa(simulatedResponse = sims, observedResponse = y, 
                              fittedPredictedResponse = fits, 
                              integerResponse = integerResponse)
    DHARMa:::plot.DHARMa(dharmaRes, quantreg = FALSE)
}

## ---- poisson_mixed_model-----------------------------------------------------
fm1 <- mixed_model(y ~ sex * time, random = ~ 1 | id, data = DF,
                   family = poisson())

## ---- poisson_mixed_model_GoF, fig.align = "center", fig.width = 8.5, fig.height = 7.5----
resids_plot(fm1, DF$y)

## ---- zi_poisson_mixed_model--------------------------------------------------
fm2 <- mixed_model(y ~ sex * time, random = ~ 1 | id, data = DF,
                   family = zi.poisson(), 
                   zi_fixed = ~ sex)

## ---- zi_poisson_mixed_model_GoF, fig.align = "center", fig.width = 8.5, fig.height = 7.5----
resids_plot(fm2, DF$y)

## ---- zi_poisson_mixed_model2-------------------------------------------------
fm3 <- mixed_model(y ~ sex * time, random = ~ 1 | id, data = DF,
                   family = zi.poisson(), 
                   zi_fixed = ~ sex, zi_random = ~ 1 | id)

## ---- zi_poisson_mixed_model2_GoF, fig.align = "center", fig.width = 8.5, fig.height = 7.5----
resids_plot(fm3, DF$y)

## ---- zi_NB_mixed_model-------------------------------------------------------
fm4 <- mixed_model(y ~ sex * time, random = ~ 1 | id, data = DF,
                   family = zi.negative.binomial(), 
                   zi_fixed = ~ sex, zi_random = ~ 1 | id)

## ---- zi_NB_mixed_model_GoF, fig.align = "center", fig.width = 8.5, fig.height = 7.5----
resids_plot(fm4, DF$y)

