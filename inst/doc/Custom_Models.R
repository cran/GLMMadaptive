## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library("GLMMadaptive")

## ------------------------------------------------------------------------
set.seed(101)
dd <- expand.grid(f1 = factor(1:3), f2 = LETTERS[1:2], g = 1:30, rep = 1:15,
                  KEEP.OUT.ATTRS = FALSE)
mu <- 5*(-4 + with(dd, as.integer(f1) + 4 * as.numeric(f2)))
dd$y <- rnbinom(nrow(dd), mu = mu, size = 0.5)

## ------------------------------------------------------------------------
my_negBinom <- function (link = "log") {
    stats <- make.link(link)
    log_dens <- function (y, eta, mu_fun, phis) {
        # the log density function
        phis <- exp(phis) # unconstrained optimization of 'phis'
        mu <- mu_fun(eta)
        log_mu_phis <- log(mu + phis)
        comp1 <- lgamma(y + phis) - lgamma(phis) - lgamma(y + 1)
        comp2 <- phis * log(phis) - phis * log_mu_phis
        comp3 <- y * log(mu) - y * log_mu_phis
        out <- comp1 + comp2 + comp3
        attr(out, "mu_y") <- mu
        out
    }
    structure(list(family = "user Neg Binom", link = stats$name, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, log_dens = log_dens),
              class = "family")
}

## ------------------------------------------------------------------------
system.time(
    fm1 <-  mixed_model(fixed = y ~ f1 * f2, random = ~ 1 | g, data = dd, 
                        family = my_negBinom(), n_phis = 1, 
                        initial_values = list("betas" = poisson()))
)

fm1

## ------------------------------------------------------------------------
my_negBinom2 <- function (link = "log") {
    stats <- make.link(link)
    log_dens <- function (y, eta, mu_fun, phis) {
        # the log density function
        phis <- exp(phis)
        mu <- mu_fun(eta)
        log_mu_phis <- log(mu + phis)
        comp1 <- lgamma(y + phis) - lgamma(phis) - lgamma(y + 1)
        comp2 <- phis * log(phis) - phis * log_mu_phis
        comp3 <- y * log(mu) - y * log_mu_phis
        out <- comp1 + comp2 + comp3
        attr(out, "mu_y") <- mu
        out
    }
    score_eta_fun <- function (y, mu, phis) {
        # the derivative of the log density w.r.t. mu
        phis <- exp(phis)
        mu_phis <- mu + phis
        comp2 <- - phis / mu_phis
        comp3 <- y / mu - y / mu_phis
        # the derivative of mu w.r.t. eta (this depends on the chosen link function)
        mu.eta <- mu
        (comp2 + comp3) * mu.eta
    }
    score_phis_fun <- function (y, mu, phis) {
        # the derivative of the log density w.r.t. phis
        phis <- exp(phis)
        mu_phis <- mu + phis
        comp1 <- digamma(y + phis) - digamma(phis)
        comp2 <- log(phis) + 1 - log(mu_phis) - phis / mu_phis
        comp3 <- - y / mu_phis
        comp1 + comp2 + comp3
    }
    structure(list(family = "user Neg Binom", link = stats$name, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, log_dens = log_dens,
                   score_eta_fun = score_eta_fun, score_phis_fun = score_phis_fun),
              class = "family")
}

## ------------------------------------------------------------------------
system.time(
    fm2 <-  mixed_model(fixed = y ~ f1 * f2, random = ~ 1 | g, data = dd, 
                        family = my_negBinom2(), n_phis = 1, 
                        initial_values = list("betas" = poisson()))
)

fm2

## ------------------------------------------------------------------------
set.seed(1234)
n <- 100 # number of subjects
K <- 8 # number of measurements per subject
t_max <- 15 # maximum follow-up time

# we constuct a data frame with the design: 
# everyone has a baseline measurment, and then measurements at random follow-up times
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects
X <- model.matrix(~ sex * time, data = DF)
Z <- model.matrix(~ time, data = DF)

betas <- c(20.1, -0.5, 0.24, -0.05) # fixed effects coefficients
sigma <- 1.2 # scale parameter for the Student's t errors
D11 <- 1.5 # variance of random intercepts
D22 <- 1.2 # variance of random slopes

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
# linear predictor
eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
# we simulate binary longitudinal data
DF$y <- eta_y + sigma * rt(n * K, df = 4)

## ------------------------------------------------------------------------
students.t <- function (df = stop("'df' must be specified"), link = "identity") {
    .df <- df
    env <- new.env(parent = .GlobalEnv)
    assign(".df", df, envir = env)
    stats <- make.link(link)
    log_dens <- function (y, eta, mu_fun, phis) {
        # the log density function
        sigma <- exp(phis)
        out <- dt(x = (y - eta) / sigma, df = .df, log = TRUE) - log(sigma)
        attr(out, "mu_y") <- eta
        out
    }
    score_eta_fun <- function (y, mu, phis) {
        # the derivative of the log density w.r.t. mu
        sigma2 <- exp(phis)^2
        y_mu <- y - mu
        (y_mu * (.df + 1) / (.df * sigma2)) / (1 + y_mu^2 / (.df * sigma2))
    }
    score_phis_fun <- function (y, mu, phis) {
        sigma <- exp(phis)
        y_mu2_df <- (y - mu)^2 / .df
        (.df + 1) * y_mu2_df * sigma^{-2} / (1 + y_mu2_df / sigma^2) - 1
    }
    environment(log_dens) <- environment(score_eta_fun) <- environment(score_phis_fun) <- env
    structure(list(family = "Student's-t", link = stats$name, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, log_dens = log_dens, 
                   score_eta_fun = score_eta_fun),
              class = "family")
}

## ------------------------------------------------------------------------
fm <-  mixed_model(y ~ sex * time, random = ~ time | id, data = DF, 
                   family = students.t(4), n_phis = 1, 
                   initial_values = list("betas" = gaussian()))

fm

