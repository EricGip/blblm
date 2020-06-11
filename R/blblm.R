#' @import purrr
#' @import stats
#' @importFrom magrittr %>%
#' @import furrr
#' @import utils
#' @aliases NULL
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("fit"))

#' Linear Regression with Bag of Little Bootstraps
#' @param formula The linear regression model you want
#' @param data The data you want to perform blblm on
#' @param m = number of observations
#' @param B = times to perform boostrap on the data, more times = more accurate; but will take longer.
#' @param parallel set to ` = TRUE` if you want to activate speed up blb if your system can handle it.
#'
#'
#' @references https://link.springer.com/chapter/10.1007/978-3-030-01418-6_53
#' @export

blblm <- function(formula, data, m = 10, B = 5000, parallel = FALSE) {

  data_list <- split_data(data, m)

  estimates <- future_map(
    data_list,

    ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))

  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"

  invisible(res)
}


#' Split data
#' @description split data into m parts of approximated equal sizes
#'
#' @param data Data you want to split
#' @param m numbers of observations
#'
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' Linear regression on each sub sample
#' compute the estimates for each subsample
#'
#' @param formula linear regression model you want
#' @param data dataset you want to perform formula on
#' @param n number of observations
#' @param B times to run bootstrap
#'
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#' Linear Regression on each bootstrap
#' compute the regression estimates for a blb dataset
#'
#' @param formula model that was specified
#' @param data dataset to perform lm on
#' @param n number of observations
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}


#' linear regression estimate
#' estimate the regression estimates based on given the number of repetitions
#'
#' @param data dataset to perform linear regression on
#' @param formula linear regression formula to perform
#' @param freqs number of repetitions, used in calcuationg regression estimates.
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

# cpp
# fastLmPure <- function(X, y) {
#   stopifnot(is.matrix(X), is.numeric(y), nrow(y) == nrow(X))
# }
#
#
# fastLm <- function(X, ...) UseMethod("fastLm") {
#   X <- as.matrix(X)
#   y <- as.numeric(y)
#
#   res <-
# }

#

#' bag of little boostraps coefficient value
#' compute the coefficients from fit
#'
#' @param fit fitted values, use to get residuals/coefficients
blbcoef <- function(fit) {
  coef(fit)
}


#' Sigma value computed from bag of litte bootstraps
#' ompute sigma from fit
#'
#' @param fit fitted values, use to get residuals/coefficients
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' Used to return coefficient estimates of our blb linear regression.
#'
#' @param x blblm model formula
#' @param ... any other argument you might want to pass through here.
#'
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' Sigma value estimate from blblm
#'
#' @param object the object/model we want to estimate
#' @param confidence set to true for confidence intervals, set to just estimate by default.
#' @param level level of significance
#' @param ... any other argument you want to pass
#'
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' Confidence interval of our blblm
#'
#' @param object object / model that we want to test
#' @param parm Parameters of regression variables we want to see. this is set to null by default and will use all variables. Set variables by using c("variable1", "variable2")
#' @param level Confidence level we want to test
#' @param ... Any other argument you might want to pass.
#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(fit$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' Model predictions on our blblm
#'
#' @param object object / model we want to test
#' @param new_data New data that we want to test our regression on
#' @param confidence Set to `TRUE` if you want a confidence
#' interval prediction.
#' @param level Confidence level we want to investigate
#' @param ... Accepts any other argument you might want to pass.
#' @references http://www.sthda.com/english/articles/40-regression-analysis/166-predict-in-r-model-predictions-and-confidence-intervals/
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (future_map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  future_map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  future_map(.x, .f, ...) %>% reduce(rbind)
}

