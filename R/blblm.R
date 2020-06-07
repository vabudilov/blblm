#' @import purrr
#' @import stats
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"
install.packages("furrr", repos = "http://cran.us.r-project.org")
library(furrr)

## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#'
#' generalized linear model
#'
#' @param formula formula
#' @param data data to be used
#' @param m split count
#' @param Bbootstrap depth
#' @param ifgeneral TRUE is generalized model to be used
#' @param family GML family parameter
#' @param wrkrs number of workers
#'
#'
#' @export
blblm <- function(formula, data, m = 10, B = 5000, ifgeneral=FALSE,family = "gaussian",parallel=FALSE, wrkrs=6) {
  if(parallel){suppressWarnings(plan(multiprocess, workers = wrkrs))}
    data_list <- split_data(data, m)
    if(parallel){estimates <- future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B, ifgeneral=ifgeneral, family=family))}
    else {estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B, ifgeneral=ifgeneral, family=family))
    }
    res <- list(estimates = estimates, formula = formula)
    class(res) <- "blblm"
    invisible(res)
    return(res)

}

#' generalized linear model
#'
#' @param formula formula (old)
#' @param data data to be used
#' @param m split count
#' @param Bbootstrap depth
#' @param ifgeneral TRUE is generalized model to be used
#' @param family GML family parameter
#' @param wrkrs number of workers
#'
#' @export
blblm1 <- function(formula, data, m = 10, B = 5000, ifgeneral=FALSE,family = "gaussian",wrkrs=6) {
  suppressWarnings(plan(multiprocess, workers = wrkrs))
  data_list <- split_data(data, m)
  estimates <- future_map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B, ifgeneral=ifgeneral, family=family))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}



#' split data into m parts of approximated equal sizes
#'
#' @param data data to be used
#' @param m splits count
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
#'
#' @param formula formula to use
#' @param data data t0 use
#' @param n number to pass (splits)
#' @param B boostrap "depth"
#' @param ifgeneral if GML to be used (defualt is FALSE)
#' @param family if GML is used, family could be splecified
lm_each_subsample <- function(formula, data, n, B, ifgeneral,family) {
  replicate(B, lm_each_boot(formula, data, n, ifgeneral,family), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
#'
#' @param formula formula to be used
#' @param data data to be used
#' @param n parameter of splits
#' @param ifgeneral if GML to be used (defulat is FALSE)
#' @param family family for GML
lm_each_boot <- function(formula, data, n, ifgeneral,family) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs,ifgeneral,family)
}


#' estimate the regression estimates based on given the number of repetitions
#'
#' @param formula formula to be used
#' @param data data to be used
#' @param freqs weigths
#' @param ifgeneral if GML to be used (default is FALSE)
#' @param family famile for GML model
lm1 <- function(formula, data, freqs,ifgeneral,family) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  if (ifgeneral)fit <- glm(formula, data,family = family, weights = freqs)
      else fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}



#' compute the coefficients from fit
#'
#' @param fit
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#'
#' @param fit
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @param x
#'
#' @param ...
#'
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' @param object
#'
#' @param confidence
#' @param level
#' @param ...
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

#' @param object
#'
#' @param ...
#'
#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(fit$formula), "term.labels")
  }
  alpha <- 1 - 0.95
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


#' Title
#'
#' @param x
#' @param level
#'
#' @return
#' @export
#'
#' @examples
mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

#' Title
#'
#' @param .x
#' @param .f
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

#' Title
#'
#' @param .x
#' @param .f
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

#' Title
#'
#' @param .x
#' @param .f
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
