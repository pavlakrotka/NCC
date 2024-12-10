#' Analysis for continuous data using the power prior approach
#'
#' @description This function performs analysis of continuous data using the power prior approach. The method allows for incorporating historical data by accounting for its likelihood with a power argument (weight parameter) to control the degree of borrowing.
#'
#' @param data Data frame with trial data, e.g. result from the `datasim_cont()` function. Must contain columns named 'treatment', 'response' and 'period'.
#' @param arm Integer. Index of the treatment arm under study to perform inference on (vector of length 1). This arm is compared to the control group.
#' @param a_0 Double. Power argument used for down-weighting the likelihood of the historical datasets (0 < a_0 < 1). Default=0.9.
#' @param alpha Double. Decision boundary (one-sided). Default=0.025
#' @param opt Integer (1 or 2). If opt==1, all former periods are used as one historical dataset; if opt==2, periods are treated as separate historical datasets. Default=2.
#' @param check Logical. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.
#' @param ... Further arguments passed by wrapper functions when running simulations.
#'
#' @importFrom BayesPPD two.grp.fixed.a0
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom stats rnorm
#' @importFrom stats rgamma
#'
#' @export
#'
#' @examples
#'
#' trial_data <- datasim_cont(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "stepwise")
#'
#' powerprior_cont(data = trial_data, arm = 3)
#'
#'
#' @return List containing the following elements regarding the results of comparing `arm` to control:
#'
#' - `p-val` - posterior probability that the difference in means is less than zero
#' - `treat_effect` - posterior mean of difference in means
#' - `lower_ci` - lower limit of the (1-2*`alpha`)*100% credible interval for difference in means
#' - `upper_ci` - upper limit of the (1-2*`alpha`)*100% credible interval for difference in means
#' - `reject_h0` - indicator of whether the null hypothesis was rejected or not (`p_val` < `alpha`)
#'
#' @author Pavla Krotka
#'

powerprior_cont <- function(data, 
                            arm, 
                            alpha = 0.025, 
                            a_0 = 0.9, 
                            opt = 2,
                            check = TRUE, ...){
  
  if (check) {
    if (!is.data.frame(data) | sum(c("treatment", "response", "period") %in% colnames(data))!=3) {
      stop("The data frame with trial data must contain the columns 'treatment', 'response' and 'period'!")
    }
    
    if(!is.numeric(arm) | length(arm)!=1){
      stop("The evaluated treatment arm (`arm`) must be one number!")
    }
    
    if(!is.numeric(alpha) | length(alpha)!=1 | alpha>=1 | alpha<=0){
      stop("The significance level (`alpha`) must be one number between 0 and 1!")
    }
    
    if(!is.numeric(a_0) | length(a_0)!=1 | a_0 < 0 | a_0 > 1){
      stop("The  power argument (`a_0`) must be one number between 0 and 1!")
    }
    
    if(!is.numeric(opt)| length(opt)!=1 | opt %in% c(1,2) == FALSE){
      stop("The parameter `opt` must be either 1 or 2!")
    }
  }
  
  min_period <- min(data[data$treatment==arm,]$period)
  max_period <- max(data[data$treatment==arm,]$period)
  
  if (min_period==1) {
    stop("There are no NCC data available for the evaluated treatment arm!")
  }
  
  data_CC <- data[data$period %in% c(min_period:max_period) & data$treatment %in% c(0,arm),]
  
  data_NCC <- data[data$period %in% c(1:(min_period-1)) & data$treatment==0,]
  
  
  y_c <- sum(data_CC[data_CC$treatment==0,]$response)
  n_c <- nrow(data_CC[data_CC$treatment==0,])
  v_c <- var(data_CC[data_CC$treatment==0,]$response)
  
  get_hist_info <- function(hist_data){
    c(sum = sum(hist_data), ss = length(hist_data), var = var(hist_data))
  }
  
  
  if (opt==1) { # treat all former periods as one data source
    historical <- matrix(c(get_hist_info(data_NCC$response), a_0), nrow = 1)
    
  }
  
  if (opt==2) { # treat former periods as separate data sources
    historical <- matrix(0, ncol=4, nrow=length(unique(data_NCC$period)))
    
    for (i in 1:max(data_NCC$period)) {
      historical[i,] <- c(get_hist_info(data_NCC[data_NCC$period==i,]$response), a_0)
    }
  }
  
  # Posterior samples for control
  
  posterior_draws <- two.grp.fixed.a0(
    data.type = "Normal",
    y.c = y_c,
    n.c = n_c,
    v.c = v_c,
    historical = historical,
    nMC = 100000,
    nBI = 250
  )
  
  # summary(posterior_draws)
  
  mu_samples_control <- posterior_draws$posterior.params$`posterior samples of mu_c`
  
  
  # Posterior samples treatment
  
  x <- data_CC[data_CC$treatment==arm,]$response # Data
  n <- length(x)
  x_bar <- mean(x)
  
  # Prior parameters
  mu_0 <- 0           # Prior mean
  kappa_0 <- 1/1000   # Prior precision (1 / variance of prior mean)
  alpha_0 <- 0.001    # Shape parameter for prior sigma^2
  beta_0 <- 0.001     # Scale parameter for prior sigma^2
  
  # Posterior parameters
  kappa_n <- kappa_0 + n
  mu_n <- (kappa_0 * mu_0 + n * x_bar) / kappa_n
  alpha_n <- alpha_0 + n / 2
  beta_n <- beta_0 + 0.5 * sum((x - x_bar)^2) + 
    (kappa_0 * n * (x_bar - mu_0)^2) / (2 * kappa_n)
  
  # Draw samples from posterior
  n_samples <- 100000
  sigma2_samples <- 1 / rgamma(n_samples, shape = alpha_n, rate = beta_n) # Draw sigma^2
  mu_samples_trt <- rnorm(n_samples, mean = mu_n, sd = sqrt(sigma2_samples / kappa_n)) # Draw mu
  
  post_diff <- mu_samples_trt - mu_samples_control
  
  # prob_superiority <- mean(mu_samples_trt > mu_samples_control)
  
  prob_superiority <- mean(post_diff > 0)
  
  reject_h0 <- prob_superiority > 1-alpha
  
  return(list(p_val = 1-prob_superiority,
              treat_effect = mean(post_diff),
              lower_ci = unname(quantile(post_diff, alpha)),
              upper_ci = unname(quantile(post_diff, 1-alpha)),
              reject_h0 = reject_h0))
  
}


