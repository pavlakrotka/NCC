#' Time machine for continuous endpoints.
#'
#' @description Performs analysis of continuous data using the Time machine approach, which uses a second-order Bayesian normal dynamic linear model (NDLM), takes into account all data until the investigated arm left the trial and includes covariate adjustment for time (separating the trial into buckets of pre-defined size) using a hierarchical model that smooths the control response rate over time.
#'
#' @param data Simulated trial data, e.g. result from the `datasim_cont()` function. Must contain columns named 'treatment' and 'response'.
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1).
#' @param alpha Type I error. Default=0.025.
#' @param prec_delta Precision of the prior regarding the treatment effect. Default=0.001.
#' @param prec_gamma Precision of the prior regarding the control response. Default=0.001.
#' @param tau_a Parameter \eqn{a} of the Gamma distribution regarding the precision of the drift parameter \eqn{\tau}. I.e., \eqn{\tau \sim Gamma(a,b)}. Default=0.1.
#' @param tau_b Parameter \eqn{b} of the Gamma distribution regarding the precision of the drift parameter \eqn{\tau}. I.e., \eqn{\tau \sim Gamma(a,b)}. Default=0.01.
#' @param prec_a Parameter \eqn{a} of the Gamma distribution regarding the precision of the responses. I.e., \eqn{\sigma \sim Gamma(a,b)}. Default=0.001.
#' @param prec_b Parameter \eqn{b} of the Gamma distribution regarding the precision of the responses. I.e., \eqn{\sigma \sim Gamma(a,b)}. Default=0.001.
#' @param bucket_size Number of patients per time bucket. Default=25.
#' @param ... Further arguments for simulation function.
#'
#' @importFrom stats aggregate
#' @importFrom stats quantile
#' @importFrom stats var
#' @importFrom rjags jags.model
#' @importFrom rjags coda.samples
#'
#' @export
#'
#' @examples
#'
#' trial_data <- datasim_cont(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "linear")
#'
#' timemachine_cont(data = trial_data, arm = 3)
#'
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval and an indicator whether the null hypothesis was rejected or not for the investigated treatment.
#' @author Dominic Magirr, Peter Jacko

timemachine_cont <- function(data,
                             arm,
                             alpha = 0.025,
                             prec_delta = 0.001,
                             prec_gamma = 0.001,
                             tau_a = 0.1,
                             tau_b = 0.01,
                             prec_a = 0.001,
                             prec_b = 0.001,
                             bucket_size = 25, ...){

  max_period <- max(data[data$treatment==arm,]$period)
  data <- data[data$period %in% c(1:max_period),]

  model_string_time_machine_cont <- "
  model {


    for (i in 1:n_trt_buckets){


       x_bar[i] ~ dnorm(mu[i], n[i] * prec)

       sigma_2[i] ~ dgamma((n[i] - 1) / 2, n[i] * prec / 2)

            ## alpha[1] is current interval, alpha[10] is most distant interval
            ## time j corresponds to alpha interval Int-(j-1)
            ## e.g. time 1 corresponds to alpha interval 10-(1-1)=10 if interim #10
            ## e.g. time 2 corresponds to alpha interval 7-(2-1)=6 if interim #7

       mu[i] = gamma0 + alpha[n_buckets - (time_bucket[i]-1)] + delta[trt[i]]

    }


    ## Priors for time, counting backwards from current time interval (k=1 for current, k=10 for most distant)
    alpha[1] = 0
    alpha[2] ~ dnorm(0, tau)
    for(k in 3:n_buckets) {
        alpha[k] ~ dnorm(2*alpha[k-1] - alpha[k-2], tau)  # 2*(most recent) - (more distant) puts trend in right direction
    }

    ## Other priors
    delta[1] <- 0
    for(i in 2:n_trts){
      delta[i] ~ dnorm(0, prec_delta)
    }

    tau ~  dgamma(tau_a, tau_b)

    ### Intercept
    gamma0 ~ dnorm(0, prec_gamma)

    ### Precision (corresponding to random error term)
    prec ~ dgamma(prec_a, prec_b)


  }
"

  ## repeat this code to take account of awkward
  ## case where n == 1
  model_string_time_machine_cont_1 <- "
  model {


    for (i in 1:n_trt_buckets){


       x_bar[i] ~ dnorm(mu[i], n[i] * prec)

       sigma_2[i] ~ dgamma((n[i] - 1) / 2, n[i] * prec / 2)

            ## alpha[1] is current interval, alpha[10] is most distant interval
            ## time j corresponds to alpha interval Int-(j-1)
            ## e.g. time 1 corresponds to alpha interval 10-(1-1)=10 if interim #10
            ## e.g. time 2 corresponds to alpha interval 7-(2-1)=6 if interim #7

       mu[i] = gamma0 + alpha[n_buckets - (time_bucket[i]-1)] + delta[trt[i]]

    }

    for (i in 1:n_trt_buckets_1){


       x_bar_1[i] ~ dnorm(mu_1[i], prec)

            ## alpha[1] is current interval, alpha[10] is most distant interval
            ## time j corresponds to alpha interval Int-(j-1)
            ## e.g. time 1 corresponds to alpha interval 10-(1-1)=10 if interim #10
            ## e.g. time 2 corresponds to alpha interval 7-(2-1)=6 if interim #7

       mu_1[i] = gamma0 + alpha[n_buckets - (time_bucket_1[i]-1)] + delta[trt_1[i]]

    }

    ## Priors for time, counting backwards from current time interval (k=1 for current, k=10 for most distant)
    alpha[1] = 0
    alpha[2] ~ dnorm(0, tau)
    for(k in 3:n_buckets) {
        alpha[k] ~ dnorm(2*alpha[k-1] - alpha[k-2], tau)  # 2*(most recent) - (more distant) puts trend in right direction
    }

    ## Other priors
    delta[1] <- 0
    for(i in 2:n_trts){
      delta[i] ~ dnorm(0, prec_delta)
    }

    tau ~  dgamma(tau_a, tau_b)

    ### Intercept
    gamma0 ~ dnorm(0, prec_gamma)

    ### Precision (corresponding to random error term)
    prec ~ dgamma(prec_a, prec_b)


  }
"

  data$time_bucket <- rep(c(1:ceiling((nrow(data)/bucket_size))), each=bucket_size)[1:nrow(data)]
  #data$time_bucket <- data$period

  trt_bucket_n <- aggregate(data$response,
                            list(treatment = data$treatment,
                                 time_bucket = data$time_bucket),
                            "length")

  trt_bucket_means <- aggregate(data$response,
                                list(treatment = data$treatment,
                                     time_bucket = data$time_bucket),
                                "mean")

  trt_bucket_vars <- aggregate(data$response,
                               list(treatment = data$treatment,
                                    time_bucket = data$time_bucket),
                               function(x) var(x) * (length(x) - 1) / length(x))


  n_trt_buckets <- dim(trt_bucket_means)[1]
  n_trts <- max(trt_bucket_means$treatment) + 1
  n_buckets <- max(trt_bucket_means$time_bucket)

  trt <- trt_bucket_means$treatment + 1 ## need index to start at 1, not 0.
  time_bucket <- trt_bucket_means$time_bucket

  x_bar <- trt_bucket_means$x
  sigma_2 <- trt_bucket_vars$x
  n <- trt_bucket_n$x

  ## check when n == 1
  size_1 <- which(n == 1)
  if (length(size_1) == 0){ ## don't split up data

    ### Arguments to pass to jags_model
    data_list = list(x_bar = x_bar,
                     sigma_2 = sigma_2,
                     n = n,
                     trt = trt,
                     time_bucket = time_bucket,
                     n_trts = n_trts,
                     n_buckets = n_buckets,
                     n_trt_buckets = n_trt_buckets,
                     prec_delta = prec_delta,
                     prec_gamma = prec_gamma,
                     tau_a = tau_a,
                     tau_b = tau_b,
                     prec_a = prec_a,
                     prec_b = prec_b)

    inits_list = list(gamma0 = 0)


    ### Fit the model
    jags_model <- jags.model(textConnection(model_string_time_machine_cont),
                             data = data_list,
                             inits = inits_list,
                             n.adapt = 1000,
                             n.chains = 3,
                             quiet = T)


  }
  else {

    ## seperate the cases where n == 1
    x_bar_1 <- x_bar[size_1]
    time_bucket_1 <- time_bucket[size_1]
    trt_1 <- trt[size_1]
    n_trt_buckets_1 <- length(x_bar_1)

    ## remove n == 1 cases
    x_bar <- x_bar[-size_1]
    time_bucket <- time_bucket[-size_1]
    trt <- trt[-size_1]
    n <- n[-size_1]
    n_trt_buckets <- length(x_bar)

    ### Arguments to pass to jags_model
    data_list = list(x_bar = x_bar,
                     x_bar_1 = x_bar_1,
                     sigma_2 = sigma_2,
                     n = n,
                     trt = trt,
                     trt_1 = trt_1,
                     time_bucket = time_bucket,
                     time_bucket_1 = time_bucket_1,
                     n_trts = n_trts,
                     n_buckets = n_buckets,
                     n_trt_buckets = n_trt_buckets,
                     n_trt_buckets_1 = n_trt_buckets_1,
                     prec_delta = prec_delta,
                     prec_gamma = prec_gamma,
                     tau_a = tau_a,
                     tau_b = tau_b,
                     prec_a = prec_a,
                     prec_b = prec_b)

    inits_list = list(gamma0 = 0)


    ### Fit the model
    jags_model <- jags.model(textConnection(model_string_time_machine_cont_1),
                             data = data_list,
                             inits = inits_list,
                             n.adapt = 1000,
                             n.chains = 3,
                             quiet = T)
  }
  ### Extract posterior samples
  mcmc_samples <- coda.samples(jags_model,
                               c("delta"),
                               n.iter = 40000)

  ### Arrange all samples together
  all_samples <- do.call(rbind.data.frame, mcmc_samples)

  ## posterior means
  post_means <- colMeans(all_samples)

  ## posterior credible intervals
  post_cis <- apply(all_samples, 2, quantile, probs = c(alpha, 1 - alpha))

  ## posterior "p-values"
  post_p <- apply(all_samples, 2, function(x) mean(x < 0))

  ## return
  list(p_val = unname(post_p[arm+1]),
       treat_effect = unname(post_means[arm+1]),
       lower_ci = post_cis[1,arm+1],
       upper_ci = post_cis[2,arm+1],
       reject_h0 = unname(post_p[arm+1] < alpha))
}



