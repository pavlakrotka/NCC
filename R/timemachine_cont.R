#' Time machine for continuous data
#'
#' @description Time machine for continuous data
#'
#' @param data Simulated trial data, e.g. result from the `datasim_cont()` function
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1)
#' @param alpha Type I error. Default=0.025
#' @param prec_delta ...
#' @param prec_gamma ...
#' @param tau_a ...
#' @param tau_b ...
#' @param prec_a ...
#' @param prec_b ...
#' @param bucket_size Number of patients per time bucket. Default=25
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
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval and an indicator whether the null hypothesis was rejected or not for the investigated treatment
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
                             bucket_size = 25){

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

  data$time_bucket <- rep(c(1:ceiling((nrow(data)/bucket_size))), each=bucket_size)[1:nrow(data)]

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

  ### Extract posterior samples
  mcmc_samples <- coda.samples(jags_model,
                               c("delta"),
                               n.iter = 4000)

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



