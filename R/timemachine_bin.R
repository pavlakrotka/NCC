#' Time machine for binary data
#'
#' @description Time machine for binary data
#'
#' @param data Simulated trial data, e.g. result from the `datasim_bin()` function
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1)
#' @param alpha Type I error. Default=0.025
#' @param prec_delta ...
#' @param prec_gamma ...
#' @param tau_a ...
#' @param tau_b ...
#'
#' @importFrom stats aggregate
#' @importFrom stats quantile
#' @importFrom rjags jags.model
#' @importFrom rjags coda.samples
#'
#' @export
#'
#' @examples
#'
#' trial_data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")
#'
#' timemachine_bin(data = trial_data, arm = 3)
#'
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval and an indicator whether the null hypothesis was rejected or not for the investigated treatment
#' @author Dominic Magirr, Peter Jacko

timemachine_bin <- function(data,
                            arm,
                            alpha = 0.025,
                            prec_delta = 0.001,
                            prec_gamma = 0.001,
                            tau_a = 0.1,
                            tau_b = 0.01){

  model_string_time_machine_bin <- "
  model {


    for (i in 1:n_trt_periods){


       x[i] ~ dbin(p[i], n[i])

            ## alpha[1] is current interval, alpha[10] is most distant interval
            ## time j corresponds to alpha interval Int-(j-1)
            ## e.g. time 1 corresponds to alpha interval 10-(1-1)=10 if interim #10
            ## e.g. time 2 corresponds to alpha interval 7-(2-1)=6 if interim #7

       logit(p[i]) = gamma0 + alpha[n_periods- (period[i]-1)] + delta[trt[i]]

    }


    ## Priors for time, counting backwards from current time interval (k=1 for current, k=10 for most distant)
    alpha[1] = 0
    alpha[2] ~ dnorm(0, tau)
    for(k in 3:n_periods) {
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

  }
"



  trt_period_n <- aggregate(data$response,
                            list(treatment = data$treatment,
                                 period = data$period),
                            "length")

  trt_period_x <- aggregate(data$response,
                            list(treatment = data$treatment,
                                 period = data$period),
                            "sum")

  n_trt_periods <- dim(trt_period_n)[1]
  n_trts <- max(trt_period_n$treatment) + 1
  n_periods <- max(trt_period_n$period)

  trt <- trt_period_n$treatment + 1 ## need index to start at 1, not 0.
  period <- trt_period_n$period

  x <- trt_period_x$x
  n <- trt_period_n$x

  ### Arguments to pass to jags_model
  data_list = list(x = x,
                   n = n,
                   trt = trt,
                   period = period,
                   n_trts = n_trts,
                   n_periods = n_periods,
                   n_trt_periods = n_trt_periods,
                   prec_delta = prec_delta,
                   prec_gamma = prec_gamma,
                   tau_a = tau_a,
                   tau_b = tau_b)

  inits_list = list(gamma0 = 0)


  ### Fit the model
  jags_model <- jags.model(textConnection(model_string_time_machine_bin),
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
  post_cis <- apply(all_samples, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2))

  ## posterior "p-values"
  post_p <- apply(all_samples, 2, function(x) mean(x < 0))

  ## return
  list(p_val = unname(post_p[arm]),
       treat_effect = unname(post_means[arm]),
       lower_ci = post_cis[1,arm],
       upper_ci = post_cis[2,arm],
       reject_h0 = unname(post_p[2] < alpha / 2))
}



