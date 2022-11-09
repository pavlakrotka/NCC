#'Remarks'
#'
#' Check input parmaters for prior specification and whether to include additional ones like robistify yes, no, specifications for robustification in the future...
#' So far: Posterior odds estimated by drawing random samples from the posterior distributions of the incidence rates of the two groups seperately
#'
#'
#' @description Borrows data from non-concurrent controls for binary endpoints for the concurrent controls arm by the use of a MAP prior
#'
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1)
#' @param data Simulated trial data, e.g. result from the `datasim_bin()` function
#' @param alpha Type I error. Default=0.025
#' @param prior_prec_tau dispersion parameter of the half normal prior, the prior for the between study heteroginity
#' @param opt binary, opt=1: all former periods are use as one source, opt=2: periods form different sources get seperately included into the final analysis
#' @param n.samples defines the number of how many random samples will get drwan for the calcualtion of the posterior mean and the CIs, default is 1000
#' @param n.chains the number of parallel chains for the rjags model
#' @param n.inter needed for coda.smaples, number of iterations to monitor of the jags.model
#' @param n.adapt needed for jags.model, defines the number of iterations for adaptation, an initial sampling phase during which the samplers adapt their behaviour to maximize their efficiency
#' @param weight weight given to the non-informative component (0 < weight < 1) for the robustification of the MAP prior according to Schmidli 2014
#' @param robustify Booelan.
#'
#' @importFrom RBesT automixfit
#' @importFrom RBesT postmix
#' @importFrom RBesT mixbeta
#' @importFrom RBesT robustify
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
#' MAP_rjags_bin(data = trial_data, arm = 3)
#'
#'
#' @return containing the p-value (one-sided), estimated treatment effect, 95% confidence interval and an indicator whether the null hypothesis was rejected or not for the investigated treatment
#' @author Katharina Hees



MAP_rjags_bin <- function(data,
                          arm,
                          alpha = 0.025,
                          opt = 2,
                          prior_prec_tau = 4,
                          n.samples = 1000,
                          n.chains = 4,
                          n.iter = 4000,
                          n.adapt = 1000,
                          robustify = TRUE,
                          weight = 0.1, ...){

  if (opt %in% c(1,2) == FALSE) stop("Wrong parameter input. Parameter opt has to be 1 or 2")

  # 1-alpha is here the decision boundary for the posterior probability,
  # in case that one uses a non-informative/flat prior instead of a MAP prior, the type one error is equal to alpha

  #Data preparation
  ##count number of patients for each treatment in each period
  tab_count <- table(data$treatment,data$period)

  ## count number of groups and number of periods
  number_of_groups <- dim(table(data$treatment,data$period))[1] # number of groups incl control
  number_of_periods <- dim(table(data$treatment,data$period))[2] #total number of periods


  ## get start and end period of each treatment
  treatment_start_period <- numeric(number_of_groups)
  treatment_end_period <- numeric(number_of_groups)
  for (i in 1:number_of_groups){
    treatment_start_period[i] <- min(which(table(data$treatment,data$period)[i,] > 0))
    treatment_end_period[i] <- max(which(table(data$treatment,data$period)[i,] > 0))
  }

  ## get concurrent and non-concurrent controls of treatment = arm
  cc <- data[data$treatment == 0 & data$period > treatment_start_period[arm + 1] - 1 & data$period < treatment_end_period[arm + 1] + 1, ]
  ncc <- data[data$treatment == 0 & data$period < treatment_start_period[arm + 1],]


  ## get patients treated with treatment=arm
  t_treatment <- data[data$treatment == arm,]

  if(dim(ncc)[1] !=0){

    if(opt==1){
      ncc$period <- 0
    }

    ## data for the model underlying for control data
    ncc_control_data_jags <- list(
      N_std = length(unique(ncc$period)),
      y = sapply(unique(ncc$period), function(x) sum(ncc[ncc$period == x, ]$response)),
      n = sapply(unique(ncc$period), function(x) length(ncc[ncc$period == x, ]$response)),
      prior_prec_tau = prior_prec_tau
    )



    model_text <- "model
      {
        for (i in 1:N_std) {
          # binomial likelihood
          y[i] ~ dbin(p[i], n[i])
          logit(p[i])  <-  logit_p[i]
          logit_p[i] ~ dnorm(mu, inv_tau2)
        }

        # prior distributions
        inv_tau2 <- pow(tau, -2)
        tau ~ dnorm(0, prior_prec_tau)I(0,) # HN(scale = 1 / sqrt(prior_prec_tau)), I(0,) means censoring on positive values
        mu  ~ dnorm(0, 0.001)

        # prediction for p in a new study
        logit_p.pred ~ dnorm(mu, inv_tau2)
        p.pred <- 1 / (1 + exp(-logit_p.pred))
      }"




    fit <- jags.model(file = textConnection(model_text),
                      data = ncc_control_data_jags,
                      #inits = inits,
                      n.chains = n.chains,
                      n.adapt = n.adapt,
                      quiet = TRUE
    )


    # Draw samples from the above fitted MCMC model
    help_mcmc_samples <- coda.samples(fit, "p.pred", n.iter = n.iter)
    help_samples <- do.call(rbind.data.frame, help_mcmc_samples)[,1]


    ## Fit Beta mixture to MCMC samples
    prior_control <- automixfit(help_samples,type="beta",Nc=3) # Nc=3: fixed number of mixture components, to speed the code up

    if (robustify== TRUE) {
      prior_control<-robustify(prior_control,weight=weight)
    }
  }
  else{

    prior_control<-mixbeta(c(1,0.5,0.5))

  }


  ## get summary data of treatment group and concurrent controls
  r.act <- sum(t_treatment$response)
  n.act <- length(t_treatment$response)

  r.pbo <- sum(cc$response)
  n.pbo <- length(cc$response)


  post_control <- postmix(prior_control,n=n.pbo,r=r.pbo)

  post_treatment <- mixbeta(c(1,0.5+r.act,0.5+n.act-r.act))



  samples_control <-rmix(post_control,n=n.samples)
  samples_treat <- rmix(post_treatment,n=n.samples)
  p1 <- samples_treat
  p2 <- samples_control

  ## calculate oddsratio and with that the estimated treatment effect
  oddsratio <- (p1 / (1 - p1)) / (p2 / (1 - p2))
  treat_effect <- mean(oddsratio)

  delta_ci <- quantile(oddsratio, probs = c(alpha, 1 - alpha))
  lower_ci <- as.numeric(delta_ci[1])
  upper_ci <- as.numeric(delta_ci[2])

  ## calculate the "p-value"
  p_value <- mean(samples_treat - samples_control < 0)
  reject_h0 <- ifelse(p_value<alpha,TRUE,FALSE)

  return(list(p_val = p_value,
              treat_effect = treat_effect,
              lower_ci = lower_ci,
              upper_ci = upper_ci,
              reject_h0 = reject_h0
  ))
}




