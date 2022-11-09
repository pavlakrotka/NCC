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
#' @param opt binary, opt=1: all former periods are used as one source, opt=2: periods form different sources get seperately included into the final analysis
#' @param n.samples defines the number of how many random samples will get drwan for the calcualtion of the posterior mean and the CIs, default is 1000
#' @param n.chains the number of parallel chains for the rjags model
#' @param n.inter needed for coda.smaples, number of iterations to monitor of the jags.model
#' @param n.adapt needed for jags.model, defines the number of iterations for adaptation, an initial sampling phase during which the samplers adapt their behaviour to maximize their efficiency
#' @param weight weight given to the non-informative component (0 < weight < 1) for the robustification of the MAP prior according to Schmidli 2014
#' @param robustify Boolean.
#'
#' @importFrom RBesT automixfit
#' @importFrom RBesT postmix
#' @importFrom RBesT mixnorm
#' @importFrom RBesT pmixdiff
#' @importFrom RBesT rmixdiff
#' @importFrom RBesT robustify
#' @importFrom stats quantile
#' @importFrom rjags jags.model
#' @importFrom rjags coda.samples
#'
#' @export
#'
#' @examples
#'
#' trial_data<-datasim_cont(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "stepwise")
#'
#' MAP_rjags_cont(data = trial_data, arm = 3)
#'
#'
#' @return containing the p-value (one-sided), estimated treatment effect, 95% confidence interval and an indicator whether the null hypothesis was rejected or not for the investigated treatment
#' @author Katharina Hees


MAP_rjags_cont <- function(data,
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
  ncc <- data[data$treatment == 0 & data$period < treatment_start_period[arm + 1], ]



  ## get patients treated with treatment=arm
  t_treatment <- data[data$treatment == arm,]

  if(dim(ncc)[1] !=0){

    if(opt==1){
      ncc$period <- 0
    }

    ## data for the model underlying for control data
    ncc_control_data_jags <- list(
      N_std = length(unique(ncc$period)),
      y = sapply(unique(ncc$period), function(x) mean(ncc[ncc$period == x, ]$response)),
      se = sapply(unique(ncc$period), function(x) sd(ncc[ncc$period == x, ]$response)/sqrt(length(ncc[ncc$period == x, ]$response))),
      prior_prec_tau = prior_prec_tau
    )



    model_text <- "model
      {
        for (i in 1:N_std) {
          y[i] ~ dnorm(theta[i],se[i])
          theta[i] ~ dnorm(mu,inv_tau2)
        }

        inv_tau2 <- pow(tau, -2)
        tau ~ dnorm(0, prior_prec_tau)I(0,) # HN(scale = 1 / sqrt(prior_prec_tau)), I(0,) means censoring on positive values
        mu  ~ dnorm(0, 0.001)

        # prediction for theta in a new study
        theta.pred ~ dnorm(mu, inv_tau2)
      }"




    fit <- jags.model(file = textConnection(model_text),
                      data = ncc_control_data_jags,
                      #inits = inits,
                      n.chains = n.chains,
                      n.adapt = n.adapt,
                      quiet = TRUE
    )


    # Draw samples from the above fitted MCMC model
    help_mcmc_samples <- coda.samples(fit, "theta.pred ", n.iter = n.iter)
    help_samples <- do.call(rbind.data.frame, help_mcmc_samples)[,1]


    ## Fit Beta mixture to MCMC samples
    prior_control <- automixfit(help_samples,type="norm",Nc=3) # Nc=3: fixed number of mixture components, to speed the code up


    if (robustify== TRUE) {
      prior_control<-robustify(prior_control,weight=weight,mean=0,sigma=2)
    }
  } else {

    prior_control<- mixnorm(c(1,0,1000), param = 'ms') # creates weak prior for the control group, in case there are no NCCs
  }

  #Create prior for treatment group
  weak_prior_treatment <- mixnorm(c(1,0,1000), param = 'ms') # creates weak prior for the treatment group
  #mixnorm creates normal mixture density, here with one component, mean 0 and sd=1000
  # param= "ms" means that parameters are given via mean and standard deviation

  ## get summary data of treatment group and concurrent controls
  y.act <- mean(t_treatment$response)
  n.act <- length(t_treatment$response)
  y.act.se <-sd(t_treatment$response)/sqrt(n.act)


  y.pbo <- mean(cc$response)
  n.pbo <- length(cc$response)
  y.pbo.se <-sd(cc$response)/sqrt(n.pbo)

  ## obtain posterior distributions
  post_act <- postmix(weak_prior_treatment, m=y.act,se=y.act.se)
  #postmix calculates the posterior distribution given a prior (here "weak_prior"),
  #where the prior is a mixture of conjugate distributions.
  #The posterior is then also a mixture of conjugate distributions.

  post_control <- postmix(prior_control,m=y.pbo,se=y.pbo.se)





  ## calculate oddsratio and with that the estimated treatment effect
  # Calculate Treatment effect and confidence interval by simulation of posterior mean
  random_gen<-as.numeric(rmixdiff(post_act,post_control,n=1000))
  #rmixdiff = random number generation for the difference of two mixture distributions
  theta_hat<-mean(random_gen)
  theta_ci<-quantile(random_gen,probs=c(alpha,1-alpha))
  lower_ci <- as.numeric(theta_ci[1])
  upper_ci <- as.numeric(theta_ci[2])



  ## calculate the "p-value"
  p_value <- 1-pmixdiff(post_control,post_act,0)
  reject_h0 <- ifelse(p_value<alpha,TRUE,FALSE)

  return(list(p_val = p_value,
              treat_effect = theta_hat,
              lower_ci = lower_ci,
              upper_ci = upper_ci,
              reject_h0 = reject_h0))
}



