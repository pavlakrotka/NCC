#' MAP Prior approach analysis for continuous data
#'
#' @description This function performs analysis of continuous data using the MAP Prior approach. The method borrows data from non-concurrent controls to obtain the prior distribution for the control response in the concurrent periods.
#'
#' @param data Simulated trial data, e.g. result from the `datasim_bin()` function. Must contain columns named 'treatment', 'response' and 'period'.
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1). This arm is compared to the control group.
#' @param alpha Type I error rate. Default=0.025
#' @param opt Binary. If opt=1, all former periods are used as one source; if opt=2, periods form different sources get separately included into the final analysis. Default=2.
#' @param prior_prec_tau Dispersion parameter of the half normal prior, the prior for the between study heterogeneity. Default=4.
#' @param n.samples How many random samples will get drawn for the calculation of the posterior mean and the CIs. Default=1000.
#' @param n.chains Number of parallel chains for the rjags model. Default=4.
#' @param n.iter Number of iterations to monitor of the jags.model. Needed for coda.samples. Default=4000.
#' @param n.adapt Number of iterations for adaptation, an initial sampling phase during which the samplers adapt their behavior to maximize their efficiency. Needed for jags.model. Default=1000.
#' @param robustify Booelan.
#' @param weight Weight given to the non-informative component (0 < weight < 1) for the robustification of the MAP prior according to Schmidli (2014). Default=0.1.
#' @param check Boolean. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.
#' @param ... Further arguments for simulation function.
#'
#' @importFrom RBesT automixfit
#' @importFrom RBesT postmix
#' @importFrom RBesT mixnorm
#' @importFrom RBesT pmixdiff
#' @importFrom RBesT rmixdiff
#' @importFrom RBesT robustify
#' @importFrom stats quantile
#' @importFrom stats sd
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
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval, and an indicator whether the null hypothesis was rejected or not (for the investigated treatment specified in the input).
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
                           weight = 0.1,
                           check = TRUE,...){

  if (check) {
    if (!is.data.frame(data) | sum(c("treatment", "response", "period") %in% colnames(data))!=3) {
      stop("The data frame with trial data must contain the columns 'treatment', 'response' and 'period'!")
    }

    if(!is.numeric(arm) | length(arm)!=1){
      stop("The evaluated treatment arm (`arm`) must be one number!")
    }

    if(!is.numeric(alpha) | length(alpha)!=1){
      stop("The significance level (`alpha`) must be one number!")
    }

    if(!is.numeric(opt)| length(opt)!=1 | opt %in% c(1,2) == FALSE){
      stop("The parameter `opt` must be either 1 or 2!")
    }

    if(!is.numeric(prior_prec_tau) | length(prior_prec_tau)!=1){
      stop("The dispersion parameter of the half normal prior, the prior for the between study heterogeneity, (`prior_prec_tau`) must be one number!")
    }

    if(!is.numeric(n.samples) | length(n.samples)!=1){
      stop("The numer of random samples (`n.samples`) must be one number!")
    }

    if(!is.numeric(n.chains) | length(n.chains)!=1){
      stop("The numer of parallel chains for the rjags model (`n.chains`) must be one number!")
    }

    if(!is.numeric(n.iter) | length(n.iter)!=1){
      stop("The number of iterations to monitor of the jags.model (`n.iter`) must be one number!")
    }

    if(!is.numeric(n.adapt) | length(n.adapt)!=1){
      stop("The number of iterations for adaptation (`n.adapt`) must be one number!")
    }

    if(!is.logical(robustify) | length(robustify)!=1){
      stop("The parameter `robustify` must be  TRUE or FALSE!")
    }

    if(!is.numeric(weight) | length(weight)!=1 | weight < 0 | weight > 1){
      stop("The weight given to the non-informative component (`weight`) must be one number between 0 and 1!")
    }
  }

  # 1-alpha is here the decision boundary for the posterior probability,
  # in case that one uses a non-informative/flat prior instead of a MAP prior, the type one error is equal to alpha

  # Data preparation
  ## count number of patients for each treatment in each period
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

  # Create prior for treatment group
  weak_prior_treatment <- mixnorm(c(1,0,1000), param = 'ms') # creates weak prior for the treatment group
  #mixnorm creates normal mixture density, here with one component, mean 0 and sd=1000
  #param= "ms" means that parameters are given via mean and standard deviation

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





  ## calculate odds ratio and with that the estimated treatment effect
  # Calculate Treatment effect and confidence interval by simulation of posterior mean
  random_gen<-as.numeric(rmixdiff(post_act,post_control,n=1000))
  #rmixdiff = random number generation for the difference of two mixture distributions
  theta_hat<-mean(random_gen)
  theta_ci<-quantile(random_gen,probs=c(alpha,1-alpha))
  lower_ci <- as.numeric(theta_ci[1])
  upper_ci <- as.numeric(theta_ci[2])



  ## calculate the "p-value"
  p_value <- 1-pmixdiff(post_control,post_act,0)
  reject_h0 <- ifelse(p_value<alpha, TRUE, FALSE)

  return(list(p_val = p_value,
              treat_effect = theta_hat,
              lower_ci = lower_ci,
              upper_ci = upper_ci,
              reject_h0 = reject_h0))
}



