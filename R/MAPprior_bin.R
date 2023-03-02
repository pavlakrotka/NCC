#' Analysis for binary data using the MAP Prior approach
#'
#' @description This function performs analysis of binary data using the Meta-Analytic-Predictive (MAP) Prior approach. The method borrows data from non-concurrent controls to obtain the prior distribution for the control response in the concurrent periods.
#'
#' @param data Data frame with trial data, e.g. result from the `datasim_bin()` function. Must contain columns named 'treatment', 'response' and 'period'.
#' @param arm Integer. Index of the treatment arm under study to perform inference on (vector of length 1). This arm is compared to the control group.
#' @param alpha Double. Decision boundary (one-sided). Default=0.025
#' @param opt Integer (1 or 2). If opt==1, all former periods are used as one source; if opt==2, periods get separately included into the final analysis. Default=2.
#' @param prior_prec_tau Double. Precision parameter (\eqn{1/\sigma^2_{\tau}}) of the half normal hyperprior, the prior for the between study heterogeneity. Default=4.
#' @param prior_prec_eta Double. Precision parameter (\eqn{1/\sigma^2_{\eta}}) of the normal hyperprior, the prior for the hyperparameter mean of the control log-odds. Default=0.001.
#' @param n_samples Integer. Number of how many random samples will get drawn for the calculation of the posterior mean, the p-value and the CI's. Default=1000.
#' @param n_chains Integer. Number of parallel chains for the rjags model. Default=4.
#' @param n_iter Integer. Number of iterations to monitor of the jags.model. Needed for coda.samples. Default=4000.
#' @param n_adapt Integer. Number of iterations for adaptation, an initial sampling phase during which the samplers adapt their behavior to maximize their efficiency. Needed for jags.model. Default=1000.
#' @param robustify Logical. Indicates whether a robust prior is to be used. If TRUE, a mixture prior is considered combining a MAP prior and a weakly non-informative component prior. Default=TRUE.
#' @param weight Double. Weight given to the non-informative component (0 < weight < 1) for the robustification of the MAP prior according to Schmidli (2014). Default=0.1.
#' @param check Logical. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.
#' @param ... Further arguments passed by wrapper functions when running simulations.
#'
#' @importFrom RBesT automixfit
#' @importFrom RBesT postmix
#' @importFrom RBesT mixbeta
#' @importFrom RBesT robustify
#' @importFrom RBesT rmix
#' @importFrom stats quantile
#' @importFrom rjags jags.model
#' @importFrom rjags coda.samples
#'
#' @export
#'
#' @details
#'
#' The MAP approach derives the prior distribution for the control response in the concurrent periods by combining the control information from the non-concurrent periods with a non-informative prior.
#'
#' The model for the binary response \eqn{y_{js}} for the control patient \eqn{j} in the non-concurrent period \eqn{s} is defined as follows:
#'
#' \deqn{g(E(y_{js})) = \eta_s}
#'
#' where \eqn{g(\cdot)} denotes the logit link function and \eqn{\eta_s} represents the control log odds in the non-concurrent period \eqn{s}.
#'
#' The log odds for the non-concurrent controls in period \eqn{s} are assumed to have a normal prior distribution with mean \eqn{\mu_{\eta}} and variance \eqn{\tau^2}:
#'
#' \deqn{\eta_s \sim \mathcal{N}(\mu_{\eta}, \tau^2)}
#'
#'
#' For the hyperparameters \eqn{\mu_{\eta}} and \eqn{\tau}, normal and half-normal hyperprior distributions are assumed, with mean 0 and variances \eqn{\sigma^2_{\eta}} and \eqn{\sigma^2_{\tau}}, respectively:
#'
#' \deqn{\mu_{\eta} \sim \mathcal{N}(0, \sigma^2_{\eta})}
#'
#' \deqn{\tau \sim HalfNormal(0, \sigma^2_{\tau})}
#'
#'
#' The MAP prior distribution \eqn{p_{MAP}(\eta_{CC})} for the control response in the concurrent periods is then obtained as the posterior distribution of the parameters \eqn{\eta_s} from the above specified model.
#'
#' If `robustify=TRUE`, the MAP prior is robustified by adding a weakly-informative mixture component \eqn{p_{\mathrm{non-inf}}}, leading to a robustified MAP prior distribution:
#'
#' \deqn{p_{rMAP}(\eta_{CC}) = (1-w) \cdot p_{MAP}(\eta_{CC}) + w \cdot p_{\mathrm{non-inf}}(\eta_{CC})}
#'
#' where \eqn{w} (parameter `weight`) may be interpreted as the degree of skepticism towards borrowing strength.
#'
#'
#'
#' In this function, the argument `alpha` corresponds to \eqn{1-\gamma}, where \eqn{\gamma} is the decision boundary. Specifically, the posterior probability of the difference distribution under the null hypothesis is such that:
#' \eqn{P(p_{treatment}-p_{control}>0) \ge 1-}`alpha`.
#' In case of a non-informative prior this coincides with the frequentist type I error.
#'
#' @examples
#'
#' trial_data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")
#'
#' MAPprior_bin(data = trial_data, arm = 3)
#'
#'
#' @return List containing the following elements regarding the results of comparing `arm` to control:
#'
#' - `p-val` - posterior probability that the log-odds ratio is less than zero
#' - `treat_effect` - posterior mean of log-odds ratio
#' - `lower_ci` - lower limit of the (1-2*`alpha`)*100% credible interval for log-odds ratio
#' - `upper_ci` - upper limit of the (1-2*`alpha`)*100% credible interval for log-odds ratio
#' - `reject_h0` - indicator of whether the null hypothesis was rejected or not (`p_val` < `alpha`)
#'
#' @author Katharina Hees
#'
#' @references
#'
#' Robust meta-analytic-predictive priors in clinical trials with historical control information. Schmidli, H., et al. Biometrics 70.4 (2014): 1023-1032.
#'
#' Applying Meta-Analytic-Predictive Priors with the R Bayesian Evidence Synthesis Tools. Weber, S., et al. Journal of Statistical Software 100.19 (2021): 1548-7660.



MAPprior_bin <- function(data,
                         arm,
                         alpha = 0.025,
                         opt = 2,
                         prior_prec_tau = 4,
                         prior_prec_eta = 0.001,
                         n_samples = 1000,
                         n_chains = 4,
                         n_iter = 4000,
                         n_adapt = 1000,
                         robustify = TRUE,
                         weight = 0.1,
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

    if(!is.numeric(opt)| length(opt)!=1 | opt %in% c(1,2) == FALSE){
      stop("The parameter `opt` must be either 1 or 2!")
    }

    if(!is.numeric(prior_prec_tau) | length(prior_prec_tau)!=1){
      stop("The dispersion parameter of the half normal prior, the prior for the between study heterogeneity, (`prior_prec_tau`) must be one number!")
    }

    if(!is.numeric(prior_prec_eta) | length(prior_prec_eta)!=1){
      stop("The dispersion parameter of the normal prior, the prior for the control log-odds, (`prior_prec_eta`) must be one number!")
    }

    if(!is.numeric(n_samples) | length(n_samples)!=1){
      stop("The numer of random samples (`n_samples`) must be one number!")
    }

    if(!is.numeric(n_chains) | length(n_chains)!=1){
      stop("The numer of parallel chains for the rjags model (`n_chains`) must be one number!")
    }

    if(!is.numeric(n_iter) | length(n_iter)!=1){
      stop("The number of iterations to monitor of the jags.model (`n_iter`) must be one number!")
    }

    if(!is.numeric(n_adapt) | length(n_adapt)!=1){
      stop("The number of iterations for adaptation (`n_adapt`) must be one number!")
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
    treatment_start_period[i] <- min(which(table(data$treatment, data$period)[i,] > 0))
    treatment_end_period[i] <- max(which(table(data$treatment, data$period)[i,] > 0))
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
      prior_prec_tau = prior_prec_tau,
      prior_prec_eta = prior_prec_eta
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
        mu  ~ dnorm(0, prior_prec_eta)

        # prediction for p in a new study
        logit_p.pred ~ dnorm(mu, inv_tau2)
        p.pred <- 1 / (1 + exp(-logit_p.pred))
      }"




    fit <- jags.model(file = textConnection(model_text),
                      data = ncc_control_data_jags,
                      #inits = inits,
                      n.chains = n_chains,
                      n.adapt = n_adapt,
                      quiet = TRUE
    )


    # Draw samples from the above fitted MCMC model
    help_mcmc_samples <- coda.samples(fit, "p.pred", n.iter = n_iter)
    help_samples <- do.call(rbind.data.frame, help_mcmc_samples)[,1]


    ## Fit Beta mixture to MCMC samples
    prior_control <- automixfit(help_samples, type="beta", Nc=3) # Nc=3: fixed number of mixture components, to speed the code up

    if (robustify==TRUE) {
      prior_control <- robustify(prior_control, weight=weight)
    }
  }
  else{

    prior_control <- mixbeta(c(1,0.5,0.5))

  }


  ## get summary data of treatment group and concurrent controls
  r.act <- sum(t_treatment$response)
  n.act <- length(t_treatment$response)

  r.pbo <- sum(cc$response)
  n.pbo <- length(cc$response)


  post_control <- postmix(prior_control, n=n.pbo, r=r.pbo)

  post_treatment <- mixbeta(c(1, 0.5+r.act, 0.5+n.act-r.act))



  samples_control <- rmix(post_control, n=n_samples)
  samples_treat <- rmix(post_treatment, n=n_samples)
  p1 <- samples_treat
  p2 <- samples_control

  ## calculate log odds ratio and with that the estimated treatment effect
  logoddsratio <- log((p1 / (1 - p1)) / (p2 / (1 - p2)))
  treat_effect <- mean(logoddsratio)

  delta_ci <- quantile(logoddsratio, probs = c(alpha, 1 - alpha))
  lower_ci <- as.numeric(delta_ci[1])
  upper_ci <- as.numeric(delta_ci[2])

  ## calculate the "p-value"
  p_value <- mean(samples_treat - samples_control < 0)
  reject_h0 <- ifelse(p_value<alpha, TRUE, FALSE)

  return(list(p_val = p_value,
              treat_effect = treat_effect,
              lower_ci = lower_ci,
              upper_ci = upper_ci,
              reject_h0 = reject_h0
  ))
}




