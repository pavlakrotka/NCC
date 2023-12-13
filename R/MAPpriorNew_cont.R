#' Analysis for continuous data using the MAP Prior approach
#'
#' @description This function performs analysis of continuous data using the Meta-Analytic-Predictive (MAP) Prior approach. The method borrows data from non-concurrent controls to obtain the prior distribution for the control response in the concurrent periods.
#'
#' @param data Data frame with trial data, e.g. result from the `datasim_bin()` function. Must contain columns named 'treatment', 'response' and 'period'.
#' @param arm Integer. Index of the treatment arm under study to perform inference on (vector of length 1). This arm is compared to the control group.
#' @param alpha Double. Decision boundary (one-sided). Default=0.025
#' @param opt Integer (1 or 2). If opt==1, all former periods are used as one source; if opt==2, periods get separately included into the final analysis. Default=2.
#' @param prior_prec_tau Double. Precision parameter (\eqn{1/\sigma^2_{\tau}}) of the half normal hyperprior, the prior for the between study heterogeneity. Default=4.
#' @param prior_prec_eta Double. Precision parameter (\eqn{1/\sigma^2_{\eta}}) of the normal hyperprior, the prior for the hyperparameter mean of the control mean. Default=0.001.
#' @param n_samples Integer. Number of how many random samples will get drawn for the calculation of the posterior mean, the p-value and the CI's. Default=1000.
#' @param robustify Logical. Indicates whether a robust prior is to be used. If TRUE, a mixture prior is considered combining a MAP prior and a weakly non-informative component prior. Default=TRUE.
#' @param weight Double. Weight given to the non-informative component (0 < weight < 1) for the robustification of the MAP prior according to Schmidli (2014). Default=0.1.
#' @param check Logical. Indicates whether the input parameters should be checked by the function. Default=TRUE, unless the function is called by a simulation function, where the default is FALSE.
#' @param ... Further arguments passed by wrapper functions when running simulations.
#'
#' @importFrom RBesT automixfit
#' @importFrom RBesT postmix
#' @importFrom RBesT mixnorm
#' @importFrom RBesT pmixdiff
#' @importFrom RBesT rmixdiff
#' @importFrom RBesT robustify
#' @importFrom RBesT gMAP
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom stats gaussian
#'
#' @export
#'
#' @details
#'
#' The MAP approach derives the prior distribution for the control response in the concurrent periods by combining the control information from the non-concurrent periods with a non-informative prior.
#'
#' The model for the continuous response \eqn{y_{js}} for the control patient \eqn{j} in the non-concurrent period \eqn{s} is defined as follows:
#'
#' \deqn{E(y_{js}) = \eta_s}
#'
#' where \eqn{\eta_s} represents the control mean in the non-concurrent period \eqn{s}.
#'
#' The means for the non-concurrent controls in period \eqn{s} are assumed to have a normal prior distribution with mean \eqn{\mu_{\eta}} and variance \eqn{\tau^2}:
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
#' \eqn{P(\mu_{treatment}-\mu_{control}>0) \ge 1-}`alpha`.
#' In case of a non-informative prior this coincides with the frequentist type I error.
#'
#' @examples
#'
#' \donttest{
#' trial_data <- datasim_cont(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' theta = rep(0.25, 3), lambda = rep(0.15, 4), sigma = 1, trend = "stepwise")
#'
#' MAPpriorNew_cont(data = trial_data, arm = 3)
#' }
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
#' @author Marta Bofill Roig, Katharina Hees
#'
#' @references
#'
#' Robust meta-analytic-predictive priors in clinical trials with historical control information. Schmidli, H., et al. Biometrics 70.4 (2014): 1023-1032.
#'
#' Applying Meta-Analytic-Predictive Priors with the R Bayesian Evidence Synthesis Tools. Weber, S., et al. Journal of Statistical Software 100.19 (2021): 1548-7660.


MAPpriorNew_cont <- function(data,
                             arm,
                             alpha = 0.025,
                             opt = 2,
                             prior_prec_tau = 4,
                             prior_prec_eta = 0.001,
                             n_samples = 1000,
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
      stop("The dispersion parameter of the normal prior, the prior for the control mean, (`prior_prec_eta`) must be one number!")
    }

    if(!is.numeric(n_samples) | length(n_samples)!=1){
      stop("The numer of random samples (`n_samples`) must be one number!")
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
  tab_count <- table(data$treatment, data$period)

  ## count number of groups and number of periods
  number_of_groups <- dim(tab_count)[1] # number of groups incl control
  number_of_periods <- dim(tab_count)[2] #total number of periods

  ## get start and end period of each treatment
  treatment_start_period <- numeric(number_of_groups)
  treatment_end_period <- numeric(number_of_groups)

  for (i in 1:number_of_groups){
    treatment_start_period[i] <- min(which(tab_count[i,] > 0))
    treatment_end_period[i] <- max(which(tab_count[i,] > 0))
  }


  ## get concurrent and non-concurrent controls of treatment = arm
  cc <- data[data$treatment == 0 & (data$period >= treatment_start_period[arm + 1]) & (data$period <= treatment_end_period[arm + 1]), ]
  ncc <- data[data$treatment == 0 & data$period < treatment_start_period[arm + 1], ]

  ## ncc sample sizes per period
  # tab_ncc_ss <- table(ncc$period)
  #
  # while (any(tab_ncc_ss<=10)) { # if less than 10 patients per period -> pool with next period
  #   if(min(which(tab_ncc_ss<=10))!=dim(tab_ncc_ss)){
  #     ncc[ncc$period==min(which(tab_ncc_ss<=10))+1,]$period <- min(which(tab_ncc_ss<=10))
  #     ncc[ncc$period>min(which(tab_ncc_ss<=10))+1,]$period <- ncc[ncc$period>min(which(tab_ncc_ss<=10)),]$period-1
  #   } else {
  #     ncc[ncc$period==min(which(tab_ncc_ss<=10)),]$period <- max(which(tab_ncc_ss>10))
  #   }
  #   tab_ncc_ss <- table(ncc$period)
  # }


  ## get patients treated with treatment=arm
  t_treatment <- data[data$treatment == arm,]

  if(dim(ncc)[1] !=0){ # in case there are NCC data available

    if(opt==1){
      ncc_data <- data.frame(period=c(0), n=c(0), y=c(0), y.se=c(0))
      y_p <- mean(ncc$response)
      n_p <- length(ncc$response)
      y_pse <- sd(ncc$response)/sqrt(n_p)
      ncc_data[1,] <- c(1,n_p,y_p,y_pse)
    }

    if(opt==2){
      # summary per period
      ncc_data <- data.frame(period=c(0), n=c(0), y=c(0), y.se=c(0))
      for(i in 1: treatment_start_period[arm + 1]-1){
        ncc_period <- ncc[ncc$period==i,]
        y_p <- mean(ncc_period$response)
        n_p <- length(ncc_period$response)
        y_pse <- sd(ncc_period$response)/sqrt(n_p)
        ncc_data[i,] <- c(i,n_p,y_p,y_pse)
      }
    }

    # prior
    options(RBesT.MC.control=list(adapt_delta=0.999))
    map_mcmc <- gMAP(cbind(y, y.se) ~ 1 | period,
                     weights=n, data=ncc_data,
                     family=gaussian,
                     beta.prior=cbind(0, sqrt(1/prior_prec_eta)), #prior_prec_eta
                     tau.dist="HalfNormal", tau.prior=cbind(0, sqrt(1/prior_prec_tau))) #prior_prec_tau


    prior_control <- automixfit(map_mcmc)

    # consider mean of the map
    # summary(map_mcmc)[3]
    if (robustify==TRUE) {
      prior_control <- suppressMessages(robustify(prior_control, weight=weight, mean=0)) # suppressMessages() to not report the default prior scale
    }
  } else {

    prior_control <- mixnorm(c(1, 0, 1000), param = 'ms') # creates weak prior for the control group, in case there are no NCCs
  }

  # Create prior for treatment group
  weak_prior_treatment <- mixnorm(c(1, 0, 1000), param = 'ms') # creates weak prior for the treatment group
  #mixnorm creates normal mixture density, here with one component, mean 0 and sd=1000
  #param= "ms" means that parameters are given via mean and standard deviation

  ## get summary data of treatment group and concurrent controls
  y.act <- mean(t_treatment$response)
  n.act <- length(t_treatment$response)
  y.act.se <- sd(t_treatment$response)/sqrt(n.act)


  y.pbo <- mean(cc$response)
  n.pbo <- length(cc$response)
  y.pbo.se <- sd(cc$response)/sqrt(n.pbo)

  ## obtain posterior distributions
  post_act <- postmix(weak_prior_treatment, m=y.act, se=y.act.se)
  # postmix calculates the posterior distribution given a prior (here "weak_prior"),
  # where the prior is a mixture of conjugate distributions.
  # The posterior is then also a mixture of conjugate distributions.

  post_control <- postmix(prior_control, m=y.pbo, se=y.pbo.se)



  ## Calculate Treatment effect and confidence interval by simulation of posterior mean
  random_gen <- as.numeric(rmixdiff(post_act, post_control, n=n_samples))
  # rmixdiff = random number generation for the difference of two mixture distributions
  theta_hat <- mean(random_gen)
  theta_ci <- quantile(random_gen, probs=c(alpha, 1-alpha))
  lower_ci <- as.numeric(theta_ci[1])
  upper_ci <- as.numeric(theta_ci[2])



  ## calculate the "p-value"
  p_value <- 1-pmixdiff(post_control, post_act, 0)
  reject_h0 <- ifelse(p_value<alpha, TRUE, FALSE)

  return(list(p_val = p_value,
              treat_effect = theta_hat,
              lower_ci = lower_ci,
              upper_ci = upper_ci,
              reject_h0 = reject_h0))
}



