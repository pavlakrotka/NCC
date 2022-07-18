#' To-DOs'
#'
#' Check input parameters for prior specification and whether to include additional ones like robistify yes, no, specifications for robustification,...
#' Check p-value and different definitions of it : posterior p-value, posterior predictive p-value, posterior probability,...
#' Posterior odds estimated by drawing random samples from the posterior distributions of the incidence rates of the two groups separately
#'
#' Check: Take into account the following warning statement?: "In total 1 divergent transitions occurred during the sampling phase.
#' Please consider increasing adapt_delta closer to 1 with the following command prior to gMAP:
#' options(RBesT.MC.control=list(adapt_delta=0.999))" ?
#'
#'
#' Use of non-concurrent controls for binary data
#'
#' @description Borrows data from non-concurrent controls for binary endpoints for the concurrent controls arm by the use of a MAP prior
#'
#' @param arm Indicator of the treatment arm under study to perform inference on (vector of length 1)
#' @param data Simulated trial data, e.g. result from the `datasim_bin()` function
#' @param alpha Type I error. Default=0.025
#' @param tau.dist Type of prior distribution for tau; supported priors are HalfNormal (default), TruncNormal, Uniform, Gamma, InvGamma, LogNormal, TruncCauchy, Exp and Fixed, default is "HalfNormal", for more details see help of the gMAP function of the RBesT package
#' @param tau.prior Parameters of prior distribution for tau, default is 1
#' @param beta.prior Mean and standard deviation for normal priors of regression coefficients, see section prior specification help of the gMAP function of the RBesT package, default is 2
#' tau.dist="HalfNormal",tau.prior=1,beta.prior=2
#' @param opt Binary, opt=1: all former periods are use as one source, opt=2: periods form different sources get separately included into the final analysis
#' @param prec.ef Defines the number of how many random samples will get drawn for the calculation of the posterior mean and the CIs, default is 1000
#'
#' @importFrom RBesT mixbeta
#' @importFrom RBesT gMAP
#' @importFrom RBesT automixfit
#' @importFrom RBesT postmix
#' @importFrom RBesT decision2S
#' @importFrom RBesT pmixdiff
#' @importFrom RBesT rmix
#' @importFrom stats quantile
#' @importFrom dplyr filter
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr n
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#'
#' trial_data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")
#'
#' MAPprior_bin(data = trial_data, arm = 3)
#'
#'
#' @return List containing the p-value (one-sided), estimated treatment effect, 95% confidence interval and an indicator whether the null hypothesis was rejected or not for the investigated treatment
#' @author Katharina Hees


MAPprior_bin <- function(data,arm,alpha=0.025,opt=1,tau.dist="HalfNormal",tau.prior=1,beta.prior=2,prec.ef=1000){

  #require("RBesT")
  #require("stats")
  #require("dplyr")

  if (opt %in% c(1,2) == FALSE) stop("Wrong parameter input. Parameter opt has to be 1 or 2")

  # 1-alpha is here the decision boundary for the posterior probability,
  # in case that one uses a non-informative/flat prior instead of a MAP prior, the type one error is equal to alpha

  #Data preparation

  ##count number of patients for each treatment in each period
  tab_count<-table(data$treatment,data$period)

  ##count number of groups and number of periods
  number_of_groups<-dim(table(data$treatment,data$period))[1] # number of groups incl control
  number_of_periods<-dim(table(data$treatment,data$period))[2] #total number of periods



  ##get start and end period of each treatment
  treatment_start_period<-numeric(number_of_groups)
  treatment_end_period<-numeric(number_of_groups)
  for (i in 1:number_of_groups){
    treatment_start_period[i]<-min(which(table(data$treatment,data$period)[i,]>0))
    treatment_end_period[i]<-max(which(table(data$treatment,data$period)[i,]>0))
  }


  ##get concurrent and non-concurrent controls of treatment=arm
  cc<-data %>% filter(treatment==0) %>% filter(period > treatment_start_period[arm+1]-1) %>%
    filter(period < treatment_end_period[arm+1]+1)
  ncc<-data %>% filter(treatment==0) %>% filter(period < treatment_start_period[arm+1])

  if(opt==1){
  ncc<-ncc %>% mutate(period=0)
  }

  ##get patients treated with treatment=arm
  t_treatment<-data %>% filter(treatment==arm)


  #Create prior for treatment group
  weak_prior <- mixbeta(c(1,0.5,0.5))

  #Create prior for controls (based on nccs)
  if (nrow(ncc)==0){ #case of no prior information == no NCCs available
    map<- mixbeta(c(1,0.5,0.5))} else {
      dat<-ncc %>% group_by(period) %>% summarise(r=sum(response), n=n())
      map_mcmc <- gMAP(cbind(r, n-r) ~ 1 | period,
                      data=dat,
                      tau.dist=tau.dist,
                      tau.prior=tau.prior,
                      beta.prior=beta.prior,
                      family=binomial)
    map <- automixfit(map_mcmc) #fits a series of mixtures of conjugate distributions
  }



  #Obtain posterior

  ## get summary data of treatment group and concurrent controls
  r.act <- as.numeric(t_treatment %>% summarise(sum(response)))
  n.act <- as.numeric(t_treatment %>% summarise(n=n()))


  r.pbo <- as.numeric(cc %>% summarise(sum(response)))
  n.pbo <- as.numeric(cc %>% summarise(n=n()))


  ## obtain posterior distributions
  post_act <- postmix(weak_prior, n=n.act,r=r.act)
  #postmix calculates the posterior distribution given a prior (here "weak_prior"),
  #where the prior is a mixture of conjugate distributions.
  #The posterior is then also a mixture of conjugate distributions.

  post_pbo <- postmix(map,n=n.pbo,r=r.pbo)




  # Define decision function
  poc <- decision2S(pc=1-alpha,qc=0,lower.tail=FALSE)
  #Creates a one-sided decision function on the basis of the difference distribution in a 2 sample situation
  #Conditions for acceptance:
  #P(x1 - x2 > 0) > 1-alpha  x_1 actice group, x2 control group

  #Get decision
  reject_h0 <- poc(post_act, post_pbo)



  # Calculate posterior probaility ("p-value"), treatment effect and confidence interval by simulation of posterior mean


  #Risk difference
  #random_gen<-as.numeric(rmixdiff(post_act,post_pbo,n=prec.ef))  #rmixdiff = random number generation for the difference of two mixture distributions
  #treat_effect<-mean(random_gen)

  #Posterior probability
  post_prob <- 1-pmixdiff(post_pbo,post_act,0)

  #Oddsratio estimation by posterior means
  p1<-as.numeric(rmix(post_act,n=prec.ef))
  p2<-as.numeric(rmix(post_pbo,n=prec.ef))
  oddsratio<-p1/(1-p1)/(p2/(1-p2))

  treat_effect<-mean(oddsratio)

  delta_ci<-quantile(oddsratio,probs=c(alpha,1-alpha))
  lower_ci<-as.numeric(delta_ci[1])
  upper_ci<-as.numeric(delta_ci[2])

  return(list(p_val = post_prob,
              treat_effect = treat_effect,
              lower_ci = lower_ci,
              upper_ci = upper_ci,
              reject_h0 = as.logical(reject_h0)))
}











