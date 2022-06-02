#' Function for visualizing the simulated trial
#'
#' @description Creates a plot visualizing the trial progress over time
#'
#' @param treatments Vector with treatment indicators ordered by time, e.g. column `treatment` from the dataframe resulting from the `datasim_bin()` or `datasim_cont()` function
#'
#' @import ggplot2
#'
#' @export
#' 
#' @examples 
#' 
#' trial_data <- datasim_bin(num_arms = 3, n_arm = 100, d = c(0, 100, 250),
#' p0 = 0.7, OR = rep(1.8, 3), lambda = rep(0.15, 4), trend="stepwise")
#' 
#' plottrial(treatments = trial_data$treatment)
#'
#'
#' @return ggplot showing trial progress over time
#' @author Pavla Krotka


plottrial <- function(treatments){
  
  ggplot() +
    geom_line(aes(x=c(1:length(treatments)), y=treatments, color=as.factor(treatments)), size=5, alpha=0.8) +
    labs(x="Time", y="Treatment") +
    theme_bw() +
    theme(legend.position="none")
  
}