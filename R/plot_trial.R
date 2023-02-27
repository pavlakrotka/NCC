#' Function for visualizing the simulated trial
#'
#' @description This function creates a plot visualizing the trial progress over time.
#'
#' @param treatments Vector with indices of assigned arms for each participant, ordered by time, e.g. column `treatment` from the dataframe resulting from the `datasim_bin()` or `datasim_cont()` function.
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
#' plot_trial(treatments = trial_data$treatment)
#'
#'
#' @return ggplot showing trial progress over time.
#' @author Pavla Krotka


plot_trial <- function(treatments){

  ggplot() +
    geom_line(aes(x=c(1:length(treatments)), y=treatments, color=as.factor(treatments)), size=5, alpha=0.8) +
    labs(x="Patients in the trial", y="Treatment") +
    theme_bw() +
    theme(legend.position="none") +
    scale_y_continuous(breaks = unique(treatments), labels = unique(treatments))

}
