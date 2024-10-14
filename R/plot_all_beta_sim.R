#' Plot All Beta Coefficients in a Simulation with HDI and True Values
#'
#' This function generates a plot for all specified beta coefficients from a given simulation,
#' displaying the mean and 95% HDI (Highest Posterior Density) intervals, along with the true
#' beta values for comparison.
#'
#' @param data Dataframe. The dataset containing the summarized statistics (mean, hdilower, hdiupper) for each beta coefficient.
#' @param B0 Dataframe. A dataset containing the true beta coefficients for comparison.
#' @param betas Character vector. The names of the beta coefficients to be plotted.
#' @param sim Integer. The simulation number to be plotted.
#' @param title Character. The title of the plot.
#'
#' @return A ggplot2 plot displaying the mean estimates, HDI intervals, and true values for all specified beta coefficients in the given simulation.
#'
#' @details
#' The function filters the data for the specified simulation and the provided beta coefficients, then combines the summarized statistics
#' (mean, hdilower, hdiupper) with the true beta values. The resulting plot shows the mean estimates as points, HDI intervals as error bars,
#' and true beta values as green points for comparison.
#'
#' @importFrom dplyr filter select
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar labs theme_minimal
#'
#' @examples
#' # Example usage:
#' data <- data.frame(sim = rep(1:5, each = 3),
#'                    metric = rep(c("beta1", "beta2", "beta3"), 5),
#'                    statistic = rep(c("mean", "hdilower", "hdiupper"), 5),
#'                    value = rnorm(15))
#' B0 <- data.frame(metric = c("beta1", "beta2", "beta3"), value = c(0.5, 0.7, -0.3))
#' betas <- c("beta1", "beta2", "beta3")
#' plot_all_beta_sim(data, B0, betas, sim = 1, title = "Simulation 1: Beta Coefficients")
#'
#' @export

plot_all_beta_sim <- function(data,
                              B0,
                              betas,
                              sim,
                              title){

  #############
  #add true to true
  #############
  B0$statistic <- "true"

  #############
  #remove iter
  #############
  d <- data %>%
    filter(!(metric == "Iter"))

  #############
  #filter for the simulation
  #############
  d <- d %>%
    filter(sim == sim) %>%
    dplyr::select(metric, value, statistic)

  #############
  #add in the truth to the simulation
  #############
  d <- rbind(d, B0)

  #############
  # Filter the data to focus on the mean, hdilower, and hdiupper statistics
  #############
  beta_plot <- d %>%
    filter(statistic %in% c("mean", "hdilower", "hdiupper", "true"))

  #############
  # Reshape the data so that mean, hdilower, and hdiupper are in separate columns
  #############
  beta_plot_wide <- beta_plot %>%
    pivot_wider(names_from = statistic, values_from = value)

  #############
  # Extract the variables of interest
  #############
  beta_plot_wide_m <- beta_plot_wide[beta_plot_wide$metric %in% betas,]
  beta_plot_wide_m$metric <- factor(beta_plot_wide_m$metric,
                                    levels = betas)

  #############
  # Create the plot
  #############
  plot <- ggplot(beta_plot_wide_m, aes(x = metric, y = mean)) +
    geom_point(size = 2, color = "blue") +  # Plot the mean
    geom_errorbar(aes(ymin = hdilower, ymax = hdiupper), width = 0.2, color = "red") +  # Plot the HDI interval
    geom_point(aes(y = true), size = 3, color = "green") +  # Add a horizontal line for the true value
    labs(title = title,
         x = "Parameter",
         y = "") +
    theme_minimal()

  print(plot)

}
