#' Plot Beta Coefficients with HDI (Highest Posterior Density) Intervals
#'
#' This function generates plots for the specified beta coefficients, displaying the mean and
#' 95% HDI (Highest Posterior Density) intervals for each simulation. The function allows
#' for comparison against true beta values, with mean estimates, HDI lower and upper bounds,
#' and a reference line for the true beta.
#'
#' @param data Dataframe. The dataset containing the summarized statistics (mean, hdilower, hdiupper) for each beta coefficient.
#' @param betas Character vector. The names of the beta coefficients to be plotted.
#' @param truth Numeric vector. The true values of the beta coefficients for comparison.
#'
#' @return A set of ggplot2 plots for each beta coefficient, displaying the mean estimates and HDI intervals.
#'
#' @details
#' The function loops through each specified beta coefficient, filtering the data to extract the mean and HDI bounds
#' for each simulation. It then generates a plot with the mean estimates as points, the HDI intervals as error bars,
#' and a dashed horizontal line indicating the true beta value. The plots are displayed separately for each beta.
#'
#' @importFrom dplyr filter
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_hline facet_wrap labs theme_minimal
#'
#' @examples
#' # Example usage:
#' data <- data.frame(sim = rep(1:5, each = 3),
#'                    metric = rep(c("beta1", "beta2", "beta3"), 5),
#'                    statistic = rep(c("mean", "hdilower", "hdiupper"), 5),
#'                    value = rnorm(15))
#' betas <- c("beta1", "beta2", "beta3")
#' truth <- c(0.5, 0.7, -0.3)
#' plot_beta_with_hdi(data, betas, truth)
#'
#' @export

plot_beta_with_hdi <- function(data, betas, truth) {

  ##############
  # Ensure that betas is a vector
  ##############
  if (!is.vector(betas)) {
    betas <- c(betas)
  }

  ##############
  # Loop through each beta and create a plot
  ##############
  for (i in seq_along(betas)) {
    beta <- betas[i]  # Get the ith beta

    ##############
    # Filter the data for the specified beta and the statistics
    ##############
    filtered_data <- data %>%
      filter(metric == beta) %>%
      filter(statistic %in% c("mean", "hdilower", "hdiupper"))

    ##############
    # Reshape the data so that mean, hdilower, and hdiupper are in separate columns
    ##############
    wide_data <- filtered_data %>%
      pivot_wider(names_from = statistic, values_from = value)

    ##############
    # Generate the plot
    ##############
    plot <- ggplot(wide_data, aes(x = sim, y = mean)) +
      geom_point(size = 2, color = "blue") +  # Plot the mean
      geom_errorbar(aes(ymin = hdilower, ymax = hdiupper), width = 0.2, color = "red") +  # Plot the HDI interval
      geom_hline(yintercept = truth[i], linetype = "dashed", color = "black") +  # Add a horizontal line at y = 1
      facet_wrap(~metric) +  # Separate plots for each beta
      labs(title = paste0(beta, ": Mean with HDI For Each Simulation"),
           x = "Simulation",
           y = "") +
      theme_minimal()

    ##############
    # Print the plot for each beta
    ##############
    print(plot)
  }
}
