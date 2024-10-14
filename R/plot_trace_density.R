#' Generate Trace Plots and Density Plots for Beta Coefficients
#'
#' This function generates trace plots and density plots for specified beta coefficients from
#' a given simulation. It removes the burn-in period from the simulation results and then
#' generates a trace plot (showing the evolution of beta values over iterations) and a
#' histogram (showing the distribution of beta values) for each beta.
#'
#' @param data Dataframe. The dataset containing the simulation output, including iteration number and beta values.
#' @param betas Character vector. The names of the beta coefficients for which trace plots and density plots should be generated.
#' @param burnin Integer. The number of burn-in iterations to remove from the simulation.
#' @param sim Integer. The simulation number to be plotted.
#'
#' @return The function generates trace plots and histograms for each beta coefficient.
#'
#' @details
#' The function first removes the burn-in period by filtering out iterations below the specified burn-in value.
#' It then filters the data to the specified simulation number and loops through each beta coefficient to generate:
#' \itemize{
#'   \item A trace plot showing the beta values over iterations.
#'   \item A histogram showing the distribution of beta values after the burn-in period.
#' }
#' The plots are displayed one after another for each beta.
#'
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes geom_line geom_histogram labs theme_minimal element_text
#' @importFrom rlang sym
#'
#' @examples
#' # Example usage:
#' data <- data.frame(Iter = rep(1:1000, 3),
#'                    sim = rep(1:3, each = 1000),
#'                    beta1 = rnorm(3000),
#'                    beta2 = rnorm(3000))
#' plot_trace_density(data, betas = c("beta1", "beta2"), burnin = 100, sim = 1)
#'
#' @export


plot_trace_density <- function(data, #output of the simulation
                               betas, #which trace plots for betas
                               burnin, #remove the burn-in (could be 0)
                               sim #which simulation to plot
){

  ##############
  # Ensure that betas is a vector
  ##############
  if (!is.vector(betas)) {
    betas <- c(betas)
  }

  ##############
  # Remove the burn-in
  ##############
  data.subset <- data %>%
    filter(Iter > burnin)

  ##############
  # Filter the data for the specified simulation
  ##############
  filtered_data <- data.subset %>%
    filter(sim == sim)

  ##############
  # Loop through each beta and create a trace plot and density plot
  ##############
  for (beta in betas) {

    ##############
    # Generate the trace plot
    ##############
    plot_trace <- ggplot(filtered_data, aes(x = Iter, y = !!sym(beta))) +
      geom_line(color = "black", linewidth = 0.25) +  # Updated from size to linewidth
      labs(title = paste0("Trace Plot for ", beta, " (Simulation ", sim, ")"),
           x = "Iteration",
           y = paste0(beta, " Value")) +
      theme_minimal() +  # Clean theme
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centered, bold title
        axis.title.x = element_text(size = 12),  # X-axis label font size
        axis.title.y = element_text(size = 12),  # Y-axis label font size
        axis.text = element_text(size = 10)  # Axis tick label font size
      )

    ##############
    # Generate the histogram plot
    ##############
    plot_hist <- ggplot(filtered_data, aes(x = !!sym(beta))) +
      geom_histogram() +
      labs(title = paste0("Histogram for ", beta, " (Simulation ", sim, ")"),
           x = beta) +
      theme_minimal() +  # Clean theme
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centered, bold title
        axis.title.x = element_text(size = 12),  # X-axis label font size
        axis.text = element_text(size = 10)  # Axis tick label font size
      )

    ##############
    # Print the plot for each beta
    ##############
    print(plot_trace)
    print(plot_hist)
  }
}
