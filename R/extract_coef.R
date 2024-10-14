#' Extract Coefficients from Hierarchical LASSO (HierLASSO) Model
#'
#' This function extracts and formats the main effects and interaction effects from a Hierarchical LASSO (HierLASSO) model into a tidy dataframe.
#' The resulting dataframe includes metric names for both the main effects and interaction terms, along with their corresponding coefficient values.
#'
#' @param ME Numeric vector. The main effects from the HierLASSO model.
#' @param IE Matrix. The interaction effects from the HierLASSO model. Only the lower triangular part of the matrix is considered for interactions.
#'
#' @return A dataframe with two columns:
#' \itemize{
#'   \item \code{metric}: The names of the effects (e.g., "C1" for main effects, "C1C2" for interaction effects).
#'   \item \code{value}: The corresponding coefficient values for the main and interaction effects.
#' }
#'
#' @details
#' The function constructs a dataframe by first extracting the main effects and labeling them with a "C" prefix (e.g., "C1" for the first main effect).
#' It then extracts interaction effects from the lower triangle of the interaction matrix (IE), assigning appropriate metric names (e.g., "C1C2" for an interaction between "C1" and "C2").
#'
#' @examples
#' # Example usage:
#' ME <- c(0.5, 0.3, 0.2)
#' IE <- matrix(c(0, 0.1, 0.2, 0, 0.3, 0, 0, 0), nrow = 3)
#' results <- extract_coef(ME, IE)
#' print(results)
#'
#' @export

extract_coef <- function(ME, #main effects
                         IE) #interaction effects
{

  #########################################
  # Set up the dataframe to return results
  #########################################
  results <- data.frame(metric = character(), value = numeric(), stringsAsFactors = FALSE)

  #########################################
  # Add main effects
  #########################################
  for (i in 1:length(ME)) {
    metric_name <- paste0("C", i)  # Main effect metric name
    results <- rbind(results, data.frame(metric = metric_name, value = ME[i]))
  }

  #########################################
  # Add interaction effects (lower triangle of the interaction matrix)
  #########################################
  for (i in 1:(nrow(IE) - 1)) {
    for (j in (i + 1):ncol(IE)) {
      metric_name <- paste0("C", i, "C", j)  # Interaction effect metric name
      interaction_value <- IE[i, j]
      results <- rbind(results, data.frame(metric = metric_name, value = interaction_value))
    }
  }

  #########################################
  # Return the dataframe
  #########################################
  return(results)
}
