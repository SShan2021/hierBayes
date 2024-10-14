#' Reformat LASSO Output
#'
#' This function reformats the output from a LASSO regression model, transforming the results into a
#' more readable dataframe format, where the coefficients are labeled by their respective variable names (metric).
#'
#' @param PostB Dataframe or matrix. The posterior estimates or coefficients from a LASSO regression model.
#'
#' @return A dataframe with two columns:
#' \itemize{
#'   \item \code{metric}: The names of the variables (e.g., DRUG or beta coefficients).
#'   \item \code{value}: The corresponding coefficient values for each variable.
#' }
#'
#' @details
#' The function transposes the input, extracts the row names (assumed to be variable names such as DRUG or beta names),
#' and assigns them to a new column labeled `metric`. It swaps the order of columns to ensure that the metric (variable names)
#' are listed first, followed by the corresponding coefficient values. The result is a tidy dataframe suitable for further analysis.
#'
#' @examples
#' # Example usage:
#' PostB <- matrix(rnorm(10), nrow = 5, dimnames = list(c("beta1", "beta2", "beta3", "beta4", "beta5"), NULL))
#' result <- lasso_reformat(PostB)
#' print(result)
#'
#' @export
lasso_reformat <- function(PostB){

  #########################################
  # Transpose the matrix or dataframe and convert to a dataframe
  #########################################
  simulation_list <- as.data.frame(t(PostB))

  #########################################
  # Extract the row names as a new column (assumed to be variable names)
  #########################################
  simulation_list$DRUG <- rownames(simulation_list)

  #########################################
  # Rename the columns to "metric" and "value"
  #########################################
  colnames(simulation_list) <- c("metric", "value")

  #########################################
  # Reset row names for the new dataframe
  #########################################
  rownames(simulation_list) <- 1:dim(simulation_list)[1]

  #########################################
  # Swap the order of the columns
  #########################################
  simulation_list[,1:2] <- simulation_list[,2:1]

  return(simulation_list)
}
