#' Generate True Parameters for Hierarchical Model with Main and Interaction Effects
#'
#' This function generates a data frame containing main effects, interaction terms, and an intercept term
#' under a hierarchical modeling assumption. Specifically, it creates interaction terms between main effects 
#' only when both main effects have non-zero values, supporting the hierarchical principle that interaction effects 
#' are present only if the corresponding main effects are non-zero.
#'
#' @param main_effects Data frame with two columns. The first column contains the names of main effects (e.g., `C1`, `C2`, etc.), 
#' and the second column contains their respective values. Main effects with values greater than 0 will contribute to non-zero interactions.
#' @param new_value Numeric. The value to assign to interaction terms when both main effect values in a pair are greater than 0.
#' @param intercept_value Numeric. The value to assign to the intercept term.
#'
#' @return A data frame with two columns:
#' \itemize{
#'   \item \code{metric}: The names of main effects, interaction terms, and the intercept term.
#'   \item \code{value}: The assigned values for each main effect, interaction term, and intercept term.
#' }
#'
#' The data frame includes:
#' \enumerate{
#'   \item Main effects with their original values from `main_effects`.
#'   \item Interaction terms between each pair of main effects, set to `new_value` if both main effects in the pair have values > 0, otherwise 0.
#'   \item An intercept term set to `intercept_value`.
#' }
#'
#' @details
#' This function is designed for use in hierarchical models, where interaction effects are present only if the corresponding main effects 
#' are non-zero. It processes the main effects provided in `main_effects`, adds an intercept term, and then computes all possible 
#' pairwise interactions. If both main effect values in a pair are greater than 0, the interaction term is assigned `new_value`; 
#' otherwise, the interaction term is set to 0, enforcing the hierarchical assumption.
#'
#' @examples
#' # Example usage:
#' main_effects <- data.frame(
#'   effect = c("C1", "C2", "C3", "C4", "C5"),
#'   value = c(1, 1, 0, 1, 0)
#' )
#' result <- set_pairwise_values(main_effects, new_value = 999, intercept_value = 1)
#' print(result)
#'
#' @importFrom utils combn
#' @export

set_pairwise_values <- function(main_effects, 
                                new_value, 
                                intercept_value) {
  
  # Extract the main effects and their values
  effects <- main_effects[,1]
  values <- main_effects[,2]
  names(values) <- effects
  
  # Initialize a result list to hold the output
  result <- list("(Intercept)" = intercept_value)
  
  # Add main effects to the result list
  for (effect in effects) {
    result[[effect]] <- values[effect]
  }
  
  # Get all pair combinations of the main effects
  pair_combinations <- combn(effects, 2, simplify = FALSE)
  
  # Loop through each pair combination
  for (pair in pair_combinations) {
    # Extract values for the current pair
    value_i <- values[pair[1]]
    value_j <- values[pair[2]]
    
    # Determine interaction value
    interaction_name <- paste(pair, collapse = "")
    if (value_i > 0 & value_j > 0) {
      result[[interaction_name]] <- new_value  # Set to new_value if both > 1
    } else {
      result[[interaction_name]] <- 0  # Set to 0 otherwise
    }
  }
  
  # Convert the result list to a data frame with columns "term" and "value"
  result_df <- data.frame(
    metric = names(result),
    value = unlist(result),
    stringsAsFactors = FALSE
  )
  rownames(result_df) <- 1:nrow(result_df)
  
  return(result_df)
}