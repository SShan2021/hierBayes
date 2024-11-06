#' Generate Parameters for Interaction-Only Model with Limited Non-Zero Interactions
#'
#' This function generates a data frame containing main effects, interaction terms, and an intercept term
#' under an interaction-only modeling assumption. All main effects are set to zero, while a specified number 
#' of interaction terms are assigned non-zero values.
#'
#' @param main_effects Data frame with one column containing the names of main effects (e.g., `C1`, `C2`, etc.). 
#' All main effects are assigned a value of 0 in the resulting data frame.
#' @param new_value Numeric. The value to assign to a specified number of interaction terms.
#' @param how_many_inter Integer. The maximum number of interaction terms to set to `new_value`.
#' @param intercept_value Numeric. The value to assign to the intercept term.
#'
#' @return A data frame with two columns:
#' \itemize{
#'   \item \code{term}: The names of main effects, interaction terms, and the intercept term.
#'   \item \code{value}: The assigned values for each term.
#' }
#'
#' The data frame includes:
#' \enumerate{
#'   \item Main effects, all set to 0.
#'   \item Interaction terms, with up to `how_many_inter` terms set to `new_value` and the rest set to 0.
#'   \item An intercept term set to `intercept_value`.
#' }
#'
#' @details
#' This function is intended for use in interaction-only models, where all main effects are constrained to zero,
#' and only a limited number of interaction terms have non-zero values. It processes the main effects provided in `main_effects`,
#' sets each main effect to 0, and generates pairwise interactions. Up to `how_many_inter` interaction terms are assigned 
#' `new_value`, with the remainder set to 0.
#'
#' @examples
#' # Example usage:
#' main_effects <- data.frame(
#'   effect = c("C1", "C2", "C3", "C4", "C5")
#' )
#' result <- set_intonly_pairwise_values(main_effects, new_value = 999, how_many_inter = 3, intercept_value = 1)
#' print(result)
#'
#' @importFrom utils combn
#' @export

set_intonly_pairwise_values <- function(main_effects, 
                                        new_value, 
                                        how_many_inter, 
                                        intercept_value) {
  
  # Extract the main effects and their values
  effects <- main_effects[,1]
  values <- rep(0, length = length(effects))
  names(values) <- effects
  
  # Initialize a result list with the intercept term
  result <- list("(Intercept)" = intercept_value)
  
  # Add main effects to the result list (set value to 0)
  for (effect in effects) {
    result[[effect]] <- values[effect]
  }
  
  # Get all pair combinations of the main effects
  pair_combinations <- combn(effects, 2, simplify = FALSE)
  
  # Intialize a counter for number of non-zero interaction terms 
  iter <- 1
  
  # Loop through each pair combination
  for (pair in pair_combinations) {
    
    interaction_name <- paste(pair, collapse = "")
    
    if (iter <= how_many_inter) {
      result[[interaction_name]] <- new_value  # Set to the new_value 
      
      iter <- iter + 1 # increment by 1
    } else{
      result[[interaction_name]] <- 0
    }
  }
  
  # Convert the result list to a data frame with columns "term" and "value"
  result_df <- data.frame(
    term = names(result),
    value = unlist(result),
    stringsAsFactors = FALSE
  )
  rownames(result_df) <- 1:nrow(result_df)
  
  return(result_df)
  
}