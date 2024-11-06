#' Generate Parameters for Anti-Hierarchical Model with Specified Interaction Terms
#'
#' This function generates a data frame containing main effects, interaction terms, and an intercept term 
#' under an anti-hierarchical modeling assumption. It assigns non-zero values to a limited number of interaction terms 
#' even when the corresponding main effects are zero, thus allowing interactions that do not depend on the main effects.
#'
#' @param main_effects Data frame with two columns. The first column contains the names of main effects (e.g., `C1`, `C2`, etc.), 
#' and the second column contains their respective values.
#' @param new_value Numeric. The value to assign to selected interaction terms.
#' @param how_many_inter Integer. The maximum number of interaction terms to set to `new_value`, even if the corresponding main effects are zero.
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
#'   \item Main effects with their original values from `main_effects`.
#'   \item Interaction terms, set to 0 if both main effects are non-zero. Up to `how_many_inter` interaction terms can be set to `new_value` 
#'   even when one or both main effects in the pair are zero.
#'   \item An intercept term set to `intercept_value`.
#' }
#'
#' @details
#' This function is intended for use in anti-hierarchical models, where interaction terms can be non-zero even if the corresponding main effects 
#' are zero. It processes the main effects provided in `main_effects`, adds an intercept term, and generates pairwise interactions. 
#' If both main effects in a pair are non-zero, the interaction term is set to 0. Otherwise, up to `how_many_inter` interaction terms 
#' are assigned `new_value`.
#'
#' @examples
#' # Example usage:
#' main_effects <- data.frame(
#'   effect = c("C1", "C2", "C3", "C4", "C5"),
#'   value = c(1, 1, 0, 1, 0)
#' )
#' result <- set_antihier_pairwise_values(main_effects, new_value = 999, how_many_inter = 3, intercept_value = 1)
#' print(result)
#'
#' @importFrom utils combn
#' @export

set_antihier_pairwise_values <- function(main_effects, 
                                         new_value, 
                                         how_many_inter, 
                                         intercept_value) {
  
  # Extract the main effects and their values
  effects <- main_effects[,1]
  values <- main_effects[,2]
  names(values) <- effects
  
  # Initialize a result list with the intercept term
  result <- list("(Intercept)" = intercept_value)
  
  # Add main effects to the result list
  for (effect in effects) {
    result[[effect]] <- values[effect]
  }
  
  # Get all pair combinations of the main effects
  pair_combinations <- combn(effects, 2, simplify = FALSE)
  
  # Intialize a counter for number of non-zero interaction terms 
  iter <- 1
  
  # Loop through each pair combination
  for (pair in pair_combinations) {
    # Extract values for the current pair
    value_i <- values[pair[1]]
    value_j <- values[pair[2]]
    
    # Determine interaction value
    interaction_name <- paste(pair, collapse = "")
    if (value_i > 0 & value_j > 0) {
      result[[interaction_name]] <- 0  # Set to 0 if both > 0
    } else if(iter <= how_many_inter){
      result[[interaction_name]] <- new_value  # Set to new_value otherwise
      
      iter <- iter + 1
    }
    else{
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