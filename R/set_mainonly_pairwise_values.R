#' Generate Parameters for Main-Effects-Only Model
#'
#' This function generates a data frame containing main effects, interaction terms, and an intercept term
#' for a main-effects-only model. Interaction terms are included but set to zero, ensuring that only the main effects
#' and the intercept term have non-zero values.
#'
#' @param main_effects Data frame with two columns. The first column contains the names of main effects (e.g., `C1`, `C2`, etc.), 
#' and the second column contains their respective values.
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
#'   \item Interaction terms between each pair of main effects, all set to 0.
#'   \item An intercept term set to `intercept_value`.
#' }
#'
#' @details
#' This function is intended for use in main-effects-only models, where interaction terms are present but set to zero.
#' It processes the main effects provided in `main_effects`, adds an intercept term with a specified value, 
#' and then generates all possible pairwise interactions between main effects, setting their values to zero.
#'
#' @examples
#' # Example usage:
#' main_effects <- data.frame(
#'   effect = c("C1", "C2", "C3", "C4", "C5"),
#'   value = c(1, 1, 0, 1, 0)
#' )
#' result <- set_mainonly_pairwise_values(main_effects, intercept_value = 1)
#' print(result)
#'
#' @importFrom utils combn
#' @export
set_mainonly_pairwise_values <- function(main_effects, 
                                         intercept_value) {
  
  # Extract the main effects and their values
  effects <- main_effects[,1]
  values <- main_effects[,2]
  names(values) <- effects
  
  # Initialize a result list with the intercept term
  result <- list("(Intercept)" = intercept_value)
  
  # Add main effects to the result list (set value to 0)
  for (effect in effects) {
    result[[effect]] <- values[effect]
  }
  
  # Get all pair combinations of the main effects
  pair_combinations <- combn(effects, 2, simplify = FALSE)
  
  
  # Loop through each pair combination
  for (pair in pair_combinations) {

    interaction_name <- paste(pair, collapse = "")
    
    result[[interaction_name]] <- 0 #Set interactions all to 0
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
