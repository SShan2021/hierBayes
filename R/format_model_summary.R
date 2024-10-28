#' Transform Data for Model Summary
#'
#' This function takes a dataframe and formats it by adding method, category, and type columns.
#' It renames the existing columns, reorders them, and removes the first row (typically used for simulation data).
#'
#' @param df Dataframe. The input dataframe that contains the model results.
#' @param model Character. The name of the model (e.g., "OG_hier", "HS_hier") to be added to the dataframe.
#' @param type Character. The type of the data (e.g., "mean", "sd") to be added to the dataframe.
#'
#' @return A transformed dataframe with columns for method, category, type, and values.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Converts the input into a dataframe (if it isn’t already).
#'   \item Adds the method and type to the dataframe.
#'   \item Renames the columns, reorders them, and removes the first row (usually containing simulation data).
#' }
#'
#' @examples
#' # Example usage:
#' df <- data.frame(Value = rnorm(10))
#' transformed_df <- data.transform(df, model = "OG_hier", type = "mean")
#' print(transformed_df)
#'
#' @export

format_model_summary <- function(df, model, type){

  #make into a data.frame
  df <- as.data.frame(df)

  #add in the category, method, type
  df$Category <- rownames(df)
  df$Method <- model
  df$Type <- type

  #rename the rownames
  rownames(df) <- 1:nrow(df)

  #rename the columns
  colnames(df)[1] <- "Value"

  #reorder the columns
  df <- df[,c("Method", "Category",
              "Type", "Value")]

  #remove the simulation row
  df <- df[-1, ]

  #return dataframe
  df
}
