#' Summarize Data with 95% HPD (Highest Posterior Density) Intervals
#'
#' This function calculates summary statistics by groups, including the 95% Highest Posterior Density (HPD) intervals.
#' It computes the minimum, mean, standard deviation, maximum, and HPD lower and upper bounds for each group.
#'
#' @param data Dataframe. The dataset to be summarized, containing the variables of interest.
#' @param group Character. The name of the grouping variable to group the data by.
#' @param num Numeric. The confidence interval for the HPD calculation. Defaults to 0.95.
#'
#' @return A long-format dataframe with summary statistics and HPD intervals for each variable, organized by group.
#'
#' @details
#' The function groups the data by the specified `group` variable and calculates the following summary statistics
#' for each variable: minimum, mean, standard deviation, maximum, and the lower and upper bounds of the HPD interval.
#' The results are then pivoted into a long format, and the variable names are split into two separate columns (metric and statistic).
#'
#' @importFrom dplyr group_by summarise across
#' @importFrom rlang sym
#' @importFrom bayestestR hdi
#' @importFrom tidyr pivot_longer separate
#'
#' @examples
#' # Example usage:
#' data <- data.frame(sim = rep(1:5, each = 10), value = rnorm(50))
#' summarize_with_hdi(data, "sim", 0.95)
#'
#' @export

summarize_with_hdi <- function(data,
                               group,
                               num = 0.95){

  ##############
  #calculate summary statistics by group
  ##############
  r1 <- data %>%
    group_by(!!sym(group)) %>% #group by simlation
    summarise(
      across(everything(),
             list(min = min,  #min
                  mean = mean, #mean
                  sd = sd,    #sd
                  max = max,  #max
                  hdilower = ~ hdi(.x, ci = num)$CI_low, #hdi lower CI
                  hdiupper = ~ hdi(.x, ci = num)$CI_high), #hdi upper CI
             .names = "{.col}_{.fn}") #names
    )


  #############
  #pivot longer
  ##############
  r1_long <- r1 %>%
    pivot_longer(
      cols = -sim,                     # Exclude the `sim` column from reshaping
      names_to = "variable",           # Name of the new column for the variable names
      values_to = "value"              # Name of the new column for the values
    )

  #############
  #split variable into two columns
  ##############
  r1_long_edited <- r1_long %>%
    separate(variable,
             into = c("metric", "statistic"), sep = "_")

  #############
  #return solution
  ##############
  r1_long_edited

}
