#' Calculate Performance Statistics for Bayesian Shrinkage Models
#'
#' This function calculates various performance metrics for a Bayesian shrinkage model, including
#' Mean Squared Error (MSE), bias, coverage probability, mean interval width (MIW), sensitivity, and specificity.
#' It compares the true regression coefficients with the posterior summary of regression coefficients and
#' computes these statistics based on the Highest Posterior Density (HPD) intervals.
#'
#' @param B0 Dataframe. A dataframe containing the true regression coefficients (`value`) and their corresponding `metric` (beta names).
#' @param PostB Dataframe. A summary of posterior regression coefficients from the model, including posterior means (`mean`), lower HDI bounds (`hdilower`), and upper HDI bounds (`hdiupper`).
#'
#' @return A list containing the following performance statistics:
#' \itemize{
#'   \item \code{MSE}: Mean Squared Error.
#'   \item \code{bias}: The average bias of the estimates.
#'   \item \code{CP}: Coverage probability, the proportion of true values falling within the HDI intervals.
#'   \item \code{MIW}: Mean interval width, the average width of the HDI intervals.
#'   \item \code{sens}: Sensitivity, the proportion of true signals correctly identified as non-zero.
#'   \item \code{spec}: Specificity, the proportion of non-signals correctly identified as zero.
#' }
#'
#' @details
#' The function removes iterations and irrelevant simulations from the posterior data (`PostB`), and joins the posterior estimates (mean and HDI intervals) with the true coefficients (`B0`).
#' It then calculates various statistics such as bias, MSE, sensitivity, specificity, coverage probability, and mean interval width based on the comparison between the true values and the posterior estimates.
#'
#' @importFrom dplyr filter left_join select rename inner_join
#'
#' @examples
#' # Example usage:
#' B0 <- data.frame(metric = c("beta1", "beta2", "beta3"), value = c(0.5, 0.7, 0))
#' PostB <- data.frame(metric = rep(c("beta1", "beta2", "beta3"), each = 3),
#'                     statistic = rep(c("mean", "hdilower", "hdiupper"), times = 3),
#'                     value = rnorm(9))
#' results <- bshrink_stats(B0, PostB)
#' print(results)
#'
#' @export
bshrink_stats <- function(B0, #true regression coefficients
                          PostB #summary of posterior regression coefficients
){

  ##############
  # Remove the iterations and simulation info
  ##############
  PostB$sim <- NULL
  PostB <- PostB %>%
    filter(!(metric == "Iter"))

  ##############
  # Extract CI bounds and Mean
  ##############
  # Extract Posterior Mean
  PMean <- PostB %>%
    filter(statistic == "mean")

  # Extract CI bounds
  CI_low <- PostB %>%
    filter(statistic == "hdilower")

  CI_high <- PostB %>%
    filter(statistic == "hdiupper")

  # Join together
  tog <- CI_low %>%
    left_join(CI_high, by = c("metric")) %>%
    left_join(PMean, by = c("metric")) %>%
    rename(hdilower = value.x, hdiupper = value.y,
           pmean = value) %>%
    dplyr::select(metric, hdilower, hdiupper, pmean)

  # Append to original matrix
  B0 <- B0 %>%
    inner_join(tog, by = c("metric"))

  ##############
  # Mean Interval Width
  ##############
  MIW <- B0$hdiupper - B0$hdilower

  ##############
  # Coverage Probability
  ##############
  CP <- (B0$value > B0$hdilower & B0$value < B0$hdiupper) * 1

  ##############
  # Bias
  ##############
  bias <- B0$value - B0$pmean

  ##############
  # Mean Squared Error
  ##############
  MSE <- bias^2

  ##############
  # Sensitivity
  ##############
  # Extract the signals
  signal <- B0 %>%
    filter(value != 0)

  # How many intervals don't cover 0
  sens <- (sign(signal$hdilower) == sign(signal$hdiupper)) * 1

  ##############
  # Specificity
  ##############
  # Extract the non-signals
  nsignal <- B0 %>%
    filter(value == 0)

  # How many intervals cover 0
  spec <- (sign(nsignal$hdilower) != sign(nsignal$hdiupper)) * 1

  ##############
  # Return the mean results
  ##############
  out <- list(
    MSE = mean(MSE),
    bias = mean(bias),
    CP = mean(CP),
    MIW = mean(MIW),
    sens = mean(sens),
    spec = mean(spec)
  )

  return(out)

}
