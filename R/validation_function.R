#' Validation Function for Penalized Regression Models
#'
#' This function validates the performance of a penalized regression model by calculating
#' the Mean Squared Error (MSE), sensitivity (SENS), and specificity (SPEC). It compares the
#' true regression coefficients (B0) to the posterior coefficients from a fitted model (PostB),
#' then refits the model using the significant coefficients and evaluates performance metrics based on this refit.
#'
#' @param input Dataframe or matrix. The covariates (X) from the dataset used to fit the model.
#' @param output Vector. The outcome variable (Y) from the dataset.
#' @param B0 Dataframe. The true regression coefficients (`metric` and `value`) for each predictor.
#' @param PostB Dataframe. A summary of the posterior regression coefficients (`metric` and `value`) from the penalized regression model.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{COEF}: A dataframe of the re-estimated significant coefficients from the refitted model.
#'   \item \code{MSE}: The Mean Squared Error between the true and estimated regression coefficients.
#'   \item \code{SENS}: The sensitivity, or the proportion of true non-zero coefficients correctly identified as non-zero in the refitted model.
#'   \item \code{SPEC}: The specificity, or the proportion of true zero coefficients correctly identified as zero in the refitted model.
#' }
#'
#' @details
#' This function takes the true coefficients (B0) and compares them with the posterior estimates (PostB) from the penalized regression model.
#' It calculates the MSE between the true and estimated coefficients and refits a new model using only the significant coefficients
#' identified. From this refitted model, sensitivity (the proportion of correctly identified non-zero coefficients)
#' and specificity (the proportion of correctly identified zero coefficients) are calculated. The function returns the re-estimated
#' coefficients, sensitivity, and specificity.
#'
#' @examples
#' # Example usage:
#' input <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' output <- rbinom(100, 1, 0.5)
#' B0 <- data.frame(metric = paste0("X", 1:10), value = rnorm(10))
#' PostB <- data.frame(metric = paste0("X", 1:10), value = rnorm(10))
#' results <- validation_function(input, output, B0, PostB)
#' print(results)
#'
#' @export

validation_function <- function(input, #X
                                output, #Y
                                B0, #true regression coefficients
                                PostB #summary of posterior regression coefficients)
){

  #########################################
  #make the parameter_list into a dataframe
  #########################################
  parameter_list <- as.data.frame(as.matrix(B0))

  #########################################
  #remove the intercept
  #########################################
  simulation_list <- PostB[!(PostB$metric == "(Intercept)"), ]

  #########################################
  #full bind the beta simulated to the beta real
  #########################################
  beta <- parameter_list %>%
    right_join(simulation_list, by = "metric")
  colnames(beta) <- c("metric", "B0", "PostB")

  beta$B0 <- as.numeric(beta$B0)
  beta$PostB <- as.numeric(beta$PostB)

  #########################################
  #calculate the MSE (beta - beta_hat)^2
  #########################################
  mse <- mean((beta[,"B0"] - beta[,"PostB"])^2)

  #########################################
  #refit the model using the beta_simulated
  #########################################
  #extract the covariates and outcomes
  x <- input
  y <- output

  #extract non-zero covariates from simulation
  keep_X <- simulation_list$metric[simulation_list$value != 0]
  x <- x[,keep_X] #extract the columns from the beta simulated
  new_data <- as.data.frame(cbind(y, x))
  colnames(new_data)[1] <- "Y" #rename the outcome column
  new_outcome <- coef(summary(glm(Y ~ ., data = new_data, family = "binomial"))) %>%
    as.data.frame()

  # Filter for significant coefficients (p-value < 0.05)
  significant_coefficients <- new_outcome %>%
    filter(`Pr(>|z|)` < 0.05) %>%
    rownames_to_column(var = "metric") %>%
    filter(!(metric == "(Intercept)")) #get rid of the intercept

  significant_coefficients <- significant_coefficients[ , c("metric", "Estimate")]
  colnames(significant_coefficients) <- c("metric", "value")

  #########################################
  #full-bind the significant beta to the real beta
  #########################################
  beta_v1 <- parameter_list %>%
    full_join(significant_coefficients, by = "metric")
  colnames(beta_v1) <- c("metric", "B0", "PostB")

  beta_v1 <- beta_v1 %>%
    mutate_all(~ifelse(is.na(.),
                       0,
                       .))

  beta_v1 <- beta_v1[!(beta_v1$metric == "(Intercept)"), ]
  beta_v1$B0 <- as.numeric(beta_v1$B0)
  beta_v1$PostB <- as.numeric(beta_v1$PostB)

  #########################################
  #Sensitivity
  #########################################
  #extract the rows with non-zero values from true parameters
  true <- beta_v1[beta_v1$B0 != 0, ]

  non_zero_count <- sum(true$PostB != 0)  # Count of non-zero elements
  total_count_nz <- length(true$B0)       # Total number of non-zero elements in the true parameters
  sens <- non_zero_count / total_count_nz  # Calculate percentage

  #########################################
  #Specificity
  #########################################
  #extract the rows with zero values from true parameters
  false <- beta_v1[beta_v1$B0 == 0, ]

  zero_count <- sum(false$PostB == 0)  # Count of non-zero elements
  total_count_z <- length(false$B0)       # Total number of zero elements in the true parameters
  spec <- zero_count / total_count_z  # Calculate percentage

  #########################################
  #Output
  #########################################
  output <- list(COEF = new_outcome,
                 MSE = mse,
                 SENS = sens,
                 SPEC = spec)
  return(output)

}
