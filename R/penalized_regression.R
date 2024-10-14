#' Penalized Regression with LASSO, Ridge, and Elastic Net
#'
#' This function performs penalized regression using LASSO, Ridge, or Elastic Net. It supports
#' both cross-validation (CV) and Bayesian Information Criterion (BIC) methods for selecting the
#' optimal lambda value. The function returns the best model and the coefficients of the fitted model.
#'
#' @param df Dataframe. The dataset containing both the outcome and predictor variables. The first column should be the outcome variable.
#' @param family.value Character. The type of regression model to fit. Default is "binomial".
#' @param model.type Character. The type of penalized regression to perform: "lasso", "ridge", or "elastic net". Default is "lasso".
#' @param max_lambda Numeric. The maximum value of lambda to try when using BIC for lambda selection. Default is 10.
#' @param method Character. The method used to select the optimal lambda: "CV" (cross-validation) or "BIC". Default is "CV".
#' @param opt_lambda Character. The lambda value to select when using cross-validation. "min" selects the lambda with minimum error, "1se" selects the lambda within one standard error of the minimum error. Default is "min".
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{COEF}: A dataframe of the coefficients for the best penalized regression model.
#'   \item \code{X}: The design matrix (covariates) used in the model.
#'   \item \code{Y}: The outcome vector.
#' }
#'
#' @details
#' This function fits a penalized regression model using `glmnet`. It supports three types of models:
#' \itemize{
#'   \item LASSO (alpha = 1)
#'   \item Ridge (alpha = 0)
#'   \item Elastic Net (alpha between 0 and 1)
#' }
#' The user can choose to perform cross-validation (CV) or use the Bayesian Information Criterion (BIC) for selecting the optimal lambda.
#' For Elastic Net, the function iterates over a sequence of alpha values to find the best alpha-lambda combination.
#'
#' @importFrom glmnet glmnet cv.glmnet
#'
#' @examples
#' # Example usage with binomial family and LASSO
#' data <- data.frame(outcome = rbinom(100, 1, 0.5), matrix(rnorm(1000), 100, 10))
#' results <- penalized_regression(data, family.value = "binomial", model.type = "lasso", method = "CV")
#' print(results$COEF)
#'
#' @export

penalized_regression <- function(df,
                                 family.value = "binomial",
                                 model.type = "lasso",
                                 max_lambda = 10,
                                 method = "CV",
                                 opt_lambda = "min"){

  #####################
  #extract the covariates and outcomes
  #####################
  x <- df[,-c(1)]
  y <- df[,1]

  #####################
  #create a data matrix from the dataframe of covariates
  #####################
  input = data.matrix(x)

  #####################
  #create a data matrix from vector of outcomes
  #####################
  output = data.matrix(y)

  #####################
  #alpha
  #####################
  if(model.type == "lasso"){
    alpha.value = 1
  } else if(model.type == "ridge"){
    alpha.value = 0
  } else if(model.type == "elastic net"){
    alpha.value = seq(0,1,length=20)^3
  }

  #maximum lambda
  # mysd <- function(y) sqrt(sum((y-mean(y))^2)/length(y))
  # sx <- scale(input, scale=apply(input, 2, mysd))
  # max_lambda <- norm(t(sx) %*% output, 'i') / nrow(input)

  #####################
  #use CV to find optimal lambda value
  #####################
  if(method == "CV"){

    #####################
    #perform k-fold cross-validation for LASSO
    #####################
    cv_model <- cv.glmnet(x = input, y = output, alpha = 1)

    #####################
    #find optimal lambda that minimizes test MSE
    #####################
    if(opt_lambda == "min"){
      optimal_lambda = cv_model$lambda.min
    }else{
      optimal_lambda = cv_model$lambda.1se
    }

    #####################
    #performs alpha selection for ELASTIC NET
    #####################
    if(model.type == "elastic net"){

      #####################
      #Initialize the deviance value
      #####################
      best_deviance <- Inf

      #####################
      #Go through all the alphas and select the one which gives best (training) deviance
      #####################
      #  s <- sample(nrow(input), 50)
      #  s_input <- input[s, ]
      #  s_output <- output[s]

      #####################
      #Iterate through all values of alpha
      #####################
      for(i in seq(alpha.value)){

        # Start with the first one
        best_alpha <- alpha.value[1]

        # Fit the model
        model <- glmnet(x = input, y = output, alpha = alpha.value[i],
                        lambda = optimal_lambda, family=family.value)

        # Calculate the deviance
        current_deviance <- deviance(model)

        #####################
        #If the model deviance is better than the previous, store it.
        #####################
        if(current_deviance < best_deviance){
          best_deviance <- current_deviance
          best_alpha <- alpha.value[i]
        }
      }
      #####################
      #Optimal alpha value
      #####################
      alpha.value <- best_alpha
    }


  }
  #####################
  #use BIC to find optimal lambda value
  #####################
  else if(method == "BIC"){

    #####################
    #set sequence of lambda
    #####################
    lambdas_to_try <- seq(0.0001, max_lambda, length.out = 100)

    #####################
    #initialized bic and alpha vectors
    #####################
    bic <- numeric(length(lambdas_to_try))
    alpha <- numeric(length(lambdas_to_try))

    #####################
    #start loop to run through all prospective lambdas
    #####################
    for (i in seq(lambdas_to_try)) {

      #set the lambda_value to the lambda we're trying out
      lambda_value <- lambdas_to_try[i]

      # Print the current index in the lambda sequence
      #  print(paste("Testing lambda index:", i, "with lambda value:", lambda_value))

      #####################
      # If we're doing elastic net, add a code chunk which iterates through a sequence of alphas
      # for each potential lambda value.
      #####################
      if(model.type == "elastic net"){

        # Initialize a variable to keep track of the best BIC for the current lambda
        best_bic <- Inf  # Start with a very large number
        best_alpha <- alpha.value[1] # Start with the first one

        #####################
        #Compute the BIC for the model for each alpha
        #####################
        for(j in seq(alpha.value)){
          model <- glmnet(x = input, y = output, alpha = alpha.value[j],
                          lambda = lambda_value, family=family.value)


          #####################
          # Compute the Total Log-Likelihood adjustment
          #####################
          tLL <- model$nulldev - deviance(model)
          k <- model$df
          n <- model$nobs
          current_bic <- log(n) * k - tLL

          # Print alpha iteration
          #      print(paste("Testing alpha index:", j, "with alpha value:", alpha.value[j]))

          #####################
          # Update best BIC if the current BIC is lower
          #####################
          if (current_bic < best_bic) {
            best_bic <- current_bic
            best_alpha <- alpha.value[j]
          }

        }

        #####################
        # Store the best BIC found for this lambda
        #####################
        bic[i] <- best_bic
        alpha[i] <- best_alpha


      }
      #####################
      #For LASSO
      #####################
      else{
        model <- glmnet(x = input, y = output, alpha = alpha.value,
                        lambda = lambda_value, family=family.value)

        #####################
        # Compute the Total Log-Likelihood adjustment
        #####################
        tLL <- model$nulldev - deviance(model)
        k <- model$df
        n <- model$nobs

        #####################
        # Compute information criteria
        #####################
        bic[i] <- log(n)*k - tLL
      }

    }

    #####################
    # Find the index of the minimum BIC
    #####################
    min_bic_index <- which.min(bic)

    #####################
    # (Elastic Net) Find the index of the minimum alpha
    #####################
    if(model.type == "elastic net"){
      alpha.value <- alpha[min_bic_index]
    }

    #####################
    # Extract the lambda corresponding to the minimum BIC
    #####################
    optimal_lambda <- lambdas_to_try[min_bic_index]

  }

  #####################
  #find the best model
  #####################
  best_model <- glmnet(x = input, y = output,
                       alpha = alpha.value,
                       lambda = optimal_lambda,
                       family = family.value)

  #####################
  #extract the coefficients of the best model
  #####################
  coef_best_model <- coef(best_model)

  #create a dataframe with the coefficients of the best model
  #along with their names
  # coef_data_frame <- data.frame(DRUG = coef_best_model@Dimnames[[1]][coef_best_model@i + 1],
  #                               COEF = coef_best_model@x)

  #output
  # coef_data_frame

  #####################
  #results
  #####################
  #store as a dataframe
  coef_df <- as.data.frame(t(as.data.frame(as.matrix(coef_best_model))))


  ###############
  #Save the list
  ###############
  out=list(COEF=coef_df,X=input,Y=output)
  return(out)
}
