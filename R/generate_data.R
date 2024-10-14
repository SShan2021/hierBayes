#' Generate Simulated Drug and Interaction Data
#'
#' This function simulates a dataset with specified main effects and interaction terms.
#' It first generates multivariate normal data for the main effects, and then computes
#' all possible two-way interaction effects between the predictors. The resulting data
#' is scaled by dividing each feature by its standard deviation (sd) to standardize the
#' values.
#'
#' @param n Integer. The number of rows (samples) in the dataset.
#' @param p Integer. The number of main effects (predictors) in the dataset.
#' @param seed Integer. A seed value for random data generation to ensure reproducibility. Default is 123.
#'
#' @return A scaled dataframe containing the simulated main effects and two-way interactions,
#' with an additional intercept column.
#'
#' @details
#' The function generates multivariate normal data for the main effects using a mean
#' vector of zeros and an identity covariance matrix. It then computes all two-way
#' interactions between these main effects and optionally scales the resulting dataset.
#' The output includes a column for the intercept, which is useful for regression models.
#'
#' @examples
#' # Generate a dataset with 100 samples and 10 main effects
#' data <- generate_data(n = 100, p = 10, seed = 123)
#'
#' @importFrom MASS mvrnorm
#' @importFrom gtools combinations
#' @export

generate_data <- function(n, p, seed = 123){

  ###############
  #set seed
  ###############
  set.seed(seed)

  ###############
  #Generate Main Effects
  ###############
  #mean vector
  mu <- rep(0, p)

  #identity covariance matrix
  sigma <- diag(p)

  #generate main effects
  X <- mvrnorm(n = n, mu = mu, Sigma = sigma)

  #make absolute value
  #  X <- abs(X)

  #add column names
  colnames(X) <- paste0("C", 1:p)

  ###############
  #Log Transforming and Scaling Main Effects
  ###############
  # log transformation
  # logX <- log(X + 1)

  # find SD
  # sd_logx <- apply(logX,2,sd)

  # scaling
  # logX_star <- apply(logX,2,function(x)(x/sd(x)))
  # apply(logX_star,2,sd)

  ###############
  #Interaction Effects
  ###############
  # Number of main effects
  p_tot <- dim(X)[2]

  # Number of interaction effects
  n_interac <- nrow(combinations(p_tot, 2))

  # Total number of individuals
  N = dim(X)[1]

  # Return all possible two-way combinations
  rr <- combinations(p_tot, 2)

  # Initialize interaction matrix
  # logZ <- array(NA,dim=c(N,n_interac))
  Z <- array(NA,dim=c(N,n_interac))

  # Fill in interaction matrix
  for(j in 1:n_interac)
  {
    Z[,j] <- X[,rr[j,1]]*X[,rr[j,2]]
  }

  # Add in colnames
  colnames(Z) <- combn(1:p, 2, function(x) paste0("C", x[1], "C", x[2]))

  # Find SD
  # sd_logz <-  apply(logZ,2,sd)

  # Scaling
  #  logZ_star <- apply(logZ,2,function(x)(x/sd(x)))
  # apply(logZ_star,2,sd)


  ###############
  #Combining Data (and Normalize Data)
  ###############
  # Combining log-transformed main & Interaction effect
  # X <- cbind(logX_star,logZ_star)
  X <- cbind(X, Z)

  # Find SD
  sd_X <- apply(X, 2, sd)

  # Scaling
  X_star <- apply(X,2,function(x)(x/sd(x)))

  # adding intercept
  X_star <- cbind(rep(1,N), X_star)

  # rename intercept
  colnames(X_star)[1] <- "Int"

  # make as a matrix
  X_star <- as.matrix(X_star)

  #return
  X_star

}
