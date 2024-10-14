#' Hierarchical Bayesian Shrinkage with Horseshoe Prior
#'
#' This function implements a hierarchical Bayesian shrinkage method with a horseshoe prior for both main and interaction effects.
#' It uses Markov Chain Monte Carlo (MCMC) to sample from the posterior distribution and updates the global and local shrinkage
#' parameters iteratively using the horseshoe prior.
#'
#' @param data Dataframe. The full dataset including the response and predictors.
#' @param main Integer. The number of main effects in the model.
#' @param demo Integer. The number of demographic effects (covariates) in the model.
#' @param MCMC Integer. The number of iterations to run MCMC.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{BETA}: A matrix of the posterior samples of regression coefficients (main and interaction effects).
#'   \item \code{ETA}: A matrix of the local shrinkage parameters for main effects.
#'   \item \code{THETA}: A matrix of the local shrinkage parameters for interaction effects.
#'   \item \code{A}: A vector of the global shrinkage parameters for main effects.
#'   \item \code{B}: A vector of the global shrinkage parameters for interaction effects.
#' }
#'
#' @details
#' This function initializes the necessary parameters for Bayesian shrinkage using the horseshoe prior, including latent variables, precision parameters, and shrinkage parameters.
#' The horseshoe prior introduces auxiliary parameters (\code{nu}, \code{delta}, \code{kappa}, and \code{alpha}) to allow for heavy-tailed shrinkage that is particularly effective for sparse models.
#' It then iteratively updates the coefficients, global, and local shrinkage parameters for both main and interaction effects using MCMC.
#'
#' @examples
#' # Example usage:
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' result <- hierBayes(data, main = 5, demo = 2, MCMC = 1000)
#' print(result$BETA)
#'
#' @note This function is adapted from the work found at \url{https://github.com/debamitakundu8/Bayesian-Inference-of-Chemical-Mixtures-in-Risk-Assessment-Incorporating-the-Hierarchical-Principle}.
#'
#' @importFrom MASS mvrnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom gtools combinations
#' @importFrom stats rgamma
#' @importFrom MCMCpack rinvgamma
#' @export

hierBayes <- function(data,  #full dataset w/ response
                      main, #number of main effects in the model
                      demo,  #the number of demographic effects
                      MCMC #number of iterations to run MCMC
){



  ###############
  #Initializes parameters in model
  ###############
  # Make as dataframe
  X <- as.data.frame(data)

  # Response
  Y <- X[,1]

  # Remove the response
  X[,1] <- NULL

  # X to be a matrix
  X <- as.matrix(X)

  # Total number of individuals
  N = dim(X)[1]

  # Number of covariates
  p_cov <- demo

  # Number of main effects
  p_tot <- main

  # Interaction combinations
  rr <- combinations(p_tot, 2)

  # Number of interaction effects
  n_interac <- nrow(rr)

  # Total number of parameters in the model
  p <- dim(X)[2]


  ###############
  #Hyperparameters and Initial Values for MCMC
  ###############
  #df (robit link w/ 7 df closely approximates the logit link)
  v <- 7

  #latent variable
  w <- rep(1,N)*0

  #precision parameter for latent variable
  lambda <- rep(.01,N)

  #coefficient vector (both main and interaction effects)
  beta <- rep(1,p)

  #local shrinkage parameter for main effects
  eta_true <- eta <- array(1,p_tot)

  #local shrinkage parameter for interaction effects
  theta_true <- theta <- array(1,n_interac)

  #global shrinkage parameter for main effects
  a <- 1

  #global shrinkage parameter for interaction effects
  b <- 1

  #shape and rate parameters for priors on eta, theta, a, b
  a1 <- a2 <- a3 <- a4 <- a5 <- 1
  b1 <- b2 <- b3 <- b4 <- b5 <- 1

  #initialize the auxiliary parameters
  nu <- replicate(p_tot, rinvgamma(n = 1, shape = 1/2, scale = 1/(b1^2)))
  delta <- replicate(n_interac, rinvgamma(n=1, shape = 1/2, scale = 1/(b2^2)))
  kappa <- rinvgamma(n=1, shape = 1/2, scale = 1/(b3^2))
  alpha <- rinvgamma(n=1, shape = 1/2, scale = 1/(b4^2))


  ###############
  #MCMC Setup
  ###############
  #beta
  BETA <- array(NA,dim=c(MCMC,p))

  #w
  W <- array(NA,dim=c(MCMC,N))

  #eta
  ETA <- array(NA,dim=c(MCMC,p_tot))

  #theta
  THETA <- array(NA,dim=c(MCMC,n_interac))

  #lambda
  LAMBDA <- array(NA,dim=c(MCMC,N))

  #a
  A <- array(MCMC)

  #b
  B <- array(MCMC)

  #start the counter
  store_count <- 0

  #demographic + main effects
  x <- X[,1:(p_tot+1+p_cov)]

  #interaction effects
  z <- X[,-(1:(p_tot+1+p_cov))]

  #main effects from beta
  beta_main <- beta[1:(p_tot+1+p_cov)]

  #interaction effects from
  beta_interac <- beta[-(1:(p_tot+1+p_cov))]

  ##############
  #Start tracking time
  ##############
  start_time <- Sys.time()


  ###############
  ## MCMC chain
  ###############
  for (m in 1:MCMC){

    #Updating latent variable W
    for (i in 1:N){
      w[i] <- ifelse(Y[i]==1,
                     #if positive response, w[i] is sampled from N(0, inf)
                     rtruncnorm(1, a=0, b=Inf, mean = X[i,]%*%beta, sd = sqrt(1/lambda[i])),
                     #if negative response, w[i] is sampled from N(-inf, 0)
                     rtruncnorm(1, a=-Inf, b=0, mean = X[i,]%*%beta, sd = sqrt(1/lambda[i])))

      #Updating precision variable lambda
      lambda[i] <- rgamma(1,0.5*(v+1),0.5*(v+(w[i]-X[i,]%*%beta)^2))

    }

    #Updating interaction shrinkage parameter
    eta_interac <- array(NA,n_interac)    ## design matrix for interaction terms
    for(i in 1:n_interac)
    {
      #product of the individual shrinkage parameters
      eta_interac[i] <- eta[rr[i,1]]*eta[rr[i,2]]

    }

    #Defining the precision matrix
    sig_inv <- diag(lambda)

    #Diagonal matrix that governs shrinkage for main effects
    sig_main_inv <- diag(c(rep(0.01,(p_cov+1)),a*eta))

    #Diagonal matrix that governs shrinkage for interaction effects
    sig_int_inv <- diag(b*eta_interac*theta)

    #Covariance matrix for main effect coefficients
    sig_main <- solve(t(x)%*%sig_inv%*%x+sig_main_inv)
    mu_main <- sig_main%*%(t(x)%*%sig_inv%*%(w-z%*%beta_interac))

    #Updating main effects coefficients
    beta_main <- mvrnorm(1,mu_main,sig_main)

    #####################################################
    #Updating interaction effects coefficients
    #####################################################
    #Getting the new mean and covariance for the interaction effects
    sig_interac <- solve(t(z)%*%sig_inv%*%z+sig_int_inv)
    mu_interac <- sig_interac%*%(t(z)%*%sig_inv%*%(w-x%*%beta_main))

    #Sample the interaction effects from its conditional posterior distribution
    beta_interac <- mvrnorm(1,mu_interac,sig_interac)

    #combine the main effects and the interactions
    beta <- c(beta_main,beta_interac)

    #####################################################
    #Update the global shrinkage parameter (main effects)
    #####################################################
    #update a^2
    a_squared <- rinvgamma(n=1, shape = (p_tot+1)/2, scale = (1/kappa) + 0.5*sum(beta[(1+p_cov+1):(p_tot+1+p_cov)]^2/(eta^2)))
    a <- sqrt(a_squared)

    #update the auxiliary parameter: kappa
    kappa <- rinvgamma(n=1, shape = 1, scale = 1 + (1/a_squared))

    #a <- rgamma(1,a3+0.5*p_tot,b3+0.5*sum(beta[(1+p_cov+1):(p_tot+1+p_cov)]^2*eta))

    #####################################################
    #Update the local shrinkage parameter (main effects)
    #####################################################
    for(i in 1:p_tot){
      eta_k <- eta[-i]
      temp  <- which(as.numeric(apply(rr,1,function(x)any(x==i)))==1)
      theta_temp <- theta[temp]
      beta_temp <- beta[temp+p_tot+1+p_cov]

      #update eta_j
      eta_squared <- rinvgamma(n = 1,
                               shape = 1,
                               scale = (1/nu[i]) + 0.5*(1/a^2)*(beta[i+1+p_cov]^2) +
                                 0.5*(1/b^2)*sum(beta_temp^2/(theta_temp*eta_k)))

      eta[i] <- sqrt(eta_squared)

      #update the auxiliary parameter
      nu[i] <- rinvgamma(n=1, shape = 1, scale = 1 + (1/eta_squared))

      #eta[i] <- rgamma(1,a1+0.5*p_tot,b1+0.5*(a*beta[i+1+p_cov]^2+b*sum(theta_temp*eta_k*beta_temp^2)))
    }



    #####################################################
    #Update the global shrinkage parameter (interaction effects)
    #####################################################
    temp2 <- array(NA,n_interac)
    for(i in 1:n_interac){
      temp2[i] <- eta[rr[i,1]]*eta[rr[i,2]]*theta[i]*beta[i+p_tot+1+p_cov]^2
    }

    #Update the local shrinkage parameter (interaction effects)
    for(i in 1:n_interac){

      theta_squared <- rinvgamma(n=1, shape = 1, scale = 1/delta[i] + 0.5*beta[i+p_tot+1+p_cov]^2/(b^2*eta[rr[i,1]]^2*eta[rr[i,2]]^2))
      # theta[i] <- rgamma(1,a2+0.5,b2+0.5*b*(eta[rr[i,1]]*eta[rr[i,2]]*beta[i+p_tot+1+p_cov]^2))

      theta[i] <- sqrt(theta_squared)

      #update the auxiliary parameters
      delta[i] <- rinvgamma(n=1, shape=1, scale=1 + 1/theta_squared)


    }



    ###############
    ## storage
    ###############
    store_count <- store_count+1
    BETA[store_count,] <- beta
    ETA[store_count,] <- eta
    THETA[store_count,] <- theta
    A [store_count] <- a
    B [store_count] <- b
    LAMBDA[store_count,] <- lambda
    W[store_count,] <- w

    # print(m)
  }

  ###############
  #Track Time
  ###############
  end_time <- Sys.time()
  end_time - start_time

  ###############
  #Save the list
  ###############
  out=list(BETA=BETA,ETA=ETA,THETA=THETA,A=A,B=B)
  return(out)
}
