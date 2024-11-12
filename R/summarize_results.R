#' Summarize Results from Bayesian Shrinkage and Penalized Regression Models
#'
#' This function processes and summarizes results from various Bayesian shrinkage models and penalized logistic regression models,
#' including hierarchical and non-hierarchical global-local shrinkage priors and horseshoe priors. It calculates metrics such as MSE, sensitivity,
#' specificity, and bias for each model across simulations, providing both mean and standard deviation values for each metric.
#'
#' @param beta.hier Dataframe. Posterior samples from the Bayesian shrinkage model with a hierarchical global-local shrinkage prior.
#' @param beta.nonhier Dataframe. Posterior samples from the Bayesian shrinkage model with a non-hierarchical global-local shrinkage prior.
#' @param beta.hier_hs Dataframe. Posterior samples from the Bayesian shrinkage model with a hierarchical horseshoe prior.
#' @param horseshoe_postmean Dataframe. Posterior mean samples from the Bayesian shrinkage model with a non-hierarchical horseshoe prior.
#' @param horseshoe_95CI Dataframe. 95% credible intervals from the Bayesian shrinkage model with a non-hierarchical horseshoe prior.
#' @param l.mse Dataframe. MSE values from logistic regression with LASSO.
#' @param l.sens Dataframe. Sensitivity values from logistic regression with LASSO.
#' @param l.spec Dataframe. Specificity values from logistic regression with LASSO.
#' @param hl.mse Dataframe. MSE values from hierarchical LASSO.
#' @param hl.sens Dataframe. Sensitivity values from hierarchical LASSO.
#' @param hl.spec Dataframe. Specificity values from hierarchical LASSO.
#' @param num.main_effects Numeric. The number of main effects included in the models.
#' @param true_param Numeric vector. The true parameter values used to generate the data.
#'
#' @return A dataframe that consolidates and summarizes the results from all models, including the mean and standard deviation for each metric.
#' The final output contains metrics such as MSE, sensitivity, specificity, and credible intervals for each model across the simulations.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Combines the main and interaction effects for each model and assigns appropriate column names.
#'   \item Removes burn-in samples from each dataset and selects every 5th sample for analysis.
#'   \item Summarizes the results by calculating MSE, sensitivity, specificity, bias, and credible intervals for each model.
#'   \item Computes both the mean and standard deviation for each metric across the simulations.
#'   \item Merges the results into a single dataframe, organized by method and metric type.
#' }
#'
#' @importFrom dplyr filter mutate select
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map_df
#'
#' @examples
#' # Example usage:
#' results <- summarize_results(beta.hier, beta.nonhier, beta.hier_hs, horseshoe_postmean, horseshoe_95CI,
#'                              l.mse, l.sens, l.spec, hl.mse, hl.sens, hl.spec, num.main_effects, true_param)
#' print(results)
#'
#' @export

summarize_results <- function(beta.hier, #bayesian shrinkage w/ hier prior
                              beta.nonhier, #bayesian shrinkage w/ nonhier prior
                              beta.hier_hs, #bayesian shrinkage w/ hier horseshoe prior
                              horseshoe_postmean, #bayesian shrinkage w/ horseshoe prior
                              horseshoe_95CI,  #bayesian shrinkage w/ horseshoe prior
                              l.mse, #LR w/ LASSO
                              l.sens, #LR w/ LASSO
                              l.spec, #LR w/ LASSO
                              hl.mse, #HierLASSO
                              hl.sens, #HierLASSO
                              hl.spec, #HierLASSO
                              num.main_effects, #number of main effects
                              true_param #the true parameters

){

  ############################
  #add in columns names
  ############################
  # Main effects
  main_effects <- paste0("C", 1:num.main_effects)

  # Interaction effects
  interaction_effects <- combn(1:num.main_effects, 2, function(x) paste0("C", x[1], "C", x[2]))

  # Combine main effects and interaction effects
  cb <- c(main_effects, interaction_effects)

  # Add to beta
  colnames(beta.hier) <- c("Iter", "Intercept", cb)
  colnames(beta.nonhier) <- c("Iter", "Intercept", cb)
  colnames(beta.hier_hs) <- c("Iter", "Intercept", cb)

  ############################
  #add simulation number
  ############################
  #how many simulation
  nsim1 <- as.numeric(table(beta.hier$Iter)[1])
  nsim2 <- as.numeric(table(beta.nonhier$Iter)[1])
  nsim3 <- as.numeric(table(beta.hier_hs$Iter)[1])
  nsim4 <- dim(horseshoe_postmean)[1]/(length(cb)+1)

  #for Bayesian
  beta.hier$sim <- rep(1:nsim1, each = 50000)
  beta.nonhier$sim <- rep(1:nsim2, each = 50000)
  beta.hier_hs$sim <- rep(1:nsim3, each = 50000)

  #find the minimum
  nsim <- min(nsim1, nsim2, nsim3)

  # Remove extra simulations to ensure consistency across models
  if(nsim1 > nsim){
    beta.hier <- beta.hier %>%
      filter(sim <= nsim)
  }

  if(nsim2 > nsim){
    beta.nonhier <- beta.nonhier %>%
      filter(sim <= nsim)
  }

  if(nsim3 > nsim){
    beta.hier_hs <- beta.hier_hs %>%
      filter(sim <= nsim)
  }

  if(nsim4 > nsim){
    horseshoe_postmean <- horseshoe_postmean[1:(nsim*(length(cb)+1)),]
    horseshoe_95CI <- horseshoe_95CI[1:(nsim*(length(cb)+1)), ]
  }

  ############################
  #Remove the burn-in and only extract every 5th sample
  ############################
  #remove first 5000 iterations from each dataset
  beta.hier.burn <- beta.hier %>%
    filter(Iter > 5000)
  beta.nonhier.burn <- beta.nonhier %>%
    filter(Iter > 5000)
  beta.hier_hs.burn <- beta.hier_hs %>%
    filter(Iter > 5000)

  #only extract every 5th sample
  beta.hier.subset <- beta.hier.burn %>%
    filter(Iter%% 5 == 0)
  beta.nonhier.subset <- beta.nonhier.burn %>%
    filter(Iter%% 5 == 0)
  beta.hier_hs.subset <- beta.hier_hs.burn %>%
    filter(Iter > 5000)

  #remove Iter b/c we don't need it anymore
  beta.hier.subset$Iter <- NULL
  beta.nonhier.subset$Iter <- NULL
  beta.hier_hs.subset$Iter <- NULL

  ############################
  #true beta
  ############################
  #add names
  tp <- as.data.frame(cbind(c("Intercept", cb), true_param))
  tp$true_param <- as.numeric(tp$true_param)
  colnames(tp) <- c("metric", "value")

  ############################
  #95% HDP CI
  ############################
  beta_hier.mm <- summarize_with_hdi(data = beta.hier.subset,
                                     group = "sim",
                                     num = 0.95)

  beta_nonhier.mm <- summarize_with_hdi(data = beta.nonhier.subset,
                                        group = "sim",
                                        num = 0.95)

  beta_hier_hs.mm <- summarize_with_hdi(data = beta.hier_hs.subset,
                                        group = "sim",
                                        num = 0.95)

  ############################
  #Calculate MSE/Sens/Spec for Bayes
  ############################

  ############
  #Hierarchical Global-Local Shrinkage Prior
  ############
  # Initialize an empty data frame to store results
  hier_results_table <- data.frame()

  # Loop through each simulation
  for (i in 1:nsim) {
    sim_1 <- beta_hier.mm %>%
      filter(sim == i)

    # Calculate bshrink_stats and store the result
    result <- bshrink_stats(B0 = tp, PostB = sim_1)

    # Add the result to the results table
    hier_results_table <- rbind(hier_results_table,
                                data.frame(simulation = i, result))
  }


  #find the mean
  OG_hier_results.mean <- format_model_summary(apply(hier_results_table, 2, mean),
                                         model = "OG_hier",
                                         type = "mean")

  #find the SD
  OG_hier_results.sd <- format_model_summary(apply(hier_results_table, 2, sd),
                                       model = "OG_hier",
                                       type = "sd")


  ############
  #Non-Hierarchical Global-Local Shrinkage Prior
  ############
  # Initialize an empty data frame to store results
  nonhier_results_table <- data.frame()

  # Loop through each simulation
  for (i in 1:nsim) {
    sim_1 <- beta_nonhier.mm %>%
      filter(sim == i)

    # Calculate bshrink_stats and store the result
    result <- bshrink_stats(B0 = tp, PostB = sim_1)

    # Add the result to the results table
    nonhier_results_table <- rbind(nonhier_results_table,
                                   data.frame(simulation = i, result))
  }


  #find the mean
  OG_nonhier_results.mean <- format_model_summary(apply(nonhier_results_table, 2, mean),
                                            model = "OG_nonhier",
                                            type = "mean")

  #find the SD
  OG_nonhier_results.sd <- format_model_summary(apply(nonhier_results_table, 2, sd),
                                          model = "OG_nonhier",
                                          type = "sd")


  ############
  #Hierarchical Horseshoe Prior
  ############
  # Initialize an empty data frame to store results
  hierhs_results_table <- data.frame()

  # Loop through each simulation
  for (i in 1:nsim) {
    sim_1 <- beta_hier_hs.mm %>%
      filter(sim == i)

    # Calculate bshrink_stats and store the result
    result <- bshrink_stats(B0 = tp, PostB = sim_1)

    # Add the result to the results table
    hierhs_results_table <- rbind(hierhs_results_table,
                                  data.frame(simulation = i, result))
  }


  #find the mean
  HS_hier_results.mean <- format_model_summary(apply(hierhs_results_table, 2, mean),
                                         model = "HS_hier",
                                         type = "mean")


  #find the SD
  HS_hier_results.sd <- format_model_summary(apply(hierhs_results_table, 2, sd),
                                       model = "HS_hier",
                                       type = "sd")


  ############
  #Horseshoe Prior
  ############
  #add the sim
  horseshoe_postmean$sim <- rep(1:nsim, each = length(cb)+1)
  horseshoe_95CI$sim <- rep(1:nsim, each = length(cb)+1)

  #add the metric (coef names )
  metric <- c(cb, "Intercept")
  horseshoe_postmean$metric <- rep(metric, nsim)
  horseshoe_95CI$metric <- rep(metric, nsim)

  #convert horseshoe_95CI to longer
  horseshoe_95CI_l <- pivot_longer(
    horseshoe_95CI,
    cols = c("X1", "X2"),
    names_to = "statistic",
    values_to = "value",
    names_repair = "minimal",   # Ensures names are left as-is before renaming
    names_transform = list(statistic = function(x) ifelse(x == "X1", "hdilower", "hdiupper"))
  )

  #add the statistic
  horseshoe_postmean$statistic <- "mean"

  #edit the colnames of horseshoe_postmean
  colnames(horseshoe_postmean)[1] <- "value"

  #bind together
  horseshoe_d <- rbind(horseshoe_95CI_l, horseshoe_postmean)

  # Initialize an empty data frame to store results
  horseshoe_results_table <- data.frame()

  # Loop through each simulation
  for (i in 1:nsim) {
    sim_1 <- horseshoe_d %>%
      filter(sim == i)

    # Calculate bshrink_stats and store the result
    result <- bshrink_stats(B0 = tp, PostB = sim_1)

    # Add the result to the results table
    horseshoe_results_table <- rbind(horseshoe_results_table,
                                     data.frame(simulation = i, result))
  }


  #find the mean
  HS_nonhier_results.mean <- format_model_summary(apply(horseshoe_results_table, 2, mean),
                                            model = "HS_nonhier",
                                            type = "mean")


  #find the SD
  HS_nonhier_results.sd <- format_model_summary(apply(horseshoe_results_table, 2, sd),
                                          model = "HS_nonhier",
                                          type = "sd")



  ############################
  #Calculate MSE/Sens/Spec for Penalized Logistic Regression
  ############################
  #LASSO
  l <- data.frame(
    Method = rep("LASSO", 6),
    Category = rep(c("MSE", "sens", "spec"), 2),
    Type = c(rep("mean", 3), rep("sd", 3)),
    Value = c(apply(l.mse, 2, mean),
              apply(l.sens, 2, mean),
              apply(l.spec, 2, mean),
              apply(l.mse, 2, sd),
              apply(l.sens, 2, sd),
              apply(l.spec, 2, sd))
  )

  #HierLASSO
  h <- data.frame(
    Method = rep("HierLASSO", 6),
    Category = rep(c("MSE", "sens", "spec"), 2),
    Type = c(rep("mean", 3), rep("sd", 3)),
    Value = c(apply(hl.mse, 2, mean),
              apply(hl.sens, 2, mean),
              apply(hl.spec, 2, mean),
              apply(hl.mse, 2, sd),
              apply(hl.sens, 2, sd),
              apply(hl.spec, 2, sd))
  )

  ############################
  #Putting all the results together
  ############################
  results <- rbind(OG_hier_results.mean,
        OG_hier_results.sd,
        OG_nonhier_results.mean,
        OG_nonhier_results.sd,
        HS_hier_results.mean,
        HS_hier_results.sd,
        HS_nonhier_results.mean,
        HS_nonhier_results.sd,
        l,
        h)

  return(results)

}


