#' create function using the Weibull, reverse Weibull ,Pareto , Frechet, gamma and lognrmal distributions to convert viral suppression estimates
#' @details probability density functions to obtain viral suppression conversion thresholds
#' obtained from Leigh et al paper
#' @param function receives viral_sup_estimate: estimate of empirical viral suppression calculated from survey
#' reported_threshold: threshold which estimate is reported from
#' adjusted_threshold: threshold which estimate is being adjusted to 
#' distribution: distribution being explored
vls_prob_threshold <- function(viral_sup_estimate, estimate_threshold, conversion_threshold, distribution, survey,
                               empirical_estimate, empirical_lower_ci, empirical_upper_ci){
  
  if(distribution == "Reverse Weibull Johnson LF et al"){
    
    adjusted_vls_estimate <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.81)
    adjusted_vls_lower_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^1.70)
    adjusted_vls_upper_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^3.92)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Reverse Weibull Johnson LF et al")
    
    return(res)
    
  }else if(distribution == "Weibull Johnson LF et al"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.85)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.43)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^1.26)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Weibull Johnson LF et al")
    
    return(res)
    
  } else if(distribution == "Reverse Weibull updated"){
   
    adjusted_vls_estimate <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.98) 
    adjusted_vls_lower_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.93) 
    adjusted_vls_upper_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^3.02) 
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Reverse Weibull updated")
    
    return(res)
    
     
  }else if(distribution == "Pareto Johnson LF et al"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.73)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^2.26)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto Johnson LF et al")
    
    return(res)
    
  }else if(distribution == "Pareto updated"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto updated")
    
    return(res)
    
  }else{
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^0.91)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^0.90) 
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^0.92) 
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Weibull updated")
    
    return(res)
    
  }
  
}

#' Reverse Weibull distribution < 400 copies/ml
rev_weib_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400, 
                                   estimate_threshold = 400, conversion_threshold = 1000, 
                                   distribution = "Reverse Weibull Johnson LF et al",
                                   survey = PHIA_vls_estimates$survey,
                                   empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

rev_weib_400$distribution <- rep("Reverse\nWeibull\nJohnson et al.\n(shape = 2.81)", dim(rev_weib_400)[1])
rev_weib_400$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                                   estimate_threshold = 200, conversion_threshold = 1000, 
                                   distribution = "Reverse Weibull Johnson LF et al",
                                   survey = PHIA_vls_estimates$survey,
                                   empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

rev_weib_200$distribution <- rep("Reverse\nWeibull\nJohnson et al.\n(shape = 2.81)", dim(rev_weib_200)[1])
rev_weib_200$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                                  estimate_threshold = 50, conversion_threshold = 1000, 
                                  distribution = "Reverse Weibull Johnson LF et al",
                                  survey = PHIA_vls_estimates$survey, 
                                  empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                  empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                  empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

rev_weib_50$distribution <- rep("Reverse\nWeibull\nJohnson et al.\n(shape = 2.81)", dim(rev_weib_50)[1])
rev_weib_50$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50)[1])


#' Weibull distribution < 400 copies/ml
weib_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                               estimate_threshold = 400, conversion_threshold = 1000, 
                               distribution = "Weibull Johnson LF et al",
                               survey = PHIA_vls_estimates$survey,
                               empirical_estimate = PHIA_vls_estimates$estimate_1000,
                               empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                               empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

weib_400$distribution <- rep("Weibull\nJohnson et al.\n(shape = 0.85)", dim(weib_400)[1])
weib_400$estimate_threshold <- rep("<400 copies/mL", dim(weib_400)[1])

#' Weibull distribution < 200 copies/ml
weib_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                               estimate_threshold = 200, conversion_threshold = 1000, 
                               distribution = "Weibull Johnson LF et al",
                               survey = PHIA_vls_estimates$survey,
                               empirical_estimate = PHIA_vls_estimates$estimate_1000,
                               empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                               empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

weib_200$distribution <- rep("Weibull\nJohnson et al.\n(shape = 0.85)", dim(weib_200)[1])
weib_200$estimate_threshold <- rep("<200 copies/mL", dim(weib_200)[1])

#' Weibull distribution < 50 copies/ml
weib_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                              estimate_threshold = 50, conversion_threshold = 1000, 
                              distribution = "Weibull Johnson LF et al",
                              survey = PHIA_vls_estimates$survey, 
                              empirical_estimate = PHIA_vls_estimates$estimate_1000,
                              empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                              empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

weib_50$distribution <- rep("Weibull\nJohnson et al.\n(shape = 0.85)", dim(weib_50)[1])
weib_50$estimate_threshold <- rep("<50 copies/mL", dim(weib_50)[1])


#' Pareto distribution < 400 copies/ml
pareto_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                               estimate_threshold = 400, conversion_threshold = 1000, 
                               distribution = "Pareto Johnson LF et al",
                               survey = PHIA_vls_estimates$survey,
                               empirical_estimate = PHIA_vls_estimates$estimate_1000,
                               empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                               empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

pareto_400$distribution <- rep("Pareto\nJohnson et al.\n(shape = 1.73)", dim(pareto_400)[1])
pareto_400$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400)[1])

#' Pareto distribution < 200 copies/ml
pareto_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                               estimate_threshold = 200, conversion_threshold = 1000, 
                               distribution = "Pareto Johnson LF et al",
                               survey = PHIA_vls_estimates$survey,
                               empirical_estimate = PHIA_vls_estimates$estimate_1000,
                               empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                               empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

pareto_200$distribution <- rep("Pareto\nJohnson et al.\n(shape = 1.73)", dim(pareto_200)[1])
pareto_200$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200)[1])

#' Pareto distribution < 50 copies/ml
pareto_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                              estimate_threshold = 50, conversion_threshold = 1000, 
                              distribution = "Pareto Johnson LF et al",
                              survey = PHIA_vls_estimates$survey, 
                              empirical_estimate = PHIA_vls_estimates$estimate_1000,
                              empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                              empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

pareto_50$distribution <- rep("Pareto\nJohnson et al.\n(shape = 1.73)", dim(pareto_50)[1])
pareto_50$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50)[1])



#' Reverse Weibull distribution with updated parameters < 400 copies/ml
new_rev_weib_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                 distribution = "Reverse Weibull updated", 
                                 survey = PHIA_vls_estimates$survey,
                                 empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                 empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                 empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

new_rev_weib_400$distribution <- rep("Reverse\nWeibull\nPHIAs\n(shape = 2.98)", dim(new_rev_weib_400)[1])
new_rev_weib_400$estimate_threshold <- rep("<400 copies/mL", dim(new_rev_weib_400)[1])

#' Reverse Weibull distribution with updated parameters < 200 copies/ml
new_rev_weib_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                 distribution = "Reverse Weibull updated",
                                 survey = PHIA_vls_estimates$survey,
                                 empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                 empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                 empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

new_rev_weib_200$distribution <- rep("Reverse\nWeibull\nPHIAs\n(shape = 2.98)", dim(new_rev_weib_200)[1])
new_rev_weib_200$estimate_threshold <- rep("<200 copies/mL", dim(new_rev_weib_200)[1])

#' Reverse Weibull distribution with updated parameters < 50 copies/ml
new_rev_weib_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                distribution = "Reverse Weibull updated", 
                                survey = PHIA_vls_estimates$survey,
                                empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

new_rev_weib_50$distribution <- rep("Reverse\nWeibull\nPHIAs\n(shape = 2.98)", dim(new_rev_weib_50)[1])
new_rev_weib_50$estimate_threshold <- rep("<50 copies/mL", dim(new_rev_weib_50)[1])


#' Weibull distribution with updated parameters < 400 copies/ml
new_weib_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                                   estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "Weibull updated",
                                   survey = PHIA_vls_estimates$survey, 
                                   empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

new_weib_400$distribution <- rep("Weibull\nPHIAs\n(shape = 0.91)", dim(new_weib_400)[1])
new_weib_400$estimate_threshold <- rep("<400 copies/mL", dim(new_weib_400)[1])

#' Weibull distribution with updated parameters < 200 copies/ml
new_weib_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                                   estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "Weibull updated", 
                                   survey = PHIA_vls_estimates$survey, empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

new_weib_200$distribution <- rep("Weibull\nPHIAs\n(shape = 0.91)", dim(new_weib_200)[1])
new_weib_200$estimate_threshold <- rep("<200 copies/mL", dim(new_weib_200)[1])

#' Weibull distribution with updated parameters < 50 copies/ml
new_weib_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                                  estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "Weibull updated",
                                  survey = PHIA_vls_estimates$survey,
                                  empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000, 
                                  empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

new_weib_50$distribution <- rep("Weibull\nPHIAs\n(shape = 0.91)", dim(new_weib_50)[1])
new_weib_50$estimate_threshold <- rep("<50 copies/mL", dim(new_weib_50)[1])


#' Pareto distribution with updated parameters < 400 copies/ml
new_pareto_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                                   estimate_threshold = 400, conversion_threshold = 1000, 
                                   distribution = "Pareto updated",
                                   survey = PHIA_vls_estimates$survey, 
                                   empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

new_pareto_400$distribution <- rep("Pareto\nPHIAs\n(shape = 1.20)", dim(new_pareto_400)[1])
new_pareto_400$estimate_threshold <- rep("<400 copies/mL", dim(new_pareto_400)[1])

#' Pareto distribution with updated parameters < 200 copies/ml
new_pareto_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                                   estimate_threshold = 200, conversion_threshold = 1000, 
                                   distribution = "Pareto updated", 
                                   survey = PHIA_vls_estimates$survey, empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

new_pareto_200$distribution <- rep("Pareto\nPHIAs\n(shape = 1.20)", dim(new_pareto_200)[1])
new_pareto_200$estimate_threshold <- rep("<200 copies/mL", dim(new_pareto_200)[1])

#' Pareto distribution with updated parameters < 50 copies/ml
new_pareto_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                                  estimate_threshold = 50, conversion_threshold = 1000, 
                                  distribution = "Pareto updated",
                                  survey = PHIA_vls_estimates$survey,
                                  empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                  empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000, 
                                  empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

new_pareto_50$distribution <- rep("Pareto\nPHIAs\n(shape = 1.20)", dim(new_pareto_50)[1])
new_pareto_50$estimate_threshold <- rep("<50 copies/mL", dim(new_pareto_50)[1])


#' calculate the root mean squared error between empirical and adjusted estimates
sqrt(mean((rev_weib_50$empirical_estimate - rev_weib_50$adjusted_vls_estimate)^2))
sqrt(mean((new_rev_weib_50$empirical_estimate - new_rev_weib_50$adjusted_vls_estimate)^2))
sqrt(mean((weib_50$empirical_estimate - weib_50$adjusted_vls_estimate)^2))
sqrt(mean((new_weib_50$empirical_estimate - new_weib_50$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50$empirical_estimate - pareto_50$adjusted_vls_estimate)^2))
sqrt(mean((new_pareto_50$empirical_estimate - new_pareto_50$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_200$empirical_estimate - rev_weib_200$adjusted_vls_estimate)^2))
sqrt(mean((new_rev_weib_200$empirical_estimate - new_rev_weib_200$adjusted_vls_estimate)^2))
sqrt(mean((weib_200$empirical_estimate - weib_200$adjusted_vls_estimate)^2))
sqrt(mean((new_weib_200$empirical_estimate - new_weib_200$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200$empirical_estimate - pareto_200$adjusted_vls_estimate)^2))
sqrt(mean((new_pareto_200$empirical_estimate - new_pareto_200$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_400$empirical_estimate - rev_weib_400$adjusted_vls_estimate)^2))
sqrt(mean((new_rev_weib_400$empirical_estimate - new_rev_weib_400$adjusted_vls_estimate)^2))
sqrt(mean((weib_400$empirical_estimate - weib_400$adjusted_vls_estimate)^2))
sqrt(mean((new_weib_400$empirical_estimate - new_weib_400$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400$empirical_estimate - pareto_400$adjusted_vls_estimate)^2))
sqrt(mean((new_pareto_400$empirical_estimate - new_pareto_400$adjusted_vls_estimate)^2))


# average bias
#' calculate the average difference between adjusted and empirical estimates
mean(rev_weib_50$adjusted_vls_estimate - rev_weib_50$empirical_estimate)
mean(new_rev_weib_50$adjusted_vls_estimate - new_rev_weib_50$empirical_estimate)
mean(weib_50$adjusted_vls_estimate - weib_50$empirical_estimate)
mean(new_weib_50$adjusted_vls_estimate - new_weib_50$empirical_estimate)
mean(pareto_50$adjusted_vls_estimate - pareto_50$empirical_estimate)
mean(new_pareto_50$adjusted_vls_estimate - new_pareto_50$empirical_estimate)

mean(rev_weib_200$adjusted_vls_estimate - rev_weib_200$empirical_estimate)
mean(new_rev_weib_200$adjusted_vls_estimate - new_rev_weib_200$empirical_estimate)
mean(weib_200$adjusted_vls_estimate - weib_200$empirical_estimate)
mean(new_weib_200$adjusted_vls_estimate - new_weib_200$empirical_estimate)
mean(pareto_200$adjusted_vls_estimate - pareto_200$empirical_estimate)
mean(new_pareto_200$adjusted_vls_estimate - new_pareto_200$empirical_estimate)

mean(rev_weib_400$adjusted_vls_estimate - rev_weib_400$empirical_estimate)
mean(new_rev_weib_400$adjusted_vls_estimate - new_rev_weib_400$empirical_estimate)
mean(weib_400$adjusted_vls_estimate - weib_400$empirical_estimate)
mean(new_weib_400$adjusted_vls_estimate - new_weib_400$empirical_estimate)
mean(pareto_400$adjusted_vls_estimate - pareto_400$empirical_estimate)
mean(new_pareto_400$adjusted_vls_estimate - new_pareto_400$empirical_estimate)


# mean absolute error
#' calculate the average difference between adjusted and empirical estimates
mean(abs(rev_weib_50$adjusted_vls_estimate - rev_weib_50$empirical_estimate))
mean(abs(new_rev_weib_50$adjusted_vls_estimate - new_rev_weib_50$empirical_estimate))
mean(abs(weib_50$adjusted_vls_estimate - weib_50$empirical_estimate))
mean(abs(new_weib_50$adjusted_vls_estimate - new_weib_50$empirical_estimate))
mean(abs(pareto_50$adjusted_vls_estimate - pareto_50$empirical_estimate))
mean(abs(new_pareto_50$adjusted_vls_estimate - new_pareto_50$empirical_estimate))

mean(abs(rev_weib_200$adjusted_vls_estimate - rev_weib_200$empirical_estimate))
mean(abs(new_rev_weib_200$adjusted_vls_estimate - new_rev_weib_200$empirical_estimate))
mean(abs(weib_200$adjusted_vls_estimate - weib_200$empirical_estimate))
mean(abs(new_weib_200$adjusted_vls_estimate - new_weib_200$empirical_estimate))
mean(abs(pareto_200$adjusted_vls_estimate - pareto_200$empirical_estimate))
mean(abs(new_pareto_200$adjusted_vls_estimate - new_pareto_200$empirical_estimate))

mean(abs(rev_weib_400$adjusted_vls_estimate - rev_weib_400$empirical_estimate))
mean(abs(new_rev_weib_400$adjusted_vls_estimate - new_rev_weib_400$empirical_estimate))
mean(abs(weib_400$adjusted_vls_estimate - weib_400$empirical_estimate))
mean(abs(new_weib_400$adjusted_vls_estimate - new_weib_400$empirical_estimate))
mean(abs(pareto_400$adjusted_vls_estimate - pareto_400$empirical_estimate))
mean(abs(new_pareto_400$adjusted_vls_estimate - new_pareto_400$empirical_estimate))

#' PLOT (code for Fgure 1 and Figure S1)
#' scatter plot showing relationship between empirical and adjusted viral suppression estimates
#' combine all estimates derived in adjusted_vls_estimates script
all_distributions_updated <- bind_rows(rev_weib_50, rev_weib_200, rev_weib_400,
                                       weib_50, weib_200, weib_400,
                                       pareto_50, pareto_200, pareto_400,
                                       new_rev_weib_50, new_rev_weib_200, new_rev_weib_400,
                                       new_weib_50, new_weib_200, new_weib_400,
                                       new_pareto_50, new_pareto_200, new_pareto_400)

#' update country names to survey names
all_distributions_updated <- all_distributions_updated %>%
  mutate(country_names = case_when(survey == "BAIS_2021" ~ "Botswana (2021)",
                                   survey == "CAMPHIA_2017" ~ "Cameroon (2017-18)",
                                   survey == "CIPHIA_2017" ~ "Côte d'Ivoire (2017-18)",
                                   survey == "EPHIA_2017" ~ "Ethiopia (2017-18)",
                                   survey == "INSIDA_2021" ~ "Mozambique (2021-22)",
                                   survey == "KENPHIA_2018" ~ "Kenya (2018-19)",
                                   survey == "LePHIA_2016" ~ "Lesotho (2016-17)",
                                   survey == "LePHIA_2020" ~ "Lesotho (2019-20)",
                                   survey == "MPHIA_2015" ~ "Malawi (2015-16)",
                                   survey == "MPHIA_2020" ~ "Malawi (2020-21)",
                                   survey == "NAIIS_2018" ~ "Nigeria (2018)",
                                   survey == "NAMPHIA_2017" ~ "Namibia (2017)",
                                   survey == "RPHIA_2018" ~ "Rwanda (2018-19)",
                                   survey == "SHIMS2_2016" ~ "Eswatini (2016-17)",
                                   survey == "SHIMS3_2021"  ~ "Eswatini (2021)",
                                   survey == "THIS_2016" ~ "Tanzania (2016-17)",
                                   survey == "UPHIA_2016" ~ "Uganda (2016-17)",
                                   survey == "ZAMPHIA_2016" ~ "Zambia (2016)",
                                   survey == "ZAMPHIA_2021" ~ "Zambia (2021)",
                                   survey == "ZIMPHIA_2015" ~ "Zimbabwe (2015-16)" ,
                                   survey == "ZIMPHIA_2020" ~ "Zimbabwe (2019-20)",
                                   TRUE ~ survey))


#' order estimates threshold from lowest to highest
all_distributions_updated$estimate_threshold <-
  factor(all_distributions_updated$estimate_threshold, levels = c("<50 copies/mL", "<200 copies/mL", "<400 copies/mL"))

all_distributions_updated$country_names <- factor(all_distributions_updated$country_names,
                                          levels = c("Côte d'Ivoire (2017-18)","Cameroon (2017-18)",
                                                     "Nigeria (2018)","Uganda (2016-17)","Zimbabwe (2015-16)",
                                                     "Tanzania (2016-17)","Ethiopia (2017-18)","Lesotho (2016-17)",
                                                     "Zambia (2016)","Mozambique (2021-22)","Rwanda (2018-19)",
                                                     "Zimbabwe (2019-20)","Kenya (2018-19)","Malawi (2015-16)",
                                                     "Namibia (2017)","Eswatini (2016-17)","Lesotho (2019-20)",
                                                     "Eswatini (2021)","Zambia (2021)","Malawi (2020-21)", "Botswana (2021)"))

## scatter plot
all_distributions_updated %>%
  filter(estimate_threshold == "<400 copies/mL") %>%
  mutate(model = case_when(distribution == "Reverse\nWeibull\nJohnson et al.\n(shape = 2.81)" ~ "Reverse Weibull",
                           distribution == "Reverse\nWeibull\nPHIAs\n(shape = 2.98)" ~ "Reverse Weibull",
                           distribution == "Weibull\nJohnson et al.\n(shape = 0.85)" ~ "Weibull",
                           distribution == "Weibull\nPHIAs\n(shape = 0.91)" ~ "Weibull",
                           distribution == "Pareto\nJohnson et al.\n(shape = 1.73)" ~ "Pareto",
                           distribution == "Pareto\nPHIAs\n(shape = 1.20)" ~ "Pareto"),
         parameters = case_when(distribution == "Reverse\nWeibull\nJohnson et al.\n(shape = 2.81)" ~ "Johnson et al\nshape parameter",
                                distribution == "Reverse\nWeibull\nPHIAs\n(shape = 2.98)" ~ "PHIA survey\nrecalibrated",
                                distribution == "Weibull\nJohnson et al.\n(shape = 0.85)" ~ "Johnson et al\nshape parameter",
                                distribution == "Weibull\nPHIAs\n(shape = 0.91)" ~ "PHIA survey\nrecalibrated",
                                distribution == "Pareto\nJohnson et al.\n(shape = 1.73)" ~ "Johnson et al\nshape parameter",
                                distribution == "Pareto\nPHIAs\n(shape = 1.20)" ~ "PHIA survey\nrecalibrated"),
         model = forcats::fct_relevel(model, c("Reverse Weibull","Weibull","Pareto"))) %>%
ggplot(aes(x = empirical_estimate, y = adjusted_vls_estimate, color = country_names)) +
  geom_point(size = 2) +
  geom_errorbar(aes(y = adjusted_vls_estimate, ymin = adjusted_vls_lower_ci, 
                    ymax = adjusted_vls_upper_ci), linewidth = 0.75, alpha = 0.3) +
  geom_errorbarh(aes(xmin = empirical_lower_ci, 
                    xmax = empirical_upper_ci), linewidth = 0.75, alpha = 0.3) +
  scale_x_continuous(limits = c(0.6, 1), labels = scales::percent) +
  scale_y_continuous(limits = c(0.6, 1), labels = scales::percent) +
  coord_fixed() +
  geom_abline(intercept = 0, slope = 1, col = "darkgrey", linewidth = 0.9, linetype = "dashed") +
  theme_bw(base_size = 13) +
  scale_color_manual(values = c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", 
                                "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70",
                                "maroon", "orchid1", "darkturquoise", "darkorange4", "brown",
                                "cyan", "red", "gold4", "orange", "blue")) +
  facet_grid(parameters~model) + 
  labs(x = "% viral suppression (PHIA survey)", 
       y = "% viral suppression (adjusted)", color = "", 
       title = "Adjustment from <400 to \u22641000 copies/mL") + 
  theme(legend.position = "bottom",
        axis.text = element_text(size = rel(1), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(0.9), family = "sans"),
        panel.grid = element_blank(), 
        panel.spacing.x = unit(1, "lines"),
        plot.title = element_text(size = rel(1.0), family = "sans", face = "bold"))


# summarise differences in predicted proportions and observed 
all_distributions_updated %>%
  filter(distribution == "Weibull\nJohnson et al.\n(shape = 0.85)" &
           estimate_threshold == "<50 copies/mL") %>%
  filter(adjusted_vls_estimate < empirical_estimate) %>%
  select(survey, estimate_threshold, adjusted_vls_estimate, empirical_estimate)


# assess surveys where CIs dont overlap
table(all_distributions_updated$distribution)

all_distributions_updated %>% filter(grepl('Johnson', distribution)) %>%
  filter(adjusted_vls_estimate < empirical_lower_ci|adjusted_vls_estimate > empirical_upper_ci) %>%
  filter(adjusted_vls_lower_ci < empirical_lower_ci|adjusted_vls_lower_ci > empirical_upper_ci) %>%
  filter(adjusted_vls_upper_ci < empirical_lower_ci|adjusted_vls_upper_ci > empirical_upper_ci) %>%
  select(country_names,estimate_threshold,distribution)

all_distributions_updated %>% filter(grepl('Johnson', distribution)) %>%
  filter(!between(adjusted_vls_estimate, empirical_lower_ci, empirical_upper_ci)) %>%
  filter(!between(adjusted_vls_lower_ci, empirical_lower_ci, empirical_upper_ci)) %>%
  filter(!between(adjusted_vls_upper_ci, empirical_lower_ci, empirical_upper_ci)) %>%
  select(country_names,estimate_threshold,distribution)




