#' create function using the Weibull, reverse Weibull and Pareto distributions to convert viral suppression estimates
#' @details probability denisity functions to obtain viral suppression conversion thresholds
#' obtained from Leigh et al paper
#' @param function receives viral_sup_estimate: estimate of empirical viral suppression calculated from survey
#' reported_threshold: threshold which estimate is reported from
#' adjusted_threshold: threshold which estimate is being adjusted to 
#' distribution: distribution being explored
vls_prob_threshold <- function(viral_sup_estimate, estimate_threshold, conversion_threshold, distribution, survey,
                               empirical_estimate, empirical_lower_ci, empirical_upper_ci){
  
  if(distribution == "Reverse Weibull Johnson"){
    
    adjusted_vls_estimate <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.81)
    adjusted_vls_lower_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^1.70)
    adjusted_vls_upper_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^3.92)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Reverse Weibull Johnson")
    
    return(res)
    
  }else if(distribution == "Weibull Johnson"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.85)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.43)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^1.26)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Weibull Johnson")
    
    return(res)
    
  } else if(distribution == "Reverse Weibull Men"){
    
    adjusted_vls_estimate <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.86) 
    adjusted_vls_lower_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.78) 
    adjusted_vls_upper_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.94) 
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Reverse Weibull Men")
    
    return(res)
    
    
  } else if(distribution == "Reverse Weibull Women"){
    
    adjusted_vls_estimate <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^3.04) 
    adjusted_vls_lower_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.94) 
    adjusted_vls_upper_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^3.13) 
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Reverse Weibull Women")
    
    return(res)
    
    
  }else if(distribution == "Pareto Johnson"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.73)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^2.26)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto Johnson")
    
    return(res)
    
  }else if(distribution == "Pareto Men"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto Men")
    
    return(res)
    
  }else if(distribution == "Pareto Women"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto Women")
    
    return(res)
    
  }else if(distribution == "Weibull Men"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^1.01)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.98)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^1.03)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Weibull Men")
    
    return(res)
    
  }else{
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.88)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.85)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.90)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Weibull Women")
    
    return(res)
    
  }
  
}


#' calculate adjusted viral suppression estimates
#' Reverse Weibull distribution < 400 copies/ml
#' males
rev_weib_400_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_400,
                                        estimate_threshold = 400, conversion_threshold = 1000, 
                                        distribution = "Reverse Weibull Men", 
                                        survey = PHIA_vls_estimates_male$survey,
                                        empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000, 
                                        empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

rev_weib_400_male_phia$distribution <- rep("Reverse Weibull Men", dim(rev_weib_400_male_phia)[1])
rev_weib_400_male_phia$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_male_phia)[1])
rev_weib_400_male_phia$sex <- rep("Men", dim(rev_weib_400_male_phia)[1])

# parametrs from Johnson
rev_weib_400_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_400,
                                        estimate_threshold = 400, conversion_threshold = 1000, 
                                        distribution = "Reverse Weibull Johnson", 
                                        survey = PHIA_vls_estimates_male$survey,
                                        empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000, 
                                        empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

rev_weib_400_male$distribution <- rep("Reverse Weibull Johnson", dim(rev_weib_400_male)[1])
rev_weib_400_male$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_male)[1])
rev_weib_400_male$sex <- rep("Men", dim(rev_weib_400_male)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_200,
                                        estimate_threshold = 200, conversion_threshold = 1000, 
                                        distribution = "Reverse Weibull Men",
                                        survey = PHIA_vls_estimates_male$survey,
                                        empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

rev_weib_200_male_phia$distribution <- rep("Reverse Weibull Men", dim(rev_weib_200_male_phia)[1])
rev_weib_200_male_phia$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_male_phia)[1])
rev_weib_200_male_phia$sex <- rep("Men", dim(rev_weib_200_male_phia)[1])

# parameters from Johnson
#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_200,
                                        estimate_threshold = 200, conversion_threshold = 1000, 
                                        distribution = "Reverse Weibull Johnson",
                                        survey = PHIA_vls_estimates_male$survey,
                                        empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

rev_weib_200_male$distribution <- rep("Reverse Weibull Johnson", dim(rev_weib_200_male)[1])
rev_weib_200_male$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_male)[1])
rev_weib_200_male$sex <- rep("Men", dim(rev_weib_200_male)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_50,
                                       estimate_threshold = 50, conversion_threshold = 1000, 
                                       distribution = "Reverse Weibull Men",
                                       survey = PHIA_vls_estimates_male$survey,
                                       empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

rev_weib_50_male_phia$distribution <- rep("Reverse Weibull Men", dim(rev_weib_50_male_phia)[1])
rev_weib_50_male_phia$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_male_phia)[1])
rev_weib_50_male_phia$sex <- rep("Men", dim(rev_weib_50_male_phia)[1])

# using parameters from Johnson et al.
rev_weib_50_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_50,
                                       estimate_threshold = 50, conversion_threshold = 1000, 
                                       distribution = "Reverse Weibull Johnson",
                                       survey = PHIA_vls_estimates_male$survey,
                                       empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

rev_weib_50_male$distribution <- rep("Reverse Weibull Johnson", dim(rev_weib_50_male)[1])
rev_weib_50_male$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_male)[1])
rev_weib_50_male$sex <- rep("Men", dim(rev_weib_50_male)[1])


#' Weibull distribution < 400 copies/ml
weib_400_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_400,
                                    estimate_threshold = 400, conversion_threshold = 1000, 
                                    distribution = "Weibull Men", 
                                    survey = PHIA_vls_estimates_male$survey,
                                    empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000, 
                                    empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

weib_400_male_phia$distribution <- rep("Weibull Men", dim(weib_400_male_phia)[1])
weib_400_male_phia$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_male_phia)[1])
weib_400_male_phia$sex <- rep("Men", dim(weib_400_male_phia)[1])

# johnson
weib_400_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_400,
                                    estimate_threshold = 400, conversion_threshold = 1000, 
                                    distribution = "Weibull Johnson", 
                                    survey = PHIA_vls_estimates_male$survey,
                                    empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000, 
                                    empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

weib_400_male$distribution <- rep("Weibull Johnson", dim(weib_400_male)[1])
weib_400_male$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_male)[1])
weib_400_male$sex <- rep("Men", dim(weib_400_male)[1])


#' Weibull distribution < 200 copies/ml
weib_200_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_200,
                                    estimate_threshold = 200, conversion_threshold = 1000, 
                                    distribution = "Weibull Men",
                                    survey = PHIA_vls_estimates_male$survey,
                                    empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

weib_200_male_phia$distribution <- rep("Weibull Men", dim(weib_200_male_phia)[1])
weib_200_male_phia$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_male_phia)[1])
weib_200_male_phia$sex <- rep("Men", dim(weib_200_male_phia)[1])

# johnson
weib_200_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_200,
                                    estimate_threshold = 200, conversion_threshold = 1000, 
                                    distribution = "Weibull Johnson",
                                    survey = PHIA_vls_estimates_male$survey,
                                    empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

weib_200_male$distribution <- rep("Weibull Johnson", dim(weib_200_male)[1])
weib_200_male$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_male)[1])
weib_200_male$sex <- rep("Men", dim(weib_200_male)[1])

#' Weibull distribution < 50 copies/ml
weib_50_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_50,
                                   estimate_threshold = 50, conversion_threshold = 1000, 
                                   distribution = "Weibull Men",
                                   survey = PHIA_vls_estimates_male$survey,
                                   empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

weib_50_male_phia$distribution <- rep("Weibull Men", dim(weib_50_male_phia)[1])
weib_50_male_phia$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_male_phia)[1])
weib_50_male_phia$sex <- rep("Men", dim(weib_50_male_phia)[1])

# johnson
weib_50_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_50,
                                   estimate_threshold = 50, conversion_threshold = 1000, 
                                   distribution = "Weibull Johnson",
                                   survey = PHIA_vls_estimates_male$survey,
                                   empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

weib_50_male$distribution <- rep("Weibull Johnson", dim(weib_50_male)[1])
weib_50_male$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_male)[1])
weib_50_male$sex <- rep("Men", dim(weib_50_male)[1])


#' Pareto distribution < 400 copies/ml
pareto_400_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_400,
                                      estimate_threshold = 400, conversion_threshold = 1000, 
                                      distribution = "Pareto Men", 
                                      survey = PHIA_vls_estimates_male$survey,
                                      empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000, 
                                      empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

pareto_400_male_phia$distribution <- rep("Pareto Men", dim(pareto_400_male_phia)[1])
pareto_400_male_phia$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_male_phia)[1])
pareto_400_male_phia$sex <- rep("Men", dim(pareto_400_male_phia)[1])

# johnson
pareto_400_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_400,
                                      estimate_threshold = 400, conversion_threshold = 1000, 
                                      distribution = "Pareto Johnson", 
                                      survey = PHIA_vls_estimates_male$survey,
                                      empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000, 
                                      empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

pareto_400_male$distribution <- rep("Pareto Johnson", dim(pareto_400_male)[1])
pareto_400_male$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_male)[1])
pareto_400_male$sex <- rep("Men", dim(pareto_400_male)[1])

#' Pareto distribution < 200 copies/ml
pareto_200_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_200,
                                      estimate_threshold = 200, conversion_threshold = 1000, 
                                      distribution = "Pareto Men",
                                      survey = PHIA_vls_estimates_male$survey,
                                      empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

pareto_200_male_phia$distribution <- rep("Pareto Men", dim(pareto_200_male_phia)[1])
pareto_200_male_phia$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_male_phia)[1])
pareto_200_male_phia$sex <- rep("Men", dim(pareto_200_male_phia)[1])

# johnson
pareto_200_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_200,
                                      estimate_threshold = 200, conversion_threshold = 1000, 
                                      distribution = "Pareto Johnson",
                                      survey = PHIA_vls_estimates_male$survey,
                                      empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

pareto_200_male$distribution <- rep("Pareto Johnson", dim(pareto_200_male)[1])
pareto_200_male$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_male)[1])
pareto_200_male$sex <- rep("Men", dim(pareto_200_male)[1])

#' Pareto distribution < 50 copies/ml
pareto_50_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_50,
                                     estimate_threshold = 50, conversion_threshold = 1000, 
                                     distribution = "Pareto Men",
                                     survey = PHIA_vls_estimates_male$survey,
                                     empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

pareto_50_male_phia$distribution <- rep("Pareto Men", dim(pareto_50_male_phia)[1])
pareto_50_male_phia$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_male_phia)[1])
pareto_50_male_phia$sex <- rep("Men", dim(pareto_50_male_phia)[1])

# johnson
pareto_50_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_50,
                                     estimate_threshold = 50, conversion_threshold = 1000, 
                                     distribution = "Pareto Johnson",
                                     survey = PHIA_vls_estimates_male$survey,
                                     empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

pareto_50_male$distribution <- rep("Pareto Johnson", dim(pareto_50_male)[1])
pareto_50_male$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_male)[1])
pareto_50_male$sex <- rep("Men", dim(pareto_50_male)[1])


#' calculate the root mean squared error between empirical and adjusted estimates
sqrt(mean((rev_weib_50_male_phia$empirical_estimate - rev_weib_50_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_male_phia$empirical_estimate - weib_50_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_male_phia$empirical_estimate - pareto_50_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_50_male$empirical_estimate - rev_weib_50_male$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_male$empirical_estimate - weib_50_male$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_male$empirical_estimate - pareto_50_male$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_200_male_phia$empirical_estimate - rev_weib_200_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_male_phia$empirical_estimate - weib_200_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_male_phia$empirical_estimate - pareto_200_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_200_male$empirical_estimate - rev_weib_200_male$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_male$empirical_estimate - weib_200_male$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_male$empirical_estimate - pareto_200_male$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_400_male_phia$empirical_estimate - rev_weib_400_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_male_phia$empirical_estimate - weib_400_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_male_phia$empirical_estimate - pareto_400_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_400_male$empirical_estimate - rev_weib_400_male$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_male$empirical_estimate - weib_400_male$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_male$empirical_estimate - pareto_400_male$adjusted_vls_estimate)^2))

# avearge bias
mean(rev_weib_50_male_phia$adjusted_vls_estimate - rev_weib_50_male_phia$empirical_estimate)
mean(weib_50_male_phia$adjusted_vls_estimate - weib_50_male_phia$empirical_estimate)
mean(pareto_50_male_phia$adjusted_vls_estimate - pareto_50_male_phia$empirical_estimate)
mean(rev_weib_50_male$adjusted_vls_estimate - rev_weib_50_male$empirical_estimate)
mean(weib_50_male$adjusted_vls_estimate - weib_50_male$empirical_estimate)
mean(pareto_50_male$adjusted_vls_estimate - pareto_50_male$empirical_estimate)

mean(rev_weib_200_male_phia$adjusted_vls_estimate - rev_weib_200_male_phia$empirical_estimate)
mean(weib_200_male_phia$adjusted_vls_estimate - weib_200_male_phia$empirical_estimate)
mean(pareto_200_male_phia$adjusted_vls_estimate - pareto_200_male_phia$empirical_estimate)
mean(rev_weib_200_male$adjusted_vls_estimate - rev_weib_200_male$empirical_estimate)
mean(weib_200_male$adjusted_vls_estimate - weib_200_male$empirical_estimate)
mean(pareto_200_male$adjusted_vls_estimate - pareto_200_male$empirical_estimate)

mean(rev_weib_400_male_phia$adjusted_vls_estimate - rev_weib_400_male_phia$empirical_estimate)
mean(weib_400_male_phia$adjusted_vls_estimate - weib_400_male_phia$empirical_estimate)
mean(pareto_400_male_phia$adjusted_vls_estimate - pareto_400_male_phia$empirical_estimate)
mean(rev_weib_400_male$adjusted_vls_estimate - rev_weib_400_male$empirical_estimate)
mean(weib_400_male$adjusted_vls_estimate - weib_400_male$empirical_estimate)
mean(pareto_400_male$adjusted_vls_estimate - pareto_400_male$empirical_estimate)


# mean absolute error
mean(abs(rev_weib_50_male_phia$adjusted_vls_estimate - rev_weib_50_male_phia$empirical_estimate))
mean(abs(weib_50_male_phia$adjusted_vls_estimate - weib_50_male_phia$empirical_estimate))
mean(abs(pareto_50_male_phia$adjusted_vls_estimate - pareto_50_male_phia$empirical_estimate))
mean(abs(rev_weib_50_male$adjusted_vls_estimate - rev_weib_50_male$empirical_estimate))
mean(abs(weib_50_male$adjusted_vls_estimate - weib_50_male$empirical_estimate))
mean(abs(pareto_50_male$adjusted_vls_estimate - pareto_50_male$empirical_estimate))

mean(abs(rev_weib_200_male_phia$adjusted_vls_estimate - rev_weib_200_male_phia$empirical_estimate))
mean(abs(weib_200_male_phia$adjusted_vls_estimate - weib_200_male_phia$empirical_estimate))
mean(abs(pareto_200_male_phia$adjusted_vls_estimate - pareto_200_male_phia$empirical_estimate))
mean(abs(rev_weib_200_male$adjusted_vls_estimate - rev_weib_200_male$empirical_estimate))
mean(abs(weib_200_male$adjusted_vls_estimate - weib_200_male$empirical_estimate))
mean(abs(pareto_200_male$adjusted_vls_estimate - pareto_200_male$empirical_estimate))

mean(abs(rev_weib_400_male_phia$adjusted_vls_estimate - rev_weib_400_male_phia$empirical_estimate))
mean(abs(weib_400_male_phia$adjusted_vls_estimate - weib_400_male_phia$empirical_estimate))
mean(abs(pareto_400_male_phia$adjusted_vls_estimate - pareto_400_male_phia$empirical_estimate))
mean(abs(rev_weib_400_male$adjusted_vls_estimate - rev_weib_400_male$empirical_estimate))
mean(abs(weib_400_male$adjusted_vls_estimate - weib_400_male$empirical_estimate))
mean(abs(pareto_400_male$adjusted_vls_estimate - pareto_400_male$empirical_estimate))

#' calculate adjusted viral suppression estimates
#' Reverse Weibull distribution < 400 copies/ml
#' females
rev_weib_400_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_400,
                                          estimate_threshold = 400, conversion_threshold = 1000, 
                                          distribution = "Reverse Weibull Women", 
                                          survey = PHIA_vls_estimates_female$survey,
                                          empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000, 
                                          empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

rev_weib_400_female_phia$distribution <- rep("Reverse Weibull Women", dim(rev_weib_400_female_phia)[1])
rev_weib_400_female_phia$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_female_phia)[1])
rev_weib_400_female_phia$sex <- rep("Women", dim(rev_weib_400_female_phia)[1])

# johnson
rev_weib_400_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_400,
                                          estimate_threshold = 400, conversion_threshold = 1000, 
                                          distribution = "Reverse Weibull Johnson", 
                                          survey = PHIA_vls_estimates_female$survey,
                                          empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000, 
                                          empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

rev_weib_400_female$distribution <- rep("Reverse Weibull Johnson", dim(rev_weib_400_female)[1])
rev_weib_400_female$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_female)[1])
rev_weib_400_female$sex <- rep("Women", dim(rev_weib_400_female)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_200,
                                          estimate_threshold = 200, conversion_threshold = 1000, 
                                          distribution = "Reverse Weibull Women",
                                          survey = PHIA_vls_estimates_female$survey,
                                          empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                          empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

rev_weib_200_female_phia$distribution <- rep("Reverse Weibull Women", dim(rev_weib_200_female_phia)[1])
rev_weib_200_female_phia$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_female_phia)[1])
rev_weib_200_female_phia$sex <- rep("Women", dim(rev_weib_200_female_phia)[1])

# johnson
rev_weib_200_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_200,
                                          estimate_threshold = 200, conversion_threshold = 1000, 
                                          distribution = "Reverse Weibull Johnson",
                                          survey = PHIA_vls_estimates_female$survey,
                                          empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                          empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

rev_weib_200_female$distribution <- rep("Reverse Weibull Johnson", dim(rev_weib_200_female)[1])
rev_weib_200_female$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_female)[1])
rev_weib_200_female$sex <- rep("Women", dim(rev_weib_200_female)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_50,
                                         estimate_threshold = 50, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull Women",
                                         survey = PHIA_vls_estimates_female$survey,
                                         empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

rev_weib_50_female_phia$distribution <- rep("Reverse Weibull Women", dim(rev_weib_50_female_phia)[1])
rev_weib_50_female_phia$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_female_phia)[1])
rev_weib_50_female_phia$sex <- rep("Women", dim(rev_weib_50_female_phia)[1])

# johnson
rev_weib_50_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_50,
                                         estimate_threshold = 50, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull Johnson",
                                         survey = PHIA_vls_estimates_female$survey,
                                         empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

rev_weib_50_female$distribution <- rep("Reverse Weibull Johnson", dim(rev_weib_50_female)[1])
rev_weib_50_female$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_female)[1])
rev_weib_50_female$sex <- rep("Women", dim(rev_weib_50_female)[1])

#' Weibull distribution < 400 copies/ml
weib_400_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_400,
                                      estimate_threshold = 400, conversion_threshold = 1000, 
                                      distribution = "Weibull Women", 
                                      survey = PHIA_vls_estimates_female$survey,
                                      empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000, 
                                      empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

weib_400_female_phia$distribution <- rep("Weibull Women", dim(weib_400_female_phia)[1])
weib_400_female_phia$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_female_phia)[1])
weib_400_female_phia$sex <- rep("Women", dim(weib_400_female_phia)[1])

# johnson
weib_400_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_400,
                                      estimate_threshold = 400, conversion_threshold = 1000, 
                                      distribution = "Weibull Johnson", 
                                      survey = PHIA_vls_estimates_female$survey,
                                      empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000, 
                                      empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

weib_400_female$distribution <- rep("Weibull Johnson", dim(weib_400_female)[1])
weib_400_female$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_female)[1])
weib_400_female$sex <- rep("Women", dim(weib_400_female)[1])

#' Weibull distribution < 200 copies/ml
weib_200_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_200,
                                      estimate_threshold = 200, conversion_threshold = 1000, 
                                      distribution = "Weibull Women",
                                      survey = PHIA_vls_estimates_female$survey,
                                      empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

weib_200_female_phia$distribution <- rep("Weibull Women", dim(weib_200_female_phia)[1])
weib_200_female_phia$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_female_phia)[1])
weib_200_female_phia$sex <- rep("Women", dim(weib_200_female_phia)[1])

# johnson
weib_200_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_200,
                                      estimate_threshold = 200, conversion_threshold = 1000, 
                                      distribution = "Weibull Johnson",
                                      survey = PHIA_vls_estimates_female$survey,
                                      empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

weib_200_female$distribution <- rep("Weibull Johnson", dim(weib_200_female)[1])
weib_200_female$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_female)[1])
weib_200_female$sex <- rep("Women", dim(weib_200_female)[1])

#' Weibull distribution < 50 copies/ml
weib_50_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_50,
                                     estimate_threshold = 50, conversion_threshold = 1000, 
                                     distribution = "Weibull Women",
                                     survey = PHIA_vls_estimates_female$survey,
                                     empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

weib_50_female_phia$distribution <- rep("Weibull Women", dim(weib_50_female_phia)[1])
weib_50_female_phia$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_female_phia)[1])
weib_50_female_phia$sex <- rep("Women", dim(weib_50_female_phia)[1])

# johnson
weib_50_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_50,
                                     estimate_threshold = 50, conversion_threshold = 1000, 
                                     distribution = "Weibull Johnson",
                                     survey = PHIA_vls_estimates_female$survey,
                                     empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

weib_50_female$distribution <- rep("Weibull Johnson", dim(weib_50_female)[1])
weib_50_female$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_female)[1])
weib_50_female$sex <- rep("Women", dim(weib_50_female)[1])


#' Pareto distribution < 400 copies/ml
pareto_400_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_400,
                                        estimate_threshold = 400, conversion_threshold = 1000, 
                                        distribution = "Pareto Women", 
                                        survey = PHIA_vls_estimates_female$survey,
                                        empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000, 
                                        empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

pareto_400_female_phia$distribution <- rep("Pareto Women", dim(pareto_400_female_phia)[1])
pareto_400_female_phia$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_female_phia)[1])
pareto_400_female_phia$sex <- rep("Women", dim(pareto_400_female_phia)[1])

# johnson
pareto_400_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_400,
                                        estimate_threshold = 400, conversion_threshold = 1000, 
                                        distribution = "Pareto Johnson", 
                                        survey = PHIA_vls_estimates_female$survey,
                                        empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000, 
                                        empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

pareto_400_female$distribution <- rep("Pareto Johnson", dim(pareto_400_female)[1])
pareto_400_female$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_female)[1])
pareto_400_female$sex <- rep("Women", dim(pareto_400_female)[1])


#' Pareto distribution < 200 copies/ml
pareto_200_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_200,
                                        estimate_threshold = 200, conversion_threshold = 1000, 
                                        distribution = "Pareto Women",
                                        survey = PHIA_vls_estimates_female$survey,
                                        empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

pareto_200_female_phia$distribution <- rep("Pareto Women", dim(pareto_200_female_phia)[1])
pareto_200_female_phia$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_female_phia)[1])
pareto_200_female_phia$sex <- rep("Women", dim(pareto_200_female_phia)[1])

# johnson
pareto_200_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_200,
                                        estimate_threshold = 200, conversion_threshold = 1000, 
                                        distribution = "Pareto Johnson",
                                        survey = PHIA_vls_estimates_female$survey,
                                        empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

pareto_200_female$distribution <- rep("Pareto Johnson", dim(pareto_200_female)[1])
pareto_200_female$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_female)[1])
pareto_200_female$sex <- rep("Women", dim(pareto_200_female)[1])

#' Pareto distribution < 50 copies/ml
pareto_50_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_50,
                                       estimate_threshold = 50, conversion_threshold = 1000, 
                                       distribution = "Pareto Women",
                                       survey = PHIA_vls_estimates_female$survey,
                                       empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

pareto_50_female_phia$distribution <- rep("Pareto Women", dim(pareto_50_female_phia)[1])
pareto_50_female_phia$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_female_phia)[1])
pareto_50_female_phia$sex <- rep("Women", dim(pareto_50_female_phia)[1])

# johnson
pareto_50_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_50,
                                       estimate_threshold = 50, conversion_threshold = 1000, 
                                       distribution = "Pareto Johnson",
                                       survey = PHIA_vls_estimates_female$survey,
                                       empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

pareto_50_female$distribution <- rep("Pareto Johnson", dim(pareto_50_female)[1])
pareto_50_female$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_female)[1])
pareto_50_female$sex <- rep("Women", dim(pareto_50_female)[1])


#' calculate the root mean squared error between empirical and adjusted estimates
sqrt(mean((rev_weib_50_female_phia$empirical_estimate - rev_weib_50_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_female_phia$empirical_estimate - weib_50_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_female_phia$empirical_estimate - pareto_50_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_50_female$empirical_estimate - rev_weib_50_female$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_female$empirical_estimate - weib_50_female$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_female$empirical_estimate - pareto_50_female$adjusted_vls_estimate)^2))


sqrt(mean((rev_weib_200_female_phia$empirical_estimate - rev_weib_200_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_female_phia$empirical_estimate - weib_200_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_female_phia$empirical_estimate - pareto_200_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_200_female$empirical_estimate - rev_weib_200_female$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_female$empirical_estimate - weib_200_female$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_female$empirical_estimate - pareto_200_female$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_400_female_phia$empirical_estimate - rev_weib_400_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_female_phia$empirical_estimate - weib_400_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_female_phia$empirical_estimate - pareto_400_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_400_female$empirical_estimate - rev_weib_400_female$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_female$empirical_estimate - weib_400_female$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_female$empirical_estimate - pareto_400_female$adjusted_vls_estimate)^2))

# average bias
mean(rev_weib_50_female_phia$adjusted_vls_estimate - rev_weib_50_female_phia$empirical_estimate)
mean(weib_50_female_phia$adjusted_vls_estimate - weib_50_female_phia$empirical_estimate)
mean(pareto_50_female_phia$adjusted_vls_estimate - pareto_50_female_phia$empirical_estimate)
mean(rev_weib_50_female$adjusted_vls_estimate - rev_weib_50_female$empirical_estimate)
mean(weib_50_female$adjusted_vls_estimate - weib_50_female$empirical_estimate)
mean(pareto_50_female$adjusted_vls_estimate - pareto_50_female$empirical_estimate)


mean(rev_weib_200_female_phia$adjusted_vls_estimate - rev_weib_200_female_phia$empirical_estimate)
mean(weib_200_female_phia$adjusted_vls_estimate - weib_200_female_phia$empirical_estimate)
mean(pareto_200_female_phia$adjusted_vls_estimate - pareto_200_female_phia$empirical_estimate)
mean(rev_weib_200_female$adjusted_vls_estimate - rev_weib_200_female$empirical_estimate)
mean(weib_200_female$adjusted_vls_estimate - weib_200_female$empirical_estimate)
mean(pareto_200_female$adjusted_vls_estimate - pareto_200_female$empirical_estimate)

mean(rev_weib_400_female_phia$adjusted_vls_estimate - rev_weib_400_female_phia$empirical_estimate)
mean(weib_400_female_phia$adjusted_vls_estimate - weib_400_female_phia$empirical_estimate)
mean(pareto_400_female_phia$adjusted_vls_estimate - pareto_400_female_phia$empirical_estimate)
mean(rev_weib_400_female$adjusted_vls_estimate - rev_weib_400_female$empirical_estimate)
mean(weib_400_female$adjusted_vls_estimate - weib_400_female$empirical_estimate)
mean(pareto_400_female$adjusted_vls_estimate - pareto_400_female$empirical_estimate)

# mean absolute error
mean(abs(rev_weib_50_female_phia$adjusted_vls_estimate - rev_weib_50_female_phia$empirical_estimate))
mean(abs(weib_50_female_phia$adjusted_vls_estimate - weib_50_female_phia$empirical_estimate))
mean(abs(pareto_50_female_phia$adjusted_vls_estimate - pareto_50_female_phia$empirical_estimate))
mean(abs(rev_weib_50_female$adjusted_vls_estimate - rev_weib_50_female$empirical_estimate))
mean(abs(weib_50_female$adjusted_vls_estimate - weib_50_female$empirical_estimate))
mean(abs(pareto_50_female$adjusted_vls_estimate - pareto_50_female$empirical_estimate))


mean(abs(rev_weib_200_female_phia$adjusted_vls_estimate - rev_weib_200_female_phia$empirical_estimate))
mean(abs(weib_200_female_phia$adjusted_vls_estimate - weib_200_female_phia$empirical_estimate))
mean(abs(pareto_200_female_phia$adjusted_vls_estimate - pareto_200_female_phia$empirical_estimate))
mean(abs(rev_weib_200_female$adjusted_vls_estimate - rev_weib_200_female$empirical_estimate))
mean(abs(weib_200_female$adjusted_vls_estimate - weib_200_female$empirical_estimate))
mean(abs(pareto_200_female$adjusted_vls_estimate - pareto_200_female$empirical_estimate))

mean(abs(rev_weib_400_female_phia$adjusted_vls_estimate - rev_weib_400_female_phia$empirical_estimate))
mean(abs(weib_400_female_phia$adjusted_vls_estimate - weib_400_female_phia$empirical_estimate))
mean(abs(pareto_400_female_phia$adjusted_vls_estimate - pareto_400_female_phia$empirical_estimate))
mean(abs(rev_weib_400_female$adjusted_vls_estimate - rev_weib_400_female$empirical_estimate))
mean(abs(weib_400_female$adjusted_vls_estimate - weib_400_female$empirical_estimate))
mean(abs(pareto_400_female$adjusted_vls_estimate - pareto_400_female$empirical_estimate))


# plot
sex_validation <- bind_rows(rev_weib_50_male, rev_weib_50_female, 
                            rev_weib_50_male_phia, rev_weib_50_female_phia,
                            rev_weib_200_male, rev_weib_200_female,
                            rev_weib_200_male_phia, rev_weib_200_female_phia,
                            rev_weib_400_male, rev_weib_400_female,
                            rev_weib_400_male_phia, rev_weib_400_female_phia,
                            weib_50_male, weib_50_female,
                            weib_50_male_phia, weib_50_female_phia,
                            weib_200_male, weib_200_female,
                            weib_200_male_phia, weib_200_female_phia,
                            weib_400_male, weib_400_female,
                            weib_400_male_phia, weib_400_female_phia,
                            pareto_50_male, pareto_50_female,
                            pareto_50_male_phia, pareto_50_female_phia,
                            pareto_200_male, pareto_200_female,
                            pareto_200_male_phia, pareto_200_female_phia,
                            pareto_400_male, pareto_400_female,
                            pareto_400_male_phia, pareto_400_female_phia,)

#' update country names to survey names
sex_validation <- sex_validation %>%
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
                                   TRUE ~ survey),
         dist_names = case_when(distribution == "Pareto Johnson" ~ "Pareto",
                                distribution == "Pareto Men" ~ "Pareto",
                                distribution == "Pareto Women" ~ "Pareto",
                                distribution == "Reverse Weibull Johnson" ~ "Reverse Weibull",
                                distribution == "Reverse Weibull Men" ~ "Reverse Weibull",
                                distribution == "Reverse Weibull Women" ~ "Reverse Weibull",
                                distribution == "Weibull Johnson" ~ "Weibull",
                                distribution == "Weibull Men" ~ "Weibull",
                                distribution == "Weibull Women" ~ "Weibull"),
         param = case_when(distribution == "Pareto Johnson" ~ "Johnson et al",
                                distribution == "Pareto Men" ~ "PHIA surveys",
                                distribution == "Pareto Women" ~ "PHIA surveys",
                                distribution == "Reverse Weibull Johnson" ~ "Johnson et al",
                                distribution == "Reverse Weibull Men" ~ "PHIA surveys",
                                distribution == "Reverse Weibull Women" ~ "PHIA surveys",
                                distribution == "Weibull Johnson" ~ "Johnson et al",
                                distribution == "Weibull Men" ~ "PHIA surveys",
                                distribution == "Weibull Women" ~ "PHIA surveys"),
         estimate_threshold = forcats::fct_relevel(estimate_threshold, 
                                                   c("<50 copies/mL", "<200 copies/mL", "<400 copies/mL")),
         dist_names = forcats::fct_relevel(dist_names,
                                             c("Reverse Weibull", "Weibull", "Pareto")),
         country_names = forcats::fct_relevel(country_names,
                                              c("Côte d'Ivoire (2017-18)","Cameroon (2017-18)",
                                                "Nigeria (2018)","Uganda (2016-17)","Zimbabwe (2015-16)",
                                                "Tanzania (2016-17)","Ethiopia (2017-18)","Lesotho (2016-17)",
                                                "Zambia (2016)","Mozambique (2021-22)","Rwanda (2018-19)",
                                                "Zimbabwe (2019-20)","Kenya (2018-19)","Malawi (2015-16)",
                                                "Namibia (2017)","Eswatini (2016-17)","Lesotho (2019-20)",
                                                "Eswatini (2021)","Zambia (2021)","Malawi (2020-21)", "Botswana (2021)")))

table(sex_validation$country_names, useNA = "always")


## scatter plot
sex_validation %>%
  #filter(estimate_threshold == "<400 copies/mL" & param == "Johnson et al") %>%
  filter(estimate_threshold == "<400 copies/mL" & param == "PHIA surveys") %>%
  mutate(sex = forcats::fct_relevel(sex, c("Women","Men"))) %>%
ggplot(aes(x = empirical_estimate, y = adjusted_vls_estimate, color = country_names)) +
  geom_point(size = 2) +
  geom_errorbar(aes(y = adjusted_vls_estimate, ymin = adjusted_vls_lower_ci, 
                    ymax = adjusted_vls_upper_ci), linewidth = 0.75, alpha = 0.3) + 
  geom_errorbarh(aes(xmin = empirical_lower_ci, 
                     xmax = empirical_upper_ci), linewidth = 0.75, alpha = 0.3) + 
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  coord_fixed(xlim = c(0.5,1), ylim = c(0.5,1)) +
  geom_abline(intercept = 0, slope = 1, col = "darkgrey", linewidth = 0.9, linetype = "dashed") +
  theme_bw(base_size = 13) +
  scale_color_manual(values = c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", 
                                "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70",
                                "maroon", "orchid1", "darkturquoise", "darkorange4", "brown",
                                "cyan", "red", "gold4", "orange","blue")) +
  facet_grid(sex~dist_names) + 
  labs(x = "% viral load suppression (PHIA survey)",
       y = "% viral load suppression (adjusted)", color = "",
       title = "Adjustment from <400 to \u22641000 copies/mL and shape parameters from calibration to PHIA surveys") + 
  theme(legend.position = "bottom",
        axis.text = element_text(size = rel(1), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(0.9), family = "sans"),
        panel.grid = element_blank(), 
        panel.spacing.x = unit(1, "lines"),
        plot.title = element_text(size = rel(1.0), family = "sans", face = "bold"))
