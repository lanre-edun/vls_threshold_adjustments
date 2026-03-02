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
    
  }else if(distribution == "Pareto Johnson"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.73)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^2.26)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto Johnson")
    
    return(res)
    
  }else if(distribution == "Reverse Weibull 15-24"){
    
    adjusted_vls_estimate <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.46) 
    adjusted_vls_lower_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.35) 
    adjusted_vls_upper_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.58) 
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Reverse Weibull 15-24")
    
    return(res)
    
    
  } else if(distribution == "Reverse Weibull 25-34"){
    
    adjusted_vls_estimate <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.52) 
    adjusted_vls_lower_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.37) 
    adjusted_vls_upper_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.66) 
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Reverse Weibull 25-34")
    
    return(res)
    
    
  }else if(distribution == "Reverse Weibull 35-44"){
    
    adjusted_vls_estimate <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.84) 
    adjusted_vls_lower_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.68) 
    adjusted_vls_upper_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^3.00) 
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Reverse Weibull 35-44")
    
    return(res)
    
    
  }else if(distribution == "Reverse Weibull 45-54"){
    
    adjusted_vls_estimate <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^3.81) 
    adjusted_vls_lower_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^3.64) 
    adjusted_vls_upper_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^3.99) 
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Reverse Weibull 45-54")
    
    return(res)
    
    
  }else if(distribution == "Reverse Weibull 55+"){
    
    adjusted_vls_estimate <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^4.32) 
    adjusted_vls_lower_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^4.10) 
    adjusted_vls_upper_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^4.55) 
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Reverse Weibull 55+")
    
    return(res)
    
    
  }else if(distribution == "Pareto 15-24"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto 15-24")
    
    return(res)
    
  }else if(distribution == "Pareto 25-34"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto 25-34")
    
    return(res)
    
  }else if(distribution == "Pareto 35-44"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto 35-44")
    
    return(res)
    
  }else if(distribution == "Pareto 45-54"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto 45-54")
    
    return(res)
    
  }else if(distribution == "Pareto 55+"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto 55+")
    
    return(res)
    
  }else if(distribution == "Weibull 15-24"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.55)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.53)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.56)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Weibull 15-24")
    
    return(res)
    
  }else if(distribution == "Weibull 25-34"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.98)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.97)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.99)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Weibull 25-34")
    
    return(res)
    
  }else if(distribution == "Weibull 35-44"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.98)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.97)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.99)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Weibull 35-44")
    
    return(res)
    
  }else if(distribution == "Weibull 45-54"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.99)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.98)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^1.02)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Weibull 45-54")
    
    return(res)
    
  }else{
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^1.00)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.98)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^1.04)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Weibull 55+")
    
    return(res)
    
  }
  
}


#' calculate adjusted viral suppression estimates
#' Reverse Weibull distribution < 400 copies/ml
#' 15-24
rev_weib_400_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull 15-24", 
                                         survey = PHIA_vls_estimates_15_24$survey,
                                         empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

rev_weib_400_15_24_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_400_15_24_phia)[1])
rev_weib_400_15_24_phia$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_15_24_phia)[1])
rev_weib_400_15_24_phia$agecat <- rep("15-24", dim(rev_weib_400_15_24_phia)[1])
rev_weib_400_15_24_phia$source <- rep("PHIA surveys", dim(rev_weib_400_15_24_phia)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull 15-24",
                                         survey = PHIA_vls_estimates_15_24$survey,
                                         empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

rev_weib_200_15_24_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_200_15_24_phia)[1])
rev_weib_200_15_24_phia$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_15_24_phia)[1])
rev_weib_200_15_24_phia$agecat <- rep("15-24", dim(rev_weib_200_15_24_phia)[1])
rev_weib_200_15_24_phia$source <- rep("PHIA surveys", dim(rev_weib_200_15_24_phia)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "Reverse Weibull 15-24",
                                        survey = PHIA_vls_estimates_15_24$survey,
                                        empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

rev_weib_50_15_24_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_50_15_24_phia)[1])
rev_weib_50_15_24_phia$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_15_24_phia)[1])
rev_weib_50_15_24_phia$agecat <- rep("15-24", dim(rev_weib_50_15_24_phia)[1])
rev_weib_50_15_24_phia$source <- rep("PHIA surveys", dim(rev_weib_50_15_24_phia)[1])

#' Weibull distribution < 400 copies/ml
weib_400_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_400,
                                     estimate_threshold = 400, conversion_threshold = 1000, 
                                     distribution = "Weibull 15-24", 
                                     survey = PHIA_vls_estimates_15_24$survey,
                                     empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000, 
                                     empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

weib_400_15_24_phia$distribution <- rep("Weibull", dim(weib_400_15_24_phia)[1])
weib_400_15_24_phia$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_15_24_phia)[1])
weib_400_15_24_phia$agecat <- rep("15-24", dim(weib_400_15_24_phia)[1])
weib_400_15_24_phia$source <- rep("PHIA surveys", dim(weib_400_15_24_phia)[1])

#' Weibull distribution < 200 copies/ml
weib_200_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_200,
                                     estimate_threshold = 200, conversion_threshold = 1000, 
                                     distribution = "Weibull 15-24",
                                     survey = PHIA_vls_estimates_15_24$survey,
                                     empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

weib_200_15_24_phia$distribution <- rep("Weibull", dim(weib_200_15_24_phia)[1])
weib_200_15_24_phia$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_15_24_phia)[1])
weib_200_15_24_phia$agecat <- rep("15-24", dim(weib_200_15_24_phia)[1])
weib_200_15_24_phia$source <- rep("PHIA surveys", dim(weib_200_15_24_phia)[1])

#' Weibull distribution < 50 copies/ml
weib_50_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_50,
                                    estimate_threshold = 50, conversion_threshold = 1000, 
                                    distribution = "Weibull 15-24",
                                    survey = PHIA_vls_estimates_15_24$survey,
                                    empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

weib_50_15_24_phia$distribution <- rep("Weibull", dim(weib_50_15_24_phia)[1])
weib_50_15_24_phia$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_15_24_phia)[1])
weib_50_15_24_phia$agecat <- rep("15-24", dim(weib_50_15_24_phia)[1])
weib_50_15_24_phia$source <- rep("PHIA surveys", dim(weib_50_15_24_phia)[1])

#' Pareto distribution < 400 copies/ml
pareto_400_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "Pareto 15-24", 
                                       survey = PHIA_vls_estimates_15_24$survey,
                                       empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000, 
                                       empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

pareto_400_15_24_phia$distribution <- rep("Pareto", dim(pareto_400_15_24_phia)[1])
pareto_400_15_24_phia$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_15_24_phia)[1])
pareto_400_15_24_phia$agecat <- rep("15-24", dim(pareto_400_15_24_phia)[1])
pareto_400_15_24_phia$source <- rep("PHIA surveys", dim(pareto_400_15_24_phia)[1])

#' Pareto distribution < 200 copies/ml
pareto_200_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "Pareto 15-24",
                                       survey = PHIA_vls_estimates_15_24$survey,
                                       empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

pareto_200_15_24_phia$distribution <- rep("Pareto", dim(pareto_200_15_24_phia)[1])
pareto_200_15_24_phia$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_15_24_phia)[1])
pareto_200_15_24_phia$agecat <- rep("15-24", dim(pareto_200_15_24_phia)[1])
pareto_200_15_24_phia$source <- rep("PHIA surveys", dim(pareto_200_15_24_phia)[1])

#' Pareto distribution < 50 copies/ml
pareto_50_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "Pareto 15-24",
                                      survey = PHIA_vls_estimates_15_24$survey,
                                      empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

pareto_50_15_24_phia$distribution <- rep("Pareto", dim(pareto_50_15_24_phia)[1])
pareto_50_15_24_phia$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_15_24_phia)[1])
pareto_50_15_24_phia$agecat <- rep("15-24", dim(pareto_50_15_24_phia)[1])
pareto_50_15_24_phia$source <- rep("PHIA surveys", dim(pareto_50_15_24_phia)[1])


#' calculate adjusted viral suppression estimates
#' Reverse Weibull distribution < 400 copies/ml
#' 25-34
rev_weib_400_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull 25-34", 
                                         survey = PHIA_vls_estimates_25_34$survey,
                                         empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

rev_weib_400_25_34_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_400_25_34_phia)[1])
rev_weib_400_25_34_phia$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_25_34_phia)[1])
rev_weib_400_25_34_phia$agecat <- rep("25-34", dim(rev_weib_400_25_34_phia)[1])
rev_weib_400_25_34_phia$source <- rep("PHIA surveys", dim(rev_weib_400_25_34_phia)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull 25-34",
                                         survey = PHIA_vls_estimates_25_34$survey,
                                         empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

rev_weib_200_25_34_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_200_25_34_phia)[1])
rev_weib_200_25_34_phia$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_25_34_phia)[1])
rev_weib_200_25_34_phia$agecat <- rep("25-34", dim(rev_weib_200_25_34_phia)[1])
rev_weib_200_25_34_phia$source <- rep("PHIA surveys", dim(rev_weib_200_25_34_phia)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "Reverse Weibull 25-34",
                                        survey = PHIA_vls_estimates_25_34$survey,
                                        empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

rev_weib_50_25_34_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_50_25_34_phia)[1])
rev_weib_50_25_34_phia$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_25_34_phia)[1])
rev_weib_50_25_34_phia$agecat <- rep("25-34", dim(rev_weib_50_25_34_phia)[1])
rev_weib_50_25_34_phia$source <- rep("PHIA surveys", dim(rev_weib_50_25_34_phia)[1])


#' Weibull distribution < 400 copies/ml
weib_400_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_400,
                                     estimate_threshold = 400, conversion_threshold = 1000, 
                                     distribution = "Weibull 25-34", 
                                     survey = PHIA_vls_estimates_25_34$survey,
                                     empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000, 
                                     empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

weib_400_25_34_phia$distribution <- rep("Weibull", dim(weib_400_25_34_phia)[1])
weib_400_25_34_phia$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_25_34_phia)[1])
weib_400_25_34_phia$agecat <- rep("25-34", dim(weib_400_25_34_phia)[1])
weib_400_25_34_phia$source <- rep("PHIA surveys", dim(weib_400_25_34_phia)[1])

#' Weibull distribution < 200 copies/ml
weib_200_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_200,
                                     estimate_threshold = 200, conversion_threshold = 1000, 
                                     distribution = "Weibull 25-34",
                                     survey = PHIA_vls_estimates_25_34$survey,
                                     empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

weib_200_25_34_phia$distribution <- rep("Weibull", dim(weib_200_25_34_phia)[1])
weib_200_25_34_phia$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_25_34_phia)[1])
weib_200_25_34_phia$agecat <- rep("25-34", dim(weib_200_25_34_phia)[1])
weib_200_25_34_phia$source <- rep("PHIA surveys", dim(weib_200_25_34_phia)[1])

#' Weibull distribution < 50 copies/ml
weib_50_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_50,
                                    estimate_threshold = 50, conversion_threshold = 1000, 
                                    distribution = "Weibull 25-34",
                                    survey = PHIA_vls_estimates_25_34$survey,
                                    empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

weib_50_25_34_phia$distribution <- rep("Weibull", dim(weib_50_25_34_phia)[1])
weib_50_25_34_phia$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_25_34_phia)[1])
weib_50_25_34_phia$agecat <- rep("25-34", dim(weib_50_25_34_phia)[1])
weib_50_25_34_phia$source <- rep("PHIA surveys", dim(weib_50_25_34_phia)[1])


#' Pareto distribution < 400 copies/ml
pareto_400_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "Pareto 25-34", 
                                       survey = PHIA_vls_estimates_25_34$survey,
                                       empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000, 
                                       empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

pareto_400_25_34_phia$distribution <- rep("Pareto", dim(pareto_400_25_34_phia)[1])
pareto_400_25_34_phia$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_25_34_phia)[1])
pareto_400_25_34_phia$agecat <- rep("25-34", dim(pareto_400_25_34_phia)[1])
pareto_400_25_34_phia$source <- rep("PHIA surveys", dim(pareto_400_25_34_phia)[1])

#' Pareto distribution < 200 copies/ml
pareto_200_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "Pareto 25-34",
                                       survey = PHIA_vls_estimates_25_34$survey,
                                       empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

pareto_200_25_34_phia$distribution <- rep("Pareto", dim(pareto_200_25_34_phia)[1])
pareto_200_25_34_phia$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_25_34_phia)[1])
pareto_200_25_34_phia$agecat <- rep("25-34", dim(pareto_200_25_34_phia)[1])
pareto_200_25_34_phia$source <- rep("PHIA surveys", dim(pareto_200_25_34_phia)[1])

#' Pareto distribution < 50 copies/ml
pareto_50_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "Pareto 25-34",
                                      survey = PHIA_vls_estimates_25_34$survey,
                                      empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

pareto_50_25_34_phia$distribution <- rep("Pareto", dim(pareto_50_25_34_phia)[1])
pareto_50_25_34_phia$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_25_34_phia)[1])
pareto_50_25_34_phia$agecat <- rep("25-34", dim(pareto_50_25_34_phia)[1])
pareto_50_25_34_phia$source <- rep("PHIA surveys", dim(pareto_50_25_34_phia)[1])


#' calculate adjusted viral suppression estimates
#' Reverse Weibull distribution < 400 copies/ml
#' 35-44
rev_weib_400_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull 35-44", 
                                         survey = PHIA_vls_estimates_35_44$survey,
                                         empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

rev_weib_400_35_44_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_400_35_44_phia)[1])
rev_weib_400_35_44_phia$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_35_44_phia)[1])
rev_weib_400_35_44_phia$agecat <- rep("35-44", dim(rev_weib_400_35_44_phia)[1])
rev_weib_400_35_44_phia$source <- rep("PHIA surveys", dim(rev_weib_400_35_44_phia)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull 35-44",
                                         survey = PHIA_vls_estimates_35_44$survey,
                                         empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

rev_weib_200_35_44_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_200_35_44_phia)[1])
rev_weib_200_35_44_phia$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_35_44_phia)[1])
rev_weib_200_35_44_phia$agecat <- rep("35-44", dim(rev_weib_200_35_44_phia)[1])
rev_weib_200_35_44_phia$source <- rep("PHIA surveys", dim(rev_weib_200_35_44_phia)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "Reverse Weibull 35-44",
                                        survey = PHIA_vls_estimates_35_44$survey,
                                        empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

rev_weib_50_35_44_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_50_35_44_phia)[1])
rev_weib_50_35_44_phia$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_35_44_phia)[1])
rev_weib_50_35_44_phia$agecat <- rep("35-44", dim(rev_weib_50_35_44_phia)[1])
rev_weib_50_35_44_phia$source <- rep("PHIA surveys", dim(rev_weib_50_35_44_phia)[1])


#' Weibull distribution < 400 copies/ml
weib_400_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_400,
                                     estimate_threshold = 400, conversion_threshold = 1000, 
                                     distribution = "Weibull 35-44", 
                                     survey = PHIA_vls_estimates_35_44$survey,
                                     empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000, 
                                     empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

weib_400_35_44_phia$distribution <- rep("Weibull", dim(weib_400_35_44_phia)[1])
weib_400_35_44_phia$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_35_44_phia)[1])
weib_400_35_44_phia$agecat <- rep("35-44", dim(weib_400_35_44_phia)[1])
weib_400_35_44_phia$source <- rep("PHIA surveys", dim(weib_400_35_44_phia)[1])

#' Weibull distribution < 200 copies/ml
weib_200_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_200,
                                     estimate_threshold = 200, conversion_threshold = 1000, 
                                     distribution = "Weibull 35-44",
                                     survey = PHIA_vls_estimates_35_44$survey,
                                     empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

weib_200_35_44_phia$distribution <- rep("Weibull", dim(weib_200_35_44_phia)[1])
weib_200_35_44_phia$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_35_44_phia)[1])
weib_200_35_44_phia$agecat <- rep("35-44", dim(weib_200_35_44_phia)[1])
weib_200_35_44_phia$source <- rep("PHIA surveys", dim(weib_200_35_44_phia)[1])

#' Weibull distribution < 50 copies/ml
weib_50_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_50,
                                    estimate_threshold = 50, conversion_threshold = 1000, 
                                    distribution = "Weibull 35-44",
                                    survey = PHIA_vls_estimates_35_44$survey,
                                    empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

weib_50_35_44_phia$distribution <- rep("Weibull", dim(weib_50_35_44_phia)[1])
weib_50_35_44_phia$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_35_44_phia)[1])
weib_50_35_44_phia$agecat <- rep("35-44", dim(weib_50_35_44_phia)[1])
weib_50_35_44_phia$source <- rep("PHIA surveys", dim(weib_50_35_44_phia)[1])


#' Pareto distribution < 400 copies/ml
pareto_400_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "Pareto 35-44", 
                                       survey = PHIA_vls_estimates_35_44$survey,
                                       empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000, 
                                       empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

pareto_400_35_44_phia$distribution <- rep("Pareto", dim(pareto_400_35_44_phia)[1])
pareto_400_35_44_phia$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_35_44_phia)[1])
pareto_400_35_44_phia$agecat <- rep("35-44", dim(pareto_400_35_44_phia)[1])
pareto_400_35_44_phia$source <- rep("PHIA surveys", dim(pareto_400_35_44_phia)[1])

#' Pareto distribution < 200 copies/ml
pareto_200_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "Pareto 35-44",
                                       survey = PHIA_vls_estimates_35_44$survey,
                                       empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

pareto_200_35_44_phia$distribution <- rep("Pareto", dim(pareto_200_35_44_phia)[1])
pareto_200_35_44_phia$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_35_44_phia)[1])
pareto_200_35_44_phia$agecat <- rep("35-44", dim(pareto_200_35_44_phia)[1])
pareto_200_35_44_phia$source <- rep("PHIA surveys", dim(pareto_200_35_44_phia)[1])

#' Pareto distribution < 50 copies/ml
pareto_50_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "Pareto 35-44",
                                      survey = PHIA_vls_estimates_35_44$survey,
                                      empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

pareto_50_35_44_phia$distribution <- rep("Pareto", dim(pareto_50_35_44_phia)[1])
pareto_50_35_44_phia$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_35_44_phia)[1])
pareto_50_35_44_phia$agecat <- rep("35-44", dim(pareto_50_35_44_phia)[1])
pareto_50_35_44_phia$source <- rep("PHIA surveys", dim(pareto_50_35_44_phia)[1])


#' calculate adjusted viral suppression estimates
#' Reverse Weibull distribution < 400 copies/ml
#' 45-54
rev_weib_400_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull 45-54", 
                                         survey = PHIA_vls_estimates_45_54$survey,
                                         empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

rev_weib_400_45_54_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_400_45_54_phia)[1])
rev_weib_400_45_54_phia$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_45_54_phia)[1])
rev_weib_400_45_54_phia$agecat <- rep("45-54", dim(rev_weib_400_45_54_phia)[1])
rev_weib_400_45_54_phia$source <- rep("PHIA surveys", dim(rev_weib_400_45_54_phia)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull 45-54",
                                         survey = PHIA_vls_estimates_45_54$survey,
                                         empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

rev_weib_200_45_54_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_200_45_54_phia)[1])
rev_weib_200_45_54_phia$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_45_54_phia)[1])
rev_weib_200_45_54_phia$agecat <- rep("45-54", dim(rev_weib_200_45_54_phia)[1])
rev_weib_200_45_54_phia$source <- rep("PHIA surveys", dim(rev_weib_200_45_54_phia)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "Reverse Weibull 45-54",
                                        survey = PHIA_vls_estimates_45_54$survey,
                                        empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

rev_weib_50_45_54_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_50_45_54_phia)[1])
rev_weib_50_45_54_phia$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_45_54_phia)[1])
rev_weib_50_45_54_phia$agecat <- rep("45-54", dim(rev_weib_50_45_54_phia)[1])
rev_weib_50_45_54_phia$source <- rep("PHIA surveys", dim(rev_weib_50_45_54_phia)[1])


#' Weibull distribution < 400 copies/ml
weib_400_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_400,
                                     estimate_threshold = 400, conversion_threshold = 1000, 
                                     distribution = "Weibull 45-54", 
                                     survey = PHIA_vls_estimates_45_54$survey,
                                     empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000, 
                                     empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

weib_400_45_54_phia$distribution <- rep("Weibull", dim(weib_400_45_54_phia)[1])
weib_400_45_54_phia$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_45_54_phia)[1])
weib_400_45_54_phia$agecat <- rep("45-54", dim(weib_400_45_54_phia)[1])
weib_400_45_54_phia$source <- rep("PHIA surveys", dim(weib_400_45_54_phia)[1])

#' Weibull distribution < 200 copies/ml
weib_200_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_200,
                                     estimate_threshold = 200, conversion_threshold = 1000, 
                                     distribution = "Weibull 45-54",
                                     survey = PHIA_vls_estimates_45_54$survey,
                                     empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

weib_200_45_54_phia$distribution <- rep("Weibull", dim(weib_200_45_54_phia)[1])
weib_200_45_54_phia$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_45_54_phia)[1])
weib_200_45_54_phia$agecat <- rep("45-54", dim(weib_200_45_54_phia)[1])
weib_200_45_54_phia$source <- rep("PHIA surveys", dim(weib_200_45_54_phia)[1])

#' Weibull distribution < 50 copies/ml
weib_50_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_50,
                                    estimate_threshold = 50, conversion_threshold = 1000, 
                                    distribution = "Weibull 45-54",
                                    survey = PHIA_vls_estimates_45_54$survey,
                                    empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

weib_50_45_54_phia$distribution <- rep("Weibull", dim(weib_50_45_54_phia)[1])
weib_50_45_54_phia$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_45_54_phia)[1])
weib_50_45_54_phia$agecat <- rep("45-54", dim(weib_50_45_54_phia)[1])
weib_50_45_54_phia$source <- rep("PHIA surveys", dim(weib_50_45_54_phia)[1])


#' Pareto distribution < 400 copies/ml
pareto_400_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "Pareto 45-54", 
                                       survey = PHIA_vls_estimates_45_54$survey,
                                       empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000, 
                                       empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

pareto_400_45_54_phia$distribution <- rep("Pareto", dim(pareto_400_45_54_phia)[1])
pareto_400_45_54_phia$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_45_54_phia)[1])
pareto_400_45_54_phia$agecat <- rep("45-54", dim(pareto_400_45_54_phia)[1])
pareto_400_45_54_phia$source <- rep("PHIA surveys", dim(pareto_400_45_54_phia)[1])

#' Pareto distribution < 200 copies/ml
pareto_200_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "Pareto 45-54",
                                       survey = PHIA_vls_estimates_45_54$survey,
                                       empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

pareto_200_45_54_phia$distribution <- rep("Pareto", dim(pareto_200_45_54_phia)[1])
pareto_200_45_54_phia$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_45_54_phia)[1])
pareto_200_45_54_phia$agecat <- rep("45-54", dim(pareto_200_45_54_phia)[1])
pareto_200_45_54_phia$source <- rep("PHIA surveys", dim(pareto_200_45_54_phia)[1])

#' Pareto distribution < 50 copies/ml
pareto_50_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "Pareto 45-54",
                                      survey = PHIA_vls_estimates_45_54$survey,
                                      empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

pareto_50_45_54_phia$distribution <- rep("Pareto", dim(pareto_50_45_54_phia)[1])
pareto_50_45_54_phia$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_45_54_phia)[1])
pareto_50_45_54_phia$agecat <- rep("45-54", dim(pareto_50_45_54_phia)[1])
pareto_50_45_54_phia$source <- rep("PHIA surveys", dim(pareto_50_45_54_phia)[1])


#' calculate adjusted viral suppression estimates
#' Reverse Weibull distribution < 400 copies/ml
#' 55+
rev_weib_400_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_400,
                                      estimate_threshold = 400, conversion_threshold = 1000, 
                                      distribution = "Reverse Weibull 55+", 
                                      survey = PHIA_vls_estimates_55$survey,
                                      empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000, 
                                      empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

rev_weib_400_55_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_400_55_phia)[1])
rev_weib_400_55_phia$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_55_phia)[1])
rev_weib_400_55_phia$agecat <- rep("55+", dim(rev_weib_400_55_phia)[1])
rev_weib_400_55_phia$source <- rep("PHIA surveys", dim(rev_weib_400_55_phia)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_200,
                                      estimate_threshold = 200, conversion_threshold = 1000, 
                                      distribution = "Reverse Weibull 55+",
                                      survey = PHIA_vls_estimates_55$survey,
                                      empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

rev_weib_200_55_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_200_55_phia)[1])
rev_weib_200_55_phia$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_55_phia)[1])
rev_weib_200_55_phia$agecat <- rep("55+", dim(rev_weib_200_55_phia)[1])
rev_weib_200_55_phia$source <- rep("PHIA surveys", dim(rev_weib_200_55_phia)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_50,
                                     estimate_threshold = 50, conversion_threshold = 1000, 
                                     distribution = "Reverse Weibull 55+",
                                     survey = PHIA_vls_estimates_55$survey,
                                     empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

rev_weib_50_55_phia$distribution <- rep("Reverse Weibull", dim(rev_weib_50_55_phia)[1])
rev_weib_50_55_phia$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_55_phia)[1])
rev_weib_50_55_phia$agecat <- rep("55+", dim(rev_weib_50_55_phia)[1])
rev_weib_50_55_phia$source <- rep("PHIA surveys", dim(rev_weib_50_55_phia)[1])


#' Weibull distribution < 400 copies/ml
weib_400_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_400,
                                  estimate_threshold = 400, conversion_threshold = 1000, 
                                  distribution = "Weibull 55+", 
                                  survey = PHIA_vls_estimates_55$survey,
                                  empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                  empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000, 
                                  empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

weib_400_55_phia$distribution <- rep("Weibull", dim(weib_400_55_phia)[1])
weib_400_55_phia$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_55_phia)[1])
weib_400_55_phia$agecat <- rep("55+", dim(weib_400_55_phia)[1])
weib_400_55_phia$source <- rep("PHIA surveys", dim(weib_400_55_phia)[1])

#' Weibull distribution < 200 copies/ml
weib_200_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_200,
                                  estimate_threshold = 200, conversion_threshold = 1000, 
                                  distribution = "Weibull 55+",
                                  survey = PHIA_vls_estimates_55$survey,
                                  empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                  empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                  empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

weib_200_55_phia$distribution <- rep("Weibull", dim(weib_200_55_phia)[1])
weib_200_55_phia$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_55_phia)[1])
weib_200_55_phia$agecat <- rep("55+", dim(weib_200_55_phia)[1])
weib_200_55_phia$source <- rep("PHIA surveys", dim(weib_200_55_phia)[1])

#' Weibull distribution < 50 copies/ml
weib_50_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_50,
                                 estimate_threshold = 50, conversion_threshold = 1000, 
                                 distribution = "Weibull 55+",
                                 survey = PHIA_vls_estimates_55$survey,
                                 empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                 empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                 empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

weib_50_55_phia$distribution <- rep("Weibull", dim(weib_50_55_phia)[1])
weib_50_55_phia$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_55_phia)[1])
weib_50_55_phia$agecat <- rep("55+", dim(weib_50_55_phia)[1])
weib_50_55_phia$source <- rep("PHIA surveys", dim(weib_50_55_phia)[1])


#' Pareto distribution < 400 copies/ml
pareto_400_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_400,
                                    estimate_threshold = 400, conversion_threshold = 1000, 
                                    distribution = "Pareto 55+", 
                                    survey = PHIA_vls_estimates_55$survey,
                                    empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000, 
                                    empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

pareto_400_55_phia$distribution <- rep("Pareto", dim(pareto_400_55_phia)[1])
pareto_400_55_phia$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_55_phia)[1])
pareto_400_55_phia$agecat <- rep("55+", dim(pareto_400_55_phia)[1])
pareto_400_55_phia$source <- rep("PHIA surveys", dim(pareto_400_55_phia)[1])

#' Pareto distribution < 200 copies/ml
pareto_200_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_200,
                                    estimate_threshold = 200, conversion_threshold = 1000, 
                                    distribution = "Pareto 55+",
                                    survey = PHIA_vls_estimates_55$survey,
                                    empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

pareto_200_55_phia$distribution <- rep("Pareto", dim(pareto_200_55_phia)[1])
pareto_200_55_phia$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_55_phia)[1])
pareto_200_55_phia$agecat <- rep("55+", dim(pareto_200_55_phia)[1])
pareto_200_55_phia$source <- rep("PHIA surveys", dim(pareto_200_55_phia)[1])

#' Pareto distribution < 50 copies/ml
pareto_50_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_50,
                                   estimate_threshold = 50, conversion_threshold = 1000, 
                                   distribution = "Pareto 55+",
                                   survey = PHIA_vls_estimates_55$survey,
                                   empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

pareto_50_55_phia$distribution <- rep("Pareto", dim(pareto_50_55_phia)[1])
pareto_50_55_phia$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_55_phia)[1])
pareto_50_55_phia$agecat <- rep("55+", dim(pareto_50_55_phia)[1])
pareto_50_55_phia$source <- rep("PHIA surveys", dim(pareto_50_55_phia)[1])


# using parameters from Johnson et al
# analyses repeated here
#' calculate adjusted viral suppression estimates
#' Reverse Weibull distribution < 400 copies/ml
#' 15-24
rev_weib_400_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull Johnson", 
                                         survey = PHIA_vls_estimates_15_24$survey,
                                         empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

rev_weib_400_15_24$distribution <- rep("Reverse Weibull", dim(rev_weib_400_15_24)[1])
rev_weib_400_15_24$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_15_24)[1])
rev_weib_400_15_24$agecat <- rep("15-24", dim(rev_weib_400_15_24)[1])
rev_weib_400_15_24$source <- rep("Johnson et al", dim(rev_weib_400_15_24)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull Johnson",
                                         survey = PHIA_vls_estimates_15_24$survey,
                                         empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

rev_weib_200_15_24$distribution <- rep("Reverse Weibull", dim(rev_weib_200_15_24)[1])
rev_weib_200_15_24$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_15_24)[1])
rev_weib_200_15_24$agecat <- rep("15-24", dim(rev_weib_200_15_24)[1])
rev_weib_200_15_24$source <- rep("Johnson et al", dim(rev_weib_200_15_24)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "Reverse Weibull Johnson",
                                        survey = PHIA_vls_estimates_15_24$survey,
                                        empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

rev_weib_50_15_24$distribution <- rep("Reverse Weibull", dim(rev_weib_50_15_24)[1])
rev_weib_50_15_24$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_15_24)[1])
rev_weib_50_15_24$agecat <- rep("15-24", dim(rev_weib_50_15_24)[1])
rev_weib_50_15_24$source <- rep("Johnson et al", dim(rev_weib_50_15_24)[1])


#' Weibull distribution < 400 copies/ml
weib_400_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_400,
                                     estimate_threshold = 400, conversion_threshold = 1000, 
                                     distribution = "Weibull Johnson", 
                                     survey = PHIA_vls_estimates_15_24$survey,
                                     empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000, 
                                     empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

weib_400_15_24$distribution <- rep("Weibull", dim(weib_400_15_24)[1])
weib_400_15_24$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_15_24)[1])
weib_400_15_24$agecat <- rep("15-24", dim(weib_400_15_24)[1])
weib_400_15_24$source <- rep("Johnson et al", dim(weib_400_15_24)[1])

#' Weibull distribution < 200 copies/ml
weib_200_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_200,
                                     estimate_threshold = 200, conversion_threshold = 1000, 
                                     distribution = "Weibull Johnson",
                                     survey = PHIA_vls_estimates_15_24$survey,
                                     empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

weib_200_15_24$distribution <- rep("Weibull", dim(weib_200_15_24)[1])
weib_200_15_24$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_15_24)[1])
weib_200_15_24$agecat <- rep("15-24", dim(weib_200_15_24)[1])
weib_200_15_24$source <- rep("Johnson et al", dim(weib_200_15_24)[1])

#' Weibull distribution < 50 copies/ml
weib_50_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_50,
                                    estimate_threshold = 50, conversion_threshold = 1000, 
                                    distribution = "Weibull Johnson",
                                    survey = PHIA_vls_estimates_15_24$survey,
                                    empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

weib_50_15_24$distribution <- rep("Weibull", dim(weib_50_15_24)[1])
weib_50_15_24$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_15_24)[1])
weib_50_15_24$agecat <- rep("15-24", dim(weib_50_15_24)[1])
weib_50_15_24$source <- rep("Johnson et al", dim(weib_50_15_24)[1])

#' Pareto distribution < 400 copies/ml
pareto_400_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "Pareto Johnson", 
                                       survey = PHIA_vls_estimates_15_24$survey,
                                       empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000, 
                                       empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

pareto_400_15_24$distribution <- rep("Pareto", dim(pareto_400_15_24)[1])
pareto_400_15_24$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_15_24)[1])
pareto_400_15_24$agecat <- rep("15-24", dim(pareto_400_15_24)[1])
pareto_400_15_24$source <- rep("Johnson et al", dim(pareto_400_15_24)[1])

#' Pareto distribution < 200 copies/ml
pareto_200_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "Pareto Johnson",
                                       survey = PHIA_vls_estimates_15_24$survey,
                                       empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

pareto_200_15_24$distribution <- rep("Pareto", dim(pareto_200_15_24)[1])
pareto_200_15_24$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_15_24)[1])
pareto_200_15_24$agecat <- rep("15-24", dim(pareto_200_15_24)[1])
pareto_200_15_24$source <- rep("Johnson et al", dim(pareto_200_15_24)[1])

#' Pareto distribution < 50 copies/ml
pareto_50_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "Pareto Johnson",
                                      survey = PHIA_vls_estimates_15_24$survey,
                                      empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

pareto_50_15_24$distribution <- rep("Pareto", dim(pareto_50_15_24)[1])
pareto_50_15_24$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_15_24)[1])
pareto_50_15_24$agecat <- rep("15-24", dim(pareto_50_15_24)[1])
pareto_50_15_24$source <- rep("Johnson et al", dim(pareto_50_15_24)[1])


#' calculate adjusted viral suppression estimates
#' Reverse Weibull distribution < 400 copies/ml
#' 25-34
rev_weib_400_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull Johnson", 
                                         survey = PHIA_vls_estimates_25_34$survey,
                                         empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

rev_weib_400_25_34$distribution <- rep("Reverse Weibull", dim(rev_weib_400_25_34)[1])
rev_weib_400_25_34$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_25_34)[1])
rev_weib_400_25_34$agecat <- rep("25-34", dim(rev_weib_400_25_34)[1])
rev_weib_400_25_34$source <- rep("Johnson et al", dim(rev_weib_400_25_34)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull Johnson",
                                         survey = PHIA_vls_estimates_25_34$survey,
                                         empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

rev_weib_200_25_34$distribution <- rep("Reverse Weibull", dim(rev_weib_200_25_34)[1])
rev_weib_200_25_34$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_25_34)[1])
rev_weib_200_25_34$agecat <- rep("25-34", dim(rev_weib_200_25_34)[1])
rev_weib_200_25_34$source <- rep("Johnson et al", dim(rev_weib_200_25_34)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "Reverse Weibull Johnson",
                                        survey = PHIA_vls_estimates_25_34$survey,
                                        empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

rev_weib_50_25_34$distribution <- rep("Reverse Weibull", dim(rev_weib_50_25_34)[1])
rev_weib_50_25_34$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_25_34)[1])
rev_weib_50_25_34$agecat <- rep("25-34", dim(rev_weib_50_25_34)[1])
rev_weib_50_25_34$source <- rep("Johnson et al", dim(rev_weib_50_25_34)[1])


#' Weibull distribution < 400 copies/ml
weib_400_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_400,
                                     estimate_threshold = 400, conversion_threshold = 1000, 
                                     distribution = "Weibull Johnson", 
                                     survey = PHIA_vls_estimates_25_34$survey,
                                     empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000, 
                                     empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

weib_400_25_34$distribution <- rep("Weibull", dim(weib_400_25_34)[1])
weib_400_25_34$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_25_34)[1])
weib_400_25_34$agecat <- rep("25-34", dim(weib_400_25_34)[1])
weib_400_25_34$source <- rep("Johnson et al", dim(weib_400_25_34)[1])

#' Weibull distribution < 200 copies/ml
weib_200_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_200,
                                     estimate_threshold = 200, conversion_threshold = 1000, 
                                     distribution = "Weibull Johnson",
                                     survey = PHIA_vls_estimates_25_34$survey,
                                     empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

weib_200_25_34$distribution <- rep("Weibull", dim(weib_200_25_34)[1])
weib_200_25_34$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_25_34)[1])
weib_200_25_34$agecat <- rep("25-34", dim(weib_200_25_34)[1])
weib_200_25_34$source <- rep("Johnson et al", dim(weib_200_25_34)[1])

#' Weibull distribution < 50 copies/ml
weib_50_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_50,
                                    estimate_threshold = 50, conversion_threshold = 1000, 
                                    distribution = "Weibull Johnson",
                                    survey = PHIA_vls_estimates_25_34$survey,
                                    empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

weib_50_25_34$distribution <- rep("Weibull", dim(weib_50_25_34)[1])
weib_50_25_34$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_25_34)[1])
weib_50_25_34$agecat <- rep("25-34", dim(weib_50_25_34)[1])
weib_50_25_34$source <- rep("Johnson et al", dim(weib_50_25_34)[1])


#' Pareto distribution < 400 copies/ml
pareto_400_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "Pareto Johnson", 
                                       survey = PHIA_vls_estimates_25_34$survey,
                                       empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000, 
                                       empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

pareto_400_25_34$distribution <- rep("Pareto", dim(pareto_400_25_34)[1])
pareto_400_25_34$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_25_34)[1])
pareto_400_25_34$agecat <- rep("25-34", dim(pareto_400_25_34)[1])
pareto_400_25_34$source <- rep("Johnson et al", dim(pareto_400_25_34)[1])

#' Pareto distribution < 200 copies/ml
pareto_200_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "Pareto Johnson",
                                       survey = PHIA_vls_estimates_25_34$survey,
                                       empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

pareto_200_25_34$distribution <- rep("Pareto", dim(pareto_200_25_34)[1])
pareto_200_25_34$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_25_34)[1])
pareto_200_25_34$agecat <- rep("25-34", dim(pareto_200_25_34)[1])
pareto_200_25_34$source <- rep("Johnson et al", dim(pareto_200_25_34)[1])

#' Pareto distribution < 50 copies/ml
pareto_50_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "Pareto Johnson",
                                      survey = PHIA_vls_estimates_25_34$survey,
                                      empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

pareto_50_25_34$distribution <- rep("Pareto", dim(pareto_50_25_34)[1])
pareto_50_25_34$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_25_34)[1])
pareto_50_25_34$agecat <- rep("25-34", dim(pareto_50_25_34)[1])
pareto_50_25_34$source <- rep("Johnson et al", dim(pareto_50_25_34)[1])


#' calculate adjusted viral suppression estimates
#' Reverse Weibull distribution < 400 copies/ml
#' 35-44
rev_weib_400_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull Johnson", 
                                         survey = PHIA_vls_estimates_35_44$survey,
                                         empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

rev_weib_400_35_44$distribution <- rep("Reverse Weibull", dim(rev_weib_400_35_44)[1])
rev_weib_400_35_44$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_35_44)[1])
rev_weib_400_35_44$agecat <- rep("35-44", dim(rev_weib_400_35_44)[1])
rev_weib_400_35_44$source <- rep("Johnson et al", dim(rev_weib_400_35_44)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull Johnson",
                                         survey = PHIA_vls_estimates_35_44$survey,
                                         empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

rev_weib_200_35_44$distribution <- rep("Reverse Weibull", dim(rev_weib_200_35_44)[1])
rev_weib_200_35_44$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_35_44)[1])
rev_weib_200_35_44$agecat <- rep("35-44", dim(rev_weib_200_35_44)[1])
rev_weib_200_35_44$source <- rep("Johnson et al", dim(rev_weib_200_35_44)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "Reverse Weibull Johnson",
                                        survey = PHIA_vls_estimates_35_44$survey,
                                        empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

rev_weib_50_35_44$distribution <- rep("Reverse Weibull", dim(rev_weib_50_35_44)[1])
rev_weib_50_35_44$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_35_44)[1])
rev_weib_50_35_44$agecat <- rep("35-44", dim(rev_weib_50_35_44)[1])
rev_weib_50_35_44$source <- rep("Johnson et al", dim(rev_weib_50_35_44)[1])


#' Weibull distribution < 400 copies/ml
weib_400_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_400,
                                     estimate_threshold = 400, conversion_threshold = 1000, 
                                     distribution = "Weibull Johnson", 
                                     survey = PHIA_vls_estimates_35_44$survey,
                                     empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000, 
                                     empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

weib_400_35_44$distribution <- rep("Weibull", dim(weib_400_35_44)[1])
weib_400_35_44$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_35_44)[1])
weib_400_35_44$agecat <- rep("35-44", dim(weib_400_35_44)[1])
weib_400_35_44$source <- rep("Johnson et al", dim(weib_400_35_44)[1])

#' Weibull distribution < 200 copies/ml
weib_200_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_200,
                                     estimate_threshold = 200, conversion_threshold = 1000, 
                                     distribution = "Weibull Johnson",
                                     survey = PHIA_vls_estimates_35_44$survey,
                                     empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

weib_200_35_44$distribution <- rep("Weibull", dim(weib_200_35_44)[1])
weib_200_35_44$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_35_44)[1])
weib_200_35_44$agecat <- rep("35-44", dim(weib_200_35_44)[1])
weib_200_35_44$source <- rep("Johnson et al", dim(weib_200_35_44)[1])

#' Weibull distribution < 50 copies/ml
weib_50_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_50,
                                    estimate_threshold = 50, conversion_threshold = 1000, 
                                    distribution = "Weibull Johnson",
                                    survey = PHIA_vls_estimates_35_44$survey,
                                    empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

weib_50_35_44$distribution <- rep("Weibull", dim(weib_50_35_44)[1])
weib_50_35_44$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_35_44)[1])
weib_50_35_44$agecat <- rep("35-44", dim(weib_50_35_44)[1])
weib_50_35_44$source <- rep("Johnson et al", dim(weib_50_35_44)[1])


#' Pareto distribution < 400 copies/ml
pareto_400_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "Pareto Johnson", 
                                       survey = PHIA_vls_estimates_35_44$survey,
                                       empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000, 
                                       empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

pareto_400_35_44$distribution <- rep("Pareto", dim(pareto_400_35_44)[1])
pareto_400_35_44$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_35_44)[1])
pareto_400_35_44$agecat <- rep("35-44", dim(pareto_400_35_44)[1])
pareto_400_35_44$source <- rep("Johnson et al", dim(pareto_400_35_44)[1])

#' Pareto distribution < 200 copies/ml
pareto_200_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "Pareto Johnson",
                                       survey = PHIA_vls_estimates_35_44$survey,
                                       empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

pareto_200_35_44$distribution <- rep("Pareto", dim(pareto_200_35_44)[1])
pareto_200_35_44$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_35_44)[1])
pareto_200_35_44$agecat <- rep("35-44", dim(pareto_200_35_44)[1])
pareto_200_35_44$source <- rep("Johnson et al", dim(pareto_200_35_44)[1])

#' Pareto distribution < 50 copies/ml
pareto_50_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "Pareto Johnson",
                                      survey = PHIA_vls_estimates_35_44$survey,
                                      empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

pareto_50_35_44$distribution <- rep("Pareto", dim(pareto_50_35_44)[1])
pareto_50_35_44$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_35_44)[1])
pareto_50_35_44$agecat <- rep("35-44", dim(pareto_50_35_44)[1])
pareto_50_35_44$source <- rep("Johnson et al", dim(pareto_50_35_44)[1])


#' calculate adjusted viral suppression estimates
#' Reverse Weibull distribution < 400 copies/ml
#' 45-54
rev_weib_400_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull Johnson", 
                                         survey = PHIA_vls_estimates_45_54$survey,
                                         empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

rev_weib_400_45_54$distribution <- rep("Reverse Weibull", dim(rev_weib_400_45_54)[1])
rev_weib_400_45_54$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_45_54)[1])
rev_weib_400_45_54$agecat <- rep("45-54", dim(rev_weib_400_45_54)[1])
rev_weib_400_45_54$source <- rep("Johnson et al", dim(rev_weib_400_45_54)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "Reverse Weibull Johnson",
                                         survey = PHIA_vls_estimates_45_54$survey,
                                         empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

rev_weib_200_45_54$distribution <- rep("Reverse Weibull", dim(rev_weib_200_45_54)[1])
rev_weib_200_45_54$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_45_54)[1])
rev_weib_200_45_54$agecat <- rep("45-54", dim(rev_weib_200_45_54)[1])
rev_weib_200_45_54$source <- rep("Johnson et al", dim(rev_weib_200_45_54)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "Reverse Weibull Johnson",
                                        survey = PHIA_vls_estimates_45_54$survey,
                                        empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

rev_weib_50_45_54$distribution <- rep("Reverse Weibull", dim(rev_weib_50_45_54)[1])
rev_weib_50_45_54$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_45_54)[1])
rev_weib_50_45_54$agecat <- rep("45-54", dim(rev_weib_50_45_54)[1])
rev_weib_50_45_54$source <- rep("Johnson et al", dim(rev_weib_50_45_54)[1])


#' Weibull distribution < 400 copies/ml
weib_400_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_400,
                                     estimate_threshold = 400, conversion_threshold = 1000, 
                                     distribution = "Weibull Johnson", 
                                     survey = PHIA_vls_estimates_45_54$survey,
                                     empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000, 
                                     empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

weib_400_45_54$distribution <- rep("Weibull", dim(weib_400_45_54)[1])
weib_400_45_54$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_45_54)[1])
weib_400_45_54$agecat <- rep("45-54", dim(weib_400_45_54)[1])
weib_400_45_54$source <- rep("Johnson et al", dim(weib_400_45_54)[1])

#' Weibull distribution < 200 copies/ml
weib_200_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_200,
                                     estimate_threshold = 200, conversion_threshold = 1000, 
                                     distribution = "Weibull Johnson",
                                     survey = PHIA_vls_estimates_45_54$survey,
                                     empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

weib_200_45_54$distribution <- rep("Weibull", dim(weib_200_45_54)[1])
weib_200_45_54$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_45_54)[1])
weib_200_45_54$agecat <- rep("45-54", dim(weib_200_45_54)[1])
weib_200_45_54$source <- rep("Johnson et al", dim(weib_200_45_54)[1])

#' Weibull distribution < 50 copies/ml
weib_50_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_50,
                                    estimate_threshold = 50, conversion_threshold = 1000, 
                                    distribution = "Weibull Johnson",
                                    survey = PHIA_vls_estimates_45_54$survey,
                                    empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

weib_50_45_54$distribution <- rep("Weibull", dim(weib_50_45_54)[1])
weib_50_45_54$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_45_54)[1])
weib_50_45_54$agecat <- rep("45-54", dim(weib_50_45_54)[1])
weib_50_45_54$source <- rep("Johnson et al", dim(weib_50_45_54)[1])


#' Pareto distribution < 400 copies/ml
pareto_400_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "Pareto Johnson", 
                                       survey = PHIA_vls_estimates_45_54$survey,
                                       empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000, 
                                       empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

pareto_400_45_54$distribution <- rep("Pareto", dim(pareto_400_45_54)[1])
pareto_400_45_54$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_45_54)[1])
pareto_400_45_54$agecat <- rep("45-54", dim(pareto_400_45_54)[1])
pareto_400_45_54$source <- rep("Johnson et al", dim(pareto_400_45_54)[1])

#' Pareto distribution < 200 copies/ml
pareto_200_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "Pareto Johnson",
                                       survey = PHIA_vls_estimates_45_54$survey,
                                       empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

pareto_200_45_54$distribution <- rep("Pareto", dim(pareto_200_45_54)[1])
pareto_200_45_54$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_45_54)[1])
pareto_200_45_54$agecat <- rep("45-54", dim(pareto_200_45_54)[1])
pareto_200_45_54$source <- rep("Johnson et al", dim(pareto_200_45_54)[1])

#' Pareto distribution < 50 copies/ml
pareto_50_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "Pareto Johnson",
                                      survey = PHIA_vls_estimates_45_54$survey,
                                      empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

pareto_50_45_54$distribution <- rep("Pareto", dim(pareto_50_45_54)[1])
pareto_50_45_54$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_45_54)[1])
pareto_50_45_54$agecat <- rep("45-54", dim(pareto_50_45_54)[1])
pareto_50_45_54$source <- rep("Johnson et al", dim(pareto_50_45_54)[1])


#' calculate adjusted viral suppression estimates
#' Reverse Weibull distribution < 400 copies/ml
#' 55+
rev_weib_400_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_400,
                                      estimate_threshold = 400, conversion_threshold = 1000, 
                                      distribution = "Reverse Weibull Johnson", 
                                      survey = PHIA_vls_estimates_55$survey,
                                      empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000, 
                                      empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

rev_weib_400_55$distribution <- rep("Reverse Weibull", dim(rev_weib_400_55)[1])
rev_weib_400_55$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400_55)[1])
rev_weib_400_55$agecat <- rep("55+", dim(rev_weib_400_55)[1])
rev_weib_400_55$source <- rep("Johnson et al", dim(rev_weib_400_55)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_200,
                                      estimate_threshold = 200, conversion_threshold = 1000, 
                                      distribution = "Reverse Weibull Johnson",
                                      survey = PHIA_vls_estimates_55$survey,
                                      empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

rev_weib_200_55$distribution <- rep("Reverse Weibull", dim(rev_weib_200_55)[1])
rev_weib_200_55$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_55)[1])
rev_weib_200_55$agecat <- rep("55+", dim(rev_weib_200_55)[1])
rev_weib_200_55$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200_55)[1])
rev_weib_200_55$source <- rep("Johnson et al", dim(rev_weib_200_55)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_50,
                                     estimate_threshold = 50, conversion_threshold = 1000, 
                                     distribution = "Reverse Weibull Johnson",
                                     survey = PHIA_vls_estimates_55$survey,
                                     empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

rev_weib_50_55$distribution <- rep("Reverse Weibull", dim(rev_weib_50_55)[1])
rev_weib_50_55$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50_55)[1])
rev_weib_50_55$agecat <- rep("55+", dim(rev_weib_50_55)[1])
rev_weib_50_55$source <- rep("Johnson et al", dim(rev_weib_50_55)[1])


#' Weibull distribution < 400 copies/ml
weib_400_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_400,
                                  estimate_threshold = 400, conversion_threshold = 1000, 
                                  distribution = "Weibull Johnson", 
                                  survey = PHIA_vls_estimates_55$survey,
                                  empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                  empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000, 
                                  empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

weib_400_55$distribution <- rep("Weibull", dim(weib_400_55)[1])
weib_400_55$estimate_threshold <- rep("<400 copies/mL", dim(weib_400_55)[1])
weib_400_55$agecat <- rep("55+", dim(weib_400_55)[1])
weib_400_55$source <- rep("Johnson et al", dim(weib_400_55)[1])

#' Weibull distribution < 200 copies/ml
weib_200_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_200,
                                  estimate_threshold = 200, conversion_threshold = 1000, 
                                  distribution = "Weibull Johnson",
                                  survey = PHIA_vls_estimates_55$survey,
                                  empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                  empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                  empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

weib_200_55$distribution <- rep("Weibull", dim(weib_200_55)[1])
weib_200_55$estimate_threshold <- rep("<200 copies/mL", dim(weib_200_55)[1])
weib_200_55$agecat <- rep("55+", dim(weib_200_55)[1])
weib_200_55$source <- rep("Johnson et al", dim(weib_200_55)[1])

#' Weibull distribution < 50 copies/ml
weib_50_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_50,
                                 estimate_threshold = 50, conversion_threshold = 1000, 
                                 distribution = "Weibull Johnson",
                                 survey = PHIA_vls_estimates_55$survey,
                                 empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                 empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                 empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

weib_50_55$distribution <- rep("Weibull", dim(weib_50_55)[1])
weib_50_55$estimate_threshold <- rep("<50 copies/mL", dim(weib_50_55)[1])
weib_50_55$agecat <- rep("55+", dim(weib_50_55)[1])
weib_50_55$source <- rep("Johnson et al", dim(weib_50_55)[1])


#' Pareto distribution < 400 copies/ml
pareto_400_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_400,
                                    estimate_threshold = 400, conversion_threshold = 1000, 
                                    distribution = "Pareto Johnson", 
                                    survey = PHIA_vls_estimates_55$survey,
                                    empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000, 
                                    empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

pareto_400_55$distribution <- rep("Pareto", dim(pareto_400_55)[1])
pareto_400_55$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400_55)[1])
pareto_400_55$agecat <- rep("55+", dim(pareto_400_55)[1])
pareto_400_55$source <- rep("Johnson et al", dim(pareto_400_55)[1])

#' Pareto distribution < 200 copies/ml
pareto_200_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_200,
                                    estimate_threshold = 200, conversion_threshold = 1000, 
                                    distribution = "Pareto Johnson",
                                    survey = PHIA_vls_estimates_55$survey,
                                    empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

pareto_200_55$distribution <- rep("Pareto", dim(pareto_200_55)[1])
pareto_200_55$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200_55)[1])
pareto_200_55$agecat <- rep("55+", dim(pareto_200_55)[1])
pareto_200_55$source <- rep("Johnson et al", dim(pareto_200_55)[1])

#' Pareto distribution < 50 copies/ml
pareto_50_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_50,
                                   estimate_threshold = 50, conversion_threshold = 1000, 
                                   distribution = "Pareto Johnson",
                                   survey = PHIA_vls_estimates_55$survey,
                                   empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

pareto_50_55$distribution <- rep("Pareto", dim(pareto_50_55)[1])
pareto_50_55$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50_55)[1])
pareto_50_55$agecat <- rep("55+", dim(pareto_50_55)[1])
pareto_50_55$source <- rep("Johnson et al", dim(pareto_50_55)[1])


#' calculate the root mean squared error between empirical and adjusted estimates
#' 15-24
sqrt(mean((rev_weib_50_15_24_phia$empirical_estimate - rev_weib_50_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_15_24_phia$empirical_estimate - weib_50_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_15_24_phia$empirical_estimate - pareto_50_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_50_15_24$empirical_estimate - rev_weib_50_15_24$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_15_24$empirical_estimate - weib_50_15_24$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_15_24$empirical_estimate - pareto_50_15_24$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_200_15_24_phia$empirical_estimate - rev_weib_200_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_15_24_phia$empirical_estimate - weib_200_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_15_24_phia$empirical_estimate - pareto_200_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_200_15_24$empirical_estimate - rev_weib_200_15_24$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_15_24$empirical_estimate - weib_200_15_24$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_15_24$empirical_estimate - pareto_200_15_24$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_400_15_24_phia$empirical_estimate - rev_weib_400_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_15_24_phia$empirical_estimate - weib_400_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_15_24_phia$empirical_estimate - pareto_400_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_400_15_24$empirical_estimate - rev_weib_400_15_24$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_15_24$empirical_estimate - weib_400_15_24$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_15_24$empirical_estimate - pareto_400_15_24$adjusted_vls_estimate)^2))

# 25-34
sqrt(mean((rev_weib_50_25_34_phia$empirical_estimate - rev_weib_50_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_25_34_phia$empirical_estimate - weib_50_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_25_34_phia$empirical_estimate - pareto_50_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_50_25_34$empirical_estimate - rev_weib_50_25_34$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_25_34$empirical_estimate - weib_50_25_34$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_25_34$empirical_estimate - pareto_50_25_34$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_200_25_34_phia$empirical_estimate - rev_weib_200_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_25_34_phia$empirical_estimate - weib_200_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_25_34_phia$empirical_estimate - pareto_200_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_200_25_34$empirical_estimate - rev_weib_200_25_34$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_25_34$empirical_estimate - weib_200_25_34$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_25_34$empirical_estimate - pareto_200_25_34$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_400_25_34_phia$empirical_estimate - rev_weib_400_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_25_34_phia$empirical_estimate - weib_400_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_25_34_phia$empirical_estimate - pareto_400_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_400_25_34$empirical_estimate - rev_weib_400_25_34$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_25_34$empirical_estimate - weib_400_25_34$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_25_34$empirical_estimate - pareto_400_25_34$adjusted_vls_estimate)^2))

# 35-44
sqrt(mean((rev_weib_50_35_44_phia$empirical_estimate - rev_weib_50_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_35_44_phia$empirical_estimate - weib_50_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_35_44_phia$empirical_estimate - pareto_50_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_50_35_44$empirical_estimate - rev_weib_50_35_44$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_35_44$empirical_estimate - weib_50_35_44$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_35_44$empirical_estimate - pareto_50_35_44$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_200_35_44_phia$empirical_estimate - rev_weib_200_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_35_44_phia$empirical_estimate - weib_200_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_35_44_phia$empirical_estimate - pareto_200_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_200_35_44$empirical_estimate - rev_weib_200_35_44$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_35_44$empirical_estimate - weib_200_35_44$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_35_44$empirical_estimate - pareto_200_35_44$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_400_35_44_phia$empirical_estimate - rev_weib_400_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_35_44_phia$empirical_estimate - weib_400_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_35_44_phia$empirical_estimate - pareto_400_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_400_35_44$empirical_estimate - rev_weib_400_35_44$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_35_44$empirical_estimate - weib_400_35_44$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_35_44$empirical_estimate - pareto_400_35_44$adjusted_vls_estimate)^2))

# 45-54
sqrt(mean((rev_weib_50_45_54_phia$empirical_estimate - rev_weib_50_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_45_54_phia$empirical_estimate - weib_50_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_45_54_phia$empirical_estimate - pareto_50_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_50_45_54$empirical_estimate - rev_weib_50_45_54$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_45_54$empirical_estimate - weib_50_45_54$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_45_54$empirical_estimate - pareto_50_45_54$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_200_45_54_phia$empirical_estimate - rev_weib_200_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_45_54_phia$empirical_estimate - weib_200_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_45_54_phia$empirical_estimate - pareto_200_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_200_45_54$empirical_estimate - rev_weib_200_45_54$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_45_54$empirical_estimate - weib_200_45_54$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_45_54$empirical_estimate - pareto_200_45_54$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_400_45_54_phia$empirical_estimate - rev_weib_400_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_45_54_phia$empirical_estimate - weib_400_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_45_54_phia$empirical_estimate - pareto_400_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_400_45_54$empirical_estimate - rev_weib_400_45_54$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_45_54$empirical_estimate - weib_400_45_54$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_45_54$empirical_estimate - pareto_400_45_54$adjusted_vls_estimate)^2))

# 55+
sqrt(mean((rev_weib_50_55_phia$empirical_estimate - rev_weib_50_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_55_phia$empirical_estimate - weib_50_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_55_phia$empirical_estimate - pareto_50_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_50_55$empirical_estimate - rev_weib_50_55$adjusted_vls_estimate)^2))
sqrt(mean((weib_50_55$empirical_estimate - weib_50_55$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50_55$empirical_estimate - pareto_50_55$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_200_55_phia$empirical_estimate - rev_weib_200_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_55_phia$empirical_estimate - weib_200_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_55_phia$empirical_estimate - pareto_200_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_200_55$empirical_estimate - rev_weib_200_55$adjusted_vls_estimate)^2))
sqrt(mean((weib_200_55$empirical_estimate - weib_200_55$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200_55$empirical_estimate - pareto_200_55$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_400_55_phia$empirical_estimate - rev_weib_400_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_55_phia$empirical_estimate - weib_400_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_55_phia$empirical_estimate - pareto_400_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((rev_weib_400_55$empirical_estimate - rev_weib_400_55$adjusted_vls_estimate)^2))
sqrt(mean((weib_400_55$empirical_estimate - weib_400_55$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400_55$empirical_estimate - pareto_400_55$adjusted_vls_estimate)^2))


# average bias
#' 15-24
mean(rev_weib_50_15_24_phia$adjusted_vls_estimate - rev_weib_50_15_24_phia$empirical_estimate)
mean(weib_50_15_24_phia$adjusted_vls_estimate - weib_50_15_24_phia$empirical_estimate)
mean(pareto_50_15_24_phia$adjusted_vls_estimate - pareto_50_15_24_phia$empirical_estimate)
mean(rev_weib_50_15_24$adjusted_vls_estimate - rev_weib_50_15_24$empirical_estimate)
mean(weib_50_15_24$adjusted_vls_estimate - weib_50_15_24$empirical_estimate)
mean(pareto_50_15_24$adjusted_vls_estimate - pareto_50_15_24$empirical_estimate)

mean(rev_weib_200_15_24_phia$adjusted_vls_estimate - rev_weib_200_15_24_phia$empirical_estimate)
mean(weib_200_15_24_phia$adjusted_vls_estimate - weib_200_15_24_phia$empirical_estimate)
mean(pareto_200_15_24_phia$adjusted_vls_estimate - pareto_200_15_24_phia$empirical_estimate)
mean(rev_weib_200_15_24$adjusted_vls_estimate - rev_weib_200_15_24$empirical_estimate)
mean(weib_200_15_24$adjusted_vls_estimate - weib_200_15_24$empirical_estimate)
mean(pareto_200_15_24$adjusted_vls_estimate - pareto_200_15_24$empirical_estimate)

mean(rev_weib_400_15_24_phia$adjusted_vls_estimate - rev_weib_400_15_24_phia$empirical_estimate)
mean(weib_400_15_24_phia$adjusted_vls_estimate - weib_400_15_24_phia$empirical_estimate)
mean(pareto_400_15_24_phia$adjusted_vls_estimate - pareto_400_15_24_phia$empirical_estimate)
mean(rev_weib_400_15_24$adjusted_vls_estimate - rev_weib_400_15_24$empirical_estimate)
mean(weib_400_15_24$adjusted_vls_estimate - weib_400_15_24$empirical_estimate)
mean(pareto_400_15_24$adjusted_vls_estimate - pareto_400_15_24$empirical_estimate)

# 25-34
mean(rev_weib_50_25_34_phia$adjusted_vls_estimate - rev_weib_50_25_34_phia$empirical_estimate)
mean(weib_50_25_34_phia$adjusted_vls_estimate - weib_50_25_34_phia$empirical_estimate)
mean(pareto_50_25_34_phia$adjusted_vls_estimate - pareto_50_25_34_phia$empirical_estimate)
mean(rev_weib_50_25_34$adjusted_vls_estimate - rev_weib_50_25_34$empirical_estimate)
mean(weib_50_25_34$adjusted_vls_estimate - weib_50_25_34$empirical_estimate)
mean(pareto_50_25_34$adjusted_vls_estimate - pareto_50_25_34$empirical_estimate)

mean(rev_weib_200_25_34_phia$adjusted_vls_estimate - rev_weib_200_25_34_phia$empirical_estimate)
mean(weib_200_25_34_phia$adjusted_vls_estimate - weib_200_25_34_phia$empirical_estimate)
mean(pareto_200_25_34_phia$adjusted_vls_estimate - pareto_200_25_34_phia$empirical_estimate)
mean(rev_weib_200_25_34$adjusted_vls_estimate - rev_weib_200_25_34$empirical_estimate)
mean(weib_200_25_34$adjusted_vls_estimate - weib_200_25_34$empirical_estimate)
mean(pareto_200_25_34$adjusted_vls_estimate - pareto_200_25_34$empirical_estimate)

mean(rev_weib_400_25_34_phia$adjusted_vls_estimate - rev_weib_400_25_34_phia$empirical_estimate)
mean(weib_400_25_34_phia$adjusted_vls_estimate - weib_400_25_34_phia$empirical_estimate)
mean(pareto_400_25_34_phia$adjusted_vls_estimate - pareto_400_25_34_phia$empirical_estimate)
mean(rev_weib_400_25_34$adjusted_vls_estimate - rev_weib_400_25_34$empirical_estimate)
mean(weib_400_25_34$adjusted_vls_estimate - weib_400_25_34$empirical_estimate)
mean(pareto_400_25_34$adjusted_vls_estimate - pareto_400_25_34$empirical_estimate)

# 35-44
mean(rev_weib_50_35_44_phia$adjusted_vls_estimate - rev_weib_50_35_44_phia$empirical_estimate)
mean(weib_50_35_44_phia$adjusted_vls_estimate - weib_50_35_44_phia$empirical_estimate)
mean(pareto_50_35_44_phia$adjusted_vls_estimate - pareto_50_35_44_phia$empirical_estimate)
mean(rev_weib_50_35_44$adjusted_vls_estimate - rev_weib_50_35_44$empirical_estimate)
mean(weib_50_35_44$adjusted_vls_estimate - weib_50_35_44$empirical_estimate)
mean(pareto_50_35_44$adjusted_vls_estimate - pareto_50_35_44$empirical_estimate)

mean(rev_weib_200_35_44_phia$adjusted_vls_estimate - rev_weib_200_35_44_phia$empirical_estimate)
mean(weib_200_35_44_phia$adjusted_vls_estimate - weib_200_35_44_phia$empirical_estimate)
mean(pareto_200_35_44_phia$adjusted_vls_estimate - pareto_200_35_44_phia$empirical_estimate)
mean(rev_weib_200_35_44$adjusted_vls_estimate - rev_weib_200_35_44$empirical_estimate)
mean(weib_200_35_44$adjusted_vls_estimate - weib_200_35_44$empirical_estimate)
mean(pareto_200_35_44$adjusted_vls_estimate - pareto_200_35_44$empirical_estimate)

mean(rev_weib_400_35_44_phia$adjusted_vls_estimate - rev_weib_400_35_44_phia$empirical_estimate)
mean(weib_400_35_44_phia$adjusted_vls_estimate - weib_400_35_44_phia$empirical_estimate)
mean(pareto_400_35_44_phia$adjusted_vls_estimate - pareto_400_35_44_phia$empirical_estimate)
mean(rev_weib_400_35_44$adjusted_vls_estimate - rev_weib_400_35_44$empirical_estimate)
mean(weib_400_35_44$adjusted_vls_estimate - weib_400_35_44$empirical_estimate)
mean(pareto_400_35_44$adjusted_vls_estimate - pareto_400_35_44$empirical_estimate)

# 45-54
mean(rev_weib_50_45_54_phia$adjusted_vls_estimate - rev_weib_50_45_54_phia$empirical_estimate)
mean(weib_50_45_54_phia$adjusted_vls_estimate - weib_50_45_54_phia$empirical_estimate)
mean(pareto_50_45_54_phia$adjusted_vls_estimate - pareto_50_45_54_phia$empirical_estimate)
mean(rev_weib_50_45_54$adjusted_vls_estimate - rev_weib_50_45_54$empirical_estimate)
mean(weib_50_45_54$adjusted_vls_estimate - weib_50_45_54$empirical_estimate)
mean(pareto_50_45_54$adjusted_vls_estimate - pareto_50_45_54$empirical_estimate)

mean(rev_weib_200_45_54_phia$adjusted_vls_estimate - rev_weib_200_45_54_phia$empirical_estimate)
mean(weib_200_45_54_phia$adjusted_vls_estimate - weib_200_45_54_phia$empirical_estimate)
mean(pareto_200_45_54_phia$adjusted_vls_estimate - pareto_200_45_54_phia$empirical_estimate)
mean(rev_weib_200_45_54$adjusted_vls_estimate - rev_weib_200_45_54$empirical_estimate)
mean(weib_200_45_54$adjusted_vls_estimate - weib_200_45_54$empirical_estimate)
mean(pareto_200_45_54$adjusted_vls_estimate - pareto_200_45_54$empirical_estimate)

mean(rev_weib_400_45_54_phia$adjusted_vls_estimate - rev_weib_400_45_54_phia$empirical_estimate)
mean(weib_400_45_54_phia$adjusted_vls_estimate - weib_400_45_54_phia$empirical_estimate)
mean(pareto_400_45_54_phia$adjusted_vls_estimate - pareto_400_45_54_phia$empirical_estimate)
mean(rev_weib_400_45_54$adjusted_vls_estimate - rev_weib_400_45_54$empirical_estimate)
mean(weib_400_45_54$adjusted_vls_estimate - weib_400_45_54$empirical_estimate)
mean(pareto_400_45_54$adjusted_vls_estimate - pareto_400_45_54$empirical_estimate)

# 55+
mean(rev_weib_50_55_phia$adjusted_vls_estimate - rev_weib_50_55_phia$empirical_estimate)
mean(weib_50_55_phia$adjusted_vls_estimate - weib_50_55_phia$empirical_estimate)
mean(pareto_50_55_phia$adjusted_vls_estimate - pareto_50_55_phia$empirical_estimate)
mean(rev_weib_50_55$adjusted_vls_estimate - rev_weib_50_55$empirical_estimate)
mean(weib_50_55$adjusted_vls_estimate - weib_50_55$empirical_estimate)
mean(pareto_50_55$adjusted_vls_estimate - pareto_50_55$empirical_estimate)

mean(rev_weib_200_55_phia$adjusted_vls_estimate - rev_weib_200_55_phia$empirical_estimate)
mean(weib_200_55_phia$adjusted_vls_estimate - weib_200_55_phia$empirical_estimate)
mean(pareto_200_55_phia$adjusted_vls_estimate - pareto_200_55_phia$empirical_estimate)
mean(rev_weib_200_55$adjusted_vls_estimate - rev_weib_200_55$empirical_estimate)
mean(weib_200_55$adjusted_vls_estimate - weib_200_55$empirical_estimate)
mean(pareto_200_55$adjusted_vls_estimate - pareto_200_55$empirical_estimate)

mean(rev_weib_400_55_phia$adjusted_vls_estimate - rev_weib_400_55_phia$empirical_estimate)
mean(weib_400_55_phia$adjusted_vls_estimate - weib_400_55_phia$empirical_estimate)
mean(pareto_400_55_phia$adjusted_vls_estimate - pareto_400_55_phia$empirical_estimate)
mean(rev_weib_400_55$adjusted_vls_estimate - rev_weib_400_55$empirical_estimate)
mean(weib_400_55$adjusted_vls_estimate - weib_400_55$empirical_estimate)
mean(pareto_400_55$adjusted_vls_estimate - pareto_400_55$empirical_estimate)


# mean absolute error
#' 15-24
mean(abs(rev_weib_50_15_24_phia$adjusted_vls_estimate - rev_weib_50_15_24_phia$empirical_estimate))
mean(abs(weib_50_15_24_phia$adjusted_vls_estimate - weib_50_15_24_phia$empirical_estimate))
mean(abs(pareto_50_15_24_phia$adjusted_vls_estimate - pareto_50_15_24_phia$empirical_estimate))
mean(abs(rev_weib_50_15_24$adjusted_vls_estimate - rev_weib_50_15_24$empirical_estimate))
mean(abs(weib_50_15_24$adjusted_vls_estimate - weib_50_15_24$empirical_estimate))
mean(abs(pareto_50_15_24$adjusted_vls_estimate - pareto_50_15_24$empirical_estimate))

mean(abs(rev_weib_200_15_24_phia$adjusted_vls_estimate - rev_weib_200_15_24_phia$empirical_estimate))
mean(abs(weib_200_15_24_phia$adjusted_vls_estimate - weib_200_15_24_phia$empirical_estimate))
mean(abs(pareto_200_15_24_phia$adjusted_vls_estimate - pareto_200_15_24_phia$empirical_estimate))
mean(abs(rev_weib_200_15_24$adjusted_vls_estimate - rev_weib_200_15_24$empirical_estimate))
mean(abs(weib_200_15_24$adjusted_vls_estimate - weib_200_15_24$empirical_estimate))
mean(abs(pareto_200_15_24$adjusted_vls_estimate - pareto_200_15_24$empirical_estimate))

mean(abs(rev_weib_400_15_24_phia$adjusted_vls_estimate - rev_weib_400_15_24_phia$empirical_estimate))
mean(abs(weib_400_15_24_phia$adjusted_vls_estimate - weib_400_15_24_phia$empirical_estimate))
mean(abs(pareto_400_15_24_phia$adjusted_vls_estimate - pareto_400_15_24_phia$empirical_estimate))
mean(abs(rev_weib_400_15_24$adjusted_vls_estimate - rev_weib_400_15_24$empirical_estimate))
mean(abs(weib_400_15_24$adjusted_vls_estimate - weib_400_15_24$empirical_estimate))
mean(abs(pareto_400_15_24$adjusted_vls_estimate - pareto_400_15_24$empirical_estimate))

# 25-34
mean(abs(rev_weib_50_25_34_phia$adjusted_vls_estimate - rev_weib_50_25_34_phia$empirical_estimate))
mean(abs(weib_50_25_34_phia$adjusted_vls_estimate - weib_50_25_34_phia$empirical_estimate))
mean(abs(pareto_50_25_34_phia$adjusted_vls_estimate - pareto_50_25_34_phia$empirical_estimate))
mean(abs(rev_weib_50_25_34$adjusted_vls_estimate - rev_weib_50_25_34$empirical_estimate))
mean(abs(weib_50_25_34$adjusted_vls_estimate - weib_50_25_34$empirical_estimate))
mean(abs(pareto_50_25_34$adjusted_vls_estimate - pareto_50_25_34$empirical_estimate))

mean(abs(rev_weib_200_25_34_phia$adjusted_vls_estimate - rev_weib_200_25_34_phia$empirical_estimate))
mean(abs(weib_200_25_34_phia$adjusted_vls_estimate - weib_200_25_34_phia$empirical_estimate))
mean(abs(pareto_200_25_34_phia$adjusted_vls_estimate - pareto_200_25_34_phia$empirical_estimate))
mean(abs(rev_weib_200_25_34$adjusted_vls_estimate - rev_weib_200_25_34$empirical_estimate))
mean(abs(weib_200_25_34$adjusted_vls_estimate - weib_200_25_34$empirical_estimate))
mean(abs(pareto_200_25_34$adjusted_vls_estimate - pareto_200_25_34$empirical_estimate))

mean(abs(rev_weib_400_25_34_phia$adjusted_vls_estimate - rev_weib_400_25_34_phia$empirical_estimate))
mean(abs(weib_400_25_34_phia$adjusted_vls_estimate - weib_400_25_34_phia$empirical_estimate))
mean(abs(pareto_400_25_34_phia$adjusted_vls_estimate - pareto_400_25_34_phia$empirical_estimate))
mean(abs(rev_weib_400_25_34$adjusted_vls_estimate - rev_weib_400_25_34$empirical_estimate))
mean(abs(weib_400_25_34$adjusted_vls_estimate - weib_400_25_34$empirical_estimate))
mean(abs(pareto_400_25_34$adjusted_vls_estimate - pareto_400_25_34$empirical_estimate))

# 35-44
mean(abs(rev_weib_50_35_44_phia$adjusted_vls_estimate - rev_weib_50_35_44_phia$empirical_estimate))
mean(abs(weib_50_35_44_phia$adjusted_vls_estimate - weib_50_35_44_phia$empirical_estimate))
mean(abs(pareto_50_35_44_phia$adjusted_vls_estimate - pareto_50_35_44_phia$empirical_estimate))
mean(abs(rev_weib_50_35_44$adjusted_vls_estimate - rev_weib_50_35_44$empirical_estimate))
mean(abs(weib_50_35_44$adjusted_vls_estimate - weib_50_35_44$empirical_estimate))
mean(abs(pareto_50_35_44$adjusted_vls_estimate - pareto_50_35_44$empirical_estimate))

mean(abs(rev_weib_200_35_44_phia$adjusted_vls_estimate - rev_weib_200_35_44_phia$empirical_estimate))
mean(abs(weib_200_35_44_phia$adjusted_vls_estimate - weib_200_35_44_phia$empirical_estimate))
mean(abs(pareto_200_35_44_phia$adjusted_vls_estimate - pareto_200_35_44_phia$empirical_estimate))
mean(abs(rev_weib_200_35_44$adjusted_vls_estimate - rev_weib_200_35_44$empirical_estimate))
mean(abs(weib_200_35_44$adjusted_vls_estimate - weib_200_35_44$empirical_estimate))
mean(abs(pareto_200_35_44$adjusted_vls_estimate - pareto_200_35_44$empirical_estimate))

mean(abs(rev_weib_400_35_44_phia$adjusted_vls_estimate - rev_weib_400_35_44_phia$empirical_estimate))
mean(abs(weib_400_35_44_phia$adjusted_vls_estimate - weib_400_35_44_phia$empirical_estimate))
mean(abs(pareto_400_35_44_phia$adjusted_vls_estimate - pareto_400_35_44_phia$empirical_estimate))
mean(abs(rev_weib_400_35_44$adjusted_vls_estimate - rev_weib_400_35_44$empirical_estimate))
mean(abs(weib_400_35_44$adjusted_vls_estimate - weib_400_35_44$empirical_estimate))
mean(abs(pareto_400_35_44$adjusted_vls_estimate - pareto_400_35_44$empirical_estimate))

# 45-54
mean(abs(rev_weib_50_45_54_phia$adjusted_vls_estimate - rev_weib_50_45_54_phia$empirical_estimate))
mean(abs(weib_50_45_54_phia$adjusted_vls_estimate - weib_50_45_54_phia$empirical_estimate))
mean(abs(pareto_50_45_54_phia$adjusted_vls_estimate - pareto_50_45_54_phia$empirical_estimate))
mean(abs(rev_weib_50_45_54$adjusted_vls_estimate - rev_weib_50_45_54$empirical_estimate))
mean(abs(weib_50_45_54$adjusted_vls_estimate - weib_50_45_54$empirical_estimate))
mean(abs(pareto_50_45_54$adjusted_vls_estimate - pareto_50_45_54$empirical_estimate))

mean(abs(rev_weib_200_45_54_phia$adjusted_vls_estimate - rev_weib_200_45_54_phia$empirical_estimate))
mean(abs(weib_200_45_54_phia$adjusted_vls_estimate - weib_200_45_54_phia$empirical_estimate))
mean(abs(pareto_200_45_54_phia$adjusted_vls_estimate - pareto_200_45_54_phia$empirical_estimate))
mean(abs(rev_weib_200_45_54$adjusted_vls_estimate - rev_weib_200_45_54$empirical_estimate))
mean(abs(weib_200_45_54$adjusted_vls_estimate - weib_200_45_54$empirical_estimate))
mean(abs(pareto_200_45_54$adjusted_vls_estimate - pareto_200_45_54$empirical_estimate))

mean(abs(rev_weib_400_45_54_phia$adjusted_vls_estimate - rev_weib_400_45_54_phia$empirical_estimate))
mean(abs(weib_400_45_54_phia$adjusted_vls_estimate - weib_400_45_54_phia$empirical_estimate))
mean(abs(pareto_400_45_54_phia$adjusted_vls_estimate - pareto_400_45_54_phia$empirical_estimate))
mean(abs(rev_weib_400_45_54$adjusted_vls_estimate - rev_weib_400_45_54$empirical_estimate))
mean(abs(weib_400_45_54$adjusted_vls_estimate - weib_400_45_54$empirical_estimate))
mean(abs(pareto_400_45_54$adjusted_vls_estimate - pareto_400_45_54$empirical_estimate))

# 55+
mean(abs(rev_weib_50_55_phia$adjusted_vls_estimate - rev_weib_50_55_phia$empirical_estimate))
mean(abs(weib_50_55_phia$adjusted_vls_estimate - weib_50_55_phia$empirical_estimate))
mean(abs(pareto_50_55_phia$adjusted_vls_estimate - pareto_50_55_phia$empirical_estimate))
mean(abs(rev_weib_50_55$adjusted_vls_estimate - rev_weib_50_55$empirical_estimate))
mean(abs(weib_50_55$adjusted_vls_estimate - weib_50_55$empirical_estimate))
mean(abs(pareto_50_55$adjusted_vls_estimate - pareto_50_55$empirical_estimate))

mean(abs(rev_weib_200_55_phia$adjusted_vls_estimate - rev_weib_200_55_phia$empirical_estimate))
mean(abs(weib_200_55_phia$adjusted_vls_estimate - weib_200_55_phia$empirical_estimate))
mean(abs(pareto_200_55_phia$adjusted_vls_estimate - pareto_200_55_phia$empirical_estimate))
mean(abs(rev_weib_200_55$adjusted_vls_estimate - rev_weib_200_55$empirical_estimate))
mean(abs(weib_200_55$adjusted_vls_estimate - weib_200_55$empirical_estimate))
mean(abs(pareto_200_55$adjusted_vls_estimate - pareto_200_55$empirical_estimate))

mean(abs(rev_weib_400_55_phia$adjusted_vls_estimate - rev_weib_400_55_phia$empirical_estimate))
mean(abs(weib_400_55_phia$adjusted_vls_estimate - weib_400_55_phia$empirical_estimate))
mean(abs(pareto_400_55_phia$adjusted_vls_estimate - pareto_400_55_phia$empirical_estimate))
mean(abs(rev_weib_400_55$adjusted_vls_estimate - rev_weib_400_55$empirical_estimate))
mean(abs(weib_400_55$adjusted_vls_estimate - weib_400_55$empirical_estimate))
mean(abs(pareto_400_55$adjusted_vls_estimate - pareto_400_55$empirical_estimate))
     
# code for Figure S8
age_validation <- bind_rows(rev_weib_50_15_24,rev_weib_50_25_34,rev_weib_50_35_44,rev_weib_50_45_54,rev_weib_50_55,
                            rev_weib_50_15_24_phia,rev_weib_50_25_34_phia,rev_weib_50_35_44_phia,rev_weib_50_45_54_phia,rev_weib_50_55_phia, 
                            rev_weib_200_15_24,rev_weib_200_25_34,rev_weib_200_35_44,rev_weib_200_45_54,rev_weib_200_55,
                            rev_weib_200_15_24_phia,rev_weib_200_25_34_phia,rev_weib_200_35_44_phia,rev_weib_200_45_54_phia,rev_weib_200_55_phia, 
                            rev_weib_400_15_24,rev_weib_400_25_34,rev_weib_400_35_44,rev_weib_400_45_54,rev_weib_400_55,
                            rev_weib_400_15_24_phia,rev_weib_400_25_34_phia,rev_weib_400_35_44_phia,rev_weib_400_45_54_phia,rev_weib_400_55_phia,
                            weib_50_15_24,weib_50_25_34,weib_50_35_44,weib_50_45_54,weib_50_55,
                            weib_50_15_24_phia,weib_50_25_34_phia,weib_50_35_44_phia,weib_50_45_54_phia,weib_50_55_phia,
                            weib_200_15_24,weib_200_25_34,weib_200_35_44,weib_200_45_54,weib_200_55,
                            weib_200_15_24_phia,weib_200_25_34_phia,weib_200_35_44_phia,weib_200_45_54_phia,weib_200_55_phia,
                            weib_400_15_24,weib_400_25_34,weib_400_35_44,weib_400_45_54,weib_400_55,
                            weib_400_15_24_phia,weib_400_25_34_phia,weib_400_35_44_phia,weib_400_45_54_phia,weib_400_55_phia,
                            pareto_50_15_24,pareto_50_25_34,pareto_50_35_44,pareto_50_45_54,pareto_50_55,
                            pareto_50_15_24_phia,pareto_50_25_34_phia,pareto_50_35_44_phia,pareto_50_45_54_phia,pareto_50_55_phia,
                            pareto_200_15_24,pareto_200_25_34,pareto_200_35_44,pareto_200_45_54,pareto_200_55,
                            pareto_200_15_24_phia,pareto_200_25_34_phia,pareto_200_35_44_phia,pareto_200_45_54_phia,pareto_200_55_phia,
                            pareto_400_15_24,pareto_400_25_34,pareto_400_35_44,pareto_400_45_54,pareto_400_55,
                            pareto_400_15_24_phia,pareto_400_25_34_phia,pareto_400_35_44_phia,pareto_400_45_54_phia,pareto_400_55_phia)


names(age_validation)
table(age_validation$source, useNA = "always")
table(age_validation$distribution, useNA = "always")

#' update country names to survey names
age_validation <- age_validation %>%
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
         estimate_threshold = forcats::fct_relevel(estimate_threshold, 
                                                   c("<50 copies/mL", "<200 copies/mL", "<400 copies/mL")),
         distribution = forcats::fct_relevel(distribution,
                                           c("Reverse Weibull", "Weibull", "Pareto")),
         country_names = forcats::fct_relevel(country_names,
                                              c("Côte d'Ivoire (2017-18)","Cameroon (2017-18)",
                                                "Nigeria (2018)","Uganda (2016-17)","Zimbabwe (2015-16)",
                                                "Tanzania (2016-17)","Ethiopia (2017-18)","Lesotho (2016-17)",
                                                "Zambia (2016)","Mozambique (2021-22)","Rwanda (2018-19)",
                                                "Zimbabwe (2019-20)","Kenya (2018-19)","Malawi (2015-16)",
                                                "Namibia (2017)","Eswatini (2016-17)","Lesotho (2019-20)",
                                                "Eswatini (2021)","Zambia (2021)","Malawi (2020-21)", "Botswana (2021)")))



## scatter plot (Fig S8)
age_validation %>%
  #filter(estimate_threshold == "<50 copies/mL" & source == "Johnson et al") %>%
  filter(estimate_threshold == "<400 copies/mL" & source == "PHIA surveys") %>%
  #filter(agecat %in% c("15-24","55+")) %>%
  ggplot(aes(x = empirical_estimate, y = adjusted_vls_estimate, color = country_names)) +
  geom_point(size = 2) +
  geom_errorbar(aes(y = adjusted_vls_estimate, ymin = adjusted_vls_lower_ci, 
                    ymax = adjusted_vls_upper_ci), linewidth = 0.75, alpha = 0.3) + 
  geom_errorbarh(aes(xmin = empirical_lower_ci, 
                     xmax = empirical_upper_ci), linewidth = 0.75, alpha = 0.3) + 
  scale_x_continuous(labels = scales::percent, breaks = seq(0,1,0.1)) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0,1,0.1)) +
  coord_fixed(xlim = c(0.4,1), ylim = c(0.4,1)) +
  geom_abline(intercept = 0, slope = 1, col = "darkgrey", linewidth = 0.9, linetype = "dashed") +
  theme_bw(base_size = 13) +
  scale_color_manual(values = c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", 
                                "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70",
                                "maroon", "orchid1", "darkturquoise", "darkorange4", "brown",
                                "cyan", "red", "gold4", "orange","blue")) +
  facet_grid(agecat~distribution) + 
  labs(x = "% viral load suppression (PHIA survey)",
       y = "% viral load suppression (adjusted)", color = "",
       title = "Adjustment from <400 to \u22641000 copies/mL and shape parameters from calibration to PHIAs") + 
  theme(legend.position = "bottom",
        axis.text = element_text(size = rel(1), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(0.9), family = "sans"),
        panel.grid = element_blank(), 
        panel.spacing.x = unit(1, "lines"),
        plot.margin = margin(0,0,0,0),
        plot.title = element_text(size = rel(1.0), family = "sans", face = "bold"))

age_validation %>% filter(is.na(agecat))
