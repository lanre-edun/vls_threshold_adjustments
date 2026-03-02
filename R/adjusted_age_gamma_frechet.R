#' create function using the Gamma, Frechet and lognormal distributions to convert viral suppression estimates
#' @details probability denisity functions to obtain viral suppression conversion thresholds
#' @param function receives viral_sup_estimate: estimate of empirical viral suppression calculated from survey
#' reported_threshold: threshold which estimate is reported from
#' adjusted_threshold: threshold which estimate is being adjusted to 
#' distribution: distribution being explored
#' function requires vectorized functions defined in adjusted_gamma_frechet script
vls_prob_threshold <- function(viral_sup_estimate, estimate_threshold, conversion_threshold, distribution, survey,
                               empirical_estimate, empirical_lower_ci, empirical_upper_ci){
  
  if(distribution == "Frechet"){
    
    adjusted_vls_estimate <- frechet_func_vec(shape = 1.86,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- frechet_func_vec(shape = 1.83,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- frechet_func_vec(shape = 1.90,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Frechet distribution")
    
    return(res)
    
  }else if(distribution == "Gamma"){
    
    adjusted_vls_estimate <- gamma_func_vec(shape = 0.81, 
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- gamma_func_vec(shape = 0.79,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- gamma_func_vec(shape = 0.84,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Gamma")
    
    return(res)
    
  }else if(distribution == "lognormal"){
    
    adjusted_vls_estimate <- lognorm_func_vec(sdlog = 0.89, 
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- lognorm_func_vec(sdlog = 0.88,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- lognorm_func_vec(sdlog = 0.90,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("lognormal")
    
    return(res)
    
  }else if(distribution == "Frechet 15-24"){
    
    adjusted_vls_estimate <- frechet_func_vec(shape = 1.52,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- frechet_func_vec(shape = 1.43,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- frechet_func_vec(shape = 1.60,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Frechet 15-24")
    
    return(res)
    
    
  } else if(distribution == "Frechet 25-34"){
    
    adjusted_vls_estimate <- frechet_func_vec(shape = 1.69,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- frechet_func_vec(shape = 1.58,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- frechet_func_vec(shape = 1.80,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Frechet 25-34")
    
    return(res)
    
    
  }else if(distribution == "Frechet 35-44"){
    
    adjusted_vls_estimate <- frechet_func_vec(shape = 1.88,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- frechet_func_vec(shape = 1.78,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- frechet_func_vec(shape = 1.98,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Frechet 35-44")
    
    return(res)
    
    
  }else if(distribution == "Frechet 45-54"){
    
    adjusted_vls_estimate <- frechet_func_vec(shape = 2.13,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- frechet_func_vec(shape = 2.02,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- frechet_func_vec(shape = 2.25,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Frechet 45-54")
    
    return(res)
    
    
  }else if(distribution == "Frechet 55+"){
    
    adjusted_vls_estimate <- frechet_func_vec(shape = 2.27,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- frechet_func_vec(shape = 2.13,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- frechet_func_vec(shape = 2.40,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Frechet 55+")
    
    return(res)
    
    
  }else if(distribution == "lognormal 15-24"){
    
    adjusted_vls_estimate <- lognorm_func_vec(sdlog = 0.99, 
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- lognorm_func_vec(sdlog = 0.94,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- lognorm_func_vec(sdlog = 1.03,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("lognormal 15-24")
    
    return(res)
    
  }else if(distribution == "lognormal 25-34"){
    
    adjusted_vls_estimate <- lognorm_func_vec(sdlog = 0.96, 
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- lognorm_func_vec(sdlog = 0.90,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- lognorm_func_vec(sdlog = 1.02,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("lognormal 25-34")
    
    return(res)
    
  }else if(distribution == "lognormal 35-44"){
    
    adjusted_vls_estimate <- lognorm_func_vec(sdlog = 0.89, 
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- lognorm_func_vec(sdlog = 0.83,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- lognorm_func_vec(sdlog = 0.94,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("lognormal 35-44")
    
    return(res)
    
  }else if(distribution == "lognormal 45-54"){
    
    adjusted_vls_estimate <- lognorm_func_vec(sdlog = 0.80, 
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- lognorm_func_vec(sdlog = 0.75,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- lognorm_func_vec(sdlog = 0.86,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("lognormal 45-54")
    
    return(res)
    
  }else if(distribution == "lognormal 55+"){
    
    adjusted_vls_estimate <- lognorm_func_vec(sdlog = 0.78, 
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- lognorm_func_vec(sdlog = 0.72,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- lognorm_func_vec(sdlog = 0.83,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("lognormal 55+")
    
    return(res)
    
  }else if(distribution == "Gamma 15-24"){
    
    adjusted_vls_estimate <- gamma_func_vec(shape = 0.82, 
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- gamma_func_vec(shape = 0.74,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- gamma_func_vec(shape = 0.90,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Gamma 15-24")
    
    return(res)
    
  }else if(distribution == "Gamma 25-34"){
    
    adjusted_vls_estimate <- gamma_func_vec(shape = 0.73, 
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- gamma_func_vec(shape = 0.65,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- gamma_func_vec(shape = 0.83,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Gamma 25-34")
    
    return(res)
    
  }else if(distribution == "Gamma 35-44"){
    
    adjusted_vls_estimate <- gamma_func_vec(shape = 0.81, 
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- gamma_func_vec(shape = 0.73,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- gamma_func_vec(shape = 0.90,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Gamma 35-44")
    
    return(res)
    
  }else if(distribution == "Gamma 45-54"){
    
    adjusted_vls_estimate <- gamma_func_vec(shape = 0.97, 
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- gamma_func_vec(shape = 0.87,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- gamma_func_vec(shape = 1.08,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Gamma 45-54")
    
    return(res)
    
  }else{
    
    adjusted_vls_estimate <- gamma_func_vec(shape = 1.01, 
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- gamma_func_vec(shape = 0.89,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- gamma_func_vec(shape = 1.14,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Gamma 55+")
    
    return(res)
    
  }
  
}


#' calculate adjusted viral suppression estimates
#' Frechet distribution < 400 copies/ml
#' 15-24
frec_400_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_400,
                                              estimate_threshold = 400, conversion_threshold = 1000, 
                                              distribution = "Frechet 15-24", 
                                              survey = PHIA_vls_estimates_15_24$survey,
                                              empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                              empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000, 
                                              empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

frec_400_15_24_phia$distribution <- rep("Frechet", dim(frec_400_15_24_phia)[1])
frec_400_15_24_phia$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_15_24_phia)[1])
frec_400_15_24_phia$agecat <- rep("15-24", dim(frec_400_15_24_phia)[1])
frec_400_15_24_phia$source <- rep("PHIA surveys", dim(frec_400_15_24_phia)[1])

#' Frechet distribution < 200 copies/ml
frec_200_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_200,
                                              estimate_threshold = 200, conversion_threshold = 1000, 
                                              distribution = "Frechet 15-24",
                                              survey = PHIA_vls_estimates_15_24$survey,
                                              empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                              empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                              empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

frec_200_15_24_phia$distribution <- rep("Frechet", dim(frec_200_15_24_phia)[1])
frec_200_15_24_phia$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_15_24_phia)[1])
frec_200_15_24_phia$agecat <- rep("15-24", dim(frec_200_15_24_phia)[1])
frec_200_15_24_phia$source <- rep("PHIA surveys", dim(frec_200_15_24_phia)[1])

#' Frechet distribution < 50 copies/ml
frec_50_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_50,
                                             estimate_threshold = 50, conversion_threshold = 1000, 
                                             distribution = "Frechet 15-24",
                                             survey = PHIA_vls_estimates_15_24$survey,
                                             empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                             empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                             empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

frec_50_15_24_phia$distribution <- rep("Frechet", dim(frec_50_15_24_phia)[1])
frec_50_15_24_phia$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_15_24_phia)[1])
frec_50_15_24_phia$agecat <- rep("15-24", dim(frec_50_15_24_phia)[1])
frec_50_15_24_phia$source <- rep("PHIA surveys", dim(frec_50_15_24_phia)[1])

#' Gamma distribution < 400 copies/ml
gamma_400_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_400,
                                          estimate_threshold = 400, conversion_threshold = 1000, 
                                          distribution = "Gamma 15-24", 
                                          survey = PHIA_vls_estimates_15_24$survey,
                                          empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000, 
                                          empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

gamma_400_15_24_phia$distribution <- rep("Gamma", dim(gamma_400_15_24_phia)[1])
gamma_400_15_24_phia$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_15_24_phia)[1])
gamma_400_15_24_phia$agecat <- rep("15-24", dim(gamma_400_15_24_phia)[1])
gamma_400_15_24_phia$source <- rep("PHIA surveys", dim(gamma_400_15_24_phia)[1])

#' Gamma distribution < 200 copies/ml
gamma_200_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_200,
                                          estimate_threshold = 200, conversion_threshold = 1000, 
                                          distribution = "Gamma 15-24",
                                          survey = PHIA_vls_estimates_15_24$survey,
                                          empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                          empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

gamma_200_15_24_phia$distribution <- rep("Gamma", dim(gamma_200_15_24_phia)[1])
gamma_200_15_24_phia$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_15_24_phia)[1])
gamma_200_15_24_phia$agecat <- rep("15-24", dim(gamma_200_15_24_phia)[1])
gamma_200_15_24_phia$source <- rep("PHIA surveys", dim(gamma_200_15_24_phia)[1])

#' Gamma distribution < 50 copies/ml
gamma_50_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_50,
                                         estimate_threshold = 50, conversion_threshold = 1000, 
                                         distribution = "Gamma 15-24",
                                         survey = PHIA_vls_estimates_15_24$survey,
                                         empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

gamma_50_15_24_phia$distribution <- rep("Gamma", dim(gamma_50_15_24_phia)[1])
gamma_50_15_24_phia$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_15_24_phia)[1])
gamma_50_15_24_phia$agecat <- rep("15-24", dim(gamma_50_15_24_phia)[1])
gamma_50_15_24_phia$source <- rep("PHIA surveys", dim(gamma_50_15_24_phia)[1])

#' lognormal distribution < 400 copies/ml
lognormal_400_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_400,
                                            estimate_threshold = 400, conversion_threshold = 1000, 
                                            distribution = "lognormal 15-24", 
                                            survey = PHIA_vls_estimates_15_24$survey,
                                            empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                            empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000, 
                                            empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

lognormal_400_15_24_phia$distribution <- rep("lognormal", dim(lognormal_400_15_24_phia)[1])
lognormal_400_15_24_phia$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_15_24_phia)[1])
lognormal_400_15_24_phia$agecat <- rep("15-24", dim(lognormal_400_15_24_phia)[1])
lognormal_400_15_24_phia$source <- rep("PHIA surveys", dim(lognormal_400_15_24_phia)[1])

#' lognormal distribution < 200 copies/ml
lognormal_200_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_200,
                                            estimate_threshold = 200, conversion_threshold = 1000, 
                                            distribution = "lognormal 15-24",
                                            survey = PHIA_vls_estimates_15_24$survey,
                                            empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                            empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                            empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

lognormal_200_15_24_phia$distribution <- rep("lognormal", dim(lognormal_200_15_24_phia)[1])
lognormal_200_15_24_phia$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_15_24_phia)[1])
lognormal_200_15_24_phia$agecat <- rep("15-24", dim(lognormal_200_15_24_phia)[1])
lognormal_200_15_24_phia$source <- rep("PHIA surveys", dim(lognormal_200_15_24_phia)[1])

#' lognormal distribution < 50 copies/ml
lognormal_50_15_24_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_50,
                                           estimate_threshold = 50, conversion_threshold = 1000, 
                                           distribution = "lognormal 15-24",
                                           survey = PHIA_vls_estimates_15_24$survey,
                                           empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                           empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                           empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

lognormal_50_15_24_phia$distribution <- rep("lognormal", dim(lognormal_50_15_24_phia)[1])
lognormal_50_15_24_phia$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_15_24_phia)[1])
lognormal_50_15_24_phia$agecat <- rep("15-24", dim(lognormal_50_15_24_phia)[1])
lognormal_50_15_24_phia$source <- rep("PHIA surveys", dim(lognormal_50_15_24_phia)[1])


#' calculate adjusted viral suppression estimates
#' Frechet distribution < 400 copies/ml
#' 25-34
frec_400_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_400,
                                              estimate_threshold = 400, conversion_threshold = 1000, 
                                              distribution = "Frechet 25-34", 
                                              survey = PHIA_vls_estimates_25_34$survey,
                                              empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                              empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000, 
                                              empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

frec_400_25_34_phia$distribution <- rep("Frechet", dim(frec_400_25_34_phia)[1])
frec_400_25_34_phia$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_25_34_phia)[1])
frec_400_25_34_phia$agecat <- rep("25-34", dim(frec_400_25_34_phia)[1])
frec_400_25_34_phia$source <- rep("PHIA surveys", dim(frec_400_25_34_phia)[1])

#' Frechet distribution < 200 copies/ml
frec_200_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_200,
                                              estimate_threshold = 200, conversion_threshold = 1000, 
                                              distribution = "Frechet 25-34",
                                              survey = PHIA_vls_estimates_25_34$survey,
                                              empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                              empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                              empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

frec_200_25_34_phia$distribution <- rep("Frechet", dim(frec_200_25_34_phia)[1])
frec_200_25_34_phia$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_25_34_phia)[1])
frec_200_25_34_phia$agecat <- rep("25-34", dim(frec_200_25_34_phia)[1])
frec_200_25_34_phia$source <- rep("PHIA surveys", dim(frec_200_25_34_phia)[1])

#' Frechet distribution < 50 copies/ml
frec_50_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_50,
                                             estimate_threshold = 50, conversion_threshold = 1000, 
                                             distribution = "Frechet 25-34",
                                             survey = PHIA_vls_estimates_25_34$survey,
                                             empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                             empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                             empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

frec_50_25_34_phia$distribution <- rep("Frechet", dim(frec_50_25_34_phia)[1])
frec_50_25_34_phia$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_25_34_phia)[1])
frec_50_25_34_phia$agecat <- rep("25-34", dim(frec_50_25_34_phia)[1])
frec_50_25_34_phia$source <- rep("PHIA surveys", dim(frec_50_25_34_phia)[1])


#' Gamma distribution < 400 copies/ml
gamma_400_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_400,
                                          estimate_threshold = 400, conversion_threshold = 1000, 
                                          distribution = "Gamma 25-34", 
                                          survey = PHIA_vls_estimates_25_34$survey,
                                          empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000, 
                                          empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

gamma_400_25_34_phia$distribution <- rep("Gamma", dim(gamma_400_25_34_phia)[1])
gamma_400_25_34_phia$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_25_34_phia)[1])
gamma_400_25_34_phia$agecat <- rep("25-34", dim(gamma_400_25_34_phia)[1])
gamma_400_25_34_phia$source <- rep("PHIA surveys", dim(gamma_400_25_34_phia)[1])

#' Gamma distribution < 200 copies/ml
gamma_200_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_200,
                                          estimate_threshold = 200, conversion_threshold = 1000, 
                                          distribution = "Gamma 25-34",
                                          survey = PHIA_vls_estimates_25_34$survey,
                                          empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                          empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

gamma_200_25_34_phia$distribution <- rep("Gamma", dim(gamma_200_25_34_phia)[1])
gamma_200_25_34_phia$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_25_34_phia)[1])
gamma_200_25_34_phia$agecat <- rep("25-34", dim(gamma_200_25_34_phia)[1])
gamma_200_25_34_phia$source <- rep("PHIA surveys", dim(gamma_200_25_34_phia)[1])

#' Gamma distribution < 50 copies/ml
gamma_50_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_50,
                                         estimate_threshold = 50, conversion_threshold = 1000, 
                                         distribution = "Gamma 25-34",
                                         survey = PHIA_vls_estimates_25_34$survey,
                                         empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

gamma_50_25_34_phia$distribution <- rep("Gamma", dim(gamma_50_25_34_phia)[1])
gamma_50_25_34_phia$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_25_34_phia)[1])
gamma_50_25_34_phia$agecat <- rep("25-34", dim(gamma_50_25_34_phia)[1])
gamma_50_25_34_phia$source <- rep("PHIA surveys", dim(gamma_50_25_34_phia)[1])


#' lognormal distribution < 400 copies/ml
lognormal_400_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_400,
                                            estimate_threshold = 400, conversion_threshold = 1000, 
                                            distribution = "lognormal 25-34", 
                                            survey = PHIA_vls_estimates_25_34$survey,
                                            empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                            empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000, 
                                            empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

lognormal_400_25_34_phia$distribution <- rep("lognormal", dim(lognormal_400_25_34_phia)[1])
lognormal_400_25_34_phia$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_25_34_phia)[1])
lognormal_400_25_34_phia$agecat <- rep("25-34", dim(lognormal_400_25_34_phia)[1])
lognormal_400_25_34_phia$source <- rep("PHIA surveys", dim(lognormal_400_25_34_phia)[1])

#' lognormal distribution < 200 copies/ml
lognormal_200_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_200,
                                            estimate_threshold = 200, conversion_threshold = 1000, 
                                            distribution = "lognormal 25-34",
                                            survey = PHIA_vls_estimates_25_34$survey,
                                            empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                            empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                            empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

lognormal_200_25_34_phia$distribution <- rep("lognormal", dim(lognormal_200_25_34_phia)[1])
lognormal_200_25_34_phia$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_25_34_phia)[1])
lognormal_200_25_34_phia$agecat <- rep("25-34", dim(lognormal_200_25_34_phia)[1])
lognormal_200_25_34_phia$source <- rep("PHIA surveys", dim(lognormal_200_25_34_phia)[1])

#' lognormal distribution < 50 copies/ml
lognormal_50_25_34_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_50,
                                           estimate_threshold = 50, conversion_threshold = 1000, 
                                           distribution = "lognormal 25-34",
                                           survey = PHIA_vls_estimates_25_34$survey,
                                           empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                           empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                           empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

lognormal_50_25_34_phia$distribution <- rep("lognormal", dim(lognormal_50_25_34_phia)[1])
lognormal_50_25_34_phia$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_25_34_phia)[1])
lognormal_50_25_34_phia$agecat <- rep("25-34", dim(lognormal_50_25_34_phia)[1])
lognormal_50_25_34_phia$source <- rep("PHIA surveys", dim(lognormal_50_25_34_phia)[1])


#' calculate adjusted viral suppression estimates
#' Frechet distribution < 400 copies/ml
#' 35-44
frec_400_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_400,
                                              estimate_threshold = 400, conversion_threshold = 1000, 
                                              distribution = "Frechet 35-44", 
                                              survey = PHIA_vls_estimates_35_44$survey,
                                              empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                              empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000, 
                                              empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

frec_400_35_44_phia$distribution <- rep("Frechet", dim(frec_400_35_44_phia)[1])
frec_400_35_44_phia$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_35_44_phia)[1])
frec_400_35_44_phia$agecat <- rep("35-44", dim(frec_400_35_44_phia)[1])
frec_400_35_44_phia$source <- rep("PHIA surveys", dim(frec_400_35_44_phia)[1])

#' Frechet distribution < 200 copies/ml
frec_200_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_200,
                                              estimate_threshold = 200, conversion_threshold = 1000, 
                                              distribution = "Frechet 35-44",
                                              survey = PHIA_vls_estimates_35_44$survey,
                                              empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                              empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                              empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

frec_200_35_44_phia$distribution <- rep("Frechet", dim(frec_200_35_44_phia)[1])
frec_200_35_44_phia$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_35_44_phia)[1])
frec_200_35_44_phia$agecat <- rep("35-44", dim(frec_200_35_44_phia)[1])
frec_200_35_44_phia$source <- rep("PHIA surveys", dim(frec_200_35_44_phia)[1])

#' Frechet distribution < 50 copies/ml
frec_50_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_50,
                                             estimate_threshold = 50, conversion_threshold = 1000, 
                                             distribution = "Frechet 35-44",
                                             survey = PHIA_vls_estimates_35_44$survey,
                                             empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                             empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                             empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

frec_50_35_44_phia$distribution <- rep("Frechet", dim(frec_50_35_44_phia)[1])
frec_50_35_44_phia$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_35_44_phia)[1])
frec_50_35_44_phia$agecat <- rep("35-44", dim(frec_50_35_44_phia)[1])
frec_50_35_44_phia$source <- rep("PHIA surveys", dim(frec_50_35_44_phia)[1])


#' Gamma distribution < 400 copies/ml
gamma_400_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_400,
                                          estimate_threshold = 400, conversion_threshold = 1000, 
                                          distribution = "Gamma 35-44", 
                                          survey = PHIA_vls_estimates_35_44$survey,
                                          empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000, 
                                          empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

gamma_400_35_44_phia$distribution <- rep("Gamma", dim(gamma_400_35_44_phia)[1])
gamma_400_35_44_phia$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_35_44_phia)[1])
gamma_400_35_44_phia$agecat <- rep("35-44", dim(gamma_400_35_44_phia)[1])
gamma_400_35_44_phia$source <- rep("PHIA surveys", dim(gamma_400_35_44_phia)[1])

#' Gamma distribution < 200 copies/ml
gamma_200_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_200,
                                          estimate_threshold = 200, conversion_threshold = 1000, 
                                          distribution = "Gamma 35-44",
                                          survey = PHIA_vls_estimates_35_44$survey,
                                          empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                          empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

gamma_200_35_44_phia$distribution <- rep("Gamma", dim(gamma_200_35_44_phia)[1])
gamma_200_35_44_phia$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_35_44_phia)[1])
gamma_200_35_44_phia$agecat <- rep("35-44", dim(gamma_200_35_44_phia)[1])
gamma_200_35_44_phia$source <- rep("PHIA surveys", dim(gamma_200_35_44_phia)[1])

#' Gamma distribution < 50 copies/ml
gamma_50_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_50,
                                         estimate_threshold = 50, conversion_threshold = 1000, 
                                         distribution = "Gamma 35-44",
                                         survey = PHIA_vls_estimates_35_44$survey,
                                         empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

gamma_50_35_44_phia$distribution <- rep("Gamma", dim(gamma_50_35_44_phia)[1])
gamma_50_35_44_phia$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_35_44_phia)[1])
gamma_50_35_44_phia$agecat <- rep("35-44", dim(gamma_50_35_44_phia)[1])
gamma_50_35_44_phia$source <- rep("PHIA surveys", dim(gamma_50_35_44_phia)[1])


#' lognormal distribution < 400 copies/ml
lognormal_400_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_400,
                                            estimate_threshold = 400, conversion_threshold = 1000, 
                                            distribution = "lognormal 35-44", 
                                            survey = PHIA_vls_estimates_35_44$survey,
                                            empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                            empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000, 
                                            empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

lognormal_400_35_44_phia$distribution <- rep("lognormal", dim(lognormal_400_35_44_phia)[1])
lognormal_400_35_44_phia$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_35_44_phia)[1])
lognormal_400_35_44_phia$agecat <- rep("35-44", dim(lognormal_400_35_44_phia)[1])
lognormal_400_35_44_phia$source <- rep("PHIA surveys", dim(lognormal_400_35_44_phia)[1])

#' lognormal distribution < 200 copies/ml
lognormal_200_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_200,
                                            estimate_threshold = 200, conversion_threshold = 1000, 
                                            distribution = "lognormal 35-44",
                                            survey = PHIA_vls_estimates_35_44$survey,
                                            empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                            empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                            empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

lognormal_200_35_44_phia$distribution <- rep("lognormal", dim(lognormal_200_35_44_phia)[1])
lognormal_200_35_44_phia$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_35_44_phia)[1])
lognormal_200_35_44_phia$agecat <- rep("35-44", dim(lognormal_200_35_44_phia)[1])
lognormal_200_35_44_phia$source <- rep("PHIA surveys", dim(lognormal_200_35_44_phia)[1])

#' lognormal distribution < 50 copies/ml
lognormal_50_35_44_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_50,
                                           estimate_threshold = 50, conversion_threshold = 1000, 
                                           distribution = "lognormal 35-44",
                                           survey = PHIA_vls_estimates_35_44$survey,
                                           empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                           empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                           empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

lognormal_50_35_44_phia$distribution <- rep("lognormal", dim(lognormal_50_35_44_phia)[1])
lognormal_50_35_44_phia$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_35_44_phia)[1])
lognormal_50_35_44_phia$agecat <- rep("35-44", dim(lognormal_50_35_44_phia)[1])
lognormal_50_35_44_phia$source <- rep("PHIA surveys", dim(lognormal_50_35_44_phia)[1])


#' calculate adjusted viral suppression estimates
#' Frechet distribution < 400 copies/ml
#' 45-54
frec_400_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_400,
                                              estimate_threshold = 400, conversion_threshold = 1000, 
                                              distribution = "Frechet 45-54", 
                                              survey = PHIA_vls_estimates_45_54$survey,
                                              empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                              empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000, 
                                              empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

frec_400_45_54_phia$distribution <- rep("Frechet", dim(frec_400_45_54_phia)[1])
frec_400_45_54_phia$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_45_54_phia)[1])
frec_400_45_54_phia$agecat <- rep("45-54", dim(frec_400_45_54_phia)[1])
frec_400_45_54_phia$source <- rep("PHIA surveys", dim(frec_400_45_54_phia)[1])

#' Frechet distribution < 200 copies/ml
frec_200_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_200,
                                              estimate_threshold = 200, conversion_threshold = 1000, 
                                              distribution = "Frechet 45-54",
                                              survey = PHIA_vls_estimates_45_54$survey,
                                              empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                              empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                              empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

frec_200_45_54_phia$distribution <- rep("Frechet", dim(frec_200_45_54_phia)[1])
frec_200_45_54_phia$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_45_54_phia)[1])
frec_200_45_54_phia$agecat <- rep("45-54", dim(frec_200_45_54_phia)[1])
frec_200_45_54_phia$source <- rep("PHIA surveys", dim(frec_200_45_54_phia)[1])

#' Frechet distribution < 50 copies/ml
frec_50_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_50,
                                             estimate_threshold = 50, conversion_threshold = 1000, 
                                             distribution = "Frechet 45-54",
                                             survey = PHIA_vls_estimates_45_54$survey,
                                             empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                             empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                             empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

frec_50_45_54_phia$distribution <- rep("Frechet", dim(frec_50_45_54_phia)[1])
frec_50_45_54_phia$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_45_54_phia)[1])
frec_50_45_54_phia$agecat <- rep("45-54", dim(frec_50_45_54_phia)[1])
frec_50_45_54_phia$source <- rep("PHIA surveys", dim(frec_50_45_54_phia)[1])


#' Gamma distribution < 400 copies/ml
gamma_400_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_400,
                                          estimate_threshold = 400, conversion_threshold = 1000, 
                                          distribution = "Gamma 45-54", 
                                          survey = PHIA_vls_estimates_45_54$survey,
                                          empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000, 
                                          empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

gamma_400_45_54_phia$distribution <- rep("Gamma", dim(gamma_400_45_54_phia)[1])
gamma_400_45_54_phia$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_45_54_phia)[1])
gamma_400_45_54_phia$agecat <- rep("45-54", dim(gamma_400_45_54_phia)[1])
gamma_400_45_54_phia$source <- rep("PHIA surveys", dim(gamma_400_45_54_phia)[1])

#' Gamma distribution < 200 copies/ml
gamma_200_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_200,
                                          estimate_threshold = 200, conversion_threshold = 1000, 
                                          distribution = "Gamma 45-54",
                                          survey = PHIA_vls_estimates_45_54$survey,
                                          empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                          empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

gamma_200_45_54_phia$distribution <- rep("Gamma", dim(gamma_200_45_54_phia)[1])
gamma_200_45_54_phia$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_45_54_phia)[1])
gamma_200_45_54_phia$agecat <- rep("45-54", dim(gamma_200_45_54_phia)[1])
gamma_200_45_54_phia$source <- rep("PHIA surveys", dim(gamma_200_45_54_phia)[1])

#' Gamma distribution < 50 copies/ml
gamma_50_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_50,
                                         estimate_threshold = 50, conversion_threshold = 1000, 
                                         distribution = "Gamma 45-54",
                                         survey = PHIA_vls_estimates_45_54$survey,
                                         empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

gamma_50_45_54_phia$distribution <- rep("Gamma", dim(gamma_50_45_54_phia)[1])
gamma_50_45_54_phia$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_45_54_phia)[1])
gamma_50_45_54_phia$agecat <- rep("45-54", dim(gamma_50_45_54_phia)[1])
gamma_50_45_54_phia$source <- rep("PHIA surveys", dim(gamma_50_45_54_phia)[1])


#' lognormal distribution < 400 copies/ml
lognormal_400_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_400,
                                            estimate_threshold = 400, conversion_threshold = 1000, 
                                            distribution = "lognormal 45-54", 
                                            survey = PHIA_vls_estimates_45_54$survey,
                                            empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                            empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000, 
                                            empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

lognormal_400_45_54_phia$distribution <- rep("lognormal", dim(lognormal_400_45_54_phia)[1])
lognormal_400_45_54_phia$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_45_54_phia)[1])
lognormal_400_45_54_phia$agecat <- rep("45-54", dim(lognormal_400_45_54_phia)[1])
lognormal_400_45_54_phia$source <- rep("PHIA surveys", dim(lognormal_400_45_54_phia)[1])

#' lognormal distribution < 200 copies/ml
lognormal_200_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_200,
                                            estimate_threshold = 200, conversion_threshold = 1000, 
                                            distribution = "lognormal 45-54",
                                            survey = PHIA_vls_estimates_45_54$survey,
                                            empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                            empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                            empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

lognormal_200_45_54_phia$distribution <- rep("lognormal", dim(lognormal_200_45_54_phia)[1])
lognormal_200_45_54_phia$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_45_54_phia)[1])
lognormal_200_45_54_phia$agecat <- rep("45-54", dim(lognormal_200_45_54_phia)[1])
lognormal_200_45_54_phia$source <- rep("PHIA surveys", dim(lognormal_200_45_54_phia)[1])

#' lognormal distribution < 50 copies/ml
lognormal_50_45_54_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_50,
                                           estimate_threshold = 50, conversion_threshold = 1000, 
                                           distribution = "lognormal 45-54",
                                           survey = PHIA_vls_estimates_45_54$survey,
                                           empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                           empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                           empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

lognormal_50_45_54_phia$distribution <- rep("lognormal", dim(lognormal_50_45_54_phia)[1])
lognormal_50_45_54_phia$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_45_54_phia)[1])
lognormal_50_45_54_phia$agecat <- rep("45-54", dim(lognormal_50_45_54_phia)[1])
lognormal_50_45_54_phia$source <- rep("PHIA surveys", dim(lognormal_50_45_54_phia)[1])


#' calculate adjusted viral suppression estimates
#' Frechet distribution < 400 copies/ml
#' 55+
frec_400_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_400,
                                           estimate_threshold = 400, conversion_threshold = 1000, 
                                           distribution = "Frechet 55+", 
                                           survey = PHIA_vls_estimates_55$survey,
                                           empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                           empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000, 
                                           empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

frec_400_55_phia$distribution <- rep("Frechet", dim(frec_400_55_phia)[1])
frec_400_55_phia$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_55_phia)[1])
frec_400_55_phia$agecat <- rep("55+", dim(frec_400_55_phia)[1])
frec_400_55_phia$source <- rep("PHIA surveys", dim(frec_400_55_phia)[1])

#' Frechet distribution < 200 copies/ml
frec_200_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_200,
                                           estimate_threshold = 200, conversion_threshold = 1000, 
                                           distribution = "Frechet 55+",
                                           survey = PHIA_vls_estimates_55$survey,
                                           empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                           empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                           empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

frec_200_55_phia$distribution <- rep("Frechet", dim(frec_200_55_phia)[1])
frec_200_55_phia$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_55_phia)[1])
frec_200_55_phia$agecat <- rep("55+", dim(frec_200_55_phia)[1])
frec_200_55_phia$source <- rep("PHIA surveys", dim(frec_200_55_phia)[1])

#' Frechet distribution < 50 copies/ml
frec_50_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_50,
                                          estimate_threshold = 50, conversion_threshold = 1000, 
                                          distribution = "Frechet 55+",
                                          survey = PHIA_vls_estimates_55$survey,
                                          empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                          empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

frec_50_55_phia$distribution <- rep("Frechet", dim(frec_50_55_phia)[1])
frec_50_55_phia$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_55_phia)[1])
frec_50_55_phia$agecat <- rep("55+", dim(frec_50_55_phia)[1])
frec_50_55_phia$source <- rep("PHIA surveys", dim(frec_50_55_phia)[1])


#' Gamma distribution < 400 copies/ml
gamma_400_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "Gamma 55+", 
                                       survey = PHIA_vls_estimates_55$survey,
                                       empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000, 
                                       empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

gamma_400_55_phia$distribution <- rep("Gamma", dim(gamma_400_55_phia)[1])
gamma_400_55_phia$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_55_phia)[1])
gamma_400_55_phia$agecat <- rep("55+", dim(gamma_400_55_phia)[1])
gamma_400_55_phia$source <- rep("PHIA surveys", dim(gamma_400_55_phia)[1])

#' Gamma distribution < 200 copies/ml
gamma_200_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "Gamma 55+",
                                       survey = PHIA_vls_estimates_55$survey,
                                       empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

gamma_200_55_phia$distribution <- rep("Gamma", dim(gamma_200_55_phia)[1])
gamma_200_55_phia$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_55_phia)[1])
gamma_200_55_phia$agecat <- rep("55+", dim(gamma_200_55_phia)[1])
gamma_200_55_phia$source <- rep("PHIA surveys", dim(gamma_200_55_phia)[1])

#' Gamma distribution < 50 copies/ml
gamma_50_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "Gamma 55+",
                                      survey = PHIA_vls_estimates_55$survey,
                                      empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

gamma_50_55_phia$distribution <- rep("Gamma", dim(gamma_50_55_phia)[1])
gamma_50_55_phia$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_55_phia)[1])
gamma_50_55_phia$agecat <- rep("55+", dim(gamma_50_55_phia)[1])
gamma_50_55_phia$source <- rep("PHIA surveys", dim(gamma_50_55_phia)[1])


#' lognormal distribution < 400 copies/ml
lognormal_400_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "lognormal 55+", 
                                         survey = PHIA_vls_estimates_55$survey,
                                         empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

lognormal_400_55_phia$distribution <- rep("lognormal", dim(lognormal_400_55_phia)[1])
lognormal_400_55_phia$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_55_phia)[1])
lognormal_400_55_phia$agecat <- rep("55+", dim(lognormal_400_55_phia)[1])
lognormal_400_55_phia$source <- rep("PHIA surveys", dim(lognormal_400_55_phia)[1])

#' lognormal distribution < 200 copies/ml
lognormal_200_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "lognormal 55+",
                                         survey = PHIA_vls_estimates_55$survey,
                                         empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

lognormal_200_55_phia$distribution <- rep("lognormal", dim(lognormal_200_55_phia)[1])
lognormal_200_55_phia$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_55_phia)[1])
lognormal_200_55_phia$agecat <- rep("55+", dim(lognormal_200_55_phia)[1])
lognormal_200_55_phia$source <- rep("PHIA surveys", dim(lognormal_200_55_phia)[1])

#' lognormal distribution < 50 copies/ml
lognormal_50_55_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "lognormal 55+",
                                        survey = PHIA_vls_estimates_55$survey,
                                        empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

lognormal_50_55_phia$distribution <- rep("lognormal", dim(lognormal_50_55_phia)[1])
lognormal_50_55_phia$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_55_phia)[1])
lognormal_50_55_phia$agecat <- rep("55+", dim(lognormal_50_55_phia)[1])
lognormal_50_55_phia$source <- rep("PHIA surveys", dim(lognormal_50_55_phia)[1])


# using parameters from Johnson et al
# analyses repeated here
#' calculate adjusted viral suppression estimates
#' Frechet distribution < 400 copies/ml
#' 15-24
frec_400_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "Frechet", 
                                         survey = PHIA_vls_estimates_15_24$survey,
                                         empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

frec_400_15_24$distribution <- rep("Frechet", dim(frec_400_15_24)[1])
frec_400_15_24$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_15_24)[1])
frec_400_15_24$agecat <- rep("15-24", dim(frec_400_15_24)[1])
frec_400_15_24$source <- rep("Johnson et al", dim(frec_400_15_24)[1])

#' Frechet distribution < 200 copies/ml
frec_200_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "Frechet",
                                         survey = PHIA_vls_estimates_15_24$survey,
                                         empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

frec_200_15_24$distribution <- rep("Frechet", dim(frec_200_15_24)[1])
frec_200_15_24$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_15_24)[1])
frec_200_15_24$agecat <- rep("15-24", dim(frec_200_15_24)[1])
frec_200_15_24$source <- rep("Johnson et al", dim(frec_200_15_24)[1])

#' Frechet distribution < 50 copies/ml
frec_50_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "Frechet",
                                        survey = PHIA_vls_estimates_15_24$survey,
                                        empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

frec_50_15_24$distribution <- rep("Frechet", dim(frec_50_15_24)[1])
frec_50_15_24$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_15_24)[1])
frec_50_15_24$agecat <- rep("15-24", dim(frec_50_15_24)[1])
frec_50_15_24$source <- rep("Johnson et al", dim(frec_50_15_24)[1])


#' Gamma distribution < 400 copies/ml
gamma_400_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_400,
                                     estimate_threshold = 400, conversion_threshold = 1000, 
                                     distribution = "Gamma", 
                                     survey = PHIA_vls_estimates_15_24$survey,
                                     empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000, 
                                     empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

gamma_400_15_24$distribution <- rep("Gamma", dim(gamma_400_15_24)[1])
gamma_400_15_24$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_15_24)[1])
gamma_400_15_24$agecat <- rep("15-24", dim(gamma_400_15_24)[1])
gamma_400_15_24$source <- rep("Johnson et al", dim(gamma_400_15_24)[1])

#' Gamma distribution < 200 copies/ml
gamma_200_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_200,
                                     estimate_threshold = 200, conversion_threshold = 1000, 
                                     distribution = "Gamma",
                                     survey = PHIA_vls_estimates_15_24$survey,
                                     empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

gamma_200_15_24$distribution <- rep("Gamma", dim(gamma_200_15_24)[1])
gamma_200_15_24$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_15_24)[1])
gamma_200_15_24$agecat <- rep("15-24", dim(gamma_200_15_24)[1])
gamma_200_15_24$source <- rep("Johnson et al", dim(gamma_200_15_24)[1])

#' Gamma distribution < 50 copies/ml
gamma_50_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_50,
                                    estimate_threshold = 50, conversion_threshold = 1000, 
                                    distribution = "Gamma",
                                    survey = PHIA_vls_estimates_15_24$survey,
                                    empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

gamma_50_15_24$distribution <- rep("Gamma", dim(gamma_50_15_24)[1])
gamma_50_15_24$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_15_24)[1])
gamma_50_15_24$agecat <- rep("15-24", dim(gamma_50_15_24)[1])
gamma_50_15_24$source <- rep("Johnson et al", dim(gamma_50_15_24)[1])

#' lognormal distribution < 400 copies/ml
lognormal_400_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "lognormal", 
                                       survey = PHIA_vls_estimates_15_24$survey,
                                       empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000, 
                                       empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

lognormal_400_15_24$distribution <- rep("lognormal", dim(lognormal_400_15_24)[1])
lognormal_400_15_24$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_15_24)[1])
lognormal_400_15_24$agecat <- rep("15-24", dim(lognormal_400_15_24)[1])
lognormal_400_15_24$source <- rep("Johnson et al", dim(lognormal_400_15_24)[1])

#' lognormal distribution < 200 copies/ml
lognormal_200_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "lognormal",
                                       survey = PHIA_vls_estimates_15_24$survey,
                                       empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

lognormal_200_15_24$distribution <- rep("lognormal", dim(lognormal_200_15_24)[1])
lognormal_200_15_24$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_15_24)[1])
lognormal_200_15_24$agecat <- rep("15-24", dim(lognormal_200_15_24)[1])
lognormal_200_15_24$source <- rep("Johnson et al", dim(lognormal_200_15_24)[1])

#' lognormal distribution < 50 copies/ml
lognormal_50_15_24 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_15_24$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "lognormal",
                                      survey = PHIA_vls_estimates_15_24$survey,
                                      empirical_estimate = PHIA_vls_estimates_15_24$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_15_24$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_15_24$upper_ci_1000)

lognormal_50_15_24$distribution <- rep("lognormal", dim(lognormal_50_15_24)[1])
lognormal_50_15_24$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_15_24)[1])
lognormal_50_15_24$agecat <- rep("15-24", dim(lognormal_50_15_24)[1])
lognormal_50_15_24$source <- rep("Johnson et al", dim(lognormal_50_15_24)[1])


#' calculate adjusted viral suppression estimates
#' Frechet distribution < 400 copies/ml
#' 25-34
frec_400_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "Frechet", 
                                         survey = PHIA_vls_estimates_25_34$survey,
                                         empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

frec_400_25_34$distribution <- rep("Frechet", dim(frec_400_25_34)[1])
frec_400_25_34$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_25_34)[1])
frec_400_25_34$agecat <- rep("25-34", dim(frec_400_25_34)[1])
frec_400_25_34$source <- rep("Johnson et al", dim(frec_400_25_34)[1])

#' Frechet distribution < 200 copies/ml
frec_200_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "Frechet",
                                         survey = PHIA_vls_estimates_25_34$survey,
                                         empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

frec_200_25_34$distribution <- rep("Frechet", dim(frec_200_25_34)[1])
frec_200_25_34$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_25_34)[1])
frec_200_25_34$agecat <- rep("25-34", dim(frec_200_25_34)[1])
frec_200_25_34$source <- rep("Johnson et al", dim(frec_200_25_34)[1])

#' Frechet distribution < 50 copies/ml
frec_50_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "Frechet",
                                        survey = PHIA_vls_estimates_25_34$survey,
                                        empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

frec_50_25_34$distribution <- rep("Frechet", dim(frec_50_25_34)[1])
frec_50_25_34$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_25_34)[1])
frec_50_25_34$agecat <- rep("25-34", dim(frec_50_25_34)[1])
frec_50_25_34$source <- rep("Johnson et al", dim(frec_50_25_34)[1])


#' Gamma distribution < 400 copies/ml
gamma_400_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_400,
                                     estimate_threshold = 400, conversion_threshold = 1000, 
                                     distribution = "Gamma", 
                                     survey = PHIA_vls_estimates_25_34$survey,
                                     empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000, 
                                     empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

gamma_400_25_34$distribution <- rep("Gamma", dim(gamma_400_25_34)[1])
gamma_400_25_34$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_25_34)[1])
gamma_400_25_34$agecat <- rep("25-34", dim(gamma_400_25_34)[1])
gamma_400_25_34$source <- rep("Johnson et al", dim(gamma_400_25_34)[1])

#' Gamma distribution < 200 copies/ml
gamma_200_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_200,
                                     estimate_threshold = 200, conversion_threshold = 1000, 
                                     distribution = "Gamma",
                                     survey = PHIA_vls_estimates_25_34$survey,
                                     empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

gamma_200_25_34$distribution <- rep("Gamma", dim(gamma_200_25_34)[1])
gamma_200_25_34$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_25_34)[1])
gamma_200_25_34$agecat <- rep("25-34", dim(gamma_200_25_34)[1])
gamma_200_25_34$source <- rep("Johnson et al", dim(gamma_200_25_34)[1])

#' Gamma distribution < 50 copies/ml
gamma_50_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_50,
                                    estimate_threshold = 50, conversion_threshold = 1000, 
                                    distribution = "Gamma",
                                    survey = PHIA_vls_estimates_25_34$survey,
                                    empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

gamma_50_25_34$distribution <- rep("Gamma", dim(gamma_50_25_34)[1])
gamma_50_25_34$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_25_34)[1])
gamma_50_25_34$agecat <- rep("25-34", dim(gamma_50_25_34)[1])
gamma_50_25_34$source <- rep("Johnson et al", dim(gamma_50_25_34)[1])


#' lognormal distribution < 400 copies/ml
lognormal_400_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "lognormal", 
                                       survey = PHIA_vls_estimates_25_34$survey,
                                       empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000, 
                                       empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

lognormal_400_25_34$distribution <- rep("lognormal", dim(lognormal_400_25_34)[1])
lognormal_400_25_34$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_25_34)[1])
lognormal_400_25_34$agecat <- rep("25-34", dim(lognormal_400_25_34)[1])
lognormal_400_25_34$source <- rep("Johnson et al", dim(lognormal_400_25_34)[1])

#' lognormal distribution < 200 copies/ml
lognormal_200_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "lognormal",
                                       survey = PHIA_vls_estimates_25_34$survey,
                                       empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

lognormal_200_25_34$distribution <- rep("lognormal", dim(lognormal_200_25_34)[1])
lognormal_200_25_34$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_25_34)[1])
lognormal_200_25_34$agecat <- rep("25-34", dim(lognormal_200_25_34)[1])
lognormal_200_25_34$source <- rep("Johnson et al", dim(lognormal_200_25_34)[1])

#' lognormal distribution < 50 copies/ml
lognormal_50_25_34 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_25_34$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "lognormal",
                                      survey = PHIA_vls_estimates_25_34$survey,
                                      empirical_estimate = PHIA_vls_estimates_25_34$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_25_34$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_25_34$upper_ci_1000)

lognormal_50_25_34$distribution <- rep("lognormal", dim(lognormal_50_25_34)[1])
lognormal_50_25_34$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_25_34)[1])
lognormal_50_25_34$agecat <- rep("25-34", dim(lognormal_50_25_34)[1])
lognormal_50_25_34$source <- rep("Johnson et al", dim(lognormal_50_25_34)[1])


#' calculate adjusted viral suppression estimates
#' Frechet distribution < 400 copies/ml
#' 35-44
frec_400_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "Frechet", 
                                         survey = PHIA_vls_estimates_35_44$survey,
                                         empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

frec_400_35_44$distribution <- rep("Frechet", dim(frec_400_35_44)[1])
frec_400_35_44$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_35_44)[1])
frec_400_35_44$agecat <- rep("35-44", dim(frec_400_35_44)[1])
frec_400_35_44$source <- rep("Johnson et al", dim(frec_400_35_44)[1])

#' Frechet distribution < 200 copies/ml
frec_200_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "Frechet",
                                         survey = PHIA_vls_estimates_35_44$survey,
                                         empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

frec_200_35_44$distribution <- rep("Frechet", dim(frec_200_35_44)[1])
frec_200_35_44$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_35_44)[1])
frec_200_35_44$agecat <- rep("35-44", dim(frec_200_35_44)[1])
frec_200_35_44$source <- rep("Johnson et al", dim(frec_200_35_44)[1])

#' Frechet distribution < 50 copies/ml
frec_50_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "Frechet",
                                        survey = PHIA_vls_estimates_35_44$survey,
                                        empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

frec_50_35_44$distribution <- rep("Frechet", dim(frec_50_35_44)[1])
frec_50_35_44$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_35_44)[1])
frec_50_35_44$agecat <- rep("35-44", dim(frec_50_35_44)[1])
frec_50_35_44$source <- rep("Johnson et al", dim(frec_50_35_44)[1])


#' Gamma distribution < 400 copies/ml
gamma_400_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_400,
                                     estimate_threshold = 400, conversion_threshold = 1000, 
                                     distribution = "Gamma", 
                                     survey = PHIA_vls_estimates_35_44$survey,
                                     empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000, 
                                     empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

gamma_400_35_44$distribution <- rep("Gamma", dim(gamma_400_35_44)[1])
gamma_400_35_44$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_35_44)[1])
gamma_400_35_44$agecat <- rep("35-44", dim(gamma_400_35_44)[1])
gamma_400_35_44$source <- rep("Johnson et al", dim(gamma_400_35_44)[1])

#' Gamma distribution < 200 copies/ml
gamma_200_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_200,
                                     estimate_threshold = 200, conversion_threshold = 1000, 
                                     distribution = "Gamma",
                                     survey = PHIA_vls_estimates_35_44$survey,
                                     empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

gamma_200_35_44$distribution <- rep("Gamma", dim(gamma_200_35_44)[1])
gamma_200_35_44$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_35_44)[1])
gamma_200_35_44$agecat <- rep("35-44", dim(gamma_200_35_44)[1])
gamma_200_35_44$source <- rep("Johnson et al", dim(gamma_200_35_44)[1])

#' Gamma distribution < 50 copies/ml
gamma_50_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_50,
                                    estimate_threshold = 50, conversion_threshold = 1000, 
                                    distribution = "Gamma",
                                    survey = PHIA_vls_estimates_35_44$survey,
                                    empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

gamma_50_35_44$distribution <- rep("Gamma", dim(gamma_50_35_44)[1])
gamma_50_35_44$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_35_44)[1])
gamma_50_35_44$agecat <- rep("35-44", dim(gamma_50_35_44)[1])
gamma_50_35_44$source <- rep("Johnson et al", dim(gamma_50_35_44)[1])


#' lognormal distribution < 400 copies/ml
lognormal_400_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "lognormal", 
                                       survey = PHIA_vls_estimates_35_44$survey,
                                       empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000, 
                                       empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

lognormal_400_35_44$distribution <- rep("lognormal", dim(lognormal_400_35_44)[1])
lognormal_400_35_44$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_35_44)[1])
lognormal_400_35_44$agecat <- rep("35-44", dim(lognormal_400_35_44)[1])
lognormal_400_35_44$source <- rep("Johnson et al", dim(lognormal_400_35_44)[1])

#' lognormal distribution < 200 copies/ml
lognormal_200_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "lognormal",
                                       survey = PHIA_vls_estimates_35_44$survey,
                                       empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

lognormal_200_35_44$distribution <- rep("lognormal", dim(lognormal_200_35_44)[1])
lognormal_200_35_44$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_35_44)[1])
lognormal_200_35_44$agecat <- rep("35-44", dim(lognormal_200_35_44)[1])
lognormal_200_35_44$source <- rep("Johnson et al", dim(lognormal_200_35_44)[1])

#' lognormal distribution < 50 copies/ml
lognormal_50_35_44 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_35_44$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "lognormal",
                                      survey = PHIA_vls_estimates_35_44$survey,
                                      empirical_estimate = PHIA_vls_estimates_35_44$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_35_44$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_35_44$upper_ci_1000)

lognormal_50_35_44$distribution <- rep("lognormal", dim(lognormal_50_35_44)[1])
lognormal_50_35_44$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_35_44)[1])
lognormal_50_35_44$agecat <- rep("35-44", dim(lognormal_50_35_44)[1])
lognormal_50_35_44$source <- rep("Johnson et al", dim(lognormal_50_35_44)[1])


#' calculate adjusted viral suppression estimates
#' Frechet distribution < 400 copies/ml
#' 45-54
frec_400_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "Frechet", 
                                         survey = PHIA_vls_estimates_45_54$survey,
                                         empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

frec_400_45_54$distribution <- rep("Frechet", dim(frec_400_45_54)[1])
frec_400_45_54$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_45_54)[1])
frec_400_45_54$agecat <- rep("45-54", dim(frec_400_45_54)[1])
frec_400_45_54$source <- rep("Johnson et al", dim(frec_400_45_54)[1])

#' Frechet distribution < 200 copies/ml
frec_200_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "Frechet",
                                         survey = PHIA_vls_estimates_45_54$survey,
                                         empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

frec_200_45_54$distribution <- rep("Frechet", dim(frec_200_45_54)[1])
frec_200_45_54$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_45_54)[1])
frec_200_45_54$agecat <- rep("45-54", dim(frec_200_45_54)[1])
frec_200_45_54$source <- rep("Johnson et al", dim(frec_200_45_54)[1])

#' Frechet distribution < 50 copies/ml
frec_50_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "Frechet",
                                        survey = PHIA_vls_estimates_45_54$survey,
                                        empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

frec_50_45_54$distribution <- rep("Frechet", dim(frec_50_45_54)[1])
frec_50_45_54$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_45_54)[1])
frec_50_45_54$agecat <- rep("45-54", dim(frec_50_45_54)[1])
frec_50_45_54$source <- rep("Johnson et al", dim(frec_50_45_54)[1])


#' Gamma distribution < 400 copies/ml
gamma_400_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_400,
                                     estimate_threshold = 400, conversion_threshold = 1000, 
                                     distribution = "Gamma", 
                                     survey = PHIA_vls_estimates_45_54$survey,
                                     empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000, 
                                     empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

gamma_400_45_54$distribution <- rep("Gamma", dim(gamma_400_45_54)[1])
gamma_400_45_54$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_45_54)[1])
gamma_400_45_54$agecat <- rep("45-54", dim(gamma_400_45_54)[1])
gamma_400_45_54$source <- rep("Johnson et al", dim(gamma_400_45_54)[1])

#' Gamma distribution < 200 copies/ml
gamma_200_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_200,
                                     estimate_threshold = 200, conversion_threshold = 1000, 
                                     distribution = "Gamma",
                                     survey = PHIA_vls_estimates_45_54$survey,
                                     empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

gamma_200_45_54$distribution <- rep("Gamma", dim(gamma_200_45_54)[1])
gamma_200_45_54$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_45_54)[1])
gamma_200_45_54$agecat <- rep("45-54", dim(gamma_200_45_54)[1])
gamma_200_45_54$source <- rep("Johnson et al", dim(gamma_200_45_54)[1])

#' Gamma distribution < 50 copies/ml
gamma_50_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_50,
                                    estimate_threshold = 50, conversion_threshold = 1000, 
                                    distribution = "Gamma",
                                    survey = PHIA_vls_estimates_45_54$survey,
                                    empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

gamma_50_45_54$distribution <- rep("Gamma", dim(gamma_50_45_54)[1])
gamma_50_45_54$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_45_54)[1])
gamma_50_45_54$agecat <- rep("45-54", dim(gamma_50_45_54)[1])
gamma_50_45_54$source <- rep("Johnson et al", dim(gamma_50_45_54)[1])


#' lognormal distribution < 400 copies/ml
lognormal_400_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "lognormal", 
                                       survey = PHIA_vls_estimates_45_54$survey,
                                       empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000, 
                                       empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

lognormal_400_45_54$distribution <- rep("lognormal", dim(lognormal_400_45_54)[1])
lognormal_400_45_54$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_45_54)[1])
lognormal_400_45_54$agecat <- rep("45-54", dim(lognormal_400_45_54)[1])
lognormal_400_45_54$source <- rep("Johnson et al", dim(lognormal_400_45_54)[1])

#' lognormal distribution < 200 copies/ml
lognormal_200_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "lognormal",
                                       survey = PHIA_vls_estimates_45_54$survey,
                                       empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

lognormal_200_45_54$distribution <- rep("lognormal", dim(lognormal_200_45_54)[1])
lognormal_200_45_54$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_45_54)[1])
lognormal_200_45_54$agecat <- rep("45-54", dim(lognormal_200_45_54)[1])
lognormal_200_45_54$source <- rep("Johnson et al", dim(lognormal_200_45_54)[1])

#' lognormal distribution < 50 copies/ml
lognormal_50_45_54 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_45_54$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "lognormal",
                                      survey = PHIA_vls_estimates_45_54$survey,
                                      empirical_estimate = PHIA_vls_estimates_45_54$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_45_54$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_45_54$upper_ci_1000)

lognormal_50_45_54$distribution <- rep("lognormal", dim(lognormal_50_45_54)[1])
lognormal_50_45_54$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_45_54)[1])
lognormal_50_45_54$agecat <- rep("45-54", dim(lognormal_50_45_54)[1])
lognormal_50_45_54$source <- rep("Johnson et al", dim(lognormal_50_45_54)[1])


#' calculate adjusted viral suppression estimates
#' Frechet distribution < 400 copies/ml
#' 55+
frec_400_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_400,
                                      estimate_threshold = 400, conversion_threshold = 1000, 
                                      distribution = "Frechet", 
                                      survey = PHIA_vls_estimates_55$survey,
                                      empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000, 
                                      empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

frec_400_55$distribution <- rep("Frechet", dim(frec_400_55)[1])
frec_400_55$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_55)[1])
frec_400_55$agecat <- rep("55+", dim(frec_400_55)[1])
frec_400_55$source <- rep("Johnson et al", dim(frec_400_55)[1])

#' Frechet distribution < 200 copies/ml
frec_200_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_200,
                                      estimate_threshold = 200, conversion_threshold = 1000, 
                                      distribution = "Frechet",
                                      survey = PHIA_vls_estimates_55$survey,
                                      empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

frec_200_55$distribution <- rep("Frechet", dim(frec_200_55)[1])
frec_200_55$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_55)[1])
frec_200_55$agecat <- rep("55+", dim(frec_200_55)[1])
frec_200_55$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_55)[1])
frec_200_55$source <- rep("Johnson et al", dim(frec_200_55)[1])

#' Frechet distribution < 50 copies/ml
frec_50_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_50,
                                     estimate_threshold = 50, conversion_threshold = 1000, 
                                     distribution = "Frechet",
                                     survey = PHIA_vls_estimates_55$survey,
                                     empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

frec_50_55$distribution <- rep("Frechet", dim(frec_50_55)[1])
frec_50_55$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_55)[1])
frec_50_55$agecat <- rep("55+", dim(frec_50_55)[1])
frec_50_55$source <- rep("Johnson et al", dim(frec_50_55)[1])


#' Gamma distribution < 400 copies/ml
gamma_400_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_400,
                                  estimate_threshold = 400, conversion_threshold = 1000, 
                                  distribution = "Gamma", 
                                  survey = PHIA_vls_estimates_55$survey,
                                  empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                  empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000, 
                                  empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

gamma_400_55$distribution <- rep("Gamma", dim(gamma_400_55)[1])
gamma_400_55$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_55)[1])
gamma_400_55$agecat <- rep("55+", dim(gamma_400_55)[1])
gamma_400_55$source <- rep("Johnson et al", dim(gamma_400_55)[1])

#' Gamma distribution < 200 copies/ml
gamma_200_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_200,
                                  estimate_threshold = 200, conversion_threshold = 1000, 
                                  distribution = "Gamma",
                                  survey = PHIA_vls_estimates_55$survey,
                                  empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                  empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                  empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

gamma_200_55$distribution <- rep("Gamma", dim(gamma_200_55)[1])
gamma_200_55$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_55)[1])
gamma_200_55$agecat <- rep("55+", dim(gamma_200_55)[1])
gamma_200_55$source <- rep("Johnson et al", dim(gamma_200_55)[1])

#' Gamma distribution < 50 copies/ml
gamma_50_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_50,
                                 estimate_threshold = 50, conversion_threshold = 1000, 
                                 distribution = "Gamma",
                                 survey = PHIA_vls_estimates_55$survey,
                                 empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                 empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                 empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

gamma_50_55$distribution <- rep("Gamma", dim(gamma_50_55)[1])
gamma_50_55$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_55)[1])
gamma_50_55$agecat <- rep("55+", dim(gamma_50_55)[1])
gamma_50_55$source <- rep("Johnson et al", dim(gamma_50_55)[1])


#' lognormal distribution < 400 copies/ml
lognormal_400_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_400,
                                    estimate_threshold = 400, conversion_threshold = 1000, 
                                    distribution = "lognormal", 
                                    survey = PHIA_vls_estimates_55$survey,
                                    empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000, 
                                    empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

lognormal_400_55$distribution <- rep("lognormal", dim(lognormal_400_55)[1])
lognormal_400_55$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_55)[1])
lognormal_400_55$agecat <- rep("55+", dim(lognormal_400_55)[1])
lognormal_400_55$source <- rep("Johnson et al", dim(lognormal_400_55)[1])

#' lognormal distribution < 200 copies/ml
lognormal_200_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_200,
                                    estimate_threshold = 200, conversion_threshold = 1000, 
                                    distribution = "lognormal",
                                    survey = PHIA_vls_estimates_55$survey,
                                    empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

lognormal_200_55$distribution <- rep("lognormal", dim(lognormal_200_55)[1])
lognormal_200_55$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_55)[1])
lognormal_200_55$agecat <- rep("55+", dim(lognormal_200_55)[1])
lognormal_200_55$source <- rep("Johnson et al", dim(lognormal_200_55)[1])

#' lognormal distribution < 50 copies/ml
lognormal_50_55 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_55$estimate_50,
                                   estimate_threshold = 50, conversion_threshold = 1000, 
                                   distribution = "lognormal",
                                   survey = PHIA_vls_estimates_55$survey,
                                   empirical_estimate = PHIA_vls_estimates_55$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates_55$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates_55$upper_ci_1000)

lognormal_50_55$distribution <- rep("lognormal", dim(lognormal_50_55)[1])
lognormal_50_55$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_55)[1])
lognormal_50_55$agecat <- rep("55+", dim(lognormal_50_55)[1])
lognormal_50_55$source <- rep("Johnson et al", dim(lognormal_50_55)[1])


#' calculate the root mean squared error between empirical and adjusted estimates
#' 15-24
sqrt(mean((frec_50_15_24_phia$empirical_estimate - frec_50_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_15_24_phia$empirical_estimate - gamma_50_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_15_24_phia$empirical_estimate - lognormal_50_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_50_15_24$empirical_estimate - frec_50_15_24$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_15_24$empirical_estimate - gamma_50_15_24$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_15_24$empirical_estimate - lognormal_50_15_24$adjusted_vls_estimate)^2))

sqrt(mean((frec_200_15_24_phia$empirical_estimate - frec_200_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_15_24_phia$empirical_estimate - gamma_200_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_15_24_phia$empirical_estimate - lognormal_200_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_200_15_24$empirical_estimate - frec_200_15_24$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_15_24$empirical_estimate - gamma_200_15_24$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_15_24$empirical_estimate - lognormal_200_15_24$adjusted_vls_estimate)^2))

sqrt(mean((frec_400_15_24_phia$empirical_estimate - frec_400_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_15_24_phia$empirical_estimate - gamma_400_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_15_24_phia$empirical_estimate - lognormal_400_15_24_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_400_15_24$empirical_estimate - frec_400_15_24$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_15_24$empirical_estimate - gamma_400_15_24$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_15_24$empirical_estimate - lognormal_400_15_24$adjusted_vls_estimate)^2))

# 25-34
sqrt(mean((frec_50_25_34_phia$empirical_estimate - frec_50_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_25_34_phia$empirical_estimate - gamma_50_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_25_34_phia$empirical_estimate - lognormal_50_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_50_25_34$empirical_estimate - frec_50_25_34$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_25_34$empirical_estimate - gamma_50_25_34$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_25_34$empirical_estimate - lognormal_50_25_34$adjusted_vls_estimate)^2))

sqrt(mean((frec_200_25_34_phia$empirical_estimate - frec_200_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_25_34_phia$empirical_estimate - gamma_200_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_25_34_phia$empirical_estimate - lognormal_200_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_200_25_34$empirical_estimate - frec_200_25_34$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_25_34$empirical_estimate - gamma_200_25_34$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_25_34$empirical_estimate - lognormal_200_25_34$adjusted_vls_estimate)^2))

sqrt(mean((frec_400_25_34_phia$empirical_estimate - frec_400_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_25_34_phia$empirical_estimate - gamma_400_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_25_34_phia$empirical_estimate - lognormal_400_25_34_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_400_25_34$empirical_estimate - frec_400_25_34$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_25_34$empirical_estimate - gamma_400_25_34$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_25_34$empirical_estimate - lognormal_400_25_34$adjusted_vls_estimate)^2))

# 35-44
sqrt(mean((frec_50_35_44_phia$empirical_estimate - frec_50_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_35_44_phia$empirical_estimate - gamma_50_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_35_44_phia$empirical_estimate - lognormal_50_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_50_35_44$empirical_estimate - frec_50_35_44$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_35_44$empirical_estimate - gamma_50_35_44$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_35_44$empirical_estimate - lognormal_50_35_44$adjusted_vls_estimate)^2))

sqrt(mean((frec_200_35_44_phia$empirical_estimate - frec_200_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_35_44_phia$empirical_estimate - gamma_200_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_35_44_phia$empirical_estimate - lognormal_200_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_200_35_44$empirical_estimate - frec_200_35_44$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_35_44$empirical_estimate - gamma_200_35_44$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_35_44$empirical_estimate - lognormal_200_35_44$adjusted_vls_estimate)^2))

sqrt(mean((frec_400_35_44_phia$empirical_estimate - frec_400_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_35_44_phia$empirical_estimate - gamma_400_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_35_44_phia$empirical_estimate - lognormal_400_35_44_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_400_35_44$empirical_estimate - frec_400_35_44$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_35_44$empirical_estimate - gamma_400_35_44$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_35_44$empirical_estimate - lognormal_400_35_44$adjusted_vls_estimate)^2))

# 45-54
sqrt(mean((frec_50_45_54_phia$empirical_estimate - frec_50_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_45_54_phia$empirical_estimate - gamma_50_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_45_54_phia$empirical_estimate - lognormal_50_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_50_45_54$empirical_estimate - frec_50_45_54$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_45_54$empirical_estimate - gamma_50_45_54$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_45_54$empirical_estimate - lognormal_50_45_54$adjusted_vls_estimate)^2))

sqrt(mean((frec_200_45_54_phia$empirical_estimate - frec_200_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_45_54_phia$empirical_estimate - gamma_200_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_45_54_phia$empirical_estimate - lognormal_200_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_200_45_54$empirical_estimate - frec_200_45_54$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_45_54$empirical_estimate - gamma_200_45_54$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_45_54$empirical_estimate - lognormal_200_45_54$adjusted_vls_estimate)^2))

sqrt(mean((frec_400_45_54_phia$empirical_estimate - frec_400_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_45_54_phia$empirical_estimate - gamma_400_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_45_54_phia$empirical_estimate - lognormal_400_45_54_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_400_45_54$empirical_estimate - frec_400_45_54$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_45_54$empirical_estimate - gamma_400_45_54$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_45_54$empirical_estimate - lognormal_400_45_54$adjusted_vls_estimate)^2))

# 55+
sqrt(mean((frec_50_55_phia$empirical_estimate - frec_50_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_55_phia$empirical_estimate - gamma_50_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_55_phia$empirical_estimate - lognormal_50_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_50_55$empirical_estimate - frec_50_55$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_55$empirical_estimate - gamma_50_55$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_55$empirical_estimate - lognormal_50_55$adjusted_vls_estimate)^2))

sqrt(mean((frec_200_55_phia$empirical_estimate - frec_200_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_55_phia$empirical_estimate - gamma_200_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_55_phia$empirical_estimate - lognormal_200_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_200_55$empirical_estimate - frec_200_55$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_55$empirical_estimate - gamma_200_55$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_55$empirical_estimate - lognormal_200_55$adjusted_vls_estimate)^2))

sqrt(mean((frec_400_55_phia$empirical_estimate - frec_400_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_55_phia$empirical_estimate - gamma_400_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_55_phia$empirical_estimate - lognormal_400_55_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_400_55$empirical_estimate - frec_400_55$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_55$empirical_estimate - gamma_400_55$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_55$empirical_estimate - lognormal_400_55$adjusted_vls_estimate)^2))


# average bias
#' 15-24
mean(frec_50_15_24_phia$adjusted_vls_estimate - frec_50_15_24_phia$empirical_estimate)
mean(gamma_50_15_24_phia$adjusted_vls_estimate - gamma_50_15_24_phia$empirical_estimate)
mean(lognormal_50_15_24_phia$adjusted_vls_estimate - lognormal_50_15_24_phia$empirical_estimate)
mean(frec_50_15_24$adjusted_vls_estimate - frec_50_15_24$empirical_estimate)
mean(gamma_50_15_24$adjusted_vls_estimate - gamma_50_15_24$empirical_estimate)
mean(lognormal_50_15_24$adjusted_vls_estimate - lognormal_50_15_24$empirical_estimate)

mean(frec_200_15_24_phia$adjusted_vls_estimate - frec_200_15_24_phia$empirical_estimate)
mean(gamma_200_15_24_phia$adjusted_vls_estimate - gamma_200_15_24_phia$empirical_estimate)
mean(lognormal_200_15_24_phia$adjusted_vls_estimate - lognormal_200_15_24_phia$empirical_estimate)
mean(frec_200_15_24$adjusted_vls_estimate - frec_200_15_24$empirical_estimate)
mean(gamma_200_15_24$adjusted_vls_estimate - gamma_200_15_24$empirical_estimate)
mean(lognormal_200_15_24$adjusted_vls_estimate - lognormal_200_15_24$empirical_estimate)

mean(frec_400_15_24_phia$adjusted_vls_estimate - frec_400_15_24_phia$empirical_estimate)
mean(gamma_400_15_24_phia$adjusted_vls_estimate - gamma_400_15_24_phia$empirical_estimate)
mean(lognormal_400_15_24_phia$adjusted_vls_estimate - lognormal_400_15_24_phia$empirical_estimate)
mean(frec_400_15_24$adjusted_vls_estimate - frec_400_15_24$empirical_estimate)
mean(gamma_400_15_24$adjusted_vls_estimate - gamma_400_15_24$empirical_estimate)
mean(lognormal_400_15_24$adjusted_vls_estimate - lognormal_400_15_24$empirical_estimate)

# 25-34
mean(frec_50_25_34_phia$adjusted_vls_estimate - frec_50_25_34_phia$empirical_estimate)
mean(gamma_50_25_34_phia$adjusted_vls_estimate - gamma_50_25_34_phia$empirical_estimate)
mean(lognormal_50_25_34_phia$adjusted_vls_estimate - lognormal_50_25_34_phia$empirical_estimate)
mean(frec_50_25_34$adjusted_vls_estimate - frec_50_25_34$empirical_estimate)
mean(gamma_50_25_34$adjusted_vls_estimate - gamma_50_25_34$empirical_estimate)
mean(lognormal_50_25_34$adjusted_vls_estimate - lognormal_50_25_34$empirical_estimate)

mean(frec_200_25_34_phia$adjusted_vls_estimate - frec_200_25_34_phia$empirical_estimate)
mean(gamma_200_25_34_phia$adjusted_vls_estimate - gamma_200_25_34_phia$empirical_estimate)
mean(lognormal_200_25_34_phia$adjusted_vls_estimate - lognormal_200_25_34_phia$empirical_estimate)
mean(frec_200_25_34$adjusted_vls_estimate - frec_200_25_34$empirical_estimate)
mean(gamma_200_25_34$adjusted_vls_estimate - gamma_200_25_34$empirical_estimate)
mean(lognormal_200_25_34$adjusted_vls_estimate - lognormal_200_25_34$empirical_estimate)

mean(frec_400_25_34_phia$adjusted_vls_estimate - frec_400_25_34_phia$empirical_estimate)
mean(gamma_400_25_34_phia$adjusted_vls_estimate - gamma_400_25_34_phia$empirical_estimate)
mean(lognormal_400_25_34_phia$adjusted_vls_estimate - lognormal_400_25_34_phia$empirical_estimate)
mean(frec_400_25_34$adjusted_vls_estimate - frec_400_25_34$empirical_estimate)
mean(gamma_400_25_34$adjusted_vls_estimate - gamma_400_25_34$empirical_estimate)
mean(lognormal_400_25_34$adjusted_vls_estimate - lognormal_400_25_34$empirical_estimate)

# 35-44
mean(frec_50_35_44_phia$adjusted_vls_estimate - frec_50_35_44_phia$empirical_estimate)
mean(gamma_50_35_44_phia$adjusted_vls_estimate - gamma_50_35_44_phia$empirical_estimate)
mean(lognormal_50_35_44_phia$adjusted_vls_estimate - lognormal_50_35_44_phia$empirical_estimate)
mean(frec_50_35_44$adjusted_vls_estimate - frec_50_35_44$empirical_estimate)
mean(gamma_50_35_44$adjusted_vls_estimate - gamma_50_35_44$empirical_estimate)
mean(lognormal_50_35_44$adjusted_vls_estimate - lognormal_50_35_44$empirical_estimate)

mean(frec_200_35_44_phia$adjusted_vls_estimate - frec_200_35_44_phia$empirical_estimate)
mean(gamma_200_35_44_phia$adjusted_vls_estimate - gamma_200_35_44_phia$empirical_estimate)
mean(lognormal_200_35_44_phia$adjusted_vls_estimate - lognormal_200_35_44_phia$empirical_estimate)
mean(frec_200_35_44$adjusted_vls_estimate - frec_200_35_44$empirical_estimate)
mean(gamma_200_35_44$adjusted_vls_estimate - gamma_200_35_44$empirical_estimate)
mean(lognormal_200_35_44$adjusted_vls_estimate - lognormal_200_35_44$empirical_estimate)

mean(frec_400_35_44_phia$adjusted_vls_estimate - frec_400_35_44_phia$empirical_estimate)
mean(gamma_400_35_44_phia$adjusted_vls_estimate - gamma_400_35_44_phia$empirical_estimate)
mean(lognormal_400_35_44_phia$adjusted_vls_estimate - lognormal_400_35_44_phia$empirical_estimate)
mean(frec_400_35_44$adjusted_vls_estimate - frec_400_35_44$empirical_estimate)
mean(gamma_400_35_44$adjusted_vls_estimate - gamma_400_35_44$empirical_estimate)
mean(lognormal_400_35_44$adjusted_vls_estimate - lognormal_400_35_44$empirical_estimate)

# 45-54
mean(frec_50_45_54_phia$adjusted_vls_estimate - frec_50_45_54_phia$empirical_estimate)
mean(gamma_50_45_54_phia$adjusted_vls_estimate - gamma_50_45_54_phia$empirical_estimate)
mean(lognormal_50_45_54_phia$adjusted_vls_estimate - lognormal_50_45_54_phia$empirical_estimate)
mean(frec_50_45_54$adjusted_vls_estimate - frec_50_45_54$empirical_estimate)
mean(gamma_50_45_54$adjusted_vls_estimate - gamma_50_45_54$empirical_estimate)
mean(lognormal_50_45_54$adjusted_vls_estimate - lognormal_50_45_54$empirical_estimate)

mean(frec_200_45_54_phia$adjusted_vls_estimate - frec_200_45_54_phia$empirical_estimate)
mean(gamma_200_45_54_phia$adjusted_vls_estimate - gamma_200_45_54_phia$empirical_estimate)
mean(lognormal_200_45_54_phia$adjusted_vls_estimate - lognormal_200_45_54_phia$empirical_estimate)
mean(frec_200_45_54$adjusted_vls_estimate - frec_200_45_54$empirical_estimate)
mean(gamma_200_45_54$adjusted_vls_estimate - gamma_200_45_54$empirical_estimate)
mean(lognormal_200_45_54$adjusted_vls_estimate - lognormal_200_45_54$empirical_estimate)

mean(frec_400_45_54_phia$adjusted_vls_estimate - frec_400_45_54_phia$empirical_estimate)
mean(gamma_400_45_54_phia$adjusted_vls_estimate - gamma_400_45_54_phia$empirical_estimate)
mean(lognormal_400_45_54_phia$adjusted_vls_estimate - lognormal_400_45_54_phia$empirical_estimate)
mean(frec_400_45_54$adjusted_vls_estimate - frec_400_45_54$empirical_estimate)
mean(gamma_400_45_54$adjusted_vls_estimate - gamma_400_45_54$empirical_estimate)
mean(lognormal_400_45_54$adjusted_vls_estimate - lognormal_400_45_54$empirical_estimate)

# 55+
mean(frec_50_55_phia$adjusted_vls_estimate - frec_50_55_phia$empirical_estimate)
mean(gamma_50_55_phia$adjusted_vls_estimate - gamma_50_55_phia$empirical_estimate)
mean(lognormal_50_55_phia$adjusted_vls_estimate - lognormal_50_55_phia$empirical_estimate)
mean(frec_50_55$adjusted_vls_estimate - frec_50_55$empirical_estimate)
mean(gamma_50_55$adjusted_vls_estimate - gamma_50_55$empirical_estimate)
mean(lognormal_50_55$adjusted_vls_estimate - lognormal_50_55$empirical_estimate)

mean(frec_200_55_phia$adjusted_vls_estimate - frec_200_55_phia$empirical_estimate)
mean(gamma_200_55_phia$adjusted_vls_estimate - gamma_200_55_phia$empirical_estimate)
mean(lognormal_200_55_phia$adjusted_vls_estimate - lognormal_200_55_phia$empirical_estimate)
mean(frec_200_55$adjusted_vls_estimate - frec_200_55$empirical_estimate)
mean(gamma_200_55$adjusted_vls_estimate - gamma_200_55$empirical_estimate)
mean(lognormal_200_55$adjusted_vls_estimate - lognormal_200_55$empirical_estimate)

mean(frec_400_55_phia$adjusted_vls_estimate - frec_400_55_phia$empirical_estimate)
mean(gamma_400_55_phia$adjusted_vls_estimate - gamma_400_55_phia$empirical_estimate)
mean(lognormal_400_55_phia$adjusted_vls_estimate - lognormal_400_55_phia$empirical_estimate)
mean(frec_400_55$adjusted_vls_estimate - frec_400_55$empirical_estimate)
mean(gamma_400_55$adjusted_vls_estimate - gamma_400_55$empirical_estimate)
mean(lognormal_400_55$adjusted_vls_estimate - lognormal_400_55$empirical_estimate)

# mean absolute error
#' 15-24
mean(abs(frec_50_15_24_phia$adjusted_vls_estimate - frec_50_15_24_phia$empirical_estimate))
mean(abs(gamma_50_15_24_phia$adjusted_vls_estimate - gamma_50_15_24_phia$empirical_estimate))
mean(abs(lognormal_50_15_24_phia$adjusted_vls_estimate - lognormal_50_15_24_phia$empirical_estimate))
mean(abs(frec_50_15_24$adjusted_vls_estimate - frec_50_15_24$empirical_estimate))
mean(abs(gamma_50_15_24$adjusted_vls_estimate - gamma_50_15_24$empirical_estimate))
mean(abs(lognormal_50_15_24$adjusted_vls_estimate - lognormal_50_15_24$empirical_estimate))

mean(abs(frec_200_15_24_phia$adjusted_vls_estimate - frec_200_15_24_phia$empirical_estimate))
mean(abs(gamma_200_15_24_phia$adjusted_vls_estimate - gamma_200_15_24_phia$empirical_estimate))
mean(abs(lognormal_200_15_24_phia$adjusted_vls_estimate - lognormal_200_15_24_phia$empirical_estimate))
mean(abs(frec_200_15_24$adjusted_vls_estimate - frec_200_15_24$empirical_estimate))
mean(abs(gamma_200_15_24$adjusted_vls_estimate - gamma_200_15_24$empirical_estimate))
mean(abs(lognormal_200_15_24$adjusted_vls_estimate - lognormal_200_15_24$empirical_estimate))

mean(abs(frec_400_15_24_phia$adjusted_vls_estimate - frec_400_15_24_phia$empirical_estimate))
mean(abs(gamma_400_15_24_phia$adjusted_vls_estimate - gamma_400_15_24_phia$empirical_estimate))
mean(abs(lognormal_400_15_24_phia$adjusted_vls_estimate - lognormal_400_15_24_phia$empirical_estimate))
mean(abs(frec_400_15_24$adjusted_vls_estimate - frec_400_15_24$empirical_estimate))
mean(abs(gamma_400_15_24$adjusted_vls_estimate - gamma_400_15_24$empirical_estimate))
mean(abs(lognormal_400_15_24$adjusted_vls_estimate - lognormal_400_15_24$empirical_estimate))

# 25-34
mean(abs(frec_50_25_34_phia$adjusted_vls_estimate - frec_50_25_34_phia$empirical_estimate))
mean(abs(gamma_50_25_34_phia$adjusted_vls_estimate - gamma_50_25_34_phia$empirical_estimate))
mean(abs(lognormal_50_25_34_phia$adjusted_vls_estimate - lognormal_50_25_34_phia$empirical_estimate))
mean(abs(frec_50_25_34$adjusted_vls_estimate - frec_50_25_34$empirical_estimate))
mean(abs(gamma_50_25_34$adjusted_vls_estimate - gamma_50_25_34$empirical_estimate))
mean(abs(lognormal_50_25_34$adjusted_vls_estimate - lognormal_50_25_34$empirical_estimate))

mean(abs(frec_200_25_34_phia$adjusted_vls_estimate - frec_200_25_34_phia$empirical_estimate))
mean(abs(gamma_200_25_34_phia$adjusted_vls_estimate - gamma_200_25_34_phia$empirical_estimate))
mean(abs(lognormal_200_25_34_phia$adjusted_vls_estimate - lognormal_200_25_34_phia$empirical_estimate))
mean(abs(frec_200_25_34$adjusted_vls_estimate - frec_200_25_34$empirical_estimate))
mean(abs(gamma_200_25_34$adjusted_vls_estimate - gamma_200_25_34$empirical_estimate))
mean(abs(lognormal_200_25_34$adjusted_vls_estimate - lognormal_200_25_34$empirical_estimate))

mean(abs(frec_400_25_34_phia$adjusted_vls_estimate - frec_400_25_34_phia$empirical_estimate))
mean(abs(gamma_400_25_34_phia$adjusted_vls_estimate - gamma_400_25_34_phia$empirical_estimate))
mean(abs(lognormal_400_25_34_phia$adjusted_vls_estimate - lognormal_400_25_34_phia$empirical_estimate))
mean(abs(frec_400_25_34$adjusted_vls_estimate - frec_400_25_34$empirical_estimate))
mean(abs(gamma_400_25_34$adjusted_vls_estimate - gamma_400_25_34$empirical_estimate))
mean(abs(lognormal_400_25_34$adjusted_vls_estimate - lognormal_400_25_34$empirical_estimate))

# 35-44
mean(abs(frec_50_35_44_phia$adjusted_vls_estimate - frec_50_35_44_phia$empirical_estimate))
mean(abs(gamma_50_35_44_phia$adjusted_vls_estimate - gamma_50_35_44_phia$empirical_estimate))
mean(abs(lognormal_50_35_44_phia$adjusted_vls_estimate - lognormal_50_35_44_phia$empirical_estimate))
mean(abs(frec_50_35_44$adjusted_vls_estimate - frec_50_35_44$empirical_estimate))
mean(abs(gamma_50_35_44$adjusted_vls_estimate - gamma_50_35_44$empirical_estimate))
mean(abs(lognormal_50_35_44$adjusted_vls_estimate - lognormal_50_35_44$empirical_estimate))

mean(abs(frec_200_35_44_phia$adjusted_vls_estimate - frec_200_35_44_phia$empirical_estimate))
mean(abs(gamma_200_35_44_phia$adjusted_vls_estimate - gamma_200_35_44_phia$empirical_estimate))
mean(abs(lognormal_200_35_44_phia$adjusted_vls_estimate - lognormal_200_35_44_phia$empirical_estimate))
mean(abs(frec_200_35_44$adjusted_vls_estimate - frec_200_35_44$empirical_estimate))
mean(abs(gamma_200_35_44$adjusted_vls_estimate - gamma_200_35_44$empirical_estimate))
mean(abs(lognormal_200_35_44$adjusted_vls_estimate - lognormal_200_35_44$empirical_estimate))

mean(abs(frec_400_35_44_phia$adjusted_vls_estimate - frec_400_35_44_phia$empirical_estimate))
mean(abs(gamma_400_35_44_phia$adjusted_vls_estimate - gamma_400_35_44_phia$empirical_estimate))
mean(abs(lognormal_400_35_44_phia$adjusted_vls_estimate - lognormal_400_35_44_phia$empirical_estimate))
mean(abs(frec_400_35_44$adjusted_vls_estimate - frec_400_35_44$empirical_estimate))
mean(abs(gamma_400_35_44$adjusted_vls_estimate - gamma_400_35_44$empirical_estimate))
mean(abs(lognormal_400_35_44$adjusted_vls_estimate - lognormal_400_35_44$empirical_estimate))

# 45-54
mean(abs(frec_50_45_54_phia$adjusted_vls_estimate - frec_50_45_54_phia$empirical_estimate))
mean(abs(gamma_50_45_54_phia$adjusted_vls_estimate - gamma_50_45_54_phia$empirical_estimate))
mean(abs(lognormal_50_45_54_phia$adjusted_vls_estimate - lognormal_50_45_54_phia$empirical_estimate))
mean(abs(frec_50_45_54$adjusted_vls_estimate - frec_50_45_54$empirical_estimate))
mean(abs(gamma_50_45_54$adjusted_vls_estimate - gamma_50_45_54$empirical_estimate))
mean(abs(lognormal_50_45_54$adjusted_vls_estimate - lognormal_50_45_54$empirical_estimate))

mean(abs(frec_200_45_54_phia$adjusted_vls_estimate - frec_200_45_54_phia$empirical_estimate))
mean(abs(gamma_200_45_54_phia$adjusted_vls_estimate - gamma_200_45_54_phia$empirical_estimate))
mean(abs(lognormal_200_45_54_phia$adjusted_vls_estimate - lognormal_200_45_54_phia$empirical_estimate))
mean(abs(frec_200_45_54$adjusted_vls_estimate - frec_200_45_54$empirical_estimate))
mean(abs(gamma_200_45_54$adjusted_vls_estimate - gamma_200_45_54$empirical_estimate))
mean(abs(lognormal_200_45_54$adjusted_vls_estimate - lognormal_200_45_54$empirical_estimate))

mean(abs(frec_400_45_54_phia$adjusted_vls_estimate - frec_400_45_54_phia$empirical_estimate))
mean(abs(gamma_400_45_54_phia$adjusted_vls_estimate - gamma_400_45_54_phia$empirical_estimate))
mean(abs(lognormal_400_45_54_phia$adjusted_vls_estimate - lognormal_400_45_54_phia$empirical_estimate))
mean(abs(frec_400_45_54$adjusted_vls_estimate - frec_400_45_54$empirical_estimate))
mean(abs(gamma_400_45_54$adjusted_vls_estimate - gamma_400_45_54$empirical_estimate))
mean(abs(lognormal_400_45_54$adjusted_vls_estimate - lognormal_400_45_54$empirical_estimate))

# 55+
mean(abs(frec_50_55_phia$adjusted_vls_estimate - frec_50_55_phia$empirical_estimate))
mean(abs(gamma_50_55_phia$adjusted_vls_estimate - gamma_50_55_phia$empirical_estimate))
mean(abs(lognormal_50_55_phia$adjusted_vls_estimate - lognormal_50_55_phia$empirical_estimate))
mean(abs(frec_50_55$adjusted_vls_estimate - frec_50_55$empirical_estimate))
mean(abs(gamma_50_55$adjusted_vls_estimate - gamma_50_55$empirical_estimate))
mean(abs(lognormal_50_55$adjusted_vls_estimate - lognormal_50_55$empirical_estimate))

mean(abs(frec_200_55_phia$adjusted_vls_estimate - frec_200_55_phia$empirical_estimate))
mean(abs(gamma_200_55_phia$adjusted_vls_estimate - gamma_200_55_phia$empirical_estimate))
mean(abs(lognormal_200_55_phia$adjusted_vls_estimate - lognormal_200_55_phia$empirical_estimate))
mean(abs(frec_200_55$adjusted_vls_estimate - frec_200_55$empirical_estimate))
mean(abs(gamma_200_55$adjusted_vls_estimate - gamma_200_55$empirical_estimate))
mean(abs(lognormal_200_55$adjusted_vls_estimate - lognormal_200_55$empirical_estimate))

mean(abs(frec_400_55_phia$adjusted_vls_estimate - frec_400_55_phia$empirical_estimate))
mean(abs(gamma_400_55_phia$adjusted_vls_estimate - gamma_400_55_phia$empirical_estimate))
mean(abs(lognormal_400_55_phia$adjusted_vls_estimate - lognormal_400_55_phia$empirical_estimate))
mean(abs(frec_400_55$adjusted_vls_estimate - frec_400_55$empirical_estimate))
mean(abs(gamma_400_55$adjusted_vls_estimate - gamma_400_55$empirical_estimate))
mean(abs(lognormal_400_55$adjusted_vls_estimate - lognormal_400_55$empirical_estimate))

# plot
age_validation <- bind_rows(frec_50_15_24,frec_50_25_34,frec_50_35_44,frec_50_45_54,frec_50_55,
                            frec_50_15_24_phia,frec_50_25_34_phia,frec_50_35_44_phia,frec_50_45_54_phia,frec_50_55_phia, 
                            frec_200_15_24,frec_200_25_34,frec_200_35_44,frec_200_45_54,frec_200_55,
                            frec_200_15_24_phia,frec_200_25_34_phia,frec_200_35_44_phia,frec_200_45_54_phia,frec_200_55_phia, 
                            frec_400_15_24,frec_400_25_34,frec_400_35_44,frec_400_45_54,frec_400_55,
                            frec_400_15_24_phia,frec_400_25_34_phia,frec_400_35_44_phia,frec_400_45_54_phia,frec_400_55_phia,
                            gamma_50_15_24,gamma_50_25_34,gamma_50_35_44,gamma_50_45_54,gamma_50_55,
                            gamma_50_15_24_phia,gamma_50_25_34_phia,gamma_50_35_44_phia,gamma_50_45_54_phia,gamma_50_55_phia,
                            gamma_200_15_24,gamma_200_25_34,gamma_200_35_44,gamma_200_45_54,gamma_200_55,
                            gamma_200_15_24_phia,gamma_200_25_34_phia,gamma_200_35_44_phia,gamma_200_45_54_phia,gamma_200_55_phia,
                            gamma_400_15_24,gamma_400_25_34,gamma_400_35_44,gamma_400_45_54,gamma_400_55,
                            gamma_400_15_24_phia,gamma_400_25_34_phia,gamma_400_35_44_phia,gamma_400_45_54_phia,gamma_400_55_phia,
                            lognormal_50_15_24,lognormal_50_25_34,lognormal_50_35_44,lognormal_50_45_54,lognormal_50_55,
                            lognormal_50_15_24_phia,lognormal_50_25_34_phia,lognormal_50_35_44_phia,lognormal_50_45_54_phia,lognormal_50_55_phia,
                            lognormal_200_15_24,lognormal_200_25_34,lognormal_200_35_44,lognormal_200_45_54,lognormal_200_55,
                            lognormal_200_15_24_phia,lognormal_200_25_34_phia,lognormal_200_35_44_phia,lognormal_200_45_54_phia,lognormal_200_55_phia,
                            lognormal_400_15_24,lognormal_400_25_34,lognormal_400_35_44,lognormal_400_45_54,lognormal_400_55,
                            lognormal_400_15_24_phia,lognormal_400_25_34_phia,lognormal_400_35_44_phia,lognormal_400_45_54_phia,lognormal_400_55_phia)


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
                                             c("Frechet", "Gamma", "lognormal")),
         country_names = forcats::fct_relevel(country_names,
                                              c("Côte d'Ivoire (2017-18)","Cameroon (2017-18)",
                                                "Nigeria (2018)","Uganda (2016-17)","Zimbabwe (2015-16)",
                                                "Tanzania (2016-17)","Ethiopia (2017-18)","Lesotho (2016-17)",
                                                "Zambia (2016)","Mozambique (2021-22)","Rwanda (2018-19)",
                                                "Zimbabwe (2019-20)","Kenya (2018-19)","Malawi (2015-16)",
                                                "Namibia (2017)","Eswatini (2016-17)","Lesotho (2019-20)",
                                                "Eswatini (2021)","Zambia (2021)","Malawi (2020-21)", "Botswana (2021)")))



## scatter plot
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
        plot.title = element_text(size = rel(1.0), family = "sans", face = "bold"))

