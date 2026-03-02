#' create function using the Gamma, Frechet and lognormal distributions to convert viral suppression estimates
#' Frechet model
#' @details probability denisity functions to obtain viral suppression conversion thresholds
#' obtained from Leigh et al paper
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
    
    print("Gamma distribution")
    
    return(res)
    
  } else if(distribution == "lognormal"){
    
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
    
    print("lognormal distribution")
    
    return(res)
    
    
  } else if(distribution == "Frechet Men"){
    
    adjusted_vls_estimate <- frechet_func_vec(shape = 1.94,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- frechet_func_vec(shape = 1.88,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- frechet_func_vec(shape = 1.99,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Frechet Men")
    
    return(res)
    
    
  }else if(distribution == "Frechet Women"){
    
    adjusted_vls_estimate <- frechet_func_vec(shape = 1.84,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- frechet_func_vec(shape = 1.77,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- frechet_func_vec(shape = 1.90,    
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Frechet Women")
    
    return(res)
    
  }else if(distribution == "Gamma Men"){
    
    adjusted_vls_estimate <- gamma_func_vec(shape = 1.02, 
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- gamma_func_vec(shape = 0.97,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- gamma_func_vec(shape = 1.07,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Gamma Men")
    
    return(res)
    
  }else if(distribution == "Gamma Women"){
    
    adjusted_vls_estimate <- gamma_func_vec(shape = 0.73, 
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- gamma_func_vec(shape = 0.69,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- gamma_func_vec(shape = 0.79,  
                                            reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Gamma Women")
    
    return(res)
    
  }else if(distribution == "lognormal Men"){
    
    adjusted_vls_estimate <- lognorm_func_vec(sdlog = 0.83, 
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- lognorm_func_vec(sdlog = 0.80,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- lognorm_func_vec(sdlog = 0.85,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("lognormal Men")
    
    return(res)
    
  }else{
    
    adjusted_vls_estimate <- lognorm_func_vec(sdlog = 0.92, 
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_lower_ci <- lognorm_func_vec(sdlog = 0.90,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    adjusted_vls_upper_ci <- lognorm_func_vec(sdlog = 0.95,  
                                              reported_threshold = estimate_threshold, reported_est = viral_sup_estimate)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("lognormal Women")
    
    return(res)
    
  }
  
}


#' calculate adjusted viral suppression estimates
#' Frechet distribution < 400 copies/ml
#' males
frec_400_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_400,
                                             estimate_threshold = 400, conversion_threshold = 1000, 
                                             distribution = "Frechet Men", 
                                             survey = PHIA_vls_estimates_male$survey,
                                             empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                             empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000, 
                                             empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

frec_400_male_phia$distribution <- rep("Frechet Men", dim(frec_400_male_phia)[1])
frec_400_male_phia$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_male_phia)[1])
frec_400_male_phia$sex <- rep("Men", dim(frec_400_male_phia)[1])

# parametrs pooled
frec_400_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_400,
                                        estimate_threshold = 400, conversion_threshold = 1000, 
                                        distribution = "Frechet", 
                                        survey = PHIA_vls_estimates_male$survey,
                                        empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000, 
                                        empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

frec_400_male$distribution <- rep("Frechet", dim(frec_400_male)[1])
frec_400_male$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_male)[1])
frec_400_male$sex <- rep("Men", dim(frec_400_male)[1])

#' Frechet distribution < 200 copies/ml
frec_200_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_200,
                                             estimate_threshold = 200, conversion_threshold = 1000, 
                                             distribution = "Frechet Men",
                                             survey = PHIA_vls_estimates_male$survey,
                                             empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                             empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                             empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

frec_200_male_phia$distribution <- rep("Frechet Men", dim(frec_200_male_phia)[1])
frec_200_male_phia$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_male_phia)[1])
frec_200_male_phia$sex <- rep("Men", dim(frec_200_male_phia)[1])

# parameters from pooled
#' Frechet distribution < 200 copies/ml
frec_200_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_200,
                                        estimate_threshold = 200, conversion_threshold = 1000, 
                                        distribution = "Frechet",
                                        survey = PHIA_vls_estimates_male$survey,
                                        empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

frec_200_male$distribution <- rep("Frechet", dim(frec_200_male)[1])
frec_200_male$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_male)[1])
frec_200_male$sex <- rep("Men", dim(frec_200_male)[1])

#' Frechet distribution < 50 copies/ml
frec_50_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_50,
                                            estimate_threshold = 50, conversion_threshold = 1000, 
                                            distribution = "Frechet Men",
                                            survey = PHIA_vls_estimates_male$survey,
                                            empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                            empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                            empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

frec_50_male_phia$distribution <- rep("Frechet Men", dim(frec_50_male_phia)[1])
frec_50_male_phia$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_male_phia)[1])
frec_50_male_phia$sex <- rep("Men", dim(frec_50_male_phia)[1])

# using parameters from pooled et al.
frec_50_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_50,
                                       estimate_threshold = 50, conversion_threshold = 1000, 
                                       distribution = "Frechet",
                                       survey = PHIA_vls_estimates_male$survey,
                                       empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

frec_50_male$distribution <- rep("Frechet", dim(frec_50_male)[1])
frec_50_male$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_male)[1])
frec_50_male$sex <- rep("Men", dim(frec_50_male)[1])


#' Gamma distribution < 400 copies/ml
gamma_400_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_400,
                                         estimate_threshold = 400, conversion_threshold = 1000, 
                                         distribution = "Gamma Men", 
                                         survey = PHIA_vls_estimates_male$survey,
                                         empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000, 
                                         empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

gamma_400_male_phia$distribution <- rep("Gamma Men", dim(gamma_400_male_phia)[1])
gamma_400_male_phia$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_male_phia)[1])
gamma_400_male_phia$sex <- rep("Men", dim(gamma_400_male_phia)[1])

# pooled
gamma_400_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_400,
                                    estimate_threshold = 400, conversion_threshold = 1000, 
                                    distribution = "Gamma", 
                                    survey = PHIA_vls_estimates_male$survey,
                                    empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000, 
                                    empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

gamma_400_male$distribution <- rep("Gamma", dim(gamma_400_male)[1])
gamma_400_male$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_male)[1])
gamma_400_male$sex <- rep("Men", dim(gamma_400_male)[1])


#' Gamma distribution < 200 copies/ml
gamma_200_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_200,
                                         estimate_threshold = 200, conversion_threshold = 1000, 
                                         distribution = "Gamma Men",
                                         survey = PHIA_vls_estimates_male$survey,
                                         empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

gamma_200_male_phia$distribution <- rep("Gamma Men", dim(gamma_200_male_phia)[1])
gamma_200_male_phia$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_male_phia)[1])
gamma_200_male_phia$sex <- rep("Men", dim(gamma_200_male_phia)[1])

# pooled
gamma_200_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_200,
                                    estimate_threshold = 200, conversion_threshold = 1000, 
                                    distribution = "Gamma",
                                    survey = PHIA_vls_estimates_male$survey,
                                    empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                    empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                    empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

gamma_200_male$distribution <- rep("Gamma", dim(gamma_200_male)[1])
gamma_200_male$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_male)[1])
gamma_200_male$sex <- rep("Men", dim(gamma_200_male)[1])

#' Gamma distribution < 50 copies/ml
gamma_50_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_50,
                                        estimate_threshold = 50, conversion_threshold = 1000, 
                                        distribution = "Gamma Men",
                                        survey = PHIA_vls_estimates_male$survey,
                                        empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

gamma_50_male_phia$distribution <- rep("Gamma Men", dim(gamma_50_male_phia)[1])
gamma_50_male_phia$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_male_phia)[1])
gamma_50_male_phia$sex <- rep("Men", dim(gamma_50_male_phia)[1])

# pooled
gamma_50_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_50,
                                   estimate_threshold = 50, conversion_threshold = 1000, 
                                   distribution = "Gamma",
                                   survey = PHIA_vls_estimates_male$survey,
                                   empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

gamma_50_male$distribution <- rep("Gamma", dim(gamma_50_male)[1])
gamma_50_male$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_male)[1])
gamma_50_male$sex <- rep("Men", dim(gamma_50_male)[1])


#' lognormal distribution < 400 copies/ml
lognormal_400_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_400,
                                           estimate_threshold = 400, conversion_threshold = 1000, 
                                           distribution = "lognormal Men", 
                                           survey = PHIA_vls_estimates_male$survey,
                                           empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                           empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000, 
                                           empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

lognormal_400_male_phia$distribution <- rep("lognormal Men", dim(lognormal_400_male_phia)[1])
lognormal_400_male_phia$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_male_phia)[1])
lognormal_400_male_phia$sex <- rep("Men", dim(lognormal_400_male_phia)[1])

# pooled
lognormal_400_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_400,
                                      estimate_threshold = 400, conversion_threshold = 1000, 
                                      distribution = "lognormal", 
                                      survey = PHIA_vls_estimates_male$survey,
                                      empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000, 
                                      empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

lognormal_400_male$distribution <- rep("lognormal", dim(lognormal_400_male)[1])
lognormal_400_male$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_male)[1])
lognormal_400_male$sex <- rep("Men", dim(lognormal_400_male)[1])

#' lognormal distribution < 200 copies/ml
lognormal_200_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_200,
                                           estimate_threshold = 200, conversion_threshold = 1000, 
                                           distribution = "lognormal Men",
                                           survey = PHIA_vls_estimates_male$survey,
                                           empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                           empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                           empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

lognormal_200_male_phia$distribution <- rep("lognormal Men", dim(lognormal_200_male_phia)[1])
lognormal_200_male_phia$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_male_phia)[1])
lognormal_200_male_phia$sex <- rep("Men", dim(lognormal_200_male_phia)[1])

# pooled
lognormal_200_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_200,
                                      estimate_threshold = 200, conversion_threshold = 1000, 
                                      distribution = "lognormal",
                                      survey = PHIA_vls_estimates_male$survey,
                                      empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

lognormal_200_male$distribution <- rep("lognormal", dim(lognormal_200_male)[1])
lognormal_200_male$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_male)[1])
lognormal_200_male$sex <- rep("Men", dim(lognormal_200_male)[1])

#' lognormal distribution < 50 copies/ml
lognormal_50_male_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_50,
                                          estimate_threshold = 50, conversion_threshold = 1000, 
                                          distribution = "lognormal Men",
                                          survey = PHIA_vls_estimates_male$survey,
                                          empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                          empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

lognormal_50_male_phia$distribution <- rep("lognormal Men", dim(lognormal_50_male_phia)[1])
lognormal_50_male_phia$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_male_phia)[1])
lognormal_50_male_phia$sex <- rep("Men", dim(lognormal_50_male_phia)[1])

# pooled
lognormal_50_male <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_male$estimate_50,
                                     estimate_threshold = 50, conversion_threshold = 1000, 
                                     distribution = "lognormal",
                                     survey = PHIA_vls_estimates_male$survey,
                                     empirical_estimate = PHIA_vls_estimates_male$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_male$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_male$upper_ci_1000)

lognormal_50_male$distribution <- rep("lognormal", dim(lognormal_50_male)[1])
lognormal_50_male$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_male)[1])
lognormal_50_male$sex <- rep("Men", dim(lognormal_50_male)[1])


#' calculate the root mean squared error between empirical and adjusted estimates
sqrt(mean((frec_50_male_phia$empirical_estimate - frec_50_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_male_phia$empirical_estimate - gamma_50_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_male_phia$empirical_estimate - lognormal_50_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_50_male$empirical_estimate - frec_50_male$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_male$empirical_estimate - gamma_50_male$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_male$empirical_estimate - lognormal_50_male$adjusted_vls_estimate)^2))

sqrt(mean((frec_200_male_phia$empirical_estimate - frec_200_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_male_phia$empirical_estimate - gamma_200_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_male_phia$empirical_estimate - lognormal_200_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_200_male$empirical_estimate - frec_200_male$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_male$empirical_estimate - gamma_200_male$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_male$empirical_estimate - lognormal_200_male$adjusted_vls_estimate)^2))

sqrt(mean((frec_400_male_phia$empirical_estimate - frec_400_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_male_phia$empirical_estimate - gamma_400_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_male_phia$empirical_estimate - lognormal_400_male_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_400_male$empirical_estimate - frec_400_male$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_male$empirical_estimate - gamma_400_male$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_male$empirical_estimate - lognormal_400_male$adjusted_vls_estimate)^2))

# avearge bias
mean(frec_50_male_phia$adjusted_vls_estimate - frec_50_male_phia$empirical_estimate)
mean(gamma_50_male_phia$adjusted_vls_estimate - gamma_50_male_phia$empirical_estimate)
mean(lognormal_50_male_phia$adjusted_vls_estimate - lognormal_50_male_phia$empirical_estimate)
mean(frec_50_male$adjusted_vls_estimate - frec_50_male$empirical_estimate)
mean(gamma_50_male$adjusted_vls_estimate - gamma_50_male$empirical_estimate)
mean(lognormal_50_male$adjusted_vls_estimate - lognormal_50_male$empirical_estimate)

mean(frec_200_male_phia$adjusted_vls_estimate - frec_200_male_phia$empirical_estimate)
mean(gamma_200_male_phia$adjusted_vls_estimate - gamma_200_male_phia$empirical_estimate)
mean(lognormal_200_male_phia$adjusted_vls_estimate - lognormal_200_male_phia$empirical_estimate)
mean(frec_200_male$adjusted_vls_estimate - frec_200_male$empirical_estimate)
mean(gamma_200_male$adjusted_vls_estimate - gamma_200_male$empirical_estimate)
mean(lognormal_200_male$adjusted_vls_estimate - lognormal_200_male$empirical_estimate)

mean(frec_400_male_phia$adjusted_vls_estimate - frec_400_male_phia$empirical_estimate)
mean(gamma_400_male_phia$adjusted_vls_estimate - gamma_400_male_phia$empirical_estimate)
mean(lognormal_400_male_phia$adjusted_vls_estimate - lognormal_400_male_phia$empirical_estimate)
mean(frec_400_male$adjusted_vls_estimate - frec_400_male$empirical_estimate)
mean(gamma_400_male$adjusted_vls_estimate - gamma_400_male$empirical_estimate)
mean(lognormal_400_male$adjusted_vls_estimate - lognormal_400_male$empirical_estimate)

# mean absolute error
mean(abs(frec_50_male_phia$adjusted_vls_estimate - frec_50_male_phia$empirical_estimate))
mean(abs(gamma_50_male_phia$adjusted_vls_estimate - gamma_50_male_phia$empirical_estimate))
mean(abs(lognormal_50_male_phia$adjusted_vls_estimate - lognormal_50_male_phia$empirical_estimate))
mean(abs(frec_50_male$adjusted_vls_estimate - frec_50_male$empirical_estimate))
mean(abs(gamma_50_male$adjusted_vls_estimate - gamma_50_male$empirical_estimate))
mean(abs(lognormal_50_male$adjusted_vls_estimate - lognormal_50_male$empirical_estimate))

mean(abs(frec_200_male_phia$adjusted_vls_estimate - frec_200_male_phia$empirical_estimate))
mean(abs(gamma_200_male_phia$adjusted_vls_estimate - gamma_200_male_phia$empirical_estimate))
mean(abs(lognormal_200_male_phia$adjusted_vls_estimate - lognormal_200_male_phia$empirical_estimate))
mean(abs(frec_200_male$adjusted_vls_estimate - frec_200_male$empirical_estimate))
mean(abs(gamma_200_male$adjusted_vls_estimate - gamma_200_male$empirical_estimate))
mean(abs(lognormal_200_male$adjusted_vls_estimate - lognormal_200_male$empirical_estimate))

mean(abs(frec_400_male_phia$adjusted_vls_estimate - frec_400_male_phia$empirical_estimate))
mean(abs(gamma_400_male_phia$adjusted_vls_estimate - gamma_400_male_phia$empirical_estimate))
mean(abs(lognormal_400_male_phia$adjusted_vls_estimate - lognormal_400_male_phia$empirical_estimate))
mean(abs(frec_400_male$adjusted_vls_estimate - frec_400_male$empirical_estimate))
mean(abs(gamma_400_male$adjusted_vls_estimate - gamma_400_male$empirical_estimate))
mean(abs(lognormal_400_male$adjusted_vls_estimate - lognormal_400_male$empirical_estimate))


#' calculate adjusted viral suppression estimates
#' Frechet distribution < 400 copies/ml
#' females
frec_400_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_400,
                                               estimate_threshold = 400, conversion_threshold = 1000, 
                                               distribution = "Frechet Women", 
                                               survey = PHIA_vls_estimates_female$survey,
                                               empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                               empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000, 
                                               empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

frec_400_female_phia$distribution <- rep("Frechet Women", dim(frec_400_female_phia)[1])
frec_400_female_phia$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_female_phia)[1])
frec_400_female_phia$sex <- rep("Women", dim(frec_400_female_phia)[1])

# pooled
frec_400_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_400,
                                          estimate_threshold = 400, conversion_threshold = 1000, 
                                          distribution = "Frechet", 
                                          survey = PHIA_vls_estimates_female$survey,
                                          empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000, 
                                          empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

frec_400_female$distribution <- rep("Frechet", dim(frec_400_female)[1])
frec_400_female$estimate_threshold <- rep("<400 copies/mL", dim(frec_400_female)[1])
frec_400_female$sex <- rep("Women", dim(frec_400_female)[1])

#' Frechet distribution < 200 copies/ml
frec_200_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_200,
                                               estimate_threshold = 200, conversion_threshold = 1000, 
                                               distribution = "Frechet Women",
                                               survey = PHIA_vls_estimates_female$survey,
                                               empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                               empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                               empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

frec_200_female_phia$distribution <- rep("Frechet Women", dim(frec_200_female_phia)[1])
frec_200_female_phia$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_female_phia)[1])
frec_200_female_phia$sex <- rep("Women", dim(frec_200_female_phia)[1])

# pooled
frec_200_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_200,
                                          estimate_threshold = 200, conversion_threshold = 1000, 
                                          distribution = "Frechet",
                                          survey = PHIA_vls_estimates_female$survey,
                                          empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                          empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

frec_200_female$distribution <- rep("Frechet", dim(frec_200_female)[1])
frec_200_female$estimate_threshold <- rep("<200 copies/mL", dim(frec_200_female)[1])
frec_200_female$sex <- rep("Women", dim(frec_200_female)[1])

#' Frechet distribution < 50 copies/ml
frec_50_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_50,
                                              estimate_threshold = 50, conversion_threshold = 1000, 
                                              distribution = "Frechet Women",
                                              survey = PHIA_vls_estimates_female$survey,
                                              empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                              empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                              empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

frec_50_female_phia$distribution <- rep("Frechet Women", dim(frec_50_female_phia)[1])
frec_50_female_phia$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_female_phia)[1])
frec_50_female_phia$sex <- rep("Women", dim(frec_50_female_phia)[1])

# pooled
frec_50_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_50,
                                         estimate_threshold = 50, conversion_threshold = 1000, 
                                         distribution = "Frechet",
                                         survey = PHIA_vls_estimates_female$survey,
                                         empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                         empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                         empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

frec_50_female$distribution <- rep("Frechet", dim(frec_50_female)[1])
frec_50_female$estimate_threshold <- rep("<50 copies/mL", dim(frec_50_female)[1])
frec_50_female$sex <- rep("Women", dim(frec_50_female)[1])

#' Gamma distribution < 400 copies/ml
gamma_400_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_400,
                                           estimate_threshold = 400, conversion_threshold = 1000, 
                                           distribution = "Gamma Women", 
                                           survey = PHIA_vls_estimates_female$survey,
                                           empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                           empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000, 
                                           empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

gamma_400_female_phia$distribution <- rep("Gamma Women", dim(gamma_400_female_phia)[1])
gamma_400_female_phia$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_female_phia)[1])
gamma_400_female_phia$sex <- rep("Women", dim(gamma_400_female_phia)[1])

# pooled
gamma_400_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_400,
                                      estimate_threshold = 400, conversion_threshold = 1000, 
                                      distribution = "Gamma", 
                                      survey = PHIA_vls_estimates_female$survey,
                                      empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000, 
                                      empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

gamma_400_female$distribution <- rep("Gamma", dim(gamma_400_female)[1])
gamma_400_female$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400_female)[1])
gamma_400_female$sex <- rep("Women", dim(gamma_400_female)[1])

#' Gamma distribution < 200 copies/ml
gamma_200_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_200,
                                           estimate_threshold = 200, conversion_threshold = 1000, 
                                           distribution = "Gamma Women",
                                           survey = PHIA_vls_estimates_female$survey,
                                           empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                           empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                           empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

gamma_200_female_phia$distribution <- rep("Gamma Women", dim(gamma_200_female_phia)[1])
gamma_200_female_phia$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_female_phia)[1])
gamma_200_female_phia$sex <- rep("Women", dim(gamma_200_female_phia)[1])

# pooled
gamma_200_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_200,
                                      estimate_threshold = 200, conversion_threshold = 1000, 
                                      distribution = "Gamma",
                                      survey = PHIA_vls_estimates_female$survey,
                                      empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

gamma_200_female$distribution <- rep("Gamma", dim(gamma_200_female)[1])
gamma_200_female$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200_female)[1])
gamma_200_female$sex <- rep("Women", dim(gamma_200_female)[1])

#' Gamma distribution < 50 copies/ml
gamma_50_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_50,
                                          estimate_threshold = 50, conversion_threshold = 1000, 
                                          distribution = "Gamma Women",
                                          survey = PHIA_vls_estimates_female$survey,
                                          empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                          empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                          empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

gamma_50_female_phia$distribution <- rep("Gamma Women", dim(gamma_50_female_phia)[1])
gamma_50_female_phia$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_female_phia)[1])
gamma_50_female_phia$sex <- rep("Women", dim(gamma_50_female_phia)[1])

# pooled
gamma_50_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_50,
                                     estimate_threshold = 50, conversion_threshold = 1000, 
                                     distribution = "Gamma",
                                     survey = PHIA_vls_estimates_female$survey,
                                     empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                     empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                     empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

gamma_50_female$distribution <- rep("Gamma", dim(gamma_50_female)[1])
gamma_50_female$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50_female)[1])
gamma_50_female$sex <- rep("Women", dim(gamma_50_female)[1])


#' lognormal distribution < 400 copies/ml
lognormal_400_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_400,
                                             estimate_threshold = 400, conversion_threshold = 1000, 
                                             distribution = "lognormal Women", 
                                             survey = PHIA_vls_estimates_female$survey,
                                             empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                             empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000, 
                                             empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

lognormal_400_female_phia$distribution <- rep("lognormal Women", dim(lognormal_400_female_phia)[1])
lognormal_400_female_phia$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_female_phia)[1])
lognormal_400_female_phia$sex <- rep("Women", dim(lognormal_400_female_phia)[1])

# pooled
lognormal_400_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_400,
                                        estimate_threshold = 400, conversion_threshold = 1000, 
                                        distribution = "lognormal", 
                                        survey = PHIA_vls_estimates_female$survey,
                                        empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000, 
                                        empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

lognormal_400_female$distribution <- rep("lognormal", dim(lognormal_400_female)[1])
lognormal_400_female$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400_female)[1])
lognormal_400_female$sex <- rep("Women", dim(lognormal_400_female)[1])


#' lognormal distribution < 200 copies/ml
lognormal_200_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_200,
                                             estimate_threshold = 200, conversion_threshold = 1000, 
                                             distribution = "lognormal Women",
                                             survey = PHIA_vls_estimates_female$survey,
                                             empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                             empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                             empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

lognormal_200_female_phia$distribution <- rep("lognormal Women", dim(lognormal_200_female_phia)[1])
lognormal_200_female_phia$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_female_phia)[1])
lognormal_200_female_phia$sex <- rep("Women", dim(lognormal_200_female_phia)[1])

# pooled
lognormal_200_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_200,
                                        estimate_threshold = 200, conversion_threshold = 1000, 
                                        distribution = "lognormal",
                                        survey = PHIA_vls_estimates_female$survey,
                                        empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                        empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                        empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

lognormal_200_female$distribution <- rep("lognormal", dim(lognormal_200_female)[1])
lognormal_200_female$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200_female)[1])
lognormal_200_female$sex <- rep("Women", dim(lognormal_200_female)[1])

#' lognormal distribution < 50 copies/ml
lognormal_50_female_phia <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_50,
                                            estimate_threshold = 50, conversion_threshold = 1000, 
                                            distribution = "lognormal Women",
                                            survey = PHIA_vls_estimates_female$survey,
                                            empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                            empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                            empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

lognormal_50_female_phia$distribution <- rep("lognormal Women", dim(lognormal_50_female_phia)[1])
lognormal_50_female_phia$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_female_phia)[1])
lognormal_50_female_phia$sex <- rep("Women", dim(lognormal_50_female_phia)[1])

# pooled
lognormal_50_female <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates_female$estimate_50,
                                       estimate_threshold = 50, conversion_threshold = 1000, 
                                       distribution = "lognormal",
                                       survey = PHIA_vls_estimates_female$survey,
                                       empirical_estimate = PHIA_vls_estimates_female$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates_female$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates_female$upper_ci_1000)

lognormal_50_female$distribution <- rep("lognormal", dim(lognormal_50_female)[1])
lognormal_50_female$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50_female)[1])
lognormal_50_female$sex <- rep("Women", dim(lognormal_50_female)[1])


#' calculate the root mean squared error between empirical and adjusted estimates
sqrt(mean((frec_50_female_phia$empirical_estimate - frec_50_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_female_phia$empirical_estimate - gamma_50_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_female_phia$empirical_estimate - lognormal_50_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_50_female$empirical_estimate - frec_50_female$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50_female$empirical_estimate - gamma_50_female$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50_female$empirical_estimate - lognormal_50_female$adjusted_vls_estimate)^2))


sqrt(mean((frec_200_female_phia$empirical_estimate - frec_200_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_female_phia$empirical_estimate - gamma_200_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_female_phia$empirical_estimate - lognormal_200_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_200_female$empirical_estimate - frec_200_female$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200_female$empirical_estimate - gamma_200_female$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200_female$empirical_estimate - lognormal_200_female$adjusted_vls_estimate)^2))

sqrt(mean((frec_400_female_phia$empirical_estimate - frec_400_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_female_phia$empirical_estimate - gamma_400_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_female_phia$empirical_estimate - lognormal_400_female_phia$adjusted_vls_estimate)^2))
sqrt(mean((frec_400_female$empirical_estimate - frec_400_female$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400_female$empirical_estimate - gamma_400_female$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400_female$empirical_estimate - lognormal_400_female$adjusted_vls_estimate)^2))

#' calculate the root mean squared error between empirical and adjusted estimates logit scale
sqrt(mean((qlogis(frec_50_female_phia$empirical_estimate) - qlogis(frec_50_female_phia$adjusted_vls_estimate))^2))
sqrt(mean((qlogis(gamma_50_female_phia$empirical_estimate) - qlogis(gamma_50_female_phia$adjusted_vls_estimate))^2))
sqrt(mean((qlogis(lognormal_50_female_phia$empirical_estimate) - qlogis(lognormal_50_female_phia$adjusted_vls_estimate))^2))

sqrt(mean((qlogis(frec_200_female_phia$empirical_estimate) - qlogis(frec_200_female_phia$adjusted_vls_estimate))^2))
sqrt(mean((qlogis(gamma_200_female_phia$empirical_estimate) - qlogis(gamma_200_female_phia$adjusted_vls_estimate))^2))
sqrt(mean((qlogis(lognormal_200_female_phia$empirical_estimate) - qlogis(lognormal_200_female_phia$adjusted_vls_estimate))^2))

sqrt(mean((qlogis(frec_400_female_phia$empirical_estimate) - qlogis(frec_400_female_phia$adjusted_vls_estimate))^2))
sqrt(mean((qlogis(gamma_400_female_phia$empirical_estimate) - qlogis(gamma_400_female_phia$adjusted_vls_estimate))^2))
sqrt(mean((qlogis(lognormal_400_female_phia$empirical_estimate) - qlogis(lognormal_400_female_phia$adjusted_vls_estimate))^2))

sqrt(mean((qlogis(frec_50_female$empirical_estimate) - qlogis(frec_50_female$adjusted_vls_estimate))^2))
sqrt(mean((qlogis(gamma_50_female$empirical_estimate) - qlogis(gamma_50_female$adjusted_vls_estimate))^2))
sqrt(mean((qlogis(lognormal_50_female$empirical_estimate) - qlogis(lognormal_50_female$adjusted_vls_estimate))^2))

sqrt(mean((qlogis(frec_200_female$empirical_estimate) - qlogis(frec_200_female$adjusted_vls_estimate))^2))
sqrt(mean((qlogis(gamma_200_female$empirical_estimate) - qlogis(gamma_200_female$adjusted_vls_estimate))^2))
sqrt(mean((qlogis(lognormal_200_female$empirical_estimate) - qlogis(lognormal_200_female$adjusted_vls_estimate))^2))

sqrt(mean((qlogis(frec_400_female$empirical_estimate) - qlogis(frec_400_female$adjusted_vls_estimate))^2))
sqrt(mean((qlogis(gamma_400_female$empirical_estimate) - qlogis(gamma_400_female$adjusted_vls_estimate))^2))
sqrt(mean((qlogis(lognormal_400_female$empirical_estimate) - qlogis(lognormal_400_female$adjusted_vls_estimate))^2))

# average bias
mean(frec_50_female_phia$adjusted_vls_estimate - frec_50_female_phia$empirical_estimate)
mean(gamma_50_female_phia$adjusted_vls_estimate - gamma_50_female_phia$empirical_estimate)
mean(lognormal_50_female_phia$adjusted_vls_estimate - lognormal_50_female_phia$empirical_estimate)
mean(frec_50_female$adjusted_vls_estimate - frec_50_female$empirical_estimate)
mean(gamma_50_female$adjusted_vls_estimate - gamma_50_female$empirical_estimate)
mean(lognormal_50_female$adjusted_vls_estimate - lognormal_50_female$empirical_estimate)


mean(frec_200_female_phia$adjusted_vls_estimate - frec_200_female_phia$empirical_estimate)
mean(gamma_200_female_phia$adjusted_vls_estimate - gamma_200_female_phia$empirical_estimate)
mean(lognormal_200_female_phia$adjusted_vls_estimate - lognormal_200_female_phia$empirical_estimate)
mean(frec_200_female$adjusted_vls_estimate - frec_200_female$empirical_estimate)
mean(gamma_200_female$adjusted_vls_estimate - gamma_200_female$empirical_estimate)
mean(lognormal_200_female$adjusted_vls_estimate - lognormal_200_female$empirical_estimate)

mean(frec_400_female_phia$adjusted_vls_estimate - frec_400_female_phia$empirical_estimate)
mean(gamma_400_female_phia$adjusted_vls_estimate - gamma_400_female_phia$empirical_estimate)
mean(lognormal_400_female_phia$adjusted_vls_estimate - lognormal_400_female_phia$empirical_estimate)
mean(frec_400_female$adjusted_vls_estimate - frec_400_female$empirical_estimate)
mean(gamma_400_female$adjusted_vls_estimate - gamma_400_female$empirical_estimate)
mean(lognormal_400_female$adjusted_vls_estimate - lognormal_400_female$empirical_estimate)

# mean absolute error
mean(abs(frec_50_female_phia$adjusted_vls_estimate - frec_50_female_phia$empirical_estimate))
mean(abs(gamma_50_female_phia$adjusted_vls_estimate - gamma_50_female_phia$empirical_estimate))
mean(abs(lognormal_50_female_phia$adjusted_vls_estimate - lognormal_50_female_phia$empirical_estimate))
mean(abs(frec_50_female$adjusted_vls_estimate - frec_50_female$empirical_estimate))
mean(abs(gamma_50_female$adjusted_vls_estimate - gamma_50_female$empirical_estimate))
mean(abs(lognormal_50_female$adjusted_vls_estimate - lognormal_50_female$empirical_estimate))

mean(abs(frec_200_female_phia$adjusted_vls_estimate - frec_200_female_phia$empirical_estimate))
mean(abs(gamma_200_female_phia$adjusted_vls_estimate - gamma_200_female_phia$empirical_estimate))
mean(abs(lognormal_200_female_phia$adjusted_vls_estimate - lognormal_200_female_phia$empirical_estimate))
mean(abs(frec_200_female$adjusted_vls_estimate - frec_200_female$empirical_estimate))
mean(abs(gamma_200_female$adjusted_vls_estimate - gamma_200_female$empirical_estimate))
mean(abs(lognormal_200_female$adjusted_vls_estimate - lognormal_200_female$empirical_estimate))

mean(abs(frec_400_female_phia$adjusted_vls_estimate - frec_400_female_phia$empirical_estimate))
mean(abs(gamma_400_female_phia$adjusted_vls_estimate - gamma_400_female_phia$empirical_estimate))
mean(abs(lognormal_400_female_phia$adjusted_vls_estimate - lognormal_400_female_phia$empirical_estimate))
mean(abs(frec_400_female$adjusted_vls_estimate - frec_400_female$empirical_estimate))
mean(abs(gamma_400_female$adjusted_vls_estimate - gamma_400_female$empirical_estimate))
mean(abs(lognormal_400_female$adjusted_vls_estimate - lognormal_400_female$empirical_estimate))


# code for Figure S7
sex_validation <- bind_rows(frec_50_male, frec_50_female, 
                            frec_50_male_phia, frec_50_female_phia,
                            frec_200_male, frec_200_female,
                            frec_200_male_phia, frec_200_female_phia,
                            frec_400_male, frec_400_female,
                            frec_400_male_phia, frec_400_female_phia,
                            gamma_50_male, gamma_50_female,
                            gamma_50_male_phia, gamma_50_female_phia,
                            gamma_200_male, gamma_200_female,
                            gamma_200_male_phia, gamma_200_female_phia,
                            gamma_400_male, gamma_400_female,
                            gamma_400_male_phia, gamma_400_female_phia,
                            lognormal_50_male, lognormal_50_female,
                            lognormal_50_male_phia, lognormal_50_female_phia,
                            lognormal_200_male, lognormal_200_female,
                            lognormal_200_male_phia, lognormal_200_female_phia,
                            lognormal_400_male, lognormal_400_female,
                            lognormal_400_male_phia, lognormal_400_female_phia,)

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
         dist_names = case_when(distribution == "lognormal" ~ "lognormal",
                                distribution == "lognormal Men" ~ "lognormal",
                                distribution == "lognormal Women" ~ "lognormal",
                                distribution == "Frechet" ~ "Frechet",
                                distribution == "Frechet Men" ~ "Frechet",
                                distribution == "Frechet Women" ~ "Frechet",
                                distribution == "Gamma" ~ "Gamma",
                                distribution == "Gamma Men" ~ "Gamma",
                                distribution == "Gamma Women" ~ "Gamma"),
         param = case_when(distribution == "lognormal" ~ "Johnson et al",
                           distribution == "lognormal Men" ~ "PHIA surveys",
                           distribution == "lognormal Women" ~ "PHIA surveys",
                           distribution == "Frechet" ~ "Johnson et al",
                           distribution == "Frechet Men" ~ "PHIA surveys",
                           distribution == "Frechet Women" ~ "PHIA surveys",
                           distribution == "Gamma" ~ "Johnson et al",
                           distribution == "Gamma Men" ~ "PHIA surveys",
                           distribution == "Gamma Women" ~ "PHIA surveys"),
         estimate_threshold = forcats::fct_relevel(estimate_threshold, 
                                                   c("<50 copies/mL", "<200 copies/mL", "<400 copies/mL")),
         dist_names = forcats::fct_relevel(dist_names,
                                           c("Frechet", "Gamma", "lognormal")),
         country_names = forcats::fct_relevel(country_names,
                                              c("Côte d'Ivoire (2017-18)","Cameroon (2017-18)",
                                                "Nigeria (2018)","Uganda (2016-17)","Zimbabwe (2015-16)",
                                                "Tanzania (2016-17)","Ethiopia (2017-18)","Lesotho (2016-17)",
                                                "Zambia (2016)","Mozambique (2021-22)","Rwanda (2018-19)",
                                                "Zimbabwe (2019-20)","Kenya (2018-19)","Malawi (2015-16)",
                                                "Namibia (2017)","Eswatini (2016-17)","Lesotho (2019-20)",
                                                "Eswatini (2021)","Zambia (2021)","Malawi (2020-21)", "Botswana (2021)")))

table(sex_validation$country_names, useNA = "always")


## scatter plot (Fig S7)
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
