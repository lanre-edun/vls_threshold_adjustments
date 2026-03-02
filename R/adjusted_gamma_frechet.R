# For a given shape, reported threshold, and reported cumulative estimate,
# this functions solves for the gamma, Frechet and lognormal scale parameter that reproduces the
# reported estimate, then returns the adjusted cumulative probability at 1000.
# Functions are vectorized to allow element-wise application over inputs.
# vectorize the functions
gamma_func_vec <- Vectorize(gamma_func <- function(shape, reported_threshold, reported_est){
  myfun <- function(scale) {
    shape <- shape
    pgamma(log10(reported_threshold), shape=shape, scale=scale) - reported_est
  }
  
  tmp <- uniroot(myfun, lower=0.1, upper=10)
  
  myscale <- tmp$root
  
  adj_vls_est <- pgamma(log10(1000), shape = shape, scale = myscale)
  return(adj_vls_est)
  
}
)

frechet_func_vec <- Vectorize(frechet_func <- function(shape, reported_threshold, reported_est){
  myfun <- function(scale) {
    shape <- shape
    pfrechet(log10(reported_threshold), shape=shape, scale=scale) - reported_est
  }
  
  tmp <- uniroot(myfun, lower=0.1, upper=10)
  
  myscale <- tmp$root
  
  adj_vls_est <- pfrechet(log10(1000), shape = shape, scale = myscale)
  return(adj_vls_est)
  
})

lognorm_func_vec <- Vectorize(lognorm_func <- function(distribution, sdlog, reported_threshold, reported_est){
  myfun <- function(meanlog) {
    sdlog <- sdlog
    plnorm(log10(reported_threshold), meanlog=meanlog, sdlog=sdlog) - reported_est
  }
  
  tmp <- uniroot(myfun, lower=-2, upper=10)
  
  mymeanlog <- tmp$root
  
  adj_vls_est <- plnorm(log10(1000), meanlog = mymeanlog, sdlog = sdlog)
  return(adj_vls_est)
  
})

#' create function using the Weibull, reverse Weibull, Pareto, Frechet, gamma and lognormal distributions to convert viral suppression estimates
#' @details probability denisity functions to obtain viral suppression conversion thresholds
#' obtained from Leigh et al paper
#' @param function receives viral_sup_estimate: estimate of empirical viral suppression calculated from survey
#' reported_threshold: threshold which estimate is reported from
#' adjusted_threshold: threshold which estimate is being adjusted to 
#' distribution: distribution being explored
vls_prob_threshold <- function(viral_sup_estimate, estimate_threshold, conversion_threshold, distribution, survey,
                               empirical_estimate, empirical_lower_ci, empirical_upper_ci){
  
  if(distribution == "Reverse Weibull"){
    
    adjusted_vls_estimate <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^2.81)
    adjusted_vls_lower_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^1.70)
    adjusted_vls_upper_ci <- (viral_sup_estimate)^(((6-log10(conversion_threshold))/(6-log10(estimate_threshold)))^3.92)
    
    res <- data.frame(survey = survey, 
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Reverse Weibull distribution")
    
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
    
    
  }else if(distribution == "Frechet"){
    
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
    
    
  }else if(distribution == "Lognormal"){
    
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
    
    print("Lognormal distribution")
    
    return(res)
    
    
}else if(distribution == "Weibull"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.85)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^0.43)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)^((((log10(conversion_threshold))/log10(estimate_threshold)))^1.26)
    
    res <- data.frame(survey = survey,
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci,empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Weibull distribution")
    
    return(res)
    
  } else{
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.73)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^2.26)
    
    res <- data.frame(survey = survey,
                      adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto distribution")
    
    return(res)
    
  }
  
}


#' calculate adjusted viral suppression estimates
#' Reverse Weibull distribution < 400 copies/ml
rev_weib_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                                   estimate_threshold = 400, conversion_threshold = 1000, 
                                   distribution = "Reverse Weibull", 
                                   survey = PHIA_vls_estimates$survey,
                                   empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

rev_weib_400$distribution <- rep("Reverse Weibull", dim(rev_weib_400)[1])
rev_weib_400$estimate_threshold <- rep("<400 copies/mL", dim(rev_weib_400)[1])

#' Reverse Weibull distribution < 200 copies/ml
rev_weib_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                                   estimate_threshold = 200, conversion_threshold = 1000, 
                                   distribution = "Reverse Weibull",
                                   survey = PHIA_vls_estimates$survey, 
                                   empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

rev_weib_200$distribution <- rep("Reverse Weibull", dim(rev_weib_200)[1])
rev_weib_200$estimate_threshold <- rep("<200 copies/mL", dim(rev_weib_200)[1])

#' Reverse Weibull distribution < 50 copies/ml
rev_weib_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                                  estimate_threshold = 50, conversion_threshold = 1000, 
                                  distribution = "Reverse Weibull", 
                                  survey = PHIA_vls_estimates$survey,
                                  empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                  empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000, 
                                  empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

rev_weib_50$distribution <- rep("Reverse Weibull", dim(rev_weib_50)[1])
rev_weib_50$estimate_threshold <- rep("<50 copies/mL", dim(rev_weib_50)[1])


#' Weibull distribution < 400 copies/ml
weib_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                               estimate_threshold = 400, 
                               conversion_threshold = 1000, 
                               distribution = "Weibull", 
                               survey = PHIA_vls_estimates$survey, 
                               empirical_estimate = PHIA_vls_estimates$estimate_1000,
                               empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                               empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

weib_400$distribution <- rep("Weibull", dim(weib_400)[1])
weib_400$estimate_threshold <- rep("<400 copies/mL", dim(weib_400)[1])

#' Weibull distribution < 200 copies/ml
weib_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200, 
                               estimate_threshold = 200, conversion_threshold = 1000, 
                               distribution = "Weibull", 
                               survey = PHIA_vls_estimates$survey, 
                               empirical_estimate = PHIA_vls_estimates$estimate_1000,
                               empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000, 
                               empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

weib_200$distribution <- rep("Weibull", dim(weib_200)[1])
weib_200$estimate_threshold <- rep("<200 copies/mL", dim(weib_200)[1])

#' Weibull distribution < 50 copies/ml
weib_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50, 
                              estimate_threshold = 50, conversion_threshold = 1000, 
                              distribution = "Weibull", 
                              survey = PHIA_vls_estimates$survey, 
                              empirical_estimate = PHIA_vls_estimates$estimate_1000,
                              empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                              empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

weib_50$distribution <- rep("Weibull", dim(weib_50)[1])
weib_50$estimate_threshold <- rep("<50 copies/mL", dim(weib_50)[1])


#' Pareto distribution < 400 copies/ml
pareto_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                                 estimate_threshold = 400, conversion_threshold = 1000, 
                                 distribution = "Pareto",
                                 survey = PHIA_vls_estimates$survey, 
                                 empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                 empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000, 
                                 empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

pareto_400$distribution <- rep("Pareto", dim(pareto_400)[1])
pareto_400$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400)[1])

#' Pareto distribution < 200 copies/ml
pareto_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                                 estimate_threshold = 200, conversion_threshold = 1000, 
                                 distribution = "Pareto",
                                 survey = PHIA_vls_estimates$survey,
                                 empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                 empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000, 
                                 empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

pareto_200$distribution <- rep("Pareto", dim(pareto_200)[1])
pareto_200$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200)[1])

#' Pareto distribution < 50 copies/ml
pareto_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                                estimate_threshold = 50, conversion_threshold = 1000, 
                                distribution = "Pareto", 
                                survey = PHIA_vls_estimates$survey, 
                                empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

pareto_50$distribution <- rep("Pareto", dim(pareto_50)[1])
pareto_50$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50)[1])

#' Gamma distribution < 400 copies/ml
gamma_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                                 estimate_threshold = 400, conversion_threshold = 1000, 
                                 distribution = "Gamma",
                                 survey = PHIA_vls_estimates$survey, 
                                 empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                 empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000, 
                                 empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

gamma_400$distribution <- rep("Gamma", dim(gamma_400)[1])
gamma_400$estimate_threshold <- rep("<400 copies/mL", dim(gamma_400)[1])

#' gamma distribution < 200 copies/ml
gamma_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                                 estimate_threshold = 200, conversion_threshold = 1000, 
                                 distribution = "Gamma",
                                 survey = PHIA_vls_estimates$survey,
                                 empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                 empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000, 
                                 empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

gamma_200$distribution <- rep("Gamma", dim(gamma_200)[1])
gamma_200$estimate_threshold <- rep("<200 copies/mL", dim(gamma_200)[1])

#' gamma distribution < 50 copies/ml
gamma_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                                estimate_threshold = 50, conversion_threshold = 1000, 
                                distribution = "Gamma", 
                                survey = PHIA_vls_estimates$survey, 
                                empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

gamma_50$distribution <- rep("Gamma", dim(gamma_50)[1])
gamma_50$estimate_threshold <- rep("<50 copies/mL", dim(gamma_50)[1])


#' frechet distribution < 400 copies/ml
frechet_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                                estimate_threshold = 400, conversion_threshold = 1000, 
                                distribution = "Frechet",
                                survey = PHIA_vls_estimates$survey, 
                                empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000, 
                                empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

frechet_400$distribution <- rep("Frechet", dim(frechet_400)[1])
frechet_400$estimate_threshold <- rep("<400 copies/mL", dim(frechet_400)[1])

#' frechet distribution < 200 copies/ml
frechet_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                                estimate_threshold = 200, conversion_threshold = 1000, 
                                distribution = "Frechet",
                                survey = PHIA_vls_estimates$survey,
                                empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000, 
                                empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

frechet_200$distribution <- rep("Frechet", dim(frechet_200)[1])
frechet_200$estimate_threshold <- rep("<200 copies/mL", dim(frechet_200)[1])

#' frechet distribution < 50 copies/ml
frechet_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                               estimate_threshold = 50, conversion_threshold = 1000, 
                               distribution = "Frechet", 
                               survey = PHIA_vls_estimates$survey, 
                               empirical_estimate = PHIA_vls_estimates$estimate_1000,
                               empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                               empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

frechet_50$distribution <- rep("Frechet", dim(frechet_50)[1])
frechet_50$estimate_threshold <- rep("<50 copies/mL", dim(frechet_50)[1])


#' lognormal distribution < 400 copies/ml
lognormal_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                                  estimate_threshold = 400, conversion_threshold = 1000, 
                                  distribution = "Lognormal",
                                  survey = PHIA_vls_estimates$survey, 
                                  empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                  empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000, 
                                  empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

lognormal_400$distribution <- rep("Lognormal", dim(lognormal_400)[1])
lognormal_400$estimate_threshold <- rep("<400 copies/mL", dim(lognormal_400)[1])

#' lognormal distribution < 200 copies/ml
lognormal_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                                  estimate_threshold = 200, conversion_threshold = 1000, 
                                  distribution = "Lognormal",
                                  survey = PHIA_vls_estimates$survey,
                                  empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                  empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000, 
                                  empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

lognormal_200$distribution <- rep("Lognormal", dim(lognormal_200)[1])
lognormal_200$estimate_threshold <- rep("<200 copies/mL", dim(lognormal_200)[1])

#' lognormal distribution < 50 copies/ml
lognormal_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                                 estimate_threshold = 50, conversion_threshold = 1000, 
                                 distribution = "Lognormal", 
                                 survey = PHIA_vls_estimates$survey, 
                                 empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                 empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                 empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

lognormal_50$distribution <- rep("Lognormal", dim(lognormal_50)[1])
lognormal_50$estimate_threshold <- rep("<50 copies/mL", dim(lognormal_50)[1])



#' calculate the root mean squared error between empirical and adjusted estimates
sqrt(mean((rev_weib_50$empirical_estimate - rev_weib_50$adjusted_vls_estimate)^2))
sqrt(mean((weib_50$empirical_estimate - weib_50$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50$empirical_estimate - pareto_50$adjusted_vls_estimate)^2))
sqrt(mean((frechet_50$empirical_estimate - frechet_50$adjusted_vls_estimate)^2))
sqrt(mean((gamma_50$empirical_estimate - gamma_50$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_50$empirical_estimate - lognormal_50$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_200$empirical_estimate - rev_weib_200$adjusted_vls_estimate)^2))
sqrt(mean((weib_200$empirical_estimate - weib_200$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200$empirical_estimate - pareto_200$adjusted_vls_estimate)^2))
sqrt(mean((frechet_200$empirical_estimate - frechet_200$adjusted_vls_estimate)^2))
sqrt(mean((gamma_200$empirical_estimate - gamma_200$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_200$empirical_estimate - lognormal_200$adjusted_vls_estimate)^2))

sqrt(mean((rev_weib_400$empirical_estimate - rev_weib_400$adjusted_vls_estimate)^2))
sqrt(mean((weib_400$empirical_estimate - weib_400$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400$empirical_estimate - pareto_400$adjusted_vls_estimate)^2))
sqrt(mean((frechet_400$empirical_estimate - frechet_400$adjusted_vls_estimate)^2))
sqrt(mean((gamma_400$empirical_estimate - gamma_400$adjusted_vls_estimate)^2))
sqrt(mean((lognormal_400$empirical_estimate - lognormal_400$adjusted_vls_estimate)^2))


# average bias
mean(rev_weib_50$adjusted_vls_estimate - rev_weib_50$empirical_estimate)
mean(weib_50$adjusted_vls_estimate - weib_50$empirical_estimate)
mean(pareto_50$adjusted_vls_estimate - pareto_50$empirical_estimate)
mean(frechet_50$adjusted_vls_estimate - frechet_50$empirical_estimate)
mean(gamma_50$adjusted_vls_estimate - gamma_50$empirical_estimate)
mean(lognormal_50$adjusted_vls_estimate - lognormal_50$empirical_estimate)

mean(rev_weib_200$adjusted_vls_estimate - rev_weib_200$empirical_estimate)
mean(weib_200$adjusted_vls_estimate - weib_200$empirical_estimate)
mean(pareto_200$adjusted_vls_estimate - pareto_200$empirical_estimate)
mean(frechet_200$adjusted_vls_estimate - frechet_200$empirical_estimate)
mean(gamma_200$adjusted_vls_estimate - gamma_200$empirical_estimate)
mean(lognormal_200$adjusted_vls_estimate - lognormal_200$empirical_estimate)

mean(rev_weib_400$adjusted_vls_estimate - rev_weib_400$empirical_estimate)
mean(weib_400$adjusted_vls_estimate - weib_400$empirical_estimate)
mean(pareto_400$adjusted_vls_estimate - pareto_400$empirical_estimate)
mean(frechet_400$adjusted_vls_estimate - frechet_400$empirical_estimate)
mean(gamma_400$adjusted_vls_estimate - gamma_400$empirical_estimate)
mean(lognormal_400$adjusted_vls_estimate - lognormal_400$empirical_estimate)


# mean absolute error
mean(abs(rev_weib_50$adjusted_vls_estimate - rev_weib_50$empirical_estimate))
mean(abs(weib_50$adjusted_vls_estimate - weib_50$empirical_estimate))
mean(abs(pareto_50$adjusted_vls_estimate - pareto_50$empirical_estimate))
mean(abs(frechet_50$adjusted_vls_estimate - frechet_50$empirical_estimate))
mean(abs(gamma_50$adjusted_vls_estimate - gamma_50$empirical_estimate))
mean(abs(lognormal_50$adjusted_vls_estimate - lognormal_50$empirical_estimate))

mean(abs(rev_weib_200$adjusted_vls_estimate - rev_weib_200$empirical_estimate))
mean(abs(weib_200$adjusted_vls_estimate - weib_200$empirical_estimate))
mean(abs(pareto_200$adjusted_vls_estimate - pareto_200$empirical_estimate))
mean(abs(frechet_200$adjusted_vls_estimate - frechet_200$empirical_estimate))
mean(abs(gamma_200$adjusted_vls_estimate - gamma_200$empirical_estimate))
mean(abs(lognormal_200$adjusted_vls_estimate - lognormal_200$empirical_estimate))

mean(abs(rev_weib_400$adjusted_vls_estimate - rev_weib_400$empirical_estimate))
mean(abs(weib_400$adjusted_vls_estimate - weib_400$empirical_estimate))
mean(abs(pareto_400$adjusted_vls_estimate - pareto_400$empirical_estimate))
mean(abs(frechet_400$adjusted_vls_estimate - frechet_400$empirical_estimate))
mean(abs(gamma_400$adjusted_vls_estimate - gamma_400$empirical_estimate))
mean(abs(lognormal_400$adjusted_vls_estimate - lognormal_400$empirical_estimate))

#' scatter plot showing relationship between empirical and adjusted viral suppression estimates
#' combine all estimates derived in adjusted_vls_estimates script
six_distributions <- bind_rows(rev_weib_50, rev_weib_200, rev_weib_400,
                                weib_50, weib_200, weib_400,
                                pareto_50, pareto_200, pareto_400,
                                gamma_50, gamma_200, gamma_400,
                                frechet_50, frechet_200, frechet_400, 
                                lognormal_50, lognormal_200, lognormal_400)

#' update country names to survey names
six_distributions <- six_distributions %>% 
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
six_distributions$estimate_threshold <- factor(six_distributions$estimate_threshold, 
                                               levels = c("<50 copies/mL", "<200 copies/mL", "<400 copies/mL"))

#' order distributions for plotting
six_distributions$distribution <- factor(six_distributions$distribution,
                                         levels = c("Frechet","Gamma","Lognormal","Reverse Weibull","Weibull","Pareto"))

six_distributions$country_names <- factor(six_distributions$country_names,
                                          levels = c("Côte d'Ivoire (2017-18)","Cameroon (2017-18)",
                                                     "Nigeria (2018)","Uganda (2016-17)","Zimbabwe (2015-16)",
                                                     "Tanzania (2016-17)","Ethiopia (2017-18)","Lesotho (2016-17)",
                                                     "Zambia (2016)","Mozambique (2021-22)","Rwanda (2018-19)",
                                                     "Zimbabwe (2019-20)","Kenya (2018-19)","Malawi (2015-16)",
                                                     "Namibia (2017)","Eswatini (2016-17)","Lesotho (2019-20)",
                                                     "Eswatini (2021)","Zambia (2021)","Malawi (2020-21)", "Botswana (2021)"))

# code for Figure 3 and S5
# scatter plot
six_distributions %>%
  filter(estimate_threshold == "<400 copies/mL") %>%
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
  #facet_grid(estimate_threshold~distribution) + 
  facet_wrap(~distribution) +
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
six_distributions %>%
  filter(distribution == "Lognormal" &
           estimate_threshold == "<400 copies/mL") %>%
  filter(adjusted_vls_estimate > empirical_estimate) %>%
  select(survey, estimate_threshold, adjusted_vls_estimate, empirical_estimate)