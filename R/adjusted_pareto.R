#' create function using the Pareto distribution to convert viral suppression estimates
#' @details probability denisity functions to obtain viral suppression conversion thresholds
#' obtained from Leigh et al paper
#' @param function receives viral_sup_estimate: estimate of empirical viral suppression calculated from survey
#' reported_threshold: threshold which estimate is reported from
#' adjusted_threshold: threshold which estimate is being adjusted to 
#' distribution: distribution being explored
vls_prob_threshold <- function(viral_sup_estimate, estimate_threshold, conversion_threshold, distribution, survey,
                               empirical_estimate, empirical_lower_ci, empirical_upper_ci){
  
  if(distribution == "Pareto Johnson LF et al"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.73)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^2.26)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto Johnson LF et al")
    
    return(res)
    
  }else if(distribution == "Pareto (1.20)"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^1.20)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto (1.20)")
    
    return(res)
    
  } else if(distribution == "Pareto (2.00)"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^2.00)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^2.00)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^2.00)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto (2.00)")
    
    return(res)
    
    
  }else if(distribution == "Pareto (2.50)"){
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^2.50)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^2.50)
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^2.50)
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto (2.50)")
    
    return(res)
    
  }else{
    
    adjusted_vls_estimate <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^3.00)
    adjusted_vls_lower_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^3.00) 
    adjusted_vls_upper_ci <- 1-(1-viral_sup_estimate)*((((log10(estimate_threshold))/log10(conversion_threshold)))^3.00) 
    
    res <- data.frame(survey = survey, adjusted_vls_estimate = adjusted_vls_estimate, adjusted_vls_lower_ci = adjusted_vls_lower_ci,
                      adjusted_vls_upper_ci = adjusted_vls_upper_ci, empirical_estimate = empirical_estimate,
                      empirical_lower_ci = empirical_lower_ci, empirical_upper_ci = empirical_upper_ci)
    
    print("Pareto (3.00)")
    
    return(res)
    
  }
  
}

#' Pareto (1.20) distribution < 400 copies/ml
par_one_two_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400, 
                                   estimate_threshold = 400, conversion_threshold = 1000, 
                                   distribution = "Pareto (1.20)",
                                   survey = PHIA_vls_estimates$survey,
                                   empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

par_one_two_400$distribution <- rep("shape = 1.20", dim(par_one_two_400)[1])
par_one_two_400$estimate_threshold <- rep("<400 copies/mL", dim(par_one_two_400)[1])

#' Pareto (1.50) distribution < 200 copies/ml
par_one_two_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                                   estimate_threshold = 200, conversion_threshold = 1000, 
                                   distribution = "Pareto (1.20)",
                                   survey = PHIA_vls_estimates$survey,
                                   empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                   empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                   empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

par_one_two_200$distribution <- rep("shape = 1.20", dim(par_one_two_200)[1])
par_one_two_200$estimate_threshold <- rep("<200 copies/mL", dim(par_one_two_200)[1])

#' Pareto (1.20) distribution < 50 copies/ml
par_one_two_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                                  estimate_threshold = 50, conversion_threshold = 1000, 
                                  distribution = "Pareto (1.20)",
                                  survey = PHIA_vls_estimates$survey, 
                                  empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                  empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                  empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

par_one_two_50$distribution <- rep("shape = 1.20", dim(par_one_two_50)[1])
par_one_two_50$estimate_threshold <- rep("<50 copies/mL", dim(par_one_two_50)[1])


#' Pareto (2.00) < 400 copies/ml
par_two_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                               estimate_threshold = 400, conversion_threshold = 1000, 
                               distribution = "Pareto (2.00)",
                               survey = PHIA_vls_estimates$survey,
                               empirical_estimate = PHIA_vls_estimates$estimate_1000,
                               empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                               empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

par_two_400$distribution <- rep("shape = 2.00", dim(par_two_400)[1])
par_two_400$estimate_threshold <- rep("<400 copies/mL", dim(par_two_400)[1])

#' Pareto (2.00) < 200 copies/ml
par_two_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                               estimate_threshold = 200, conversion_threshold = 1000, 
                               distribution = "Pareto (2.00)",
                               survey = PHIA_vls_estimates$survey,
                               empirical_estimate = PHIA_vls_estimates$estimate_1000,
                               empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                               empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

par_two_200$distribution <- rep("shape = 2.00", dim(par_two_200)[1])
par_two_200$estimate_threshold <- rep("<200 copies/mL", dim(par_two_200)[1])

#' Pareto (2.00) < 50 copies/ml
par_two_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                              estimate_threshold = 50, conversion_threshold = 1000, 
                              distribution = "Pareto (2.00)",
                              survey = PHIA_vls_estimates$survey, 
                              empirical_estimate = PHIA_vls_estimates$estimate_1000,
                              empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                              empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

par_two_50$distribution <- rep("shape = 2.00", dim(par_two_50)[1])
par_two_50$estimate_threshold <- rep("<50 copies/mL", dim(par_two_50)[1])


#' Pareto distribution < 400 copies/ml
pareto_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                                 estimate_threshold = 400, conversion_threshold = 1000, 
                                 distribution = "Pareto Johnson LF et al",
                                 survey = PHIA_vls_estimates$survey,
                                 empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                 empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                 empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

pareto_400$distribution <- rep("shape = 1.73\n(Johnson et al.)", dim(pareto_400)[1])
pareto_400$estimate_threshold <- rep("<400 copies/mL", dim(pareto_400)[1])

#' Pareto distribution < 200 copies/ml
pareto_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                                 estimate_threshold = 200, conversion_threshold = 1000, 
                                 distribution = "Pareto Johnson LF et al",
                                 survey = PHIA_vls_estimates$survey,
                                 empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                 empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                 empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

pareto_200$distribution <- rep("shape = 1.73\n(Johnson et al.)", dim(pareto_200)[1])
pareto_200$estimate_threshold <- rep("<200 copies/mL", dim(pareto_200)[1])

#' Pareto distribution < 50 copies/ml
pareto_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                                estimate_threshold = 50, conversion_threshold = 1000, 
                                distribution = "Pareto Johnson LF et al",
                                survey = PHIA_vls_estimates$survey, 
                                empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

pareto_50$distribution <- rep("shape = 1.73\n(Johnson et al.)", dim(pareto_50)[1])
pareto_50$estimate_threshold <- rep("<50 copies/mL", dim(pareto_50)[1])

#' Pareto (2.50) < 400 copies/ml
pareto_two_five_400 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_400,
                                       estimate_threshold = 400, conversion_threshold = 1000, 
                                       distribution = "Pareto (2.50)",
                                       survey = PHIA_vls_estimates$survey,
                                       empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

pareto_two_five_400$distribution <- rep("shape = 2.50", dim(pareto_two_five_400)[1])
pareto_two_five_400$estimate_threshold <- rep("<400 copies/mL", dim(pareto_two_five_400)[1])

#' Reverse Weibull distribution with updated parameters < 200 copies/ml
pareto_two_five_200 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_200,
                                       estimate_threshold = 200, conversion_threshold = 1000, 
                                       distribution = "Pareto (2.50)",
                                       survey = PHIA_vls_estimates$survey,
                                       empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                       empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                       empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

pareto_two_five_200$distribution <- rep("shape = 2.50", dim(pareto_two_five_200)[1])
pareto_two_five_200$estimate_threshold <- rep("<200 copies/mL", dim(pareto_two_five_200)[1])

#' Reverse Weibull distribution with updated parameters < 50 copies/ml
pareto_two_five_50 <- vls_prob_threshold(viral_sup_estimate = PHIA_vls_estimates$estimate_50,
                                      estimate_threshold = 50, conversion_threshold = 1000, 
                                      distribution = "Pareto (2.50)",
                                      survey = PHIA_vls_estimates$survey,
                                      empirical_estimate = PHIA_vls_estimates$estimate_1000,
                                      empirical_lower_ci = PHIA_vls_estimates$lower_ci_1000,
                                      empirical_upper_ci = PHIA_vls_estimates$upper_ci_1000)

pareto_two_five_50$distribution <- rep("shape = 2.50", dim(pareto_two_five_50)[1])
pareto_two_five_50$estimate_threshold <- rep("<50 copies/mL", dim(pareto_two_five_50)[1])

#' calculate the root mean squared error between empirical and adjusted estimates
sqrt(mean((par_one_two_50$empirical_estimate - par_one_two_50$adjusted_vls_estimate)^2))
sqrt(mean((pareto_50$empirical_estimate - pareto_50$adjusted_vls_estimate)^2))
sqrt(mean((par_two_50$empirical_estimate - par_two_50$adjusted_vls_estimate)^2))
sqrt(mean((pareto_two_five_50$empirical_estimate - pareto_two_five_50$adjusted_vls_estimate)^2))
      
sqrt(mean((par_one_two_200$empirical_estimate - par_one_two_200$adjusted_vls_estimate)^2))
sqrt(mean((pareto_200$empirical_estimate - pareto_200$adjusted_vls_estimate)^2))
sqrt(mean((par_two_200$empirical_estimate - par_two_200$adjusted_vls_estimate)^2))
sqrt(mean((pareto_two_five_200$empirical_estimate - pareto_two_five_200$adjusted_vls_estimate)^2))

sqrt(mean((par_one_two_400$empirical_estimate - par_one_two_400$adjusted_vls_estimate)^2))
sqrt(mean((pareto_400$empirical_estimate - pareto_400$adjusted_vls_estimate)^2))
sqrt(mean((par_two_400$empirical_estimate - par_two_400$adjusted_vls_estimate)^2))
sqrt(mean((pareto_two_five_400$empirical_estimate - pareto_two_five_400$adjusted_vls_estimate)^2))


#' PLOT
#' scatter plot showing relationship between empirical and adjusted viral suppression estimates
#' combine all estimates derived in adjusted_vls_estimates script
all_distributions_pareto <- bind_rows(par_one_two_50, pareto_50, par_two_50, pareto_two_five_50,
                                      par_one_two_200, pareto_200, par_two_200, pareto_two_five_200,
                                      par_one_two_400, pareto_400, par_two_400, pareto_two_five_400)


#' update country names to survey names
all_distributions_pareto <- all_distributions_pareto %>%
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
         estimate_threshold = forcats::fct_relevel(estimate_threshold, c("<50 copies/mL", "<200 copies/mL", "<400 copies/mL")),
         distribution = forcats::fct_relevel(distribution, c("shape = 1.73\n(Johnson et al.)", "shape = 1.20", "shape = 2.00", "shape = 2.50")),
         country_names = forcats::fct_relevel(country_names, 
                                              c("Côte d'Ivoire (2017-18)","Cameroon (2017-18)",
                                                "Nigeria (2018)","Uganda (2016-17)","Zimbabwe (2015-16)",
                                                "Tanzania (2016-17)","Ethiopia (2017-18)","Lesotho (2016-17)",
                                                "Zambia (2016)","Mozambique (2021-22)","Rwanda (2018-19)",
                                                "Zimbabwe (2019-20)","Kenya (2018-19)","Malawi (2015-16)",
                                                "Namibia (2017)","Eswatini (2016-17)","Lesotho (2019-20)",
                                                "Eswatini (2021)","Zambia (2021)","Malawi (2020-21)", "Botswana (2021)")))

names(all_distributions_pareto)

# Figure S3 code
## scatter plot
all_distributions_pareto %>%
  filter(distribution %in% c("shape = 1.73\n(Johnson et al.)", "shape = 1.20")) %>%
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
  facet_grid(distribution~estimate_threshold) + 
  labs(x = "% viral suppression (PHIA survey)", 
       y = "% viral suppression (adjusted)", color = "",
       title = "Using the Pareto model",) + 
  theme(legend.position = "bottom",
        axis.text = element_text(size = rel(1), family = "sans"),
        axis.title = element_text(size = rel(1.0), family = "sans"),
        strip.text = element_text(size = rel(1.0), family = "sans", face = "bold"),
        legend.text = element_text(size = rel(0.9), family = "sans"),
        panel.grid = element_blank(), 
        panel.spacing.x = unit(1, "lines"),
        plot.title = element_text(size = rel(1.0), family = "sans", face = "bold"))
