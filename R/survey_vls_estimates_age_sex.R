#' calculate survey weighted VLS estimates for each survey at the defined thresholds
#' create function to create survey objects for all surveys
create_svy_object <- function(dataset){
  svy_object <- svydesign(ids = ~varunit, weights = ~btwt0, strata = ~varstrat,
                          nest = TRUE, data = dataset)
  return(svy_object)
}

svy_object_list <- lapply(X = dataset_subset, FUN = create_svy_object)

#' create function to calculate viral suppression estimates from survey data for all surveys
survey_vls_estimates_males <- function(svy_object){
  #' create subset of individuals sero-positive, on ART and with available viral load result
  onART_male <- svy_object[svy_object$variables$hivstatusfinal == 1 & svy_object$variables$gender == 1 &
                        svy_object$variables$art == 1 & !is.na(svy_object$variables$viral_sup1000)]
  
  res_1000 <- svyciprop(formula = ~I(viral_sup1000 == "Suppressed"),
                        design = onART_male, method = "mean",level = 0.95, df= 25)
  
  res_400 <- svyciprop(formula = ~I(viral_sup400 == "Suppressed"),
                       design = onART_male, method = "mean",level = 0.95, df= 25)
  
  res_200 <- svyciprop(formula = ~I(viral_sup200 == "Suppressed"),
                       design = onART_male, method = "mean",level = 0.95, df= 25)
  
  res_50 <- svyciprop(formula = ~I(viral_sup50 == "Suppressed"),
                      design = onART_male, method = "mean",level = 0.95, df= 25)
  
  res_df <- data.frame(country = unique(svy_object$variables$country), 
                       estimate_1000 = res_1000[[1]], lower_ci_1000 = attr(res_1000, "ci")[["2.5%"]], upper_ci_1000 = attr(res_1000, "ci")[["97.5%"]],
                       estimate_400 = res_400[[1]], lower_ci_400 = attr(res_400, "ci")[["2.5%"]], upper_ci_400 = attr(res_400, "ci")[["97.5%"]],
                       estimate_200 = res_200[[1]], lower_ci_200 = attr(res_200, "ci")[["2.5%"]], upper_ci_200 = attr(res_200, "ci")[["97.5%"]], 
                       estimate_50 = res_50[[1]], lower_ci_50 = attr(res_50, "ci")[["2.5%"]], upper_ci_50 = attr(res_50, "ci")[["97.5%"]])
  
  return(res_df)
} 


#' females
survey_vls_estimates_females <- function(svy_object){
  #' create subset of individuals sero-positive, on ART and with available viral load result
  onART_female <- svy_object[svy_object$variables$hivstatusfinal == 1 & svy_object$variables$gender == 2 &
                        svy_object$variables$art == 1 & !is.na(svy_object$variables$viral_sup1000)]
  
  res_1000 <- svyciprop(formula = ~I(viral_sup1000 == "Suppressed"),
                        design = onART_female, method = "mean",level = 0.95, df= 25)
  
  res_400 <- svyciprop(formula = ~I(viral_sup400 == "Suppressed"),
                       design = onART_female, method = "mean",level = 0.95, df= 25)
  
  res_200 <- svyciprop(formula = ~I(viral_sup200 == "Suppressed"),
                       design = onART_female, method = "mean",level = 0.95, df= 25)
  
  res_50 <- svyciprop(formula = ~I(viral_sup50 == "Suppressed"),
                      design = onART_female, method = "mean",level = 0.95, df= 25)
  
  res_df <- data.frame(country = unique(svy_object$variables$country), 
                       estimate_1000 = res_1000[[1]], lower_ci_1000 = attr(res_1000, "ci")[["2.5%"]], upper_ci_1000 = attr(res_1000, "ci")[["97.5%"]],
                       estimate_400 = res_400[[1]], lower_ci_400 = attr(res_400, "ci")[["2.5%"]], upper_ci_400 = attr(res_400, "ci")[["97.5%"]],
                       estimate_200 = res_200[[1]], lower_ci_200 = attr(res_200, "ci")[["2.5%"]], upper_ci_200 = attr(res_200, "ci")[["97.5%"]], 
                       estimate_50 = res_50[[1]], lower_ci_50 = attr(res_50, "ci")[["2.5%"]], upper_ci_50 = attr(res_50, "ci")[["97.5%"]])
  
  return(res_df)
} 

#' create dataframe with all estimates for all surveys
PHIA_vls_estimates_male <- purrr::map_df(lapply(svy_object_list, survey_vls_estimates_males), ~as.data.frame(.x), .id="survey")
PHIA_vls_estimates_female <- purrr::map_df(lapply(svy_object_list, survey_vls_estimates_females), ~as.data.frame(.x), .id="survey")


#' create function to calculate viral suppression estimates from survey data for all surveys
survey_vls_estimates_15_24 <- function(svy_object){
  #' create subset of individuals sero-positive, on ART and with available viral load result
  onART_15_24 <- svy_object[svy_object$variables$hivstatusfinal == 1 & svy_object$variables$age < 25 &
                             svy_object$variables$art == 1 & !is.na(svy_object$variables$viral_sup1000)]
  
  res_1000 <- svyciprop(formula = ~I(viral_sup1000 == "Suppressed"),
                        design = onART_15_24, method = "mean",level = 0.95, df= 25)
  
  res_400 <- svyciprop(formula = ~I(viral_sup400 == "Suppressed"),
                       design = onART_15_24, method = "mean",level = 0.95, df= 25)
  
  res_200 <- svyciprop(formula = ~I(viral_sup200 == "Suppressed"),
                       design = onART_15_24, method = "mean",level = 0.95, df= 25)
  
  res_50 <- svyciprop(formula = ~I(viral_sup50 == "Suppressed"),
                      design = onART_15_24, method = "mean",level = 0.95, df= 25)
  
  res_df <- data.frame(country = unique(svy_object$variables$country), 
                       estimate_1000 = res_1000[[1]], lower_ci_1000 = attr(res_1000, "ci")[["2.5%"]], upper_ci_1000 = attr(res_1000, "ci")[["97.5%"]],
                       estimate_400 = res_400[[1]], lower_ci_400 = attr(res_400, "ci")[["2.5%"]], upper_ci_400 = attr(res_400, "ci")[["97.5%"]],
                       estimate_200 = res_200[[1]], lower_ci_200 = attr(res_200, "ci")[["2.5%"]], upper_ci_200 = attr(res_200, "ci")[["97.5%"]], 
                       estimate_50 = res_50[[1]], lower_ci_50 = attr(res_50, "ci")[["2.5%"]], upper_ci_50 = attr(res_50, "ci")[["97.5%"]])
  
  return(res_df)
} 

#' create function to calculate viral suppression estimates from survey data for all surveys
survey_vls_estimates_25_34 <- function(svy_object){
  #' create subset of individuals sero-positive, on ART and with available viral load result
  onART_25_34 <- svy_object[svy_object$variables$hivstatusfinal == 1 & svy_object$variables$age > 24 &
                              svy_object$variables$age < 35 &
                              svy_object$variables$art == 1 & !is.na(svy_object$variables$viral_sup1000)]
  
  res_1000 <- svyciprop(formula = ~I(viral_sup1000 == "Suppressed"),
                        design = onART_25_34, method = "mean",level = 0.95, df= 25)
  
  res_400 <- svyciprop(formula = ~I(viral_sup400 == "Suppressed"),
                       design = onART_25_34, method = "mean",level = 0.95, df= 25)
  
  res_200 <- svyciprop(formula = ~I(viral_sup200 == "Suppressed"),
                       design = onART_25_34, method = "mean",level = 0.95, df= 25)
  
  res_50 <- svyciprop(formula = ~I(viral_sup50 == "Suppressed"),
                      design = onART_25_34, method = "mean",level = 0.95, df= 25)
  
  res_df <- data.frame(country = unique(svy_object$variables$country), 
                       estimate_1000 = res_1000[[1]], lower_ci_1000 = attr(res_1000, "ci")[["2.5%"]], upper_ci_1000 = attr(res_1000, "ci")[["97.5%"]],
                       estimate_400 = res_400[[1]], lower_ci_400 = attr(res_400, "ci")[["2.5%"]], upper_ci_400 = attr(res_400, "ci")[["97.5%"]],
                       estimate_200 = res_200[[1]], lower_ci_200 = attr(res_200, "ci")[["2.5%"]], upper_ci_200 = attr(res_200, "ci")[["97.5%"]], 
                       estimate_50 = res_50[[1]], lower_ci_50 = attr(res_50, "ci")[["2.5%"]], upper_ci_50 = attr(res_50, "ci")[["97.5%"]])
  
  return(res_df)
}

#' create function to calculate viral suppression estimates from survey data for all surveys
survey_vls_estimates_35_44 <- function(svy_object){
  #' create subset of individuals sero-positive, on ART and with available viral load result
  onART_35_44 <- svy_object[svy_object$variables$hivstatusfinal == 1 & svy_object$variables$age > 34 &
                              svy_object$variables$age < 45 &
                              svy_object$variables$art == 1 & !is.na(svy_object$variables$viral_sup1000)]
  
  res_1000 <- svyciprop(formula = ~I(viral_sup1000 == "Suppressed"),
                        design = onART_35_44, method = "mean",level = 0.95, df= 25)
  
  res_400 <- svyciprop(formula = ~I(viral_sup400 == "Suppressed"),
                       design = onART_35_44, method = "mean",level = 0.95, df= 25)
  
  res_200 <- svyciprop(formula = ~I(viral_sup200 == "Suppressed"),
                       design = onART_35_44, method = "mean",level = 0.95, df= 25)
  
  res_50 <- svyciprop(formula = ~I(viral_sup50 == "Suppressed"),
                      design = onART_35_44, method = "mean",level = 0.95, df= 25)
  
  res_df <- data.frame(country = unique(svy_object$variables$country), 
                       estimate_1000 = res_1000[[1]], lower_ci_1000 = attr(res_1000, "ci")[["2.5%"]], upper_ci_1000 = attr(res_1000, "ci")[["97.5%"]],
                       estimate_400 = res_400[[1]], lower_ci_400 = attr(res_400, "ci")[["2.5%"]], upper_ci_400 = attr(res_400, "ci")[["97.5%"]],
                       estimate_200 = res_200[[1]], lower_ci_200 = attr(res_200, "ci")[["2.5%"]], upper_ci_200 = attr(res_200, "ci")[["97.5%"]], 
                       estimate_50 = res_50[[1]], lower_ci_50 = attr(res_50, "ci")[["2.5%"]], upper_ci_50 = attr(res_50, "ci")[["97.5%"]])
  
  return(res_df)
}


#' create function to calculate viral suppression estimates from survey data for all surveys
survey_vls_estimates_45_54 <- function(svy_object){
  #' create subset of individuals sero-positive, on ART and with available viral load result
  onART_45_54 <- svy_object[svy_object$variables$hivstatusfinal == 1 & svy_object$variables$age > 44 &
                              svy_object$variables$age < 55 &
                              svy_object$variables$art == 1 & !is.na(svy_object$variables$viral_sup1000)]
  
  res_1000 <- svyciprop(formula = ~I(viral_sup1000 == "Suppressed"),
                        design = onART_45_54, method = "mean",level = 0.95, df= 25)
  
  res_400 <- svyciprop(formula = ~I(viral_sup400 == "Suppressed"),
                       design = onART_45_54, method = "mean",level = 0.95, df= 25)
  
  res_200 <- svyciprop(formula = ~I(viral_sup200 == "Suppressed"),
                       design = onART_45_54, method = "mean",level = 0.95, df= 25)
  
  res_50 <- svyciprop(formula = ~I(viral_sup50 == "Suppressed"),
                      design = onART_45_54, method = "mean",level = 0.95, df= 25)
  
  res_df <- data.frame(country = unique(svy_object$variables$country), 
                       estimate_1000 = res_1000[[1]], lower_ci_1000 = attr(res_1000, "ci")[["2.5%"]], upper_ci_1000 = attr(res_1000, "ci")[["97.5%"]],
                       estimate_400 = res_400[[1]], lower_ci_400 = attr(res_400, "ci")[["2.5%"]], upper_ci_400 = attr(res_400, "ci")[["97.5%"]],
                       estimate_200 = res_200[[1]], lower_ci_200 = attr(res_200, "ci")[["2.5%"]], upper_ci_200 = attr(res_200, "ci")[["97.5%"]], 
                       estimate_50 = res_50[[1]], lower_ci_50 = attr(res_50, "ci")[["2.5%"]], upper_ci_50 = attr(res_50, "ci")[["97.5%"]])
  
  return(res_df)
}


#' create function to calculate viral suppression estimates from survey data for all surveys
survey_vls_estimates_55 <- function(svy_object){
  #' create subset of individuals sero-positive, on ART and with available viral load result
  onART_55 <- svy_object[svy_object$variables$hivstatusfinal == 1 & svy_object$variables$age > 54 &
                              svy_object$variables$art == 1 & !is.na(svy_object$variables$viral_sup1000)]
  
  res_1000 <- svyciprop(formula = ~I(viral_sup1000 == "Suppressed"),
                        design = onART_55, method = "mean",level = 0.95, df= 25)
  
  res_400 <- svyciprop(formula = ~I(viral_sup400 == "Suppressed"),
                       design = onART_55, method = "mean",level = 0.95, df= 25)
  
  res_200 <- svyciprop(formula = ~I(viral_sup200 == "Suppressed"),
                       design = onART_55, method = "mean",level = 0.95, df= 25)
  
  res_50 <- svyciprop(formula = ~I(viral_sup50 == "Suppressed"),
                      design = onART_55, method = "mean",level = 0.95, df= 25)
  
  res_df <- data.frame(country = unique(svy_object$variables$country), 
                       estimate_1000 = res_1000[[1]], lower_ci_1000 = attr(res_1000, "ci")[["2.5%"]], upper_ci_1000 = attr(res_1000, "ci")[["97.5%"]],
                       estimate_400 = res_400[[1]], lower_ci_400 = attr(res_400, "ci")[["2.5%"]], upper_ci_400 = attr(res_400, "ci")[["97.5%"]],
                       estimate_200 = res_200[[1]], lower_ci_200 = attr(res_200, "ci")[["2.5%"]], upper_ci_200 = attr(res_200, "ci")[["97.5%"]], 
                       estimate_50 = res_50[[1]], lower_ci_50 = attr(res_50, "ci")[["2.5%"]], upper_ci_50 = attr(res_50, "ci")[["97.5%"]])
  
  return(res_df)
}

#' create dataframe with all estimates for all surveys
PHIA_vls_estimates_15_24 <- purrr::map_df(lapply(svy_object_list, survey_vls_estimates_15_24), ~as.data.frame(.x), .id="survey")
PHIA_vls_estimates_25_34 <- purrr::map_df(lapply(svy_object_list, survey_vls_estimates_25_34), ~as.data.frame(.x), .id="survey")
PHIA_vls_estimates_35_44 <- purrr::map_df(lapply(svy_object_list, survey_vls_estimates_35_44), ~as.data.frame(.x), .id="survey")
PHIA_vls_estimates_45_54 <- purrr::map_df(lapply(svy_object_list, survey_vls_estimates_45_54), ~as.data.frame(.x), .id="survey")
PHIA_vls_estimates_55 <- purrr::map_df(lapply(svy_object_list, survey_vls_estimates_55), ~as.data.frame(.x), .id="survey")

# correct viral load suppression estimates greater than 100% error
PHIA_vls_estimates_55$upper_ci_1000 <- ifelse(PHIA_vls_estimates_55$upper_ci_1000 > 1, 1, 
                                              PHIA_vls_estimates_55$upper_ci_1000)
summary(PHIA_vls_estimates_55$upper_ci_1000)
