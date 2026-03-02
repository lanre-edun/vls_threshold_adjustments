#' calculate survey weighted VLS estimates for each survey at the defined thresholds
#' create function to create survey objects for all surveys
create_svy_object <- function(dataset){
  svy_object <- svydesign(ids = ~varunit, weights = ~btwt0, strata = ~varstrat,
                          nest = TRUE, data = dataset)
  return(svy_object)
}

svy_object_list <- lapply(X = dataset_subset, FUN = create_svy_object)

#' create function to calculate viral suppression estimates from survey data for all surveys
survey_vls_estimates <- function(svy_object){
  #' create subset of individuals sero-positive, on ART and with available viral load result
  onART <- svy_object[svy_object$variables$hivstatusfinal == 1 &
                                svy_object$variables$art == 1 & !is.na(svy_object$variables$viral_sup1000)]
  
  res_1000 <- svyciprop(formula = ~I(viral_sup1000 == "Suppressed"),
                 design = onART, method = "mean",level = 0.95, df= 25)
  
  res_400 <- svyciprop(formula = ~I(viral_sup400 == "Suppressed"),
                        design = onART, method = "mean",level = 0.95, df= 25)
  
  res_200 <- svyciprop(formula = ~I(viral_sup200 == "Suppressed"),
                       design = onART, method = "mean",level = 0.95, df= 25)
  
  res_50 <- svyciprop(formula = ~I(viral_sup50 == "Suppressed"),
                       design = onART, method = "mean",level = 0.95, df= 25)
  
  res_df <- data.frame(country = unique(svy_object$variables$country), 
                       estimate_1000 = res_1000[[1]],
                       lower_ci_1000 = attr(res_1000, "ci")[["2.5%"]],
                       upper_ci_1000 = attr(res_1000, "ci")[["97.5%"]],
                       estimate_400 = res_400[[1]],
                       lower_ci_400 = attr(res_400, "ci")[["2.5%"]], 
                       upper_ci_400 = attr(res_400, "ci")[["97.5%"]],
                       estimate_200 = res_200[[1]],
                       lower_ci_200 = attr(res_200, "ci")[["2.5%"]], 
                       upper_ci_200 = attr(res_200, "ci")[["97.5%"]], 
                       estimate_50 = res_50[[1]], 
                       lower_ci_50 = attr(res_50, "ci")[["2.5%"]],
                       upper_ci_50 = attr(res_50, "ci")[["97.5%"]])
  
  return(res_df)
} 

#' create dataframe with all estimates for all surveys
PHIA_vls_estimates <- purrr::map_df(lapply(svy_object_list, survey_vls_estimates), ~as.data.frame(.x), .id="survey")

# summary stats of survey weighted vls estimates
PHIA_vls_estimates %>% 
  summarise(mean(estimate_1000), median(estimate_1000), 
            quantile(estimate_1000, 1/4), quantile(estimate_1000, 3/4))

