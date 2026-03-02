# install cmdstan
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

#set_cmdstan_path("~/.cmdstan/cmdstan-2.34.1")

# subset to only those on ART > 12 months
# indicator definitions
# 1 - On ART 24 months or more
# 2 - On ART 12-23 months
# 3 - On ART <12 months
# 4 - Not on ART
# 99 - Missing
ssa_surveys_art_experienced <- ssa_surveys_art %>%
  filter(artduration %in% 1:2)

ssa_surveys_art_201519 <- ssa_surveys_art %>%
  filter(survey_year <2020)

# create parameter to censor viral load results not recorded
ssa_surveys_art_experienced <- ssa_surveys_art_experienced %>% 
  mutate(censored = case_when(resultvlc == "< LLOD" ~ "left",
                              resultvlc == "< LLOQ: 40" ~ "left",
                              resultvlc == "< LLOQ: 182" ~ "left",
                              resultvlc == "< LLOQ: 20" ~ "left",
                              resultvlc == "< LLOQ: 80" ~ "left",
                              resultvlc == "TND" ~ "left",
                              resultvlc == "< LLOQ: 400" ~ "left",
                              resultvlc == "< LLOQ: 839" ~ "left",
                              resultvlc == "< LLOQ: 550" ~ "left",
                              resultvlc == "< 20" ~ "left",
                              resultvlc == "< 40" ~ "left",
                              resultvlc == "< 400" ~ "left",
                              resultvlc == "<40" ~ "left",
                              resultvlc == "Target Not Detected" ~ "left",
                              resultvlc == "< Titer min" ~ "left",
                              TRUE ~ "none"))
table(ssa_surveys_art_experienced$censored, useNA = "always")

table(ssa_surveys_art_experienced$country)
table(ssa_surveys_art_experienced$countryname)

# change column titled country to survey as it includes countrynames with survey year to separate surveys
ssa_surveys_art_experienced <- ssa_surveys_art_experienced %>% rename(survey = country)
table(ssa_surveys_art_experienced$survey, useNA = "always")

ssa_surveys_exp_revwb <- ssa_surveys_art_experienced %>% 
  mutate(resultvlc_numeric = case_when(resultvlc_numeric > 999998 ~ 999998,
                                       TRUE ~ resultvlc_numeric))

options (mc.cores=parallel::detectCores())
parallel::detectCores()

table(ssa_surveys_art_experienced$country, useNA = "always")

#' random effetcs weibull model using brms
#' left censored
re_model <- bf(log10(resultvlc_numeric)|cens(censored) ~ 1 + 
                 (1|survey) + (1|countryname))

weibull_model_re <- brm(
  re_model,
  family = weibull(link = "log", link_shape = "identity"),
  #data = ssa_surveys_art_experienced,
  data = ssa_surveys_art,
  warmup = 100,
  iter   = 500,
  chains = 1,
  cores  = 4,
  #backend = "cmdstanr",
  init  = "random",
  prior = set_prior(prior = "uniform(0.43,1.26)", lb = 0.43, ub = 1.26, class = "shape")
)

summary(weibull_model_re) # shape: 0.89 ; 95% CI:0.87 - 0.90


#' random effetcs reverse weibull model using brms
re_model_rev_weib <- bf((6-log10(resultvlc_numeric))|cens(censored) ~ 1 +
                          (1|survey) + (1|countryname))
rev_weib_model_re <- brm(re_model_rev_weib,
                         family = weibull(link = "log", link_shape = "identity"),
                         data = ssa_surveys_exp_revwb,
                         warmup = 100, 
                         iter   = 500, 
                         chains = 1, 
                         init  = "random",
                         #prior = set_prior(prior = "uniform(1.70,3.92)", lb = 1.70, ub = 3.92, class = "shape"), 
                         cores  = 4)
summary(rev_weib_model_re)


# random effects gamma model using brms
#' set initial values
# Initial kicks
init_func <- function(chain_id=1) {
  list ( Intercept  =  0.23 ,
         beta  =  0.3 ,
         shape  =  1.64)
}
# See an example output
init_func(chain_id = 1)
# In this case just four chains are needed.
init_list <- list(
  init_func(chain_id = 1),
  init_func(chain_id = 2),
  init_func(chain_id = 3),
  init_func(chain_id = 4)
)

# for only one chain
init_list <- list(
  init_func(chain_id = 1)
)

gamma_model_re <- brm(re_model,
                      family = Gamma(link = "log"),
                      data = ssa_surveys_art_experienced,
                      warmup = 150, 
                      iter   = 600, 
                      chains = 1,
                      init = init_list,
                      #init  = "random",
                      cores  = 4)


summary(gamma_model_re)

# frechet model
frechet_model_re <- brm(re_model,
                        family = frechet(link = "log"),
                        data = ssa_surveys_art_experienced,
                        warmup = 150, 
                        iter   = 600, 
                        chains = 1,
                        init = init_list,
                        #init  = "random",
                        cores  = 4)

summary(frechet_model_re) # shape: 1.66; 95%CI: 1.62-1.70
# 21 surveys shape: 1.80; 1.77-1.83

# lognormal model
lnorm_model_re <- brm(re_model,
                      family = lognormal(link = "identity", link_sigma = "identity"),
                      data = ssa_surveys_art,
                      warmup = 150, 
                      iter   = 600, 
                      chains = 1,
                      #init = init_list,
                      init  = "random",
                      cores  = 4)

summary(lnorm_model_re) 


# subset those on ART >12 months
library(dplyr)

filtered_list <- dataset_subset %>%
  # keep only dfs with artduration
  purrr::keep(~ "artduration" %in% names(.x)) %>%
  # filter rows
  purrr::map(~ filter(.x, artduration %in% c(1, 2)))

#' calculate survey weighted VLS estimates for each survey at the defined thresholds

#' create function to create survey objects for all surveys
options(survey.lonely.psu = "adjust")

create_svy_object <- function(dataset){
  svy_object <- svydesign(ids = ~varunit, weights = ~btwt0, strata = ~varstrat,
                          nest = TRUE, data = dataset)
  return(svy_object)
}

svy_object_list <- lapply(X = filtered_list, FUN = create_svy_object)

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
#PHIA_vls_estimates <- do.call(rbind, lapply(svy_object_list, survey_vls_estimates))
PHIA_vls_estimates_experienced <- purrr::map_df(lapply(svy_object_list, survey_vls_estimates), ~as.data.frame(.x), .id="survey")

# 21 surveys all PLHIV
PHIA_vls_estimates %>% 
  #filter(!grepl("2020$|2021$", survey)) %>%
  summarise(mean(estimate_1000), median(estimate_1000), 
            quantile(estimate_1000, 1/4), quantile(estimate_1000, 3/4))

# 14 surveys PLHIV on ART >12 months
PHIA_vls_estimates_experienced %>% 
  summarise(mean(estimate_1000), median(estimate_1000), 
            quantile(estimate_1000, 1/4), quantile(estimate_1000, 3/4))

# 14 surveys all PLHIV
PHIA_vls_estimates %>% 
  filter(!grepl("2020$|2021$", survey)) %>%
  summarise(mean(estimate_1000), median(estimate_1000), 
            quantile(estimate_1000, 1/4), quantile(estimate_1000, 3/4))

# review sample on ART >12 months
ssa_surveys_art %>% group_by(survey, artduration) %>% summarise(n()) %>% print(n =Inf)

# review numbers of PLHIV on ART with VL <50, <200, <400 and <1000
ssa_surveys_art %>% 
  group_by(country) %>%
  summarise(less_50 = sum(resultvlc_numeric < 50), 
            less_200 = sum(resultvlc_numeric < 200),
            less_400 = sum(resultvlc_numeric < 400),
            less1000 = sum(resultvlc_numeric < 1000)) %>%
  print(n = Inf) %>%
  write.csv(file = "C:/Users/oedun/OneDrive - Imperial College London/R projects/viremia_thresholds/jaids_supp_data.csv")

# JAIDS review
PHIA_vls_estimates %>%
  mutate(vls_less50 = sprintf("%.1f (%.1f, %.1f)", estimate_50*100, lower_ci_50*100, upper_ci_50*100),
         vls_less200 = sprintf("%.1f (%.1f, %.1f)", estimate_200*100, lower_ci_200*100, upper_ci_200*100),
         vls_less400 = sprintf("%.1f (%.1f, %.1f)", estimate_400*100, lower_ci_400*100, upper_ci_400*100)) %>%
  select(survey, country, vls_less50, vls_less200, vls_less400) %>%
  arrange(country) %>%
  write.csv(file = "C:/Users/oedun/OneDrive - Imperial College London/R projects/viremia_thresholds/jaids_vls_est.csv")


# temporal hold out regressions
# create parameter to censor viral load results not recorded
ssa_surveys_art_201519 <- ssa_surveys_art_201519 %>% 
  mutate(censored = case_when(resultvlc == "< LLOD" ~ "left",
                              resultvlc == "< LLOQ: 40" ~ "left",
                              resultvlc == "< LLOQ: 182" ~ "left",
                              resultvlc == "< LLOQ: 20" ~ "left",
                              resultvlc == "< LLOQ: 80" ~ "left",
                              resultvlc == "TND" ~ "left",
                              resultvlc == "< LLOQ: 400" ~ "left",
                              resultvlc == "< LLOQ: 839" ~ "left",
                              resultvlc == "< LLOQ: 550" ~ "left",
                              resultvlc == "< 20" ~ "left",
                              resultvlc == "< 40" ~ "left",
                              resultvlc == "< 400" ~ "left",
                              resultvlc == "<40" ~ "left",
                              resultvlc == "Target Not Detected" ~ "left",
                              resultvlc == "< Titer min" ~ "left",
                              TRUE ~ "none"))
table(ssa_surveys_art_201519$censored, useNA = "always")

table(ssa_surveys_art_201519$country)
table(ssa_surveys_art_201519$countryname)

# change column titled country to survey as it includes countrynames with survey year to separate surveys
ssa_surveys_art_201519 <- ssa_surveys_art_201519 %>% rename(survey = country)
table(ssa_surveys_art_201519$survey, useNA = "always")

ssa_surveys_201519_revwb <- ssa_surveys_art_201519 %>% 
  mutate(resultvlc_numeric = case_when(resultvlc_numeric > 999998 ~ 999998,
                                       TRUE ~ resultvlc_numeric))

options (mc.cores=parallel::detectCores())
parallel::detectCores()

table(ssa_surveys_201519_revwb$country, useNA = "always")

#' random effetcs weibull model using brms
#' left censored
re_model <- bf(log10(resultvlc_numeric)|cens(censored) ~ 1 + 
                 (1|survey) + (1|countryname))

weibull_model_re <- brm(
  re_model,
  family = weibull(link = "log", link_shape = "identity"),
  data = ssa_surveys_art_201519,
  warmup = 100,
  iter   = 500,
  chains = 1,
  cores  = 4,
  #backend = "cmdstanr",
  init  = "random",
  prior = set_prior(prior = "uniform(0.43,1.26)", lb = 0.43, ub = 1.26, class = "shape")
)

summary(weibull_model_re) # shape: 0.89 ; 95% CI:0.87 - 0.90


#' random effetcs reverse weibull model using brms
re_model_rev_weib <- bf((6-log10(resultvlc_numeric))|cens(censored) ~ 1 +
                          (1|survey) + (1|countryname))
rev_weib_model_re <- brm(re_model_rev_weib,
                         family = weibull(link = "log", link_shape = "identity"),
                         data = ssa_surveys_201519_revwb,
                         warmup = 100, 
                         iter   = 500, 
                         chains = 1, 
                         init  = "random",
                         #prior = set_prior(prior = "uniform(1.70,3.92)", lb = 1.70, ub = 3.92, class = "shape"), 
                         cores  = 4)
summary(rev_weib_model_re)


# random effects gamma model using brms
#' set initial values
# Initial kicks
init_func <- function(chain_id=1) {
  list ( Intercept  =  0.23 ,
         beta  =  0.3 ,
         shape  =  1.64)
}
# See an example output
init_func(chain_id = 1)
# In this case just four chains are needed.
init_list <- list(
  init_func(chain_id = 1),
  init_func(chain_id = 2),
  init_func(chain_id = 3),
  init_func(chain_id = 4)
)

# for only one chain
init_list <- list(
  init_func(chain_id = 1)
)

gamma_model_re <- brm(re_model,
                      family = Gamma(link = "log"),
                      data = ssa_surveys_art_201519,
                      warmup = 200, 
                      iter   = 1000, 
                      chains = 1,
                      init = init_list,
                      #init  = "random",
                      cores  = 4)


summary(gamma_model_re)

# frechet model
frechet_model_re <- brm(re_model,
                        family = frechet(link = "log"),
                        data = ssa_surveys_art_201519,
                        warmup = 150, 
                        iter   = 600, 
                        chains = 1,
                        init = init_list,
                        #init  = "random",
                        cores  = 4)

summary(frechet_model_re) # shape: 1.66; 95%CI: 1.62-1.70
# 21 surveys shape: 1.80; 1.77-1.83

# lognormal model
lnorm_model_re <- brm(re_model,
                      family = lognormal(link = "identity", link_sigma = "identity"),
                      data = ssa_surveys_art_201519,
                      warmup = 300, 
                      iter   = 1500, 
                      chains = 1,
                      #init = init_list,
                      init  = "random",
                      cores  = 4)


summary(lnorm_model_re) 
