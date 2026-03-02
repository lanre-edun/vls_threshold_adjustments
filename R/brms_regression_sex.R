# define age categories
ssa_surveys_art <- ssa_surveys_art %>% 
  mutate(agecat = case_when(age < 25 ~ "15-24",
                            age < 35 ~ "25-34",
                            age < 45 ~ "35-44",
                            age < 55 ~ "45-54",
                            TRUE ~ "55+"))

table(ssa_surveys_art$agecat, ssa_surveys_art$age, useNA = "always")
options (mc.cores=parallel::detectCores ())

# AGE #####
# parameters by sex and ages
shape_est <- function(new_est, ll, ul, intercept_shape){
  est <- intercept_shape + new_est
  l_est <- intercept_shape + ll
  u_est <- intercept_shape + ul
  return(data.frame(est = est, l_est = l_est, u_est = u_est))
}

# explore different shapes by age
dist_weib_age <- bf(log10(resultvlc_numeric)|cens(censored) ~ 1 + 
                      (1|survey) + (1|countryname) + (1|agecat), shape ~ agecat)
weibull_model_age <- brm(dist_weib_age,
                         family = weibull(link = "log", link_shape = "identity"),
                         data = ssa_surveys_art,
                         warmup = 1500, 
                         iter   = 4000, 
                         chains = 4, 
                         init  = "random",
                         control=list(adapt_delta=0.9999, max_treedepth=100),
                         cores  = 4)

summary(weibull_model_age)  


#' reverse weibull model
dist_rwb_age <- bf((6-log10(resultvlc_numeric))|cens(censored) ~ 1 +
                     (1|survey) + (1|countryname) + (1|agecat), shape ~ agecat)

rev_weib_model_age <- brm(dist_rwb_age,
                          family = weibull(link = "log", link_shape = "identity"),
                          data = ssa_surveys_art,
                          warmup = 1500, 
                          iter   = 4000, 
                          chains = 4, 
                          init  = "random",
                          control=list(adapt_delta=0.9999, max_treedepth=100),  
                          cores  = 4)
summary(rev_weib_model_age)  

#' gamma model
dist_gamma_age <- bf(log10(resultvlc_numeric)|cens(censored) ~ 1 + 
                       (1|survey) + (1|countryname) + (1|agecat), shape ~ agecat)

gamma_model_age <- brm(dist_gamma_age,
                       family = Gamma(link = identity),
                       data = ssa_surveys_art,
                       warmup = 1500, 
                       iter   = 4000, 
                       chains = 4, 
                       init  = "random",
                       control=list(adapt_delta=0.9999, max_treedepth=100),
                       cores  = 4)
summary(gamma_model_age)  

#' frechet model
dist_frechet_age <- bf(log10(resultvlc_numeric)|cens(censored) ~ 1 + 
                         (1|survey) + (1|countryname) + (1|agecat), nu ~ agecat)

frechet_model_age <- brm(dist_frechet_age,
                         family = frechet(link = "log", link_nu = "identity"),
                         data = ssa_surveys_art,
                         warmup = 1500, 
                         iter   = 4000, 
                         chains = 4,
                         init  = "random",
                         cores  = 4,
                         control=list(adapt_delta=0.9999, max_treedepth=100))

summary(frechet_model_age)  

#' lognormal model
dist_lnorm_age <- bf(log10(resultvlc_numeric)|cens(censored) ~ 1 + 
                         (1|survey) + (1|countryname) + (1|agecat), sigma ~ agecat)

lnorm_model_age <- brm(dist_lnorm_age,
                       family = lognormal(link = "identity", link_sigma = "identity"),
                         data = ssa_surveys_art,
                         warmup = 1500, 
                         iter   = 4000, 
                         chains = 4, 
                         init  = "random",
                         cores  = 4,
                         control=list(adapt_delta=0.9999, max_treedepth=100))

summary(lnorm_model_age)  

# Sex
# parameters by sex 
ssa_surveys_art <- ssa_surveys_art %>% mutate(sex = case_when(gender == 1 ~ "men",
                                                              TRUE ~ "women"))

ssa_surveys_art_revwb <- ssa_surveys_art_revwb %>% mutate(sex = case_when(gender == 1 ~ "men",
                                                                          TRUE ~ "women"))

table(ssa_surveys_art$gender, ssa_surveys_art$sex)
table(ssa_surveys_art_revwb$gender, ssa_surveys_art_revwb$sex)

dist_weib_sex <- bf(log10(resultvlc_numeric)|cens(censored) ~ 1 + 
                      (1|survey) + (1|countryname) + (1|sex), shape ~ sex)
weibull_model_sex <- brm(dist_weib_sex,
                         family = weibull(link = "log", link_shape = "identity"),
                         data = ssa_surveys_art,
                         warmup = 1500, 
                         iter   = 4000, 
                         chains = 4, 
                         init  = "random",
                         control=list(adapt_delta=0.9999, max_treedepth=100), 
                         cores  = 4)
summary(weibull_model_sex)  

#' reverse weibull model
dist_rwb_sex <- bf((6-log10(resultvlc_numeric))|cens(censored) ~ 1 +
                     (1|survey) + (1|countryname) + (1|sex), shape ~ sex)

rev_weib_model_sex <- brm(dist_rwb_sex,
                          family = weibull(link = "log", link_shape = "identity"),
                          data = ssa_surveys_art_revwb,
                          warmup = 1500, 
                          iter   = 4000, 
                          chains = 4, 
                          init  = "random",
                          control=list(adapt_delta=0.9999, max_treedepth=100),
                          cores  = 4)
summary(rev_weib_model_sex)  

#' gamma model
dist_gamma_sex <- bf(log10(resultvlc_numeric)|cens(censored) ~ 1 + 
                       (1|survey) + (1|countryname) + (1|sex), shape ~ sex)

gamma_model_sex <- brm(dist_gamma_sex,
                       family = Gamma(link = identity),
                       data = ssa_surveys_art,
                       warmup = 1500, 
                       iter   = 4000, 
                       chains = 4, 
                       init  = "random",
                       control=list(adapt_delta=0.9999, max_treedepth=100),
                       cores  = 4)
summary(gamma_model_sex)  

#' frechet model
dist_frechet_sex <- bf(log10(resultvlc_numeric)|cens(censored) ~ 1 + 
                       (1|survey) + (1|countryname)  + (1|sex), nu ~ sex)

frechet_model_sex <- brm(dist_frechet_sex,
                       family = frechet(link = "log", link_nu = "identity"),
                       data = ssa_surveys_art,
                       warmup = 400, 
                       iter   = 1500, 
                       chains = 4, 
                       init  = "random",
                       cores  = 4,
                       control=list(adapt_delta=0.9999, max_treedepth=100))

summary(frechet_model_sex)  

#' lognormal model
dist_lnorm_sex <- bf(log10(resultvlc_numeric)|cens(censored) ~ 1 + 
                         (1|survey) + (1|countryname) + (1|sex), sigma ~ sex)

lnorm_model_sex <- brm(dist_lnorm_sex,
                         family = lognormal(link = "identity", link_sigma = "identity"),
                         data = ssa_surveys_art,
                         warmup = 1500, 
                         iter   = 4000, 
                         chains = 4, 
                         init  = "random",
                         cores  = 4,
                         control=list(adapt_delta=0.9999, max_treedepth=100))

summary(lnorm_model_sex)  
save(lnorm_model_sex = lnorm_model_sex,
     file = "C:/Users/oae19/OneDrive - Imperial College London/R projects/viremia_thresholds/regression-outputs/lnorm_model_sex")


