# create parameter to censor viral load results reported at lower limit of detection
ssa_surveys_art <- ssa_surveys_art %>% 
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
table(ssa_surveys_art$censored, useNA = "always")

table(ssa_surveys_art$country)
table(ssa_surveys_art$countryname)

# change column titled country to survey as it includes countrynames with survey year to separate surveys
ssa_surveys_art <- ssa_surveys_art %>% rename(survey = country)
table(ssa_surveys_art$survey, useNA = "always")

# for reverse Weibull fitting set VL values greater than or equal to 6 log10 to 5.99 log10 to enable fitting
ssa_surveys_art_revwb <- ssa_surveys_art %>% mutate(resultvlc_numeric = case_when(resultvlc_numeric > 999998 ~ 999998,
                                                                                  TRUE ~ resultvlc_numeric))

options (mc.cores=parallel::detectCores())
parallel::detectCores()

#' random effetcs weibull model using brms
#' left censored
re_model <- bf(log10(resultvlc_numeric)|cens(censored) ~ 1 + 
                 (1|survey) + (1|countryname))

weibull_model_re <- brm(re_model,
                        family = weibull(link = "log", link_shape = "identity"),
                        data = ssa_surveys_art,
                        warmup = 1500, 
                        iter   = 4000, 
                        chains = 4, 
                        init  = "random",
                        cores  = 4)

summary(weibull_model_re) 

#' random effetcs reverse weibull model using brms
re_model_rev_weib <- bf((6 - log10(resultvlc_numeric))|cens(censored) ~ 1 +
                          (1|survey) + (1|countryname))
rev_weib_model_re <- brm(re_model_rev_weib,
                        family = weibull(link = "log", link_shape = "identity"),
                        data = ssa_surveys_art_revwb,
                        warmup = 1500, 
                        iter   = 4000, 
                        chains = 4, 
                        init  = "random",
                        cores  = 4)
summary(rev_weib_model_re)

gamma_model_re <- brm(re_model,
                      family = Gamma(link = "log"),
                      data = ssa_surveys_art,
                      warmup = 1500, 
                      iter   = 4000, 
                      chains = 4,
                      init  = "random",
                      cores  = 4)

summary(gamma_model_re) 

# frechet model
frechet_model_re <- brm(re_model,
                        family = frechet(link = "log"),
                        data = ssa_surveys_art,
                        warmup = 1500, 
                        iter   = 4000, 
                        chains = 4,
                        init  = "random",
                        cores  = 4)

summary(frechet_model_re) 

# lognormal model
lnorm_model_re <- brm(re_model,
                        family = lognormal(link = "identity", link_sigma = "identity"),
                        data = ssa_surveys_art,
                        warmup = 1500, 
                        iter   = 4000, 
                        chains = 4,
                        init  = "random",
                        cores  = 4)

summary(lnorm_model_re) 