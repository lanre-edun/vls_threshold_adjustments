#' load packages
library(haven) 
library(dplyr)
library(ggplot2)
library(survey)
library(brms)
library(readxl)
library(scales)

#' read in datasets
#' prepare file names for all datasets used
 filenames <- c("STATA/MPHIA/MPHIA2015adultbio.dta",
                "STATA/ZIMPHIA/ZIMPHIA2015adultbio.dta",
                "STATA/LePHIA/LePHIA2016adultbio.dta",
                "STATA/SHIMS/SHIMS22016adultbio.dta",
                "STATA/THIS/THIS2016adultbio.dta",
                "STATA/UPHIA/UPHIA2016adultbio.dta",
                "STATA/ZAMPHIA/ZAMPHIA2016adultbio.dta",
                "STATA/CAMPHIA/CAMPHIA2017adultbio.dta",
                "STATA/CIPHIA/CIPHIA2017adultbio.dta",
                "STATA/EPHIA/EPHIA2017adultbio.dta",
                "STATA/NAMPHIA/NAMPHIA2017adultbio.dta",
                "STATA/NAIIS/NAIIS2018adultbio.dta",
                "STATA/RPHIA/RPHIA2018adultbio.dta",
                "STATA/KENPHIA/KENPHIA2018adultbio.dta",
                "STATA/BAIS/BAISV2021adultbio.dta",
                "STATA/INSIDA/INSIDA2021adultbio.dta",
                "STATA/LePHIA/LePHIA2020adultbio.dta",
                "STATA/MPHIA/MPHIA2020adultbio.dta",
                "STATA/SHIMS/SHIMS32021adultbio.dta",
                "STATA/ZIMPHIA/ZIMPHIA2020adultbio.dta",
                "STATA/ZAMPHIA/ZAMPHIA2021adultbio.dta")

#' prepare names for list using survey name
listnames <- c("MPHIA_2015", "ZIMPHIA_2015", "LePHIA_2016", "SHIMS2_2016",
               "THIS_2016", "UPHIA_2016", "ZAMPHIA_2016", "CAMPHIA_2017",
               "CIPHIA_2017", "EPHIA_2017", "NAMPHIA_2017", "NAIIS_2018",
               "RPHIA_2018", "KENPHIA_2018", "BAIS_2021", "INSIDA_2021",
               "LePHIA_2020", "MPHIA_2020", "SHIMS3_2021", "ZIMPHIA_2020",
               "ZAMPHIA_2021")

#' function to read in data set
#' @details it takes a file name and reads in the data set
#' @param filename which is a string of texts with the file names and location
read_dataset <- function(filename){
  #' read in dataset
  dataset <- read_dta(file = filename)
  
  #' select only participants with definite blood test result
  #' 1 = lab test has definite results and defacto participant
  dataset_eligible <- dataset %>% filter(bt_status == 1)
  
  return(data.frame(dataset_eligible))
}

#' create a list with all survey datasets
datasets_list <- lapply(X = filenames, FUN = read_dataset)

#' name the list with survey names
names(datasets_list) <- listnames
names(datasets_list)

#' ZAMPHIA 2016 lacks an art variable, it needs to be created also using
#' aware, arvstatus and artselfreported variables
datasets_list[["ZAMPHIA_2016"]]$art <-
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 1, 1,                                   # on ART
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 2 & datasets_list[["ZAMPHIA_2016"]]$aware == 1 &  datasets_list[["ZAMPHIA_2016"]]$artselfreported == 1, 1,
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 2 & datasets_list[["ZAMPHIA_2016"]]$aware == 99 & datasets_list[["ZAMPHIA_2016"]]$artselfreported == 1, 1,
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 99 & datasets_list[["ZAMPHIA_2016"]]$aware == 1 & datasets_list[["ZAMPHIA_2016"]]$artselfreported == 1, 1,
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 99 & datasets_list[["ZAMPHIA_2016"]]$aware == 99 & datasets_list[["ZAMPHIA_2016"]]$artselfreported == 1, 1,
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 2 & datasets_list[["ZAMPHIA_2016"]]$aware == 2, 2, # not on ART
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 99 & datasets_list[["ZAMPHIA_2016"]]$aware == 2, 2,
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 2 & datasets_list[["ZAMPHIA_2016"]]$aware == 1 & datasets_list[["ZAMPHIA_2016"]]$artselfreported == 2, 2,
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 99 & datasets_list[["ZAMPHIA_2016"]]$aware == 1 & datasets_list[["ZAMPHIA_2016"]]$artselfreported == 2, 2,
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 2 & datasets_list[["ZAMPHIA_2016"]]$aware == 99 & datasets_list[["ZAMPHIA_2016"]]$artselfreported == 2, 2,
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 99 & datasets_list[["ZAMPHIA_2016"]]$aware == 99 & datasets_list[["ZAMPHIA_2016"]]$artselfreported == 2, 2,
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 2 & datasets_list[["ZAMPHIA_2016"]]$aware == 1 & datasets_list[["ZAMPHIA_2016"]]$artselfreported == 99, 99, 
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 99 & datasets_list[["ZAMPHIA_2016"]]$aware == 1 & datasets_list[["ZAMPHIA_2016"]]$artselfreported == 99, 99,   ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 2 & datasets_list[["ZAMPHIA_2016"]]$aware == 99 & datasets_list[["ZAMPHIA_2016"]]$artselfreported == 99, 99,
  ifelse(datasets_list[["ZAMPHIA_2016"]]$arvstatus == 99 & datasets_list[["ZAMPHIA_2016"]]$aware == 99 & datasets_list[["ZAMPHIA_2016"]]$artselfreported == 99, 99, NA)))))))))))))))

table(datasets_list[["ZAMPHIA_2016"]]$art, useNA = "always")

#' exclude 3 participants with missing blood test weights 
datasets_list[["NAIIS_2018"]] <- datasets_list[["NAIIS_2018"]] %>% filter(!is.na(btwt0))

#' exclude missing strata for THIS (Tanzania  PHIA) survey
datasets_list[["THIS_2016"]] <- datasets_list[["THIS_2016"]] %>% filter(varstrat != 36) %>%
  filter(varstrat != 37) %>% filter(varstrat != 58)

#' add a region_id variable for second level of random effect
datasets_list[["MPHIA_2015"]]$survey_year <- 2015
datasets_list[["ZIMPHIA_2015"]]$survey_year <- 2015
datasets_list[["LePHIA_2016"]]$survey_year <- 2016
datasets_list[["SHIMS2_2016"]]$survey_year <- 2016
datasets_list[["THIS_2016"]]$survey_year <- 2016
datasets_list[["UPHIA_2016"]]$survey_year <- 2016
datasets_list[["ZAMPHIA_2016"]]$survey_year <- 2016
datasets_list[["CAMPHIA_2017"]]$survey_year <- 2017
datasets_list[["CIPHIA_2017"]]$survey_year <- 2017
datasets_list[["EPHIA_2017"]]$survey_year <- 2017
datasets_list[["NAMPHIA_2017"]]$survey_year <- 2017
datasets_list[["NAIIS_2018"]]$survey_year <- 2018
datasets_list[["RPHIA_2018"]]$survey_year <- 2018
datasets_list[["KENPHIA_2018"]]$survey_year <- 2018
datasets_list[["BAIS_2021"]]$survey_year <- 2021
datasets_list[["INSIDA_2021"]]$survey_year <- 2021
datasets_list[["LePHIA_2020"]]$survey_year <- 2020
datasets_list[["MPHIA_2020"]]$survey_year <- 2020
datasets_list[["SHIMS3_2021"]]$survey_year <- 2021
datasets_list[["ZIMPHIA_2020"]]$survey_year <- 2020
datasets_list[["ZAMPHIA_2021"]]$survey_year <- 2021

table(datasets_list[["ZAMPHIA_2021"]]$resultvlc, useNA = "always")

#' create a function to convert viral load results to all numeric values
#' also create viral suppression variable based on different cut-offs
#' @details function converts character viral load results to numeric values
#'  then creates viral suppression variable using different thresholds
#'  function also specifies the lower numeric value for surveys reporting lower limits of detection
#'  lower limits varies by surveys
#' @param receives dataset which is survey dataset
numeric_vl_sup <- function(dataset){

    dataset$resultvlc_numeric <- ifelse(dataset$resultvlc == "< LLOD" & dataset$country %in% c("Malawi","Cameroon") , 40, 
                             ifelse(dataset$resultvlc == "< LLOD" & !dataset$country %in% c("Malawi","Cameroon") , 20,          
                             ifelse(dataset$resultvlc == "< LLOQ: 40", 40,
                             ifelse(dataset$resultvlc == "< LLOQ: 182", 182,
                             ifelse(dataset$resultvlc == "< LLOQ: 20", 20,
                             ifelse(dataset$resultvlc == "< LLOQ: 80", 80,
                             ifelse(dataset$resultvlc == "> ULOQ: 10000000", 999998,
                             ifelse(dataset$resultvlc == "TND" & dataset$country %in% c("Malawi","Cameroon"), 40,
                             ifelse(dataset$resultvlc == "TND" & !dataset$country %in% c("Malawi","Cameroon"), 20,
                             ifelse(dataset$resultvlc == "< LLOQ: 400", 400,
                             ifelse(dataset$resultvlc == "< LLOQ: 839", 839,
                             ifelse(dataset$resultvlc == "< LLOQ: 550", 550,
                             ifelse(dataset$resultvlc == "Failed", NA,
                             ifelse(dataset$resultvlc == "< 20", 20,
                             ifelse(dataset$resultvlc == "< 40", 40,
                             ifelse(dataset$resultvlc == "< 400", 400,
                             ifelse(dataset$resultvlc == "<40", 40,
                             ifelse(dataset$resultvlc == "Target Not Detected", 20,
                             ifelse(dataset$resultvlc == "< Titer min", 20,
                             ifelse(dataset$resultvlc == "Rejected", NA,
                             ifelse(dataset$resultvlc == "Result pending", NA, dataset$resultvlc)))))))))))))))))))))
  
  dataset$resultvlc_numeric <- as.numeric(dataset$resultvlc_numeric)
  
  dataset$viral_sup1000 <- ifelse(dataset$resultvlc_numeric < 1000, "Suppressed",
                         ifelse(dataset$resultvlc_numeric > 999, "Non-suppressed", NA))
  
  dataset$viral_sup400 <- ifelse(dataset$resultvlc_numeric < 400, "Suppressed",
                        ifelse(dataset$resultvlc_numeric > 399, "Non-suppressed", NA))
  
  dataset$viral_sup200 <- ifelse(dataset$resultvlc_numeric < 200, "Suppressed",
                        ifelse(dataset$resultvlc_numeric > 199, "Non-suppressed", NA))
  
  dataset$viral_sup50 <- ifelse(dataset$resultvlc_numeric < 50, "Suppressed",
                       ifelse(dataset$resultvlc_numeric > 49, "Non-suppressed", NA))
  
  # Define the columns you want
  cols <- c("country", "gender", "age", "resultvlc_numeric", "varstrat", "varunit", "btwt0",
            "hivstatusfinal", "art", "survey_year", "resultvlc", "arvstatus", "artselfreported",
            "artduration", "viral_sup1000", "viral_sup400", "viral_sup200", "viral_sup50")
  
  # Select only those that exist in the dataset
  dataset <- dataset %>%
    dplyr::select(dplyr::any_of(cols))
 
  dataset$resultvlc_numeric <- as.numeric(dataset$resultvlc_numeric)
  return(data.frame(dataset))
}

#' create a list with all survey datasets including the viral suppression variable and numeric viral load results
#' this is done by applying the function above to the survey datasets in the list using the lapply function
dataset_subset <- lapply(X = datasets_list, FUN = numeric_vl_sup)

#' unlist data set 
list2env(dataset_subset, .GlobalEnv)

#' assign country name and year of survey only for plotting purpose
MPHIA_2015$country <- rep("Malawi (2015-16)", dim(MPHIA_2015)[1])
ZIMPHIA_2015$country <- rep("Zimbabwe (2015-16)", dim(ZIMPHIA_2015)[1])
LePHIA_2016$country <- rep("Lesotho (2016-17)", dim(LePHIA_2016)[1])
SHIMS2_2016$country <- rep("Eswatini (2016-17)", dim(SHIMS2_2016)[1])
THIS_2016$country <- rep("Tanzania (2016-17)", dim(THIS_2016)[1])
UPHIA_2016$country <- rep("Uganda (2016-17)", dim(UPHIA_2016)[1])
ZAMPHIA_2016$country <- rep("Zambia (2016)", dim(ZAMPHIA_2016)[1])
CAMPHIA_2017$country <- rep("Cameroon (2017-18)", dim(CAMPHIA_2017)[1])
CIPHIA_2017$country <- rep("Côte d'Ivoire (2017-18)", dim(CIPHIA_2017)[1])
EPHIA_2017$country <- rep("Ethiopia (2017-18)", dim(EPHIA_2017)[1])
NAMPHIA_2017$country <- rep("Namibia (2017)", dim(NAMPHIA_2017)[1])
NAIIS_2018$country <- rep("Nigeria (2018)", dim(NAIIS_2018)[1])
RPHIA_2018$country <- rep("Rwanda (2018-19)", dim(RPHIA_2018)[1])
KENPHIA_2018$country <- rep("Kenya (2018-19)", dim(KENPHIA_2018)[1])
BAIS_2021$country <- rep("Botswana (2021)", dim(BAIS_2021)[1])
INSIDA_2021$country <- rep("Mozambique (2021-22)", dim(INSIDA_2021)[1])
LePHIA_2020$country <- rep("Lesotho (2019-20)", dim(LePHIA_2020)[1])
MPHIA_2020$country <- rep("Malawi (2020-21)", dim(MPHIA_2020)[1])
SHIMS3_2021$country <- rep("Eswatini (2021)", dim(SHIMS3_2021)[1])
ZIMPHIA_2020$country <- rep("Zimbabwe (2019-20)", dim(ZIMPHIA_2020)[1])
ZAMPHIA_2021$country <- rep("Zambia (2021)", dim(ZAMPHIA_2021)[1])

# create normalized country weights
CAMPHIA_2017$norm_wt <- CAMPHIA_2017$btwt0/mean(CAMPHIA_2017$btwt0)
CIPHIA_2017$norm_wt <- CIPHIA_2017$btwt0/mean(CIPHIA_2017$btwt0)
NAIIS_2018$norm_wt <- NAIIS_2018$btwt0/mean(NAIIS_2018$btwt0)
EPHIA_2017$norm_wt <- EPHIA_2017$btwt0/mean(EPHIA_2017$btwt0)
MPHIA_2015$norm_wt <- MPHIA_2015$btwt0/mean(MPHIA_2015$btwt0)
RPHIA_2018$norm_wt <- RPHIA_2018$btwt0/mean(RPHIA_2018$btwt0)
KENPHIA_2018$norm_wt <- KENPHIA_2018$btwt0/mean(KENPHIA_2018$btwt0)
ZAMPHIA_2016$norm_wt <- ZAMPHIA_2016$btwt0/mean(ZAMPHIA_2016$btwt0)
THIS_2016$norm_wt <- THIS_2016$btwt0/mean(THIS_2016$btwt0)
UPHIA_2016$norm_wt <- UPHIA_2016$btwt0/mean(UPHIA_2016$btwt0)
ZIMPHIA_2015$norm_wt <- ZIMPHIA_2015$btwt0/mean(ZIMPHIA_2015$btwt0)
LePHIA_2016$norm_wt <- LePHIA_2016$btwt0/mean(LePHIA_2016$btwt0)
SHIMS2_2016$norm_wt <- SHIMS2_2016$btwt0/mean(SHIMS2_2016$btwt0)
NAMPHIA_2017$norm_wt <- NAMPHIA_2017$btwt0/mean(NAMPHIA_2017$btwt0)
BAIS_2021$norm_wt <- BAIS_2021$btwt0/mean(BAIS_2021$btwt0)
INSIDA_2021$norm_wt <- INSIDA_2021$btwt0/mean(INSIDA_2021$btwt0)
LePHIA_2020$norm_wt <- LePHIA_2020$btwt0/mean(LePHIA_2020$btwt0)
MPHIA_2020$norm_wt <- MPHIA_2020$btwt0/mean(MPHIA_2020$btwt0)
SHIMS3_2021$norm_wt <- SHIMS3_2021$btwt0/mean(SHIMS3_2021$btwt0)
ZIMPHIA_2020$norm_wt <- ZIMPHIA_2020$btwt0/mean(ZIMPHIA_2020$btwt0)
ZAMPHIA_2021$norm_wt <- ZAMPHIA_2021$btwt0/mean(ZAMPHIA_2021$btwt0)

#' assign country name alone
MPHIA_2015$countryname <- "Malawi"
ZIMPHIA_2015$countryname <- "Zimbabwe"
LePHIA_2016$countryname <- "Lesotho"
SHIMS2_2016$countryname <- "Eswatini"
THIS_2016$countryname <- "Tanzania"
UPHIA_2016$countryname <- "Uganda"
ZAMPHIA_2016$countryname <- "Zambia"
CAMPHIA_2017$countryname <- "Cameroon"
CIPHIA_2017$countryname <- "Côte d'Ivoire"
EPHIA_2017$countryname <- "Ethiopia"
NAMPHIA_2017$countryname <- "Namibia"
NAIIS_2018$countryname <- "Nigeria"
RPHIA_2018$countryname <- "Rwanda"
KENPHIA_2018$countryname <- "Kenya"
BAIS_2021$countryname <- "Botswana"
INSIDA_2021$countryname <- "Mozambique"
LePHIA_2020$countryname <- "Lesotho"
MPHIA_2020$countryname <- "Malawi"
SHIMS3_2021$countryname <- "Eswatini"
ZIMPHIA_2020$countryname <- "Zimbabwe"
ZAMPHIA_2021$countryname <- "Zambia"

# merge surveys
ssa_surveys <- bind_rows(MPHIA_2015, ZIMPHIA_2015, LePHIA_2016, SHIMS2_2016,
                     THIS_2016, UPHIA_2016, ZAMPHIA_2016, CAMPHIA_2017,
                     CIPHIA_2017, EPHIA_2017, NAMPHIA_2017, NAIIS_2018,
                     RPHIA_2018, KENPHIA_2018, BAIS_2021, INSIDA_2021,
                     LePHIA_2020, MPHIA_2020, SHIMS3_2021, ZIMPHIA_2020,
                     ZAMPHIA_2021)

#' subset to those on ART
ssa_surveys_art <- ssa_surveys %>% filter(art == 1)

table(ssa_surveys_art$resultvlc_numeric, useNA = "always")
# 60 NAs (missing VLS)
tail(table(ssa_surveys_art$resultvlc_numeric, useNA = "always"))

table(ssa_surveys_art$country, useNA = "always")

ssa_surveys_art <- ssa_surveys_art %>% filter(!is.na(resultvlc_numeric))

ssa_surveys_art %>% filter(log10(resultvlc_numeric) > 6) %>% select(resultvlc_numeric, countryname)
ssa_surveys_art %>% filter(resultvlc == "Target Not Detected") %>% distinct(countryname)
  
table(ssa_surveys_art$artduration, useNA = "always")

# indicator definitions for duration on ART variable in PHIAs
# 1 - On ART 24 months or more
# 2 - On ART 12-23 months
# 3 - On ART <12 months
# 4 - Not on ART
# 99 - Missing

table(ssa_surveys_art$country, useNA = "always")

ssa_surveys_art %>%
  group_by(country, artduration) %>%
  summarise(n()) %>% print(n = Inf)

