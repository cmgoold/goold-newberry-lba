############################################################################
# R analysis script for Goold & Newberry (2021):
#     - Predicting individual shelter dog behaviour after adoption using 
#       longitudinal behavioural assessment: a hierarchical Bayesian approach

# Copyright(C) Conor Goold  2021
# goold.conor@gmail.com
############################################################################

path <- "~/Dropbox/PhD/PhD_NMBU/PaperIV/goold-newberry-lba/"
figure_path <- paste0(path, "paper/figures/")
source(paste0(path, "analysis-code/helper_functions.R"))
setwd(path)

# install packages required
get_needed_packages()

# load the data
d_s <- read.table(paste0(path, "data/d_shelter_behaviour.txt"), header = T)
d_a <- read.table(paste0(path, "data/d_post_adoption_behaviour.txt"), header=T)
d_dem <- read.csv(paste0(path, "data/d_demographic_details.txt"), header = T, na.strings = "")

# // contexts:
#       1 = Interactions with dogs
#       2 = Sociability with unfamiliar people
#       3 = Sociability with familiar people
#       4 = Handling
#       5 = In kennel
#       6 = Outside kennel
#       7 = Interactions with toys
#       8 = Eating food

# Confirm data variable types ----------------------------------------------------------------------------------------

numeric_dem_vars <- c("length_of_stay", "number_of_total_observations", "latest_weight", "estimated_age_at_departure")
factor_dem_vars <- c("sex", "neuter_status", "source_type", "adoption_site")

d_dem[, numeric_dem_vars] <- apply(d_dem[,numeric_dem_vars], 2, function(x) as_numeric_safe(x, remove_NAs = FALSE))
d_dem[, factor_dem_vars]  <- apply(d_dem[,factor_dem_vars], 2, function(x) as.factor(x))

# TABLE 1 ----------------------------------------------------------------------------------------------------------
dog_demographics <- c(numeric_dem_vars, factor_dem_vars)

print(table_1 <- sapply(dog_demographics, 
       function(x){
         if(x %in% c("sex", "neuter_status", "source_type", "adoption_site")){
           table( d_dem[, x])
         }
         else{
           list(mean(d_dem[,x], na.rm = T), sd(d_dem[,x], na.rm = T))
         }
       })
)

# TABLE 2 ----------------------------------------------------------------------------------------------------------

shelter_contexts <- unique(d_s$context_id)

print(table_2 <- data.frame(
  contexts = shelter_contexts, 
  median =  sapply(shelter_contexts, 
                   function(x) {
                     median(unlist(
                       lapply(
                         split(d_s, d_s$dog_id), 
                         FUN = function(z) {
                           length( 
                             z[z$context_id %in% x & !is.na(z$behaviour_code_colour), "behaviour_code_colour"])
                         }
                       )), na.rm=T)
                   }),
  IQR =  sapply(shelter_contexts, 
                function(x) {
                  IQR(unlist(
                    lapply(
                      split(d_s, d_s$dog_id), 
                      FUN = function(z) {
                        length( 
                          z[z$context_id %in% x & !is.na(z$behaviour_code_colour), "behaviour_code_colour"])
                      }
                    )), na.rm=T)
                }),
  avg_missing = sapply(shelter_contexts, 
                       function(x){
                         mean(
                           unlist(lapply(split(d_s, d_s$dog_id), 
                                         function(z) {
                                           sum( is.na(z[ z$context_id %in% x, "behaviour_code_colour"]) )/nrow(z[z$context_id %in% x, ])
                                         })
                           )
                         )
                       } 
  ),
  sd_missing = sapply(shelter_contexts, 
                      function(x){
                        sd(
                          unlist(lapply(split(d_s, d_s$dog_id), 
                                        function(z) {
                                          sum( is.na(z[ z$context_id %in% x, "behaviour_code_colour"]) )/nrow(z[z$context_id %in% x, ])
                                        })
                          )
                        )
                      } 
  )
))

mean(table_2[,"avg_missing"])
sd(table_2[,"avg_missing"])

# TABLE 3 ----------------------------------------------------------------------------------------------------------

# first need to create number of surveys complete variable
only_first <- lapply( split(d_a, d_a$dog_id), function(x) all( is.na(unique(x[x$survey_number=="Second","behaviour_code_number"]))))
print(num_with_only_first <- sum(unlist(only_first)))

# add a new variable saying how many surveys there are for each dog
d_dem$total_number_of_surveys <- unlist(lapply(split(d_a, d_a$dog_id),
                                       function(x) {
                                         ifelse( unique(x$dog_id) %in% which(unlist(only_first)==1), 1, 2) 
                                       }
                                  )
                            )

adoption_contexts <- unique(d_a$context_id)

table_3 <- data.frame(
  contexts = adoption_contexts,
  n_two_surveys = sapply(adoption_contexts, 
                         function(x){
                           sum(unlist(
                             lapply(split(d_a, d_a$dog_id), 
                                    function(z) ifelse( any(is.na(z[ z$context_id == x, "behaviour_code_colour"])), 0, 1))
                           ))
                         })
)


table_3[,"percent"] <- round(table_3$n_two_surveys/length(unique(d_a$dog_id))*100, 1)
print(table_3)

days_between_surveys <- unlist(lapply(split(d_a, d_a$dog_id, drop = T), 
                                      function(x) {
                                        if( all(is.na(x[x$survey_number %in% "Second", "behaviour_code_number"]))){# if one survey
                                          NA
                                        }
                                        else{
                                          diff(unique(x$days_after_adoption))
                                        }
                                        })
                               )

mean(days_between_surveys, na.rm = T)
sd(days_between_surveys, na.rm = T)
range(days_between_surveys, na.rm = T)


####################################################################################
# Prepare the predictor variables for the model ------------------------------------
####################################################################################

# missing demographic data
apply(d_dem, 2, function(x) sum( x %in% "NA" | is.na(x)))

# One dog is missing neuter status -- we enter the most likely value (ON SITE) for this small case
dog_missing_neuter_status <- which( is.na( d_dem[ , "neuter_status"]) )
d_dem[ dog_missing_neuter_status, "neuter_status"] <- "On Site"

# two dogs are missing weight and age at departure values -- mean impute because the size is small
dogs_missing_weights <- unique(d_dem[ is.na(d_dem$latest_weight), "dog_id"])
d_dem[ d_dem$dog_id %in% dogs_missing_weights, "latest_weight"] <- table_1$latest_weight[[1]]

# age at departure missing
dogs_missing_age <- unique(d_dem[ is.na(d_dem$estimated_age_at_departure), "dog_id"])
d_dem[ d_dem$dog_id %in% dogs_missing_age, "estimated_age_at_departure"] <- table_1$estimated_age_at_departure[[1]]

# we mean center and standardise continuous variables
# categorical variables are coded as sum-to-zero deflections from the mean

# shelter observation day variable
d_s$day_since_arrival_Z <- center_scale(x = d_s$day_since_arrival, 
                                        # center around the average length of stay ~ 30 days
                                        center = table_1$length_of_stay[[1]], 
                                        scale = sd(d_s$day_since_arrival)
                                        )

# post-adoption days_after_adoption variable
dog_means_days_after_adoption <- 
  as_numeric_safe(
    unlist(
  lapply(split(d_a, d_a$dog_id), 
         function(x) {
           # dogs with two surveys
           if( all(is.na(x[x$survey_number %in% "Second", "behaviour_code_number"] ))){
             mean(unique(x$days_after_adoption))
           }
           else{
             # single survey
             unique(x$days_after_adoption)
           }
           }
         )
))

d_a$days_after_adoption_Z <- center_scale(x = d_a$days_after_adoption, 
                                          center = mean(dog_means_days_after_adoption),
                                          scale = sd(d_a$days_after_adoption)
                                          )

# length_of_stay_Z
d_dem$length_of_stay_Z <- center_scale(x = d_dem$length_of_stay, 
                                       center = table_1$length_of_stay[[1]],
                                       scale = table_1$length_of_stay[[2]]
                                       )

# latest_weight_z
d_dem$latest_weight_Z <- center_scale(x = d_dem$latest_weight, 
                                      center = table_1$latest_weight[[1]],
                                      scale = table_1$latest_weight[[2]]
                                      )

# estimated_age_at_departure_Z
d_dem$estimated_age_at_departure_Z <- center_scale(x = d_dem$estimated_age_at_departure, 
                                                   center = table_1$estimated_age_at_departure[[1]],
                                                   scale = table_1$estimated_age_at_departure[[2]]
                                                   )

# sex_MALE = sum to zero contrasts
d_dem$sex_MALE <- ifelse(d_dem$sex %in% "Female", -1, 1)

# neuter_status sum to zero contrasts (YES as reference)
d_dem$neutered_NO <- ifelse(d_dem$neuter_status %in% "No", 1, 
                                ifelse(d_dem$neuter_status %in% "On Site", 0, 
                                       -1)
                            )

d_dem$neutered_ONSITE <- ifelse(d_dem$neuter_status %in% "On Site", 1, 
                                    ifelse(d_dem$neuter_status %in% "No", 0, 
                                           -1)
                                )

# source_type sum to zero contrasts (STRAY as reference)
d_dem$source_GIFT <- ifelse(d_dem$source_type == "Gift", 1, 
                                ifelse(d_dem$source_type == "Return", 0, 
                                       -1)
                            )

d_dem$source_RETURN <- ifelse(d_dem$source_type == "Return", 1, 
                                  ifelse(d_dem$source_type == "Gift", 0, 
                                         -1)
                              )

d_dem$total_number_of_surveys_TWO <- ifelse(d_dem$total_number_of_surveys %in% 2, 1, -1) 

# order the data sets by dog_id
d_s_order_by_id <- d_s[ order(d_s$dog_id), ]
d_a_order_by_id <- d_a[ order(d_a$dog_id), ]

# create matrices of dog-level predictor variables for the shelter data
dog_level_predictors_s <- c("length_of_stay_Z", "latest_weight_Z", "estimated_age_at_departure_Z",
                          "sex_MALE", "neutered_NO", "neutered_ONSITE", "source_GIFT", "source_RETURN" 
                          )

dog_level_predictors_a <- c(dog_level_predictors_s, "total_number_of_surveys_TWO")

# observation-length dog-level predictor matrices
X_s <- as.matrix(d_dem[ d_s_order_by_id$dog_id, dog_level_predictors_s])
X_a <- as.matrix(d_dem[ d_a_order_by_id$dog_id, dog_level_predictors_a])

# finally, we need to add a dog x context interaction variable
d_s_order_by_id$dog_context_id <- make_cc_interaction(jj = d_s_order_by_id$dog_id, d_s_order_by_id$context_id)
d_a_order_by_id$dog_context_id <- make_cc_interaction(jj = d_a_order_by_id$dog_id, d_a_order_by_id$context_id)


########################################################################################################
# package the data needed for the Stan model
########################################################################################################

# the Stan model is run outside R using cmdstan -- the command line interface for Stan
# alternatively, use cmdstanr or Rstan

d_s_stan <- d_s_order_by_id 
d_a_stan <- d_a_order_by_id 
X_s_stan <- X_s 
X_a_stan <- X_a 

stan_data <- list(
  N = c(nrow(d_s_stan), nrow(d_a_stan)), 
  N_j = length(unique(d_s_stan$dog_id)),
  N_g = length(unique(d_s_stan$context_id)),
  N_jg = length(unique(d_s_stan$dog_context_id)),
  N_x = c(ncol(X_s), ncol(X_a)), 
  x = d_s_stan$day_since_arrival_Z,
  x_joint = d_a_stan$days_after_adoption_Z,
  X_1 = X_s_stan, 
  X_2 = X_a_stan,
  y = ifelse( is.na(d_s_stan$behaviour_code_colour), 0, d_s_stan$behaviour_code_colour),
  y_joint = ifelse( is.na(d_a_stan$behaviour_code_colour), 0, d_a_stan$behaviour_code_colour),
  K = 3, 
  theta = c(1.5, 2.5),
  jj_1 = d_s_stan$dog_id, 
  jj_2 = d_a_stan$dog_id, 
  gg_1 = d_s_stan$context_id, 
  gg_2 = d_a_stan$context_id, 
  jjgg_1 = d_s_stan$dog_context_id, 
  jjgg_2 = d_a_stan$dog_context_id,
  jj_levels = 6, gg_levels = 6, jjgg_levels = 6, 
  use_reduce_sum = 1, 
  grainsize = 20
)

rstan::stan_rdump( list = ls(stan_data), file = "~/Documents/cmdstan-2.23.0/jhb-shelterdogs/real_analysis/dog_data_for_cmdstan.R", envir = list2env(stan_data))

# OR USING CMDSTANR
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# cmdstanr::set_cmdstan_path("~/.cmdstan/cmdstan-2.26.0")
# jhb_model <- cmdstanr::cmdstan_model(paste0(path, "analysis-code/jhb-ordinal-hurdle-model.stan"), cpp_options = list(stan_threads = TRUE))
# 
# 
# start <- Sys.time()
# fit <- jhb_model$sample(
#   data = stan_data,
#   chains = 2, paralel_chains=2, threads_per_chain=2,
#   init = 0, refresh = 100,
#   iter_warmup = 1000, iter_sampling = 5000
# )
# print(end <- Sys.time() - start)


