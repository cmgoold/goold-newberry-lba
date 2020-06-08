############################################################################
# R script for Goold & Newberry (2020):
#     - Longitudinal behavioural assessment of shelter dogs predicts 
#       behaviour post-adoption

# Copyright(C) Conor Goold  2020
# c.goold@leeds.ac.uk
############################################################################

source("~/Dropbox/PhD/PhD_NMBU/PaperIV/GooldNewberry2020-lba/analysis_code/helper_functions.R")

# install packages required
get_needed_packages()

# load the data
d_s <- read.table("~/Dropbox/PhD/PhD_NMBU/PaperIV/GooldNewberry2020-lba/data/d_shelter_behaviour.txt", header = T)
d_a <- read.table("~/Dropbox/PhD/PhD_NMBU/PaperIV/GooldNewberry2020-lba/data/d_post_adoption_behaviour.txt", header=T)
d_dem <- read.csv("~/Dropbox/PhD/PhD_NMBU/PaperIV/GooldNewberry2020-lba/data/d_demographic_details.txt", header = T)

# TABLE 1 ----------------------------------------------------------------------------------------------------------
shelter_demographics <- c(
  "length_of_stay", "number_of_total_observations", "sex", "latest_weight", "estimated_age_at_departure",
  "neuter_status", "source_type", "adoption_site"
)

print(table_1 <- sapply(shelter_demographics, 
       function(x){
         if(x %in% c("sex", "neuter_status", "source_type", "adoption_site")){
           table( d_dem[, x])
         }
         else{
           list(mean(d_dem[,x], na.rm = T), sd(d_dem[,x], na.rm = T))
         }
       })
)
