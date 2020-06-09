############################################################################
# R functions for Goold & Newberry (2020):
#     - Longitudinal behavioural assessment of shelter dogs predicts 
#       behaviour post-adoption

# Copyright Conor Goold (2020)
# c.goold@leeds.ac.uk
############################################################################

get_needed_packages <- function(x){
  needed <- c("rstan","coda")
  if(sum(installed.packages() %in% needed) < length(needed)){
    install.packages( needed[which( !(needed %in% installed.packages()))] )
  } 
  else{
    "All packages are already installed"
  }
}

# transform factors safely to numerics -----------------------------------------------
as_numeric_safe <- function(x, remove_NAs=TRUE) {
  message("Converting to character, then to numeric.")
  message("Non-number formats will be converted to NAs.")
  as_numeric <- suppressWarnings( as.numeric(as.character(x) ) )
  if(remove_NAs) return( as_numeric[!is.na(as_numeric)] )
  if(!remove_NAs) return( as_numeric )
}

#-----------------------------------------------------------------------------------------------
# inverse logistic function
inv_logit <- function(x) exp(x) / (1 + exp(x))

#-----------------------------------------------------------------------------------------------
# my highest density interval function, set with probability = 0.9, i.e. .5% - 95%
my_HDI <- function(x) HPDinterval( obj = as.mcmc(x), prob = .9 )

#-----------------------------------------------------------------------------------------------
# clean scale
center_scale <- function(x, mu, sd){
  return(
    (x - mu)/sd
  )
}

#-----------------------------------------------------------------------------------------------
# make cross-classified interaction
make_cc_interaction <- function(jj, gg){
  N_g <- length(unique(gg))
  jj*gg - ( (jj - 1)*(gg - 1) ) + (N_g-1)*(jj-1)
}

