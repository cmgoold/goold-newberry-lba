############################################################################
# R helper functions for Goold & Newberry (2021):
#     - Predicting individual shelter dog behaviour after adoption using 
#       longitudinal behavioural assessment: a hierarchical Bayesian approach

# Copyright(C) Conor Goold  2021
# goold.conor@gmail.com
############################################################################

get_needed_packages <- function(x){
  needed <- c("rstan", "posterior", "coda", "scales")
  which_not_installed <- needed[which( !(needed %in% installed.packages()))]
  if(length(which_not_installed) > 0){
    for(i in which_not_installed){
      # posterior isn't on CRAN yet
      if(i == "posterior"){
        install.packages(
          "posterior", 
          repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
        )
      } else{
        install.packages(which_not_installed)
      }
    }
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
# clean scale
center_scale <- function(x, center=NULL, scale=NULL){
  
  if(is.null(center)){
    center = mean(x, na.rm = TRUE)
  }
  
  if(is.null(scale)){
    scale = sd(x, na.rm = TRUE)
  }
  
  return(
    (x - center) / scale
  )
}

#-----------------------------------------------------------------------------------------------
# make cross-classified interaction
make_cc_interaction <- function(jj, gg){
  N_g <- 8
  jj*gg - ( (jj - 1)*(gg - 1) ) + (N_g-1)*(jj-1)
}

#-----------------------------------------------------------------------------------------------
# hdi function
hdi <- function(x, prob=0.95, type="lower_upper"){
  
  types = c("lower_upper", "lower", "upper")
  
  if( !sum(types %in% type) ){
    stop( past0("type must be one of ", types) )
  }
  
  if(all(is.na(x))){
    return(NA)
  }
  
  hdi_ <- coda::HPDinterval(coda::as.mcmc(x), prob = prob)
  
  if(type == "lower"){
    hdi_ <- hdi_[1]
  }
  if(type == "upper"){
    hdi_ <- hdi_[2] 
  }
  
  return(hdi_)
}
#-----------------------------------------------------------------------------------------------
# get_ordinal_probs
#   : obtain ordinal category probabilities from latent scale

get_ordinal_probs <- function(latent_mu, latent_sigma)
{
  #if( dim(latent_mu)[1] != length(latent_sigma) ) stop("mu and sigma are not the same length")
  
  green <- pnorm( (1.5 - latent_mu)/latent_sigma )
  amber <- pnorm( (2.5 - latent_mu)/latent_sigma ) - pnorm( (1.5 - latent_mu)/latent_sigma )
  red <- 1 - pnorm( (2.5 - latent_mu)/latent_sigma )
  
  colour_df <- data.frame(
    green, amber, red
  )
  
  return(colour_df)
}

#-----------------------------------------------------------------------------------------------
# get_marginal_samples
#   : obtain samples of parameters marginal of the random effects
get_marginal_samples <- function(){
  
  out <- rep(list(list()), 18)
  
  draws_j <- draws[, grep("\\bz_j\\b", colnames(draws))]
  draws_g <- draws[, grep("\\bz_g\\b", colnames(draws))]
  draws_jg <- draws[, grep("\\bz_jg\\b", colnames(draws))]
  draws_sigmas <- as.matrix(draws[ , grep("sigma_", colnames(draws))])
  
  n_dogs <- length(unique(d_s_stan$dog_id))
  n_contexts <- length(unique(d_s_stan$context_id))
  n_dogs_contexts <- length(unique(d_s_stan$dog_context_id))
  
  for(i in 1:6){
      out[[i]] <- draws_j[((n_dogs - n_dogs + 1) + n_dogs*(i-1)):(n_dogs*i) ] * draws_sigmas[,i]
      out[[i+6]] <- draws_g[((n_contexts - n_contexts + 1) + n_contexts*(i-1)):(n_contexts*i) ] * draws_sigmas[,i+6]
      out[[i+6*2]] <- draws_jg[((n_dogs_contexts - n_dogs_contexts + 1) + n_dogs_contexts*(i-1)):(n_dogs_contexts*i) ] * draws_sigmas[,i+6*2]
  }
  
  names(out) <- c("r_s_j_intercept", "r_s_j_slope", "r_s_j_missing", 
                  "r_a_j_intercept", "r_a_j_slope", "r_a_j_missing",
                  "r_s_g_intercept", "r_s_g_slope", "r_s_g_missing", 
                  "r_a_g_intercept", "r_a_g_slope", "r_a_g_missing",
                  "r_s_jg_intercept", "r_s_jg_slope", "r_s_jg_missing", 
                  "r_a_jg_intercept", "r_a_jg_slope", "r_a_jg_missing"
                  )
  
  out
}

#---------------------------------------------------------------------------------------------------------
# get true positive, true negative, false positive, false negative values
get_diagnostic <- function(d1, d2, code_colour, type)
{
  return(
    unlist(
      lapply( split(d1, d1$dog_id), 
            function(x){
              any_red_s <- ifelse( any(x$behaviour_code_colour == code_colour, na.rm = T), 1, 0)
              any_red_a <- ifelse( any(d2[d2$dog_id %in% unique(x$dog_id), "behaviour_code_colour"] == code_colour, na.rm = T), 1, 0)
              
              if(type=="fp") if(any_red_s & !any_red_a) return(1) else return(0)
              if(type=="fn") if(!any_red_s & any_red_a) return(1) else return(0)
              if(type=="tp") if(any_red_s & any_red_a) return(1) else return(0)
              if(type=="tn") if(!any_red_s & !any_red_a) return(1) else return(0)
              
            })
    )
  )
}

#---------------------------------------------------------------------------------------------------------
# predictive values
ppv <- function(tp, fp){  
  if(!tp & fp > 0) 0 
  if(!tp & !fp) NA
  else tp/(tp + fp) 
  }
npv <- function(tn, fn){  
  if(!tn & fn > 0) 0 
  if(!tn & !fn) NA
  else tn/(tn + fn) 
  }
#end
#---------------------------------------------------------------------------------------------------------
