###############################################################################
# Simulation example of the hierarchical Bayesian ordinal hurdle model with category
#   one inflated responses, accompanying:
#   Goold & Newberry (2021).
###############################################################################

# requires cmdstan, cmdstanr and rstan installed on your local machine
# alternatively, you can just use rstan to fit the models below (not cmdstanr)

path <- "" # set to the goold-newberry-lba folder
figure_path <- paste0(path, "paper/figures/")
source(paste0(path, "simulation-code/joint-hierarchical-bayesian-model_simulation-example_functions.R"))
setwd(path)

# DESCRIPTION:
#   The model being estimated in the paper is a hierarchical Bayesian joint likelihood model for the shelter and post-adoption data
#   For the first and second (referred to as 'joint') likelihoods, we fit a mixture model comprised of
#     1) a hurdle component for the missing data,
#     2) an inflation component for the hypothesised extra probability mass at category 1 (i.e. green codes)
#     3) an ordinal probit model for the non-missing and non-inflated categories
#   The full text and supplementary materials has more infomation

# set the seed
set.seed(2020)
# number of individuals
N_j = 20;
# number of contexts
N_g = 8
# number of data points for the first likelihood
N_i = 3;
# number of data points for the second (joint) likelihood
N_i_joint = 2;
# correlation matrices for the individual, context and interaction effects
Rho_j <- create_rand_corr_matrix(n = 6)
Rho_g <- create_rand_corr_matrix(n = 6)
Rho_jg <- create_rand_corr_matrix(n = 6)
# standard deviations of the random effects
sigma_j = runif(6, 0.5, 3);
sigma_g = runif(6, 0.5, 3);
sigma_jg = runif(6, 0.5, 3)
# number of ordinal categories
K = 3
# first likelihood intercept, slope and residual SD
alpha = 1; beta = -1; sigma = 1.5;
# missing data model's intercept and slope for the first likelihood
alpha_miss = logit(0.3); beta_miss = log(1.2)
# cutpoints for the ordinal model
cut_points = seq(1.5, 2.5, 1)
# category one inflation probability
kappa = c(0.2, 0.3)
# intercept, slope and residual SD for the joint likelihood ordinal model
delta = 0.5; gamma = 0.5; epsilon = 1
# intercept and slope for the joint likelihood missing data (Bernoulli) model
delta_miss = logit(0.05); gamma_miss = log(5)

# generate the data
d <- sim_ord_data(N_i = N_i, N_j = N_j, N_g = N_g,
                  N_i_joint = N_i_joint,
                  K = K,
                  alpha = alpha,
                  beta = beta, sigma = sigma,
                  alpha_miss = alpha_miss, beta_miss = beta_miss,
                  sigma_j = sigma_j, sigma_g = sigma_g, sigma_jg = sigma_jg,
                  Rho_j = Rho_j, Rho_g = Rho_g, Rho_jg = Rho_jg,
                  cut_points = cut_points,
                  kappa = kappa,
                  delta = delta, epsilon = epsilon,
                  gamma = gamma, delta_miss = delta_miss, gamma_miss = gamma_miss
)

# set the path to your installation of cmdstan
cmdstanr::set_cmdstan_path("YOUR PATH HERE")
ordinal_model <- cmdstanr::cmdstan_model("YOUR PATH HERE",cpp_options = list(stan_threads = TRUE))

# fit the model
#   :- times for the running the model vary! Try 20 individuals, 8 contexts and < 10 time points to start
start <- Sys.time()
fit <- ordinal_model$sample(
  # data list for Stan
  data = list(use_reduce_sum = 1, grainsize=1,
              y = d$d_1$y,
              y_joint = d$d_2$y_joint,
              x = d$d_1$x, x_joint = d$d_2$x_joint,
              N = c(length(d$d_1$y), length(d$d_2$y_joint)),
              jj_levels = 6, gg_levels = 6, jjgg_levels = 6,
              N_j = length(unique(d$d_1$jj_1)),
              N_g = length(unique(d$d_1$gg_1)),
              N_jg = length(unique(d$d_1$jjgg_1)),
              jj_1 = d$d_1$jj_1, jj_2 = d$d_2$jj_2,
              gg_1 = d$d_1$gg_1, gg_2 = d$d_2$gg_2,
              jjgg_1 = d$d_1$jjgg_1, jjgg_2 = d$d_2$jjgg_2,
              K = 3,
              theta=c(1.5,2.5)
  ),
  # chains and cores
  chains = 1, parallel_chains = 2, threads_per_chain=2,
  # inital values
  init = 0,
  refresh = 100
)
print(end <- Sys.time() - start)

## divergent transitions may occur when there is not a lot of data
## the model parameters are recovered however, even when divergent transitions are present.

# inspect the results
print( rstan::read_stan_csv(fit$output_files()), pars = c("alpha", "beta", "sigma", 
                                                          "alpha_miss", "beta_miss",
                                                          "sigma_j", "sigma_g", "sigma_jg",
                                                          "kappa",
                                                          "delta","epsilon","gamma",
                                                          "delta_miss", "gamma_miss",
                                                          "Rho_j", "Rho_g", "Rho_jg"
                                                          )
)
