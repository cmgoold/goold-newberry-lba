###############################################################################
# Simulation functions for the example of the 
#   hierarchical Bayesian ordinal hurdle model with category one inflated responses, accompanying:
#   Goold & Newberry (2020). Longitudinal assessment of shelter dogs predicts beahviour post-adoption
###############################################################################

# Helper functions------------------------------------------------------------
# inverse logit function
inv_logit <- function(x) exp(x)/(1+exp(x))

# logit function
logit <- function(x) log(x/(1-x))

# creates a random correlation matrix
#   : setting reg = 1 simulates the correlations from a Normal(0, 0.4) distribution
#   : setting reg = 0 (default) simulates the correlations from a Uniform(-1, 1) distribution
create_rand_corr_matrix <- function(n, reg=0){
  r <- matrix(NA, nrow = n, ncol = n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(reg == 1){
        val <- rnorm(1, 0, 0.4)
        r[i,j] <- ifelse(abs(val) > 1, 0.9*sign(val), val)
      }
      else{
        r[i,j] <- runif(1,-1,1)
      }
      
      r[j,i] <- r[i,j]
    }
  }
  
  for(i in 1:n) r[i,i] <- 0
  
  return(r)
}

# creates a covariance matrix from standard deviations and correlation matrix
create_cov_matrix <- function(sigmas, Rho){
  symm_dim <- ncol(Rho)
  u <- matrix(NA, nrow = symm_dim, ncol = symm_dim)
  for(i in 1:(symm_dim-1)){
    for(j in (i+1):(symm_dim)){
      u[i,j] <- sigmas[i]*Rho[i,j]*sigmas[j]
      u[j,i] <- u[i,j]
    }
  }
  for(i in 1:(nrow(Rho))) u[i,i] <- sigmas[i]^2
  
  return(u)
}

# the main simulation function
sim_ord_data <- function(
  N_i, N_j, N_g,
  K, 
  N_i_joint,
  alpha, beta, sigma, 
  sigma_j, Rho_j,
  sigma_g, Rho_g,
  sigma_jg, Rho_jg,
  kappa,
  cut_points,
  alpha_miss, beta_miss,
  delta, gamma, epsilon, 
  delta_miss, gamma_miss)
{
  if(length(cut_points) != (K-1) ){
    stop("Number of cut-points does not equal K-1")
  }
  
  N <- N_i * N_j * N_g
  N_jg <- N_j * N_g
  U_cov_j <- create_cov_matrix(sigmas = sigma_j, Rho = Rho_j)
  U_cov_g <- create_cov_matrix(sigmas = sigma_g, Rho = Rho_g)
  U_cov_jg <- create_cov_matrix(sigmas = sigma_jg, Rho = Rho_jg)
  zeros <- as.vector(rep(0, 6))
  U_j <- MASS::mvrnorm(n = N_j, 
                       mu = zeros, 
                       Sigma = U_cov_j, 
                       tol = 1e3, 
                       empirical = T
  )
  jj_1 <- rep(1:N_j, each = N_i*N_g) # vector of ids
  
  U_g <- MASS::mvrnorm(n = N_g,
                       mu = zeros,
                       Sigma = U_cov_g,
                       tol = 1e3,
                       empirical = T
  )
  gg_1 <- rep(1:N_g, N_j*N_i)
  
  U_jg <- MASS::mvrnorm(n = N_jg,
                        mu = zeros,
                        Sigma = U_cov_jg,
                        tol = 1e3,
                        empirical = T
  )
  jjgg_1 <- jj_1*gg_1 - ( (jj_1-1)*(gg_1-1) ) + (N_g-1)*(jj_1-1)
  
  x <- rnorm(N, 0, 1)
  
  # linear predictor
  mu <- alpha + U_j[jj_1, 1] + U_g[gg_1, 1] + U_jg[jjgg_1, 1] +
    (beta + U_j[jj_1, 2] + U_g[gg_1, 2] + U_jg[jjgg_1, 2]) * x
  
  # missing model linear predictor
  eta <- alpha_miss + U_j[jj_1, 3] + U_g[gg_1, 3] + U_jg[jjgg_1, 3] + beta_miss * x
  is_miss <- rbinom(n = N, size = 1, prob = inv_logit(eta))
  
  # inflated component
  is_inflated <- rbinom(n = N, size = 1, prob = kappa[1])
  
  pi_ <- array(NA, dim=c(N,K))
  
  pi_[,1] <- pnorm( (cut_points[1] - mu) / sigma )
  for(j in 2:(K-1)){
    pi_[,j] <- pnorm( (cut_points[j] - mu)/sigma) - pnorm( (cut_points[j-1]-mu)/sigma)
  }
  pi_[,K] <- 1 - pnorm( (cut_points[K-1] - mu)/sigma )
  
  y <- numeric(N)
  for(i in 1:N){
    y[i] <- ifelse(is_inflated[i] == 1, (1 - is_miss[i]) * 1, (1 - is_miss[i]) * sample(1:K, 1, replace = TRUE, prob = pi_[i,]) )
  }
  
  # joint model
  N_joint <- N_i_joint * N_j * N_g
  x_joint_raw <- rep(1:N_i_joint, N_j*N_g) 
  x_joint <- (x_joint_raw - mean(x_joint_raw))/sd(x_joint_raw)
  
  jj_2 <- rep(1:N_j, each = N_i_joint*N_g)
  gg_2 <- rep(1:N_g, N_i_joint*N_j)
  jjgg_2 <- jj_2*gg_2 - ( (jj_2-1)*(gg_2-1) ) + (N_g-1)*(jj_2-1)
  
  zeta <- delta + U_j[jj_2, 4] + U_g[gg_2, 4] + U_g[gg_2, 4] + 
    (gamma + U_j[jj_2, 5] + U_g[gg_2, 5] + U_g[gg_2, 5]) * x_joint
  
  nu <- delta_miss + U_j[jj_2, 6] + U_g[gg_2, 6] + U_g[gg_2, 6] +
    gamma_miss * x_joint
  
  is_miss_joint <- rbinom(N_joint, 1, prob=inv_logit(nu))
  is_inflated_joint <- rbinom(N_joint, 1, prob=kappa[2])
  
  pi_joint_ <- array(NA, dim=c(N_joint,K))
  pi_joint_[,1] <- pnorm( (cut_points[1] - zeta) / epsilon )
  for(j in 2:(K-1)){
    pi_joint_[,j] <- pnorm( (cut_points[j] - zeta)/epsilon) - pnorm( (cut_points[j-1]-zeta)/epsilon)
  }
  pi_joint_[,K] <- 1 - pnorm( (cut_points[K-1] - zeta)/epsilon )
  
  y_joint <- numeric(N_joint)
  for(i in 1:N_joint){
    y_joint[i] <- ifelse(is_inflated_joint[i]==1, 
                         (1 - is_miss_joint[i]) * 1, 
                         (1 - is_miss_joint[i]) * sample(1:K, 1, replace=TRUE, prob = pi_joint_[i,])
    )
  }
  
  d_1 <- data.frame(
    y = y, 
    x = x, 
    jj_1 = jj_1, 
    gg_1 = gg_1, 
    jjgg_1 = jjgg_1
  )
  
  d_2 <- data.frame(
    y_joint = y_joint, 
    x_joint = x_joint, 
    jj_2 = jj_2, 
    gg_2 = gg_2, 
    jjgg_2 = jjgg_2
  )
  
  return(list(
    probs = pi_, 
    probs_joint = pi_joint_,
    d_1 = d_1, 
    d_2 = d_2
  ))
}