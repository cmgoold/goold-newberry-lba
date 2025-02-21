// custom joint Bayesian hierarchical ordinal hurdle model with one-inflation
// Copyright (C) 2021 Conor Goold

functions{

  /****************************************************/
  /* log-PDF for a single response from the custom ordinal-hurdle model with one-inflation
   *    y:   ordinal data point
   *    pi: linear predictor for the hurdle model
   *    pi: vector of ordinal category probabilities
   *    kappa: probability of one-inflation
  */
  real ordinal_hurdle_one_inflated_lpmf(int y, real eta, vector pi, real kappa){

    real lp = 0;

    if(y == 0){
      lp += bernoulli_logit_lpmf(1 | eta);
    }

    if(y == 1){
      vector[2] lp1;
      lp1[1] = bernoulli_logit_lpmf(0 | eta) + log(kappa);
      lp1[2] = bernoulli_logit_lpmf(0 | eta) + log1m(kappa) + categorical_lpmf(y | pi);
      lp += log_sum_exp(lp1);
    }

    if(y > 1){
      lp += bernoulli_logit_lpmf(0 | eta) + log1m(kappa) + categorical_lpmf(y | pi);
    }

    return lp;
  }
   /****************************************************/
   /* num_unique_values
    *
    */
   int num_unique_values(int[] x){
     int counter = 1;

     for(i in 1:(size(x)-1)){
       if( (x[i+1] - x[i]) != 0){
         counter += 1;
       }
     }
     return counter;
   }
  /****************************************************/
  /* find_id_rows_start
   * Find the row indices for the first observation for each hierarchical unit
   *    Only applies to sorted values
   */
  int[] find_id_rows_start( int[] id){
    // size of id
    int n_id_vec = size(id);
    // number of unique id units
    int n_id = num_unique_values(id);
    // position of first elemtent for each id
    int rows_start[n_id];

    for(i in 1:n_id){// loop through unique ids
      int pos = 1;
      for(j in 1:n_id_vec){// loop through id vector
        if(id[j] == i){
          // add the current position to rows_start
          rows_start[i] = pos;
          // break the loop
          break;
        }//if
        else{
          // if not the current id, increment position
          pos += 1;
        }//else
      }//for j
    }//for i

    return rows_start;
  }//find_id_rows_start
  /****************************************************/
  /* get_obs_end
   ** return the index of the last observation in the slice
   */
   int get_obs_end(int start, int end, int[] id_variable, int[] id_rows_start){
     //int n_id = size(unique_values(id_variable));
     if(end != num_unique_values(id_variable)){// if last unit in slice is not the last unit in the data
       // the last index
       return id_rows_start[end + 1] - 1;
     }//if
     else{// if we are at the last unit in the data
       return size(id_variable);
     }//else
   }//get_obs_end
  /****************************************************/
  /* make_slice_jj
   * Returns a new id vector for each slice from 1:n_obs_in_slice
   *    assumes that the vector is ordered
   */
   int[] make_slice_id(int[] id){
     // number of id units in slice
     //int n_id_in_slice = end - start + 1;
     // starting index of the current slice
     //int obs_start = j_rows_start[start];
     // end index of the current slice
     //int obs_end = get_obs_end(start, end, jj, j_rows_start);
     // number of observations in slice
     int n_obs_in_slice = size(id); //obs_end - obs_start + 1;
     // new slice_id
     int slice_id[n_obs_in_slice];

     slice_id[1] = 1;

     for(i in 1:(n_obs_in_slice-1)){
       if( (id[i+1] - id[i]) != 0){
         slice_id[i+1] = slice_id[i] + 1;
       }//if
       else
        slice_id[i+1] = slice_id[i];
     }//for
     return slice_id[];
   }//get_slice_id
  /****************************************************/
  /* make_interaction
   *
   */
  int[] make_interaction(int[] jj, int[] gg, int N_g){
    int n = size(jj);
    int res[n];

    for(i in 1:n){
      res[i] = jj[i] * gg[i] - ( (jj[i]-1)*(gg[i]-1)) + (N_g-1)*(jj[i]-1);
    }

    return res;
  }
  /****************************************************/
  /* partial_sum
   */
  real partial_sum(
    vector[] slice_z_j, // slice over id units (e.g. schools)
    int start, int end,
    int[] y, int[] y_joint,
    vector x, vector x_joint,
    real alpha, real beta, real sigma,
    real alpha_miss, real beta_miss,
    vector Beta, vector Beta_miss, matrix X_1,
    vector sigma_j, vector sigma_g, vector sigma_jg,
    vector[] z_g, vector[] z_jg,
    matrix L_j, matrix L_g, matrix L_jg,
    int jj_levels, int gg_levels, int jjgg_levels,
    int[] jj_1, int[] j_1_rows_start,
    int[] jj_2, int[] j_2_rows_start,
    int[] gg_1, int[] gg_2,
    int[] jjgg_1, int[] jjgg_2,
    int K, vector theta,
    vector kappa,
    real delta, real epsilon, real gamma,
    real delta_miss, real gamma_miss,
    vector Gamma, vector Gamma_miss, matrix X_2,
    int N_g, int N_jg, int[] N_x
  ){
    // log probability accumulator
    real lp = 0;

    /* get numbers and indices of different units in the slice */
    // number of j units in slice
    int n_j_in_slice = end - start + 1;
    // starting index of the current slice for the first dataset
    int obs_1_start = j_1_rows_start[start];
    // end index of the current slice for the first dataset
    int obs_1_end = get_obs_end(start, end, jj_1, j_1_rows_start);
    // starting index of the current slice for the second dataset
    int obs_2_start = j_2_rows_start[start];
    // end index of the current slice
    int obs_2_end = get_obs_end(start, end, jj_2, j_2_rows_start);
    // number of observations in slice in first dataset
    int n_1_obs_in_slice = obs_1_end - obs_1_start + 1;
    // number of observations in slice in second dataset
    int n_2_obs_in_slice = obs_2_end - obs_2_start + 1;
    // number of g units in slice
    int n_g_in_slice = N_g;
    // number of jg units in slice
    int n_jg_in_slice = n_j_in_slice * n_g_in_slice;

    /* create new id vectors needed for the slice */
    // new j vectors
    int slice_jj_1[n_1_obs_in_slice] = make_slice_id(jj_1[obs_1_start:obs_1_end]);
    int slice_jj_2[n_2_obs_in_slice] = make_slice_id(jj_2[obs_2_start:obs_2_end]);
    // new g vectors
    int slice_gg_1[n_1_obs_in_slice] = gg_1[obs_1_start:obs_1_end];
    int slice_gg_2[n_2_obs_in_slice] = gg_2[obs_2_start:obs_2_end];
    // new jg vectors
    int slice_jjgg_1[n_1_obs_in_slice] = make_interaction(slice_jj_1, slice_gg_1, N_g);
    int slice_jjgg_2[n_2_obs_in_slice] = make_interaction(slice_jj_2, slice_gg_2, N_g);


    /* slice up the data */
    matrix[n_j_in_slice, jj_levels] slice_r_j;
    matrix[n_g_in_slice, gg_levels] slice_r_g;
    vector[jjgg_levels] slice_z_jg[n_jg_in_slice] = z_jg[(jjgg_1[obs_1_start]):(jjgg_1[obs_1_end]),];
    matrix[n_jg_in_slice, jjgg_levels] slice_r_jg;
    // new y variable in slice
    int slice_y[n_1_obs_in_slice] = y[obs_1_start:obs_1_end];
    // new y_joint variable in slice
    int slice_y_joint[n_2_obs_in_slice] = y_joint[obs_2_start:obs_2_end];
    // new x variable
    vector[n_1_obs_in_slice] slice_x = x[obs_1_start:obs_1_end];
    // new x_joint variable
    vector[n_2_obs_in_slice] slice_x_joint = x_joint[obs_2_start:obs_2_end];
    // new X_1 matrix of j-level predictors
    matrix[n_1_obs_in_slice, N_x[1]] slice_X_1 = X_1[obs_1_start:obs_1_end, ];
    // new X_2 matrix of j-level predictors
    matrix[n_2_obs_in_slice, N_x[2]] slice_X_2 = X_2[obs_2_start:obs_2_end, ];
    // linear predictor
    vector[n_1_obs_in_slice] slice_mu = alpha + beta * slice_x + slice_X_1 * Beta;
    // linear predictor for missing model
    vector[n_1_obs_in_slice] slice_eta = alpha_miss + beta_miss * slice_x + slice_X_1 * Beta_miss;
    // linear predictor for the joint model
    vector[n_2_obs_in_slice] slice_zeta = delta + gamma * slice_x_joint + slice_X_2 * Gamma;
    // linear predictor for joint missing model
    vector[n_2_obs_in_slice] slice_nu = delta_miss + gamma_miss * slice_x_joint + slice_X_2 * Gamma_miss;


    /* additional variables */
    // ordinal category probabilities
    vector[K] pi[2];

    for(j in 1:n_j_in_slice) slice_r_j[j, ] = (diag_pre_multiply(sigma_j, L_j) * slice_z_j[j,])';
    for(g in 1:n_g_in_slice) slice_r_g[g, ] = (diag_pre_multiply(sigma_g, L_g) * z_g[g,])';
    for(jg in 1:n_jg_in_slice) slice_r_jg[jg,] = (diag_pre_multiply(sigma_jg, L_jg) * slice_z_jg[jg,])';

    for(n in 1:n_1_obs_in_slice){
      slice_mu[n] += slice_r_j[slice_jj_1[n], 1] + slice_r_g[slice_gg_1[n], 1] + slice_r_jg[slice_jjgg_1[n], 1] +
                     (slice_r_j[slice_jj_1[n], 2] + slice_r_g[slice_gg_1[n], 2] + slice_r_jg[slice_jjgg_1[n], 2]) * slice_x[n];
      slice_eta[n] += slice_r_j[slice_jj_1[n], 3] + slice_r_g[slice_gg_1[n], 3] + slice_r_jg[slice_jjgg_1[n], 3];

      pi[1,1] = Phi( (theta[1] - slice_mu[n])/sigma);
      pi[1,2] = Phi( (theta[2] - slice_mu[n])/sigma) -  Phi( (theta[1] - slice_mu[n])/sigma);
      pi[1,3] = 1 - Phi( (theta[2] - slice_mu[n])/sigma);

      lp += ordinal_hurdle_one_inflated_lpmf(slice_y[n] | slice_eta[n], pi[1,], kappa[1]);
    }//for

    for(n in 1:n_2_obs_in_slice){
      slice_zeta[n] += slice_r_j[slice_jj_2[n], 4] + slice_r_g[slice_gg_2[n], 4] + slice_r_jg[slice_jjgg_2[n], 4] +
                      (slice_r_j[slice_jj_2[n], 5] + slice_r_g[slice_gg_2[n], 5] + slice_r_jg[slice_jjgg_2[n], 5]) * slice_x_joint[n];
      slice_nu[n] += slice_r_j[slice_jj_2[n], 6] + slice_r_g[slice_gg_2[n], 6] + slice_r_jg[slice_jjgg_2[n], 6];

      pi[2,1] = Phi( (theta[1] - slice_zeta[n])/epsilon);
      pi[2,2] = Phi( (theta[2] - slice_zeta[n])/epsilon) -  Phi( (theta[1] - slice_zeta[n])/epsilon);
      pi[2,3] = 1 - Phi( (theta[2] - slice_zeta[n])/epsilon);

      lp += ordinal_hurdle_one_inflated_lpmf(slice_y_joint[n] | slice_nu[n], pi[2,], kappa[2]);
    }//for joint

    return lp;
  }//partial_sum
/****************************************************/
}//functions

data{
  int<lower=1> N[2];
	int<lower=1> N_j;
  int<lower=1> N_g;
  int<lower=1> N_jg;
  int<lower=1> N_x[2];
  vector[N[1]] x;
  vector[N[2]] x_joint;
  matrix[N[1], N_x[1]] X_1;
  matrix[N[2], N_x[2]] X_2;
  int y[N[1]];
  int y_joint[N[2]];
  int K;
	vector[K-1] theta;
	int jj_1[N[1]];
  int jj_2[N[2]];
  int jj_levels;
  int gg_1[N[1]];
  int gg_2[N[2]];
  int gg_levels;
  int jjgg_1[N[1]];
  int jjgg_2[N[2]];
  int jjgg_levels;
  int use_reduce_sum;
  int grainsize;
}

transformed data{
  int j_1_rows_start[N_j] = find_id_rows_start(jj_1);
  int j_2_rows_start[N_j] = find_id_rows_start(jj_2);
}

parameters{
  real alpha;
  real beta;
  real<lower=0> sigma;
  real alpha_miss;
  real beta_miss;
  real delta;
  real gamma;
  real delta_miss;
  real gamma_miss;
  real<lower=0> epsilon;
  vector[N_x[1]] Beta;
  vector[N_x[1]] Beta_miss;
  vector[N_x[2]] Gamma;
  vector[N_x[2]] Gamma_miss;
  vector<lower=0, upper=1>[2] kappa;
	vector<lower=0>[jj_levels] sigma_j;
  vector<lower=0>[gg_levels] sigma_g;
  vector<lower=0>[jjgg_levels] sigma_jg;
  cholesky_factor_corr[jj_levels] L_j;
  cholesky_factor_corr[gg_levels] L_g;
  cholesky_factor_corr[jjgg_levels] L_jg;
	vector[jj_levels] z_j[N_j];
  vector[gg_levels] z_g[N_g];
  vector[jjgg_levels] z_jg[N_jg];
}

model{
  // priors
  alpha ~ normal(0, 5);
  beta ~ normal(0, 1);
  Beta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  epsilon ~ normal(0, 1);
  alpha_miss ~ normal(0, 5);
  beta_miss ~ normal(0, 1);
  Beta_miss ~ normal(0, 1);
	for(i in 1:jj_levels) {
    z_j[,i] ~ normal(0, 1);
    z_g[,i] ~ normal(0, 1);
    z_jg[,i] ~ normal(0, 1);
  }
	sigma_j ~ normal(0, 1);
  sigma_g ~ normal(0, 1);
  sigma_jg ~ normal(0, 1);
  kappa ~ beta(2,2);
  delta ~ normal(0, 5);
  gamma ~ normal(0, 1);
  Gamma ~ normal(0, 1);
  delta_miss ~ normal(0, 5);
  gamma_miss ~ normal(0, 1);
  Gamma_miss ~ normal(0, 1);
  L_j ~ lkj_corr_cholesky(3.0);
  L_g ~ lkj_corr_cholesky(3.0);
  L_jg ~ lkj_corr_cholesky(3.0);

  // reduce_sum likelihood
  if(use_reduce_sum == 1){
    target += reduce_sum(
      partial_sum,
      z_j,
      grainsize,
      y, y_joint,
      x, x_joint,
      alpha, beta, sigma,
      alpha_miss, beta_miss,
      Beta, Beta_miss, X_1,
      sigma_j, sigma_g, sigma_jg,
      z_g, z_jg,
      L_j, L_g, L_jg,
      jj_levels, gg_levels, jjgg_levels,
      jj_1, j_1_rows_start,
      jj_2, j_2_rows_start,
      gg_1, gg_2,
      jjgg_1, jjgg_2,
      K, theta,
      kappa,
      delta, epsilon, gamma,
      delta_miss, gamma_miss,
      Gamma, Gamma_miss, X_2,
      N_g, N_jg, N_x
    );
  }
  else{
    vector[K] pi[2];
    // linear predictors (vectorized)
    vector[N[1]] mu = alpha + beta * x + X_1 * Beta;
    vector[N[1]] eta = alpha_miss + beta_miss * x + X_1 * Beta_miss;
    vector[N[2]] zeta = delta + gamma * x_joint + X_2 * Gamma;
    vector[N[2]] nu = delta_miss + gamma_miss * x_joint + X_2 * Gamma_miss;
    // random effect containers
    matrix[N_j, jj_levels] r_j;
    matrix[N_g, gg_levels] r_g;
    matrix[N_jg, jjgg_levels] r_jg;

    // fill the random effects
    for(j in 1:N_j) r_j[j,] = (diag_pre_multiply(sigma_j, L_j) * z_j[j,])';
    for(g in 1:N_g) r_g[g,] = (diag_pre_multiply(sigma_g, L_g) * z_g[g,])';
    for(jg in 1:N_jg) r_jg[jg,] = (diag_pre_multiply(sigma_jg, L_jg) * z_jg[jg,])';

    // likelihood

    for(n in 1:N[1]){
  			mu[n] += r_j[jj_1[n], 1] + r_g[gg_1[n], 1] + r_jg[jjgg_1[n], 1] +
                 (r_j[jj_1[n], 2] + r_g[gg_1[n], 2] + r_jg[jjgg_1[n], 2]) * x[n];
        eta[n] += r_j[jj_1[n], 3] + r_g[gg_1[n], 3] + r_jg[jjgg_1[n], 3];

        pi[1,1] = Phi( (theta[1] - mu[n])/sigma);
        pi[1,2] = Phi( (theta[2] - mu[n])/sigma) -  Phi( (theta[1] - mu[n])/sigma);
        pi[1,3] = 1 - Phi( (theta[2] - mu[n])/sigma);

        target += ordinal_hurdle_one_inflated_lpmf(y[n] | eta[n], pi[1,], kappa[1]);
      }//for

    for(n in 1:N[2]){
      zeta[n] += r_j[jj_2[n], 4] + r_g[gg_2[n], 4] + r_jg[jjgg_2[n], 4] +
                (r_j[jj_2[n], 5] + r_g[gg_2[n], 5] + r_jg[jjgg_2[n], 5]) * x_joint[n];
      nu[n] += r_j[jj_2[n], 6] + r_g[gg_2[n], 6] + r_jg[jjgg_2[n], 6];

      pi[2,1] = Phi( (theta[1] - zeta[n])/epsilon);
      pi[2,2] = Phi( (theta[2] - zeta[n])/epsilon) -  Phi( (theta[1] - zeta[n])/epsilon);
      pi[2,3] = 1 - Phi( (theta[2] - zeta[n])/epsilon);

      target += ordinal_hurdle_one_inflated_lpmf(y_joint[n] | nu[n], pi[2,], kappa[2]);
    }
  }//else
}//model

generated quantities{
  matrix[jj_levels,jj_levels] Rho_j;
  matrix[gg_levels, gg_levels] Rho_g;
  matrix[jjgg_levels, jjgg_levels] Rho_jg;
  Rho_j = L_j * L_j';
  Rho_g = L_g * L_g';
  Rho_jg = L_jg * L_jg';
}
