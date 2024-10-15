data {
  int N;
  int J;
  int K;
  int K_hunt;
  array[N, J-2] int y;  // Response variables for the complexity dimensions
  array[N] int storage;  // Food storage data
  array[N] int hunting;  // Hunting data
}

parameters {
  array[J-2] ordered[K-1] c;  // Cutpoints for ordered logit for response variables
  ordered[K_hunt-1] c_hunt;  // Cutpoints for hunting variable
  real a_storage;  // Intercept for food storage
  vector<lower=0>[6] lambda_TSD;  // Factor loadings for TSD
  vector<lower=0>[5] lambda_RI;  // Factor loadings for RI
  vector<lower=0>[2] lambda_CS;  // Factor loadings for CS (Community Size)
  matrix[3, N] res_z;  // Uncorrelated and unscaled latent factors for RI, TSD, and CS
  vector<lower=0>[3] sigma_eta;  // Standard deviations for the latent factors
  cholesky_factor_corr[3] L_Rho_res;  // Cholesky factor of the correlation matrix for the latent factors
}

transformed parameters {
  vector[J] l_TSD;  // Full vector of loadings for TSD
  vector[J] l_RI;  // Full vector of loadings for RI
  vector[J] l_CS;  // Full vector of loadings for CS
  matrix[N, 3] res_v;  // Correlated latent factors

  // Fill in factor loadings for TSD, RI, and CS
  l_TSD[1:4] = rep_vector(0.0, 4);  // No loadings for first 4 variables for TSD
  l_TSD[5:J-2] = lambda_TSD;  // Loadings for variables 5 to J-2 for TSD
  l_TSD[(J-1):J] = rep_vector(0.0, 2);  // No loadings for last 2 variables for TSD

  l_RI[1] = 1;  // Fix the first factor loading to 1 for RI
  l_RI[2:4] = lambda_RI[1:3];
  l_RI[5:(J-2)] = rep_vector(0.0, 6);  // No loadings for middle variables for RI
  l_RI[(J-1):J] = lambda_RI[4:5];  // Loadings for last 2 variables for RI

  l_CS[1:2] = lambda_CS;  // Loadings for first two variables (Urbanization and Density of Population) for CS
  l_CS[3:J] = rep_vector(0.0, J-2);  // No loadings for other variables for CS

  // Scale latent factors with Cholesky factor of the correlation matrix
  res_v = (diag_pre_multiply(sigma_eta, L_Rho_res) * res_z)';
}

model {
  // Priors for cutpoints, intercept, and factor loadings
  for (j in 1:J-2) {
    c[j] ~ normal(0, 2);  // Prior for cutpoints for response variables
  }
  c_hunt ~ normal(0, 2);  // Prior for cutpoints for hunting variable
  a_storage ~ normal(0, 2);  // Prior for intercept of food storage

  to_vector(res_z) ~ std_normal();  // Prior for the uncorrelated latent factors
  sigma_eta ~ normal(0, 1);  // Stronger prior for the standard deviations of the latent factors
  L_Rho_res ~ lkj_corr_cholesky(3);  // LKJ prior on the Cholesky factor of the correlation matrix with stronger regularization

  lambda_TSD ~ std_normal();  // Prior for the loadings of TSD
  lambda_RI ~ std_normal();  // Prior for the loadings of RI
  lambda_CS ~ std_normal();  // Prior for the loadings of CS

  // Likelihood for the complexity variables
  for (j in 1:J-2) {
    for (i in 1:N) {
      y[i,j] ~ ordered_logistic( res_v[i,1]*l_RI[j] + res_v[i,2]*l_TSD[j] + res_v[i,3]*l_CS[j], c[j] );
    }
  }

  // Likelihood for the hunting variable
  for (i in 1:N) {
    hunting[i] ~ ordered_logistic( res_v[i,1]*l_RI[J-1] + res_v[i,2]*l_TSD[J-1] + res_v[i,3]*l_CS[J-1], c_hunt );

    // Likelihood for the food storage variable (mixture model for missing data)
    if (storage[i] == -99) {
      target += log_mix( inv_logit( a_storage + res_v[i,1]*l_RI[J] + res_v[i,2]*l_TSD[J] + res_v[i,3]*l_CS[J] ),
                        bernoulli_logit_lpmf( 1 | a_storage + res_v[i,1]*l_RI[J] + res_v[i,2]*l_TSD[J] + res_v[i,3]*l_CS[J] ),
                        bernoulli_logit_lpmf( 0 | a_storage + res_v[i,1]*l_RI[J] + res_v[i,2]*l_TSD[J] + res_v[i,3]*l_CS[J] ));
    } else {
      storage[i] ~ bernoulli_logit( a_storage + res_v[i,1]*l_RI[J] + res_v[i,2]*l_TSD[J] + res_v[i,3]*l_CS[J] );
    }
  }
}

generated quantities {
  array[N, J] real log_lik;  // Log likelihood
  array[N] int store_hat;  // Food storage predictions for missing values

  for (i in 1:N) {
    store_hat[i] = storage[i];
    
    if (storage[i] == -99) {
      store_hat[i] = bernoulli_logit_rng(a_storage + res_v[i,1]*l_RI[J] + res_v[i,2]*l_TSD[J] + res_v[i,3]*l_CS[J]);
    }
    
    for (j in 1:J-2) {
      log_lik[i,j] = ordered_logistic_lpmf( y[i,j] | res_v[i,1]*l_RI[j] + res_v[i,2]*l_TSD[j] + res_v[i,3]*l_CS[j], c[j] );
    }
    log_lik[i,J-1] = ordered_logistic_lpmf( hunting[i] | res_v[i,1]*l_RI[J-1] + res_v[i,2]*l_TSD[J-1] + res_v[i,3]*l_CS[J-1], c_hunt );
    log_lik[i,J] = bernoulli_logit_lpmf( store_hat[i] | a_storage + res_v[i,1]*l_RI[J] + res_v[i,2]*l_TSD[J] + res_v[i,3]*l_CS[J] );
  }
}
