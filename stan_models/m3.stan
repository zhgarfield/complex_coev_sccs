data {
  int N;
  int J;
  int K;
  int K_hunt;
  int y[N,J-2];
  int storage[N];
  int hunting[N];
}

parameters {
  ordered[K-1] c[J-2];  // cutpoints for ordered logit
  ordered[K_hunt-1] c_hunt;  // cutpoints for hunting
  real a_storage;  // intercept for food storage
  vector<lower=0>[6] lambda_TSD;  // factor loadings for TSD
  vector<lower=0>[5] lambda_RI;  // factor loadings for RI
  vector<lower=0>[2] lambda_CS;  // factor loadings for CS (only for Urbanization and Population Density)
  matrix[3,N] res_z;  // uncorrelated and unscaled latent factors
  vector<lower=0>[3] sigma_eta;
  real<lower=0,upper=1> res_cor[3,3];  // correlation matrix for three latent factors
}

transformed parameters {
  vector[J] l_TSD;  // vector of loadings for TSD
  vector[J] l_RI;  // vector of loadings for RI
  vector[J] l_CS;  // vector of loadings for CS (Community Size)
  matrix[N,3] res_v;  // correlated latent factors
  matrix[3,3] Rho_res;  // correlation matrix for the latent factors

  // Fill in positive-constrained correlation matrix for latent factors
  Rho_res = diag_matrix(rep_vector(1,3));
  Rho_res[1,2] = res_cor[1,2];
  Rho_res[2,1] = res_cor[1,2];
  Rho_res[1,3] = res_cor[1,3];
  Rho_res[3,1] = res_cor[1,3];
  Rho_res[2,3] = res_cor[2,3];
  Rho_res[3,2] = res_cor[2,3];

  // Assign loadings to RI, TSD, and CS factors
  l_RI[1] = 1;  // Fix first loading for RI to 1
  l_RI[2:4] = lambda_RI[1:3];
  l_RI[5:J] = rep_vector(0.0, J - 5);  // Set remaining to 0

  l_TSD[5:J-2] = lambda_TSD;
  l_TSD[1:4] = rep_vector(0.0, 4);  // No loadings for first 4 variables in TSD

  l_CS[3:4] = lambda_CS;  // CS loads on Urbanization and Density of Population
  l_CS[1:2] = rep_vector(0.0, 2);  // No loadings on other variables

  // Generate correlated latent factors
  res_v = (diag_pre_multiply(sigma_eta, cholesky_decompose(Rho_res)) * res_z)';
}

model {
  // Priors for cutpoints and intercepts
  for (j in 1:J-2) {
    c[j,] ~ normal(0,2);
  }
  c_hunt ~ normal(0,2);
  a_storage ~ normal(0,2);

  // Priors for factor loadings
  lambda_RI ~ std_normal();
  lambda_TSD ~ std_normal();
  lambda_CS ~ std_normal();

  // Priors for latent factors
  to_vector(res_z) ~ std_normal();
  sigma_eta ~ exponential(1);

  // Correlation constraints
  res_cor ~ beta(2,2);

  // Likelihood for each variable
  for (j in 1:J-2) {
    for (i in 1:N) {
      y[i,j] ~ ordered_logistic( res_v[i,1]*l_RI[j] + res_v[i,2]*l_TSD[j] + res_v[i,3]*l_CS[j], c[j,] );
    }
  }

  // Likelihood for hunting and storage
  for (i in 1:N) {
    hunting[i] ~ ordered_logistic( res_v[i,1]*l_RI[J-1] + res_v[i,2]*l_TSD[J-1] + res_v[i,3]*l_CS[J-1], c_hunt );
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
  matrix[N,J] log_lik;
  int store_hat[N];  // predictions for food storage

  for (i in 1:N) {
    store_hat[i] = storage[i];
    if (storage[i] == -99) {
      store_hat[i] = bernoulli_logit_rng(a_storage + res_v[i,1]*l_RI[J] + res_v[i,2]*l_TSD[J] + res_v[i,3]*l_CS[J]);
    }

    // Log likelihood for each variable
    for (j in 1:J-2) {
      log_lik[i,j] = ordered_logistic_lpmf( y[i,j] | res_v[i,1]*l_RI[j] + res_v[i,2]*l_TSD[j] + res_v[i,3]*l_CS[j], c[j,] );
    }

    log_lik[i,J-1] = ordered_logistic_lpmf( hunting[i] | res_v[i,1]*l_RI[J-1] + res_v[i,2]*l_TSD[J-1] + res_v[i,3]*l_CS[J-1], c_hunt );
    log_lik[i,J] = bernoulli_logit_lpmf( store_hat[i] | a_storage + res_v[i,1]*l_RI[J] + res_v[i,2]*l_TSD[J] + res_v[i,3]*l_CS[J] );
  }
}
