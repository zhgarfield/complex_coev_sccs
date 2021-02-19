data{
  int N;
  int J;
  int K;
  int K_hunt;
  int y[N,J-2];
  int storage[N];
  int hunting[N];
  matrix[N,N] dist_mat;
}

parameters{
  ordered[K-1] c[J-2];  // cutpoints for ordered logit
  ordered[K_hunt-1] c_hunt; // cutpoints for hunting
  real a_storage; // intercept for food storage
  matrix[N,J] phy_z;  // phylogenetic random effects, unscaled and uncorrelated
  vector<lower=0>[J] eta;      // maximum covariance for GP functions
  vector<lower=0>[J] rho;      // rate of decay for GP functions
  vector<lower=0>[N] eta_C; // complexity latent scores, contrained > 0 for identifiability
  vector[J] lambda_C; // factor loadings
}

transformed parameters{
  matrix[N,J] phy_v;  // phylogenetic random effects, scaled and correlated
  
  for (j in 1:J) {
    matrix[N,N] phy_cov; // phylogenetic correlation matrices
    matrix[N,N] L_phy_cov; // cholesky decomposition matrix
    
    for ( i in 1:(N-1) )
      for ( m in (i+1):N ) {
        phy_cov[i,m] = eta[j]*exp(-(rho[j]*dist_mat[i,m]));
        phy_cov[m,i] = phy_cov[i,m];
      }
    for ( q in 1:N )
      phy_cov[q,q] = eta[j] + 0.01;
    L_phy_cov = cholesky_decompose(phy_cov);
    phy_v[,j] = L_phy_cov * phy_z[,j];
  }
}

model{
for (j in 1:J-2) {
c[j,] ~ normal(0,1);
}
c_hunt ~ normal(0,1);
a_storage ~ normal(0,1);

eta_C ~ lognormal(0,1);
lambda_C ~ normal(0,1);

to_vector(phy_z) ~ normal(0,1);
rho ~ exponential(1);
eta ~ exponential(1);

for (j in 1:J-2)
for (i in 1:N) {
y[i,j] ~ ordered_logistic( phy_v[i,j] + eta_C[i]*lambda_C[j], c[j,] );
}

for (i in 1:N) {
  hunting[i] ~ ordered_logistic( phy_v[i,J-1] + eta_C[i]*lambda_C[J-1], c_hunt );
  
  if (storage[i] == -99) {
  target += log_mix( inv_logit( a_storage + phy_v[i,J] + eta_C[i]*lambda_C[J] ),
  bernoulli_logit_lpmf( 1 | a_storage + phy_v[i,J] + eta_C[i]*lambda_C[J] ),
  bernoulli_logit_lpmf( 0 | a_storage + phy_v[i,J] + eta_C[i]*lambda_C[J] )
  );
  }
  
  else
  storage[i] ~ bernoulli_logit( a_storage + phy_v[i,J] + eta_C[i]*lambda_C[J] );
}
}

generated quantities{
matrix[N,J] log_lik;
int store_hat[N]; // food storage, rng for missing

for (i in 1:N) {
  store_hat[i] = storage[i];
  
  if (storage[i] == -99) {
    store_hat[i] = bernoulli_logit_rng(a_storage + phy_v[i,J] + eta_C[i]*lambda_C[J]);
  }
  
  for (j in 1:J-2) {
    log_lik[i,j] = ordered_logistic_lpmf( y[i,j] | phy_v[i,j] + eta_C[i]*lambda_C[j], c[j,] );
  }
  log_lik[i,J-1] = ordered_logistic_lpmf( hunting[i] | phy_v[i,J-1] + eta_C[i]*lambda_C[J-1], c_hunt );
  log_lik[i,J] = bernoulli_logit_lpmf( store_hat[i] | a_storage + phy_v[i,J] + eta_C[i]*lambda_C[J] );
}

}
