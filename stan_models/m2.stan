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
  vector[N] phy_z;  // phylogenetic random effects, unscaled and uncorrelated
  real<lower=0> eta;      // maximum covariance for GP functions
  real<lower=0> rho;      // rate of decay for GP functions
  vector[N] eta_z; // complexity latent scores, unscaled
  real<lower=0> sigma_eta; // scale of latent variable
  vector[J-1] loadings; // factor loadings, -1 for agriculture which is set to 1 for identifiability
}

transformed parameters{
  vector[N] phy_v;  // phylogenetic random effects, scaled and correlated
  vector[J] lambda_C; // factor loadings for all observed variables
  vector[N] eta_C;
  
  lambda_C[1] = 1;
  lambda_C[2:J] = loadings;
  eta_C = eta_z * sigma_eta;
  
  
  {
    matrix[N,N] phy_cov; // phylogenetic correlation matrices
    matrix[N,N] L_phy_cov; // cholesky decomposition matrix
    
    for ( i in 1:(N-1) )
      for ( m in (i+1):N ) {
        phy_cov[i,m] = eta*exp(-(rho*dist_mat[i,m]));
        phy_cov[m,i] = phy_cov[i,m];
      }
    for ( q in 1:N )
      phy_cov[q,q] = eta + 0.01;
    L_phy_cov = cholesky_decompose(phy_cov);
    phy_v = L_phy_cov * phy_z;
  }
}

model{
for (j in 1:J-2) {
c[j,] ~ normal(0,2);
}
c_hunt ~ normal(0,2);
a_storage ~ normal(0,2);

eta_z ~ std_normal();
sigma_eta ~ exponential(1);
loadings ~ std_normal();

to_vector(phy_z) ~ std_normal();
rho ~ exponential(1);
eta ~ exponential(1);

for (j in 1:J-2)
for (i in 1:N) {
y[i,j] ~ ordered_logistic( (eta_C[i] + phy_v[i])*lambda_C[j], c[j,] );
}

for (i in 1:N) {
  hunting[i] ~ ordered_logistic( (eta_C[i] + phy_v[i])*lambda_C[J-1], c_hunt );
  
  if (storage[i] == -99) {
  target += log_mix( inv_logit( a_storage + (eta_C[i] + phy_v[i])*lambda_C[J] ),
  bernoulli_logit_lpmf( 1 | a_storage + (eta_C[i] + phy_v[i])*lambda_C[J] ),
  bernoulli_logit_lpmf( 0 | a_storage + (eta_C[i] + phy_v[i])*lambda_C[J] )
  );
  }
  
  else
  storage[i] ~ bernoulli_logit( a_storage + (eta_C[i] + phy_v[i])*lambda_C[J] );
}
}

generated quantities{
matrix[N,J] log_lik;
int store_hat[N]; // food storage, rng for missing

for (i in 1:N) {
  store_hat[i] = storage[i];
  
  if (storage[i] == -99) {
    store_hat[i] = bernoulli_logit_rng(a_storage + (eta_C[i] + phy_v[i])*lambda_C[J]);
  }
  
  for (j in 1:J-2) {
    log_lik[i,j] = ordered_logistic_lpmf( y[i,j] | (eta_C[i] + phy_v[i])*lambda_C[j], c[j,] );
  }
  log_lik[i,J-1] = ordered_logistic_lpmf( hunting[i] | (eta_C[i] + phy_v[i])*lambda_C[J-1], c_hunt );
  log_lik[i,J] = bernoulli_logit_lpmf( store_hat[i] | a_storage + (eta_C[i] + phy_v[i])*lambda_C[J] );
}

}
