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
  
  matrix[N,J] phy_z;  // phylogeN_socetic raN_socdom effects, uN_socscaled aN_socd uN_soccorrelated
  vector<lower=0>[J] sigma_phy; // how quickly covariaN_socce decliN_soces with phylo distaN_socce
  vector<lower=0>[J] rho_phy; // how quickly covariaN_socce decliN_soces with phylo distaN_socce
  cholesky_factor_corr[J] L_phy;  // correlatioN_socs betweeN_soc phy effects
  
  matrix[J,N] eta_z;  // residual random effects, unscaled and uncorrelated
  vector<lower=0>[J] sigma_eta; // scaling of OLV
  cholesky_factor_corr[J] L_eta;  // residual correlation matrix
}

transformed parameters{
  matrix[N,J] eta_v;  // residual random effects, scaled and correlated (note the transpose here)
  matrix[N,J] phy_v;  // phylogenetic random effects, scaled and correlated
  
    // Phylogenetic covariance function
    for (f in 1:J) {
    matrix[N,N] phy_cov; // phylogenetic/geographic/temporal correlation matrices
    matrix[N,N] L_phy_cov; // cholesky decomposition matrix
    
    for ( i in 1:(N-1) )
      for ( m in (i+1):N ) {
        phy_cov[i,m] = exp( -(rho_phy[f]*dist_mat[i,m]) );
        // Make upper tri and lower tri equal
        phy_cov[m,i] = phy_cov[i,m];
      }
    for ( q in 1:N )
    phy_cov[q,q] = 1 + 0.001; // adding small constant to diagonal to keep positive semidefinite
    
    L_phy_cov = cholesky_decompose(phy_cov);
    phy_v[,f] = L_phy_cov * phy_z[,f];
    }
    /// Also correlate phy effects across factors, matrix normal sampling
    phy_v = phy_v * diag_pre_multiply(sigma_phy, L_phy)';
    
    
  // Scale and correlate residual variance
  eta_v = (diag_pre_multiply(sigma_eta, L_eta) * eta_z)';
  
  // Update latent variable based on phylogeny
  for (n in 1:N) {
  for (f in 1:J) {
    eta_v[n,f] = eta_v[n,f] + phy_v[n,f];
  }
  }
}

model{
for (j in 1:J-2) {
c[j,] ~ normal(0,2);
}
c_hunt ~ normal(0,2);
a_storage ~ normal(0,2);

// Latent-scale covariance
to_vector(eta_z) ~ std_normal();
L_eta ~ lkj_corr_cholesky(2);
sigma_eta ~ std_normal();

to_vector(phy_z) ~ std_normal();
rho_phy  ~ std_normal(); // scaling for phylogenetic distance
sigma_phy  ~ std_normal();
L_phy ~ lkj_corr_cholesky(2);


for (j in 1:J-2)
for (n in 1:N) {
y[n,j] ~ ordered_logistic( eta_v[n,j] + phy_v[n,j] , c[j,] );
}

for (n in 1:N) {
  hunting[n] ~ ordered_logistic( eta_v[n,J-1] + phy_v[n,J-1], c_hunt );
  
  if (storage[n] == -99) {
  target += log_mix( inv_logit( a_storage + eta_v[n,J] + phy_v[n,J] ),
  bernoulli_logit_lpmf( 1 | a_storage + eta_v[n,J] + phy_v[n,J] ),
  bernoulli_logit_lpmf( 0 | a_storage + eta_v[n,J] + phy_v[n,J] )
  );
  }
  
  else
  storage[n] ~ bernoulli( inv_logit( a_storage + eta_v[n,J] + phy_v[n,J] ));
}
}

generated quantities{
matrix[N,J] log_lik;
matrix[J,J] Rho_eta;
matrix[J,J] Rho_phy;
int store_hat[N]; // food storage, rng for missing

for (i in 1:N) {
  store_hat[i] = storage[i];
  
  if (storage[i] == -99) store_hat[i] = bernoulli_logit_rng(a_storage + eta_v[i,J] + phy_v[i,J]);
  
  for (j in 1:J-2) {
    log_lik[i,j] = ordered_logistic_lpmf( y[i,j] | eta_v[i,j] + phy_v[i,j], c[j,] );
  }
  log_lik[i,J-1] = ordered_logistic_lpmf( hunting[i] | eta_v[i,J-1] + phy_v[i,J-1], c_hunt );
  log_lik[i,J] = bernoulli_logit_lpmf( store_hat[i] | a_storage + eta_v[i,J] + phy_v[i,J] );
}

Rho_phy = L_phy * L_phy';
Rho_eta = L_eta * L_eta';
}
