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
  matrix[2,N] phy_z;  // phylogenetic random effects, unscaled and uncorrelated
  vector<lower=0>[2] eta;      // maximum covariance for GP functions
  vector<lower=0>[2] rho;      // rate of decay for GP functions
  vector<lower=0>[6] lambda_TSD; // factor loadings
  vector<lower=0>[5] lambda_RI; // -1 for agr which is fixed to 1 for identifiability
  matrix[2,N] res_z; // uncorrelated and unscaled latent factors
  vector<lower=0>[2] sigma_eta;
  real<lower=0,upper=1> res_cor; // we'll constrain the correlation to be positive for identifiability, and use transformed beta distribution rather than lkj
  real<lower=0,upper=1> phy_cor; //  same but for phylogenetic correlation
}

transformed parameters{
  vector[J] l_TSD; // full vector of loadings that we'll fill with 0's where appropriate
  vector[J] l_RI; // same
  matrix[N,2] phy_v;  // phylogenetic random effects, scaled and correlated
  matrix[N,2] res_v; // correlated latent factors
  matrix[2,2] Rho_res;
  matrix[2,2] Rho_phy;

  // Fill in positive-constrained correlation matrix for latent factors
  Rho_res = diag_matrix(rep_vector(1,2));
  Rho_res[1,2] = res_cor;
  Rho_res[2,1] = Rho_res[1,2];
  
  Rho_phy = diag_matrix(rep_vector(1,2));
  Rho_phy[1,2] = phy_cor;
  Rho_phy[2,1] = Rho_phy[1,2];
  
  l_TSD[1:4] = rep_vector(0.0, 4);
  l_TSD[5:J-2] = lambda_TSD;
  l_TSD[(J-1):J] = rep_vector(0.0, 2);
  
  l_RI[1] = 1;
  l_RI[2:4] = lambda_RI[1:3];
  l_RI[5:(J-2)] = rep_vector(0.0,6);
  l_RI[(J-1):J] = lambda_RI[4:5];
  
  res_v = (diag_pre_multiply(sigma_eta, cholesky_decompose(Rho_res)) * res_z)';
  
  // Between variable phylogenetic correlations
  phy_v = (cholesky_decompose(Rho_phy) * phy_z)';
  
  // Within variable phylo correlations (e.g., between society covariance)
  for (j in 1:2) {
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
    phy_v[,j] = L_phy_cov * phy_v[,j];
  }
  
}

model{
for (j in 1:J-2) {
c[j,] ~ normal(0,2);
}
c_hunt ~ normal(0,2);
a_storage ~ normal(0,2);

to_vector(res_z) ~ std_normal();
sigma_eta ~ exponential(1);

res_cor ~ beta(2,2);
phy_cor ~ beta(2,2);

lambda_TSD ~ std_normal();
lambda_RI ~ std_normal();

to_vector(phy_z) ~ std_normal();
rho ~ exponential(1);
eta ~ exponential(1);

for (j in 1:J-2)
for (i in 1:N) {
y[i,j] ~ ordered_logistic( (res_v[i,1] + phy_v[i,1])*l_RI[j] + (res_v[i,2] + phy_v[i,2])*l_TSD[j], c[j,] );
}

for (i in 1:N) {
  hunting[i] ~ ordered_logistic( (res_v[i,1] + phy_v[i,1])*l_RI[J-1] + (res_v[i,2] + phy_v[i,2])*l_TSD[J-1], c_hunt );
  
  if (storage[i] == -99) {
  target += log_mix( inv_logit( a_storage + (res_v[i,1] + phy_v[i,1])*l_RI[J] + (res_v[i,2] + phy_v[i,2])*l_TSD[J] ),
  bernoulli_logit_lpmf( 1 | a_storage + (res_v[i,1] + phy_v[i,1])*l_RI[J] + (res_v[i,2] + phy_v[i,2])*l_TSD[J]  ),
  bernoulli_logit_lpmf( 0 | a_storage + (res_v[i,1] + phy_v[i,1])*l_RI[J] + (res_v[i,2] + phy_v[i,2])*l_TSD[J]  )
  );
  }
  
  else storage[i] ~ bernoulli_logit( a_storage + (res_v[i,1] + phy_v[i,1])*l_RI[J] + (res_v[i,2] + phy_v[i,2])*l_TSD[J]  );
}
}

generated quantities{
matrix[N,J] log_lik;
int store_hat[N]; // food storage, rng for missing

for (i in 1:N) {
  store_hat[i] = storage[i];
  
  if (storage[i] == -99) {
    store_hat[i] = bernoulli_logit_rng(a_storage + (res_v[i,1] + phy_v[i,1])*l_RI[J] + (res_v[i,2] + phy_v[i,2])*l_TSD[J]);
  }
  
  for (j in 1:J-2) {
    log_lik[i,j] = ordered_logistic_lpmf( y[i,j] | (res_v[i,1] + phy_v[i,1])*l_RI[j] + (res_v[i,2] + phy_v[i,2])*l_TSD[j], c[j,] );
  }
  log_lik[i,J-1] = ordered_logistic_lpmf( hunting[i] | (res_v[i,1] + phy_v[i,1])*l_RI[J-1] + (res_v[i,2] + phy_v[i,2])*l_TSD[J-1], c_hunt );
  log_lik[i,J] = bernoulli_logit_lpmf( store_hat[i] | a_storage + (res_v[i,1] + phy_v[i,1])*l_RI[J] + (res_v[i,2] + phy_v[i,2])*l_TSD[J]  );
}

}
