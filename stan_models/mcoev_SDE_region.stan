functions { 
matrix kronecker_prod(matrix A, matrix B) { // returns the Kronecker Product
  matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
  int m;
  int n;
  int p;
  int q;
  m = rows(A);
  n = cols(A);
  p = rows(B);
  q = cols(B);
  for (i in 1:m)
    for (j in 1:n)
      for (k in 1:p)
        for (l in 1:q)
      C[p*(i-1)+k,q*(j-1)+l] = A[i,j]*B[k,l];
  return C;
}

matrix A_dt(matrix A, real t) {  // expected auto and cross effects over a discrete time t
 return( matrix_exp(A * t) );
}

matrix A_sharp(matrix A) {
  matrix[rows(A) * rows(A), cols(A) * cols(A)] A_temp;
  matrix[rows(A),cols(A)] I; // identity matrix
  
  I = diag_matrix(rep_vector(1,rows(A)));

  A_temp = kronecker_prod(A,I) + kronecker_prod(I,A);
  return(A_temp);
}

matrix cov_drift(matrix A, matrix Q, real ts) {
  matrix[rows(A) * rows(A), cols(A) * cols(A)] A_sharp_temp;
  matrix[rows(A) * rows(A), cols(A) * cols(A)] I; // identity matrix
  vector[rows(Q)*cols(Q)] row_Q;
  vector[rows(A)*cols(A)] irow_vec;
  matrix[rows(A),cols(A)] irow_mat;
  
  I = diag_matrix(rep_vector(1,rows(A_sharp_temp)));
  
  A_sharp_temp = A_sharp(A);
  
  // row operation takes elements of a matrix rowwise and puts them into a column vector
  for (i in 1:rows(Q))
  for (j in 1:cols(Q)) {
    row_Q[i + (j-1)*rows(Q)] = Q[j,i];
  }
  
  irow_vec = inverse(A_sharp_temp) * (matrix_exp(A_sharp_temp * ts) - I) * row_Q;
  
  // irow takes elements of a column vector and puts them in a matrix rowwise
  {
  int row_size = rows(A);
  int row_ticker = 1;
  int col_ticker = 0;
  
  for (i in 1:num_elements(irow_vec)) {
    col_ticker += 1;
    if (col_ticker > row_size) {
      row_ticker += 1;
      col_ticker = 1;
    }
    irow_mat[row_ticker,col_ticker] = irow_vec[i];
  }
  }
  return(irow_mat);
}

} // end of functions block


data{
  int N;
  int J;
  int N_region;
  int N_soc[N_region];
  int N_seg[N_region];
  int node_seq[sum(N_seg)];
  int parent[sum(N_seg)];
  vector[sum(N_seg)] ts;
  int soc_order[N];
  int y[N,J-2];
  int storage[N];
  int hunting[N];
  int K;
  int K_hunt;
}

parameters{
// selection coefficients (2 for auto-regressive, 2 for cross-effects)
vector[2] alpha_auto_mu; // <0, assumes mean reverting process
vector[2] alpha_cross_mu;

vector[2] sigma_mu; // drift, log-scale

vector[2] b_mu; // SDE intercepts
vector[2] eta_anc_mu; // ancestral states
matrix[sum(N_seg),2] z_drift; // stochastic drift, unscaled and uncorrelated

matrix[10,N_region] region_z; // region-specific random effects, unscaled and uncorrelated
vector<lower=0>[10] sigma_region; // scale parameter for random effects
cholesky_factor_corr[10] L_region; // correlations between parameters

ordered[K-1] c[J-2];  // cutpoints for ordered logit
ordered[K_hunt-1] c_hunt; // cutpoints for hunting
real a_storage; // intercept for food storage
vector<lower=0>[J-1] loadings; // loadings of latent OU variables onto observed, except agriculture which is fixed = 1
}

transformed parameters{
matrix[sum(N_seg),2] eta;
vector[J] lambda; // factor loadings
matrix[N_region,10] region_v;
matrix[2,2] I; 

region_v = (diag_pre_multiply(sigma_region, L_region) * region_z)';
lambda = to_vector( append_col(1, to_row_vector(loadings)) );

// identity matrix
I = diag_matrix(rep_vector(1,2));

for (r in 1:N_region) {
matrix[2,2] Q;
matrix[2,2] A;
vector[2] sigmasq;
vector[2] b;
int start_i;
int end_i;
matrix[N_seg[r],2] eta_temp;

// selection matrix, col 1 is RI, col 2 is TPC
A[1,1] = -exp(alpha_auto_mu[1] + region_v[r,1]); // autoregressive effect of RI on itself
A[2,2] = -exp(alpha_auto_mu[2] + region_v[r,2]); // autoregressive effect of TPC on itself

A[2,1] = alpha_cross_mu[1] + region_v[r,3];
A[1,2] = alpha_cross_mu[2] + region_v[r,4];

// drift matrix
sigmasq[1] = exp( sigma_mu[1] + region_v[r,5] );
sigmasq[2] = exp( sigma_mu[2] + region_v[r,6] );
Q = diag_matrix( square(sigmasq) );

// SDE intercepts
b[1] = b_mu[1] + region_v[r,7];
b[2] = b_mu[2] + region_v[r,8];

// Where to start the index, depending on region
if (r == 1) {
  start_i = 1;
  end_i = N_seg[r];
}

if (r > 1) {
  start_i = sum(N_seg[1:(r - 1)]) + 1;
  end_i = sum(N_seg[1:(r - 1)]) + N_seg[r];
}

// setting ancestral states
eta_temp[node_seq[start_i],1] = eta_anc_mu[1] + region_v[r,9]; // RI ancestral state
eta_temp[node_seq[start_i],2] = eta_anc_mu[2] + region_v[r,10]; // TPC ancesrtal state

  for (i in 1:(N_seg[r]-1)) {
  matrix[2,2] A_delta; // amount of deterministic change ("selection")
  matrix[2,2] VCV; // variance-covariance matrix of stochastic change ("drift")
  vector[2] drift_seg; // accumulated drift over the segment
  
  A_delta = A_dt(A, ts[start_i + i]);
  VCV = cov_drift(A, Q, ts[start_i + i]);
  
  // No drift on the interaction, bc its simply a product of RI and TPC
  drift_seg = cholesky_decompose(VCV) * to_vector( z_drift[start_i + i,] );
  
  eta_temp[node_seq[start_i + i],] = to_row_vector(
    A_delta * to_vector(eta_temp[parent[start_i + i],]) + (inverse(A) * (A_delta - I) * b) + drift_seg
  );
  }
  
  eta[start_i:end_i,] = eta_temp;
}

} // end transformed parameters block

model{
  alpha_auto_mu ~ std_normal();
  alpha_cross_mu ~ std_normal();
  to_vector(region_z) ~ std_normal();
  b_mu ~ std_normal();
  sigma_mu ~ std_normal();
  sigma_region ~ std_normal();
  L_region ~ lkj_corr_cholesky(2);
  eta_anc_mu ~ std_normal();
  to_vector(z_drift) ~ std_normal();
  loadings ~ std_normal();
  
  // Priors when latent variables = 0
  for (n in 1:J-2) c[n,] ~ normal(0,2);
  c_hunt ~ normal(0,2);
  a_storage ~ normal(0,2);
  
  // Model likelihood /////////////
  for (i in 1:N) {
    // Murdock complexity variables
  for (j in 1:J-2) {
    // RI variables
      if (j <= 4) y[i,j] ~ ordered_logistic( eta[soc_order[i],1]*lambda[j], c[j,] );
    // TPC variables
     else y[i,j] ~ ordered_logistic( eta[soc_order[i],2]*lambda[j], c[j,] );
  }
  // Hunting and food storage
    hunting[i] ~ ordered_logistic( eta[soc_order[i],1]*lambda[J-1], c_hunt );
  
  // Mix over missingness in food storage
  if (storage[i] == -99) {
    target += log_mix( inv_logit( a_storage + eta[soc_order[i],1]*lambda[J] ),
    bernoulli_logit_lpmf( 1 | a_storage + eta[soc_order[i],1]*lambda[J] ),
    bernoulli_logit_lpmf( 0 | a_storage + eta[soc_order[i],1]*lambda[J] )
    );
  }
    else
    storage[i] ~ bernoulli_logit( a_storage + eta[soc_order[i],1]*lambda[J] );
}
}
