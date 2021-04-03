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
  int N_seg;
  int node_seq[N_seg];
  int parent[N_seg];
  vector[N_seg] ts;
  int y[N,J-2];
  int storage[N];
  int hunting[N];
  int K;
  int K_hunt;
  matrix[N,N] geo_dist;
}

parameters{
// selection coefficients (2 for auto-regressive, 2 for cross-effects)
vector[2] alpha_auto; // <0, assumes mean reverting process
vector[2] alpha_cross; 

vector<lower=0>[2] sigma; // drift scale
vector[2] b; // SDE intercepts
vector[2] eta_anc; // ancestral states
matrix[N_seg - 1,2] z_drift; // stochastic drift, unscaled and uncorrelated

ordered[K-1] c[J-2];  // cutpoints for ordered logit
ordered[K_hunt-1] c_hunt; // cutpoints for hunting
real a_storage; // intercept for food storage
vector<lower=0>[J-1] loadings; // loadings of latent OU variables onto observed, except agriculture which is fixed = 1

matrix[N,2] geo_z; // geographic covariance random effects, unscaled and uncorrelated
vector<lower=0>[2]rho_geo; // how quickly does geographic covariance decline with distance
vector<lower=0>[2]sigma_geo; // maximum covariance due to geographic loc
}

transformed parameters{
matrix[N_seg,2] eta;
matrix[2,2] Q;
matrix[2,2] I; 
matrix[2,2] A;
vector[J] lambda; // factor loadings
matrix[N,2] geo_v; // geographic covariance random effects

lambda = to_vector( append_col(1, to_row_vector(loadings)) );

// selection matrix, col 1 is RI, col 2 is TSD
A[1,1] = -exp(alpha_auto[1]); // autoregressive effect of RI on itself
A[2,2] = -exp(alpha_auto[2]); // autoregressive effect of TSD on itself

A[2,1] = alpha_cross[1];
A[1,2] = alpha_cross[2];

// drift matrix
Q = diag_matrix(square(sigma));

// identity matrix
I = diag_matrix(rep_vector(1,2));

// setting ancestral states
eta[node_seq[1],1] = eta_anc[1]; // RI ancestral state
eta[node_seq[1],2] = eta_anc[2]; // TSD ancesrtal state

  for (i in 2:N_seg) {
  matrix[2,2] A_delta; // amount of deterministic change ("selection")
  matrix[2,2] VCV; // variance-covariance matrix of stochastic change ("drift")
  vector[2] drift_seg; // accumulated drift over the segment
  
  A_delta = A_dt(A, ts[i]);
  VCV = cov_drift(A, Q, ts[i]);
  
  // No drift on the interaction, bc its simply a product of RI and TSD
  drift_seg = cholesky_decompose(VCV) * to_vector( z_drift[i-1,] );
  
  eta[node_seq[i],] = to_row_vector(
    A_delta * to_vector(eta[parent[i],]) + (inverse(A) * (A_delta - I) * b) + drift_seg
  );
  }
  
// geographic covariance functions
for (j in 1:2) {
  matrix[N,N] geo_cov;
  matrix[N,N] L_geo_cov;
  
      for ( i in 1:(N-1) )
      for ( m in (i+1):N ) {
        geo_cov[i,m] = sigma_geo[j]*exp(-(rho_geo[j]*geo_dist[i,m]));
        geo_cov[m,i] = geo_cov[i,m];
      }
    for ( q in 1:N )
    geo_cov[q,q] = sigma_geo[j] + 0.01;
    
    L_geo_cov = cholesky_decompose(geo_cov);
    geo_v[,j] = L_geo_cov * geo_z[,j];
}
  

} // end transformed parameters block

model{
  alpha_auto ~ std_normal();
  alpha_cross ~ std_normal();
  b ~ std_normal();
  sigma ~ std_normal();
  eta_anc ~ std_normal();
  to_vector(z_drift) ~ std_normal();
  loadings ~ std_normal();
  
  to_vector(geo_z) ~ std_normal();
  sigma_geo ~ exponential(1);
  rho_geo ~ exponential(1);
  
  // Priors when latent variables = 0, based on sample means
  for (n in 1:J-2) c[n,] ~ normal(0,2);
  c_hunt ~ normal(0,2); // vague prior on hunting
  a_storage ~ normal(0,2);
  
  for (i in 1:N) {
    // Murdock complexity variables
  for (j in 1:J-2) {
    // RI variables
    if (j <= 4) y[i,j] ~ ordered_logistic( (eta[i,1] + geo_v[i,1])*lambda[j], c[j,] );
    // TSD variables
    else y[i,j] ~ ordered_logistic( (eta[i,2] + geo_v[i,2])*lambda[j], c[j,] );
  }
  // Hunting and food storage
  hunting[i] ~ ordered_logistic( (eta[i,1] + geo_v[i,1])*lambda[J-1], c_hunt );
  
  // Mix over missingness in food storage
  if (storage[i] == -99) {
  target += log_mix( inv_logit( a_storage + (eta[i,1] + geo_v[i,1])*lambda[J] ),
  bernoulli_logit_lpmf( 1 | a_storage + (eta[i,1] + geo_v[i,1])*lambda[J] ),
  bernoulli_logit_lpmf( 0 | a_storage + (eta[i,1] + geo_v[i,1])*lambda[J] )
  );
  }
  
  else
  storage[i] ~ bernoulli_logit( a_storage + (eta[i,1] + geo_v[i,1])*lambda[J] );
}
}
