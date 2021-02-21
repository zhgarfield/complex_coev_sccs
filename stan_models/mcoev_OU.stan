functions {
  real[] dz_dt(real t,       // time
               real[] z,     // system state 
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {  // unused data
                 
    real y1 = z[1];
    real y2 = z[2];
    
    // System of differential equations
    real dy1_dt = theta[13]*(theta[1] + theta[2]*y1 + theta[3]*square(y1) + theta[4]*y2 + theta[5]*square(y2) + theta[6]*y1*y2 - y1);
    
    real dy2_dt = theta[14]*(theta[7] + theta[8]*y2 + theta[9]*square(y2) + theta[10]*y1 + theta[11]*square(y1) + theta[12]*y1*y2 - y2);
    
    return { dy1_dt, dy2_dt };
  }
}

data{
  int N;
  int J;
  int N_seg;
  int node_seq[N_seg];
  int parent[N_seg];
  real ts[N_seg,2];
  int y[N,J-2];
  int storage[N];
  int hunting[N];
  real date[N_seg];
  int K;
  int K_hunt;
  real prior_y[J-2,K-1];
  real prior_storage;
  real prior_hunting[K_hunt-1];
}

parameters{
real<upper=0> z_anc[2]; // ancestral states, constrained to be lower than sample mean (0)
vector[12] B_theta; // interaction of the two latent variables
matrix[N_seg-1,2] z_drift;

real<lower=0> alpha[2]; // strength of selection
real<lower=0> sigma[2]; // strength of drift
ordered[K-1] c[J-2];  // cutpoints for ordered logit
ordered[K_hunt-1] c_hunt; // cutpoints for hunting
real a_storage; // intercept for food storage
vector<lower=0>[J-1] lambda; // loadings of latent OU variables onto observed, except agriculture which is fixed = 1
}

transformed parameters{
real z[N_seg,2,2];
real theta[14];

for (p in 1:12) theta[p] = B_theta[p];
theta[13] = alpha[1];
theta[14] = alpha[2];

// setting ancestral states
for (n in 1:2) {
  z[node_seq[1],1,n] = z_anc[1];
  z[node_seq[1],2,n] = z_anc[2];
}

// running ODE over tree sequence
  for (i in 2:N_seg) {
  real starts[2];

  for (n in 1:2) {
  starts[n] = z[parent[i],2,n];
  }
  
  // Solve the system at appropriate time points
  z[node_seq[i],,] = integrate_ode_rk45(dz_dt , starts , 0 , ts[i,] , theta ,  rep_array(0.0, 0), rep_array(0, 0));
  
  for (n in 1:2) {
  // adding gaussian noise, scaled by expected sd of the OU model
  z[node_seq[i],2,n] = z[node_seq[i],2,n] + z_drift[i-1,n]*sigma[n]*ts[i,2];
  }

}
} // end transformed parameters block

model{
  B_theta ~ std_normal();
  to_vector(z_drift) ~ std_normal();
  alpha ~ normal(0,1);
  sigma ~ normal(0,1);
  lambda ~ std_normal();
  
  // priors for ancestral state
  z_anc ~ normal(0,3);
  
  // Priors when latent variables = 0, based on sample means
  for (n in 1:J-2) c[n,] ~ normal(prior_y[n,],1);
  c_hunt ~ normal(prior_hunting,1); // vague prior on hunting
  a_storage ~ normal(prior_storage,1);
  
  for (i in 1:N) {
    // Murdock complexity variables
  for (j in 1:J-2) {
    if (j == 1) y[i,j] ~ ordered_logistic( z[i,2,1], c[j,] );
    else if (j <= 4) y[i,j] ~ ordered_logistic( z[i,2,1]*lambda[j-1], c[j,] );
    else y[i,j] ~ ordered_logistic( z[i,2,2]*lambda[j-1], c[j,] );
  }
  // Hunting and food storage
  hunting[i] ~ ordered_logistic( z[i,2,1]*lambda[J-2], c_hunt );
  
  if (storage[i] == -99) {
  target += log_mix( inv_logit( a_storage + z[i,2,1]*lambda[J-1] ),
  bernoulli_logit_lpmf( 1 | a_storage + z[i,2,1]*lambda[J-1] ),
  bernoulli_logit_lpmf( 0 | a_storage + z[i,2,1]*lambda[J-1] )
  );
  }
  
  else
  storage[i] ~ bernoulli_logit( a_storage + z[i,2,1]*lambda[J-1] );
}
}
