functions {
  real[] dz_dt(real t,       // time
               real[] z,     // system state 
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {  // unused data
    
    real dydt[2];             
    real y1 = z[1];
    real y2 = z[2];
    
    // System of differential equations
    dydt[1] = theta[3]*(theta[1] + theta[5]*(y1 - theta[1]) + theta[6]*square(y1 - theta[1]) + theta[7]*(y2 - theta[2]) + theta[8]*square(y2 - theta[2]) + theta[9]*(y1 - theta[1])*(y2 - theta[2]) - y1);
    
    dydt[2] = theta[4]*(theta[2] + theta[10]*(y2 - theta[2]) + theta[11]*square(y2 - theta[2]) + theta[12]*(y1 - theta[1]) + theta[13]*square(y1 - theta[1]) + theta[14]*(y1 - theta[1])*(y2 - theta[2]) - y2);
    
    return dydt;
  }
}

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
  real date[N_seg];
  int K;
  int K_hunt;
  real prior_y[J-2,K-1];
  real prior_storage;
  real prior_hunting[K_hunt-1];
}

parameters{
vector<upper=0>[2] theta0; // ancestral equilibrium, putting an informative prior on <0
vector[10] B_theta; // interaction of the two latent variables
matrix[N_seg-1,2] z_drift;

real<lower=0> alpha[2]; // strength of selection
real<lower=0> sigma[2]; // strength of drift
ordered[K-1] c[J-2];  // cutpoints for ordered logit
ordered[K_hunt-1] c_hunt; // cutpoints for hunting
real a_storage; // intercept for food storage
vector<lower=0>[J-1] lambda; // loadings of latent OU variables onto observed, except agriculture which is fixed = 1
}

transformed parameters{
real z[N_seg,2];
real theta[14];

// Pack in parameters for ode solver into a single vector theta
for (p in 1:2) theta[p] = theta0[p];
theta[3] = alpha[1];
theta[4] = alpha[2];
for (p in 5:14) theta[p] = B_theta[p-4];

// setting ancestral states
z[node_seq[1],1] = theta[1];
z[node_seq[1],2] = theta[2];

// running ODE over tree sequence
  for (i in 2:N_seg) {
  real starts[2];
  real z_temp[1,2];
  vector[1] time_length;

  time_length[1] = ts[i];

  for (n in 1:2) {
  starts[n] = z[parent[i],n];
  }
  
  // Solve the system at appropriate time points
  z_temp = integrate_ode_rk45(dz_dt, starts , 0.0 , to_array_1d(time_length), theta ,  rep_array(0.0, 0), rep_array(0, 0));
  z[node_seq[i],] = z_temp[1,];
  
  for (n in 1:2) {
  // adding gaussian noise, scaled by expected sd of the OU process
  real X_SD = sqrt( (square(sigma[n])/(2*alpha[n])) * (1 - exp(-2*alpha[n]*ts[i])) );
  z[node_seq[i],n] = z[node_seq[i],n] + z_drift[i-1,n]*X_SD;
  }

}
} // end transformed parameters block

model{
  B_theta ~ std_normal(); // std_normal() is a more efficient version of normal(0,1)
  theta0 ~ std_normal();
  to_vector(z_drift) ~ std_normal();
  alpha ~ std_normal();
  sigma ~ std_normal();
  lambda ~ std_normal();
  
  // Priors when latent variables = 0, based on sample means
  for (n in 1:J-2) c[n,] ~ normal(prior_y[n,],1);
  c_hunt ~ normal(prior_hunting,1); // vague prior on hunting
  a_storage ~ normal(prior_storage,1);
  
  for (i in 1:N) {
    // Murdock complexity variables
  for (j in 1:J-2) {
    if (j == 1) y[i,j] ~ ordered_logistic( z[i,1], c[j,] );
    else if (j <= 4) y[i,j] ~ ordered_logistic( z[i,1]*lambda[j-1], c[j,] );
    else y[i,j] ~ ordered_logistic( z[i,2]*lambda[j-1], c[j,] );
  }
  // Hunting and food storage
  hunting[i] ~ ordered_logistic( z[i,1]*lambda[J-2], c_hunt );
  
  if (storage[i] == -99) {
  target += log_mix( inv_logit( a_storage + z[i,1]*lambda[J-1] ),
  bernoulli_logit_lpmf( 1 | a_storage + z[i,1]*lambda[J-1] ),
  bernoulli_logit_lpmf( 0 | a_storage + z[i,1]*lambda[J-1] )
  );
  }
  
  else
  storage[i] ~ bernoulli_logit( a_storage + z[i,1]*lambda[J-1] );
}
}
