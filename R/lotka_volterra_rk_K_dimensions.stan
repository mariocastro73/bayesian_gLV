data {
  int<lower=1> K;  // Number of species
  int<lower=0> N;          // number of data points
  array[N,K] real<lower=0> Y;   // observed data
  array[N] real<lower=0> T;      // time points
}

parameters {
  array [K] real<lower=0> r;               // intrinsic growth rates
  array[K,K]real b;            // interaction matrix
  array[K] real<lower=0> sigma;
}

model {
  for (n in 2:N) {
    real dt = T[n] - T[n-1];
    array[K] real k1;
    array[K] real k2;
    array[K] real k3;
    array[K] real k4;
    array[K] real y_temp;
    
    // Priors
    sigma ~ exponential(1); // Signed vague prior
    r ~ exponential(0.1); // Signed vague prior
    for (i in 1:K) {
      for (j in 1:K) {
        b[i, j] ~ normal(0,2); // Vague prior
      }
    }
    // Calculate k1
    for (i in 1:K) {
      k1[i] = r[i] * Y[n-1, i];
      for (j in 1:K) {
        k1[i] += b[i, j] * Y[n-1, i] * Y[n-1, j];
      }
    }
    
    // Calculate k2
    for (i in 1:K) {
      y_temp[i] = Y[n-1, i] + 0.5 * dt * k1[i];
    }
    for (i in 1:K) {
      k2[i] = r[i] * y_temp[i];
      for (j in 1:K) {
        k2[i] += b[i, j] * y_temp[i] * y_temp[j];
      }
    }
    
    // Calculate k3
    for (i in 1:K) {
      y_temp[i] = Y[n-1, i] + 0.5 * dt * k2[i];
    }
    for (i in 1:K) {
      k3[i] = r[i] * y_temp[i];
      for (j in 1:K) {
        k3[i] += b[i, j] * y_temp[i] * y_temp[j];
      }
    }
    
    // Calculate k4
    for (i in 1:K) {
      y_temp[i] = Y[n-1, i] + dt * k3[i];
    }
    for (i in 1:K) {
      k4[i] = r[i] * y_temp[i];
      for (j in 1:K) {
        k4[i] += b[i, j] * y_temp[i] * y_temp[j];
      }
    }
    
    // Update Y using RK4
    for (i in 1:K) {
      real dPop = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
      // Y[n, i] ~ normal(Y[n-1, i] + dPop * dt, sigma[i]); // Assuming a noise level
      Y[n, i] ~ lognormal(log(Y[n-1, i] + dPop * dt), sigma[i]); // Assuming a noise level
    }
    
  }
}

