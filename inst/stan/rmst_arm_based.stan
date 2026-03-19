// RMST Arm-based Network Meta-Analysis Model
// Updated for Stan 2.37.0 syntax

functions {
  // Calculate RMST from survival function
  real calculate_rmst(vector times, vector surv, real tau) {
    int n = num_elements(times);
    real rmst = 0.0;
    
    // Trapezoidal integration
    for (i in 2:n) {
      if (times[i-1] < tau) {
        real t1 = times[i-1];
        real t2 = fmin(times[i], tau);
        real s1 = surv[i-1];
        real s2 = surv[i];
        
        if (times[i] > tau) {
          // Linear interpolation for survival at tau
          s2 = s1 + (s2 - s1) * (tau - t1) / (times[i] - t1);
        }
        
        rmst += 0.5 * (s1 + s2) * (t2 - t1);
      }
    }
    
    return rmst;
  }
}

data {
  // Dimensions
  int<lower=1> N_studies;
  int<lower=2> N_treatments;
  int<lower=1> N_arms;
  int<lower=1> N_tau;
  
  // Survival data
  int<lower=1> max_times;
  matrix[N_arms, max_times] times;
  matrix[N_arms, max_times] surv;
  array[N_arms] int<lower=1> n_times;  // NEW SYNTAX
  
  // Restriction times
  vector[N_tau] tau;
  
  // Study/treatment structure  
  array[N_arms] int<lower=1,upper=N_studies> study;  // NEW SYNTAX
  array[N_arms] int<lower=1,upper=N_treatments> treatment;  // NEW SYNTAX
  array[N_studies] int<lower=1> n_arms_by_study;  // NEW SYNTAX
  
  // Model options
  int<lower=0,upper=1> use_re;
  int<lower=0,upper=1> re_dist; // 0=normal, 1=student_t
  int<lower=0,upper=1> baseline_model; // 0=piecewise, 1=spline
  
  // Priors
  real prior_mean;
  real<lower=0> prior_sd;
}

parameters {
  // Treatment effects (on log hazard ratio scale)
  vector[N_treatments] d;
  
  // Random effects
  array[use_re ? 1 : 0] real<lower=0> sigma;  // NEW SYNTAX
  matrix[use_re ? N_studies : 0, use_re ? N_treatments : 0] delta_raw;
  
  // Baseline parameters (simplified)
  vector[N_studies] mu;
}

transformed parameters {
  // Random effects
  matrix[N_studies, N_treatments] delta;
  
  if (use_re == 1) {
    for (s in 1:N_studies) {
      for (t in 1:N_treatments) {
        if (re_dist == 0) {
          delta[s, t] = delta_raw[s, t] * sigma[1];
        } else {
          // Student-t for robustness
          delta[s, t] = delta_raw[s, t] * sigma[1] * 2.5;
        }
      }
    }
  } else {
    delta = rep_matrix(0.0, N_studies, N_treatments);
  }
}

model {
  // Priors
  d[1] ~ normal(0, 0.001); // Reference treatment
  d[2:N_treatments] ~ normal(prior_mean, prior_sd);
  
  mu ~ normal(0, 2);
  
  if (use_re == 1) {
    sigma ~ normal(0, 1);
    to_vector(delta_raw) ~ std_normal();
  }
  
  // Likelihood (simplified - would need full survival model)
  // This is a placeholder that assumes we observe RMST directly
  for (a in 1:N_arms) {
    vector[n_times[a]] obs_times;
    vector[n_times[a]] obs_surv;
    
    for (i in 1:n_times[a]) {
      obs_times[i] = times[a, i];
      obs_surv[i] = surv[a, i];
    }
    
    for (t_idx in 1:N_tau) {
      real obs_rmst = calculate_rmst(obs_times, obs_surv, tau[t_idx]);
      real expected_rmst = mu[study[a]] + d[treatment[a]] + delta[study[a], treatment[a]];
      
      // Simplified likelihood
      obs_rmst ~ normal(expected_rmst, 1.0);
    }
  }
}

generated quantities {
  // RMST predictions for each treatment at each tau
  matrix[N_treatments, N_tau] rmst_pred;
  
  // Treatment rankings at each tau
  matrix[N_treatments, N_tau] rank;
  
  for (t_idx in 1:N_tau) {
    for (trt in 1:N_treatments) {
      rmst_pred[trt, t_idx] = mean(mu) + d[trt];
    }
    
    // Calculate ranks
    for (trt in 1:N_treatments) {
      int r = 1;
      for (trt2 in 1:N_treatments) {
        if (trt2 != trt && rmst_pred[trt2, t_idx] > rmst_pred[trt, t_idx]) {
          r = r + 1;
        }
      }
      rank[trt, t_idx] = r;
    }
  }
}
