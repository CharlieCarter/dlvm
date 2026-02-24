// dlvm_student_t.stan — Dynamic Latent Variable Model (Student-t observation family)
//
// Observation model: Student-t likelihood for continuous data.
//   y[i] ~ student_t(nu_obs, theta[s(i)], sigma_eff[i])
//
// This is the default observation family for continuous, unbounded data.
// The Student-t likelihood provides robustness to outlier observations that
// are far from the latent theta, preventing single extreme values from
// unduly influencing the latent state.  The degrees of freedom parameter
// nu_obs (passed as data) controls the tail heaviness:
//   - nu_obs = 4 (default): strong outlier robustness
//   - nu_obs >= 30: effectively Gaussian
//
// NOTE on the two Student-t distributions in this model:
//   nu_obs   → controls the OBSERVATION distribution's tails (likelihood).
//              Determines how much a single extreme y[i] can pull theta.
//   nu_state → controls the INNOVATION distribution's tails (state process).
//              Determines how large structural breaks in the latent trajectory
//              can be.  This is the "robust" property from Reuning et al.
//   These are orthogonal: nu_obs is about measurement noise robustness,
//   nu_state is about structural break accommodation.

functions {
  real partial_sum_lpmf(array[] int slice_obs,
                        int start, int end,
                        vector y,
                        vector mu,
                        vector sigma_eff,
                        real nu_obs) {
    return student_t_lpdf(y[start:end] | nu_obs, mu[start:end], sigma_eff[start:end]);
  }
}

data {
  int<lower=1> n_states;
  int<lower=1> n_obs;
  int<lower=1> n_units;

  vector[n_obs] y;
  array[n_obs] int<lower=1,upper=n_states> obs_to_state;

  int<lower=0,upper=1> has_se;
  vector<lower=0>[n_obs] se;

  array[n_obs] int<lower=1> obs_n;

  array[n_states] int<lower=0,upper=n_states> state_prev;
  vector<lower=0>[n_states] delta_t;

  real<lower=1> nu_obs;
  real<lower=1> nu_state;
  real<lower=0> scale_state;
  int<lower=1> grainsize;
  int<lower=0,upper=1> compute_gq;
}

transformed data {
  vector[n_states] state_nu;
  array[n_obs] int seq_obs;
  vector[n_obs] sqrt_obs_n;

  for (s in 1:n_states) {
    state_nu[s] = (state_prev[s] == 0) ? 4.0 : nu_state;
  }

  for (i in 1:n_obs) {
    seq_obs[i] = i;
    sqrt_obs_n[i] = sqrt(obs_n[i]);
  }

  int n_gq = compute_gq ? n_obs : 0;
}

parameters {
  vector[n_states] theta_raw;
  real<lower=0> innov;
  real<lower=0> sigma;
}

transformed parameters {
  vector[n_states] theta;

  for (s in 1:n_states) {
    if (state_prev[s] == 0) {
      theta[s] = theta_raw[s];
    } else {
      theta[s] = theta[state_prev[s]] + theta_raw[s] * innov * sqrt(delta_t[s]);
    }
  }
}

model {
  theta_raw ~ student_t(state_nu, 0, scale_state);
  innov ~ normal(0, 1);
  sigma ~ normal(0, 1);

  vector[n_obs] mu;
  vector[n_obs] sigma_eff;

  for (i in 1:n_obs) {
    mu[i] = theta[obs_to_state[i]];
  }

  if (has_se) {
    for (i in 1:n_obs) {
      sigma_eff[i] = sqrt(square(sigma) / obs_n[i] + square(se[i]));
    }
  } else {
    for (i in 1:n_obs) {
      sigma_eff[i] = sigma / sqrt_obs_n[i];
    }
  }

  target += reduce_sum(partial_sum_lpmf,
                       seq_obs,
                       grainsize,
                       y, mu, sigma_eff, nu_obs);
}

generated quantities {
  vector[n_gq] log_lik;
  vector[n_gq] y_rep;

  if (compute_gq) {
    for (i in 1:n_obs) {
      real mu_i = theta[obs_to_state[i]];
      real sigma_eff_i;

      if (has_se) {
        sigma_eff_i = sqrt(square(sigma) / obs_n[i] + square(se[i]));
      } else {
        sigma_eff_i = sigma / sqrt_obs_n[i];
      }

      log_lik[i] = student_t_lpdf(y[i] | nu_obs, mu_i, sigma_eff_i);
      y_rep[i] = student_t_rng(nu_obs, mu_i, sigma_eff_i);
    }
  }
}
