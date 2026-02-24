// dlvm_bernoulli.stan — Dynamic Latent Variable Model (Bernoulli observation family)
//
// Observation model: Bernoulli logit for binary (0/1) data.
//   y[i] ~ Bernoulli(inv_logit(theta[s(i)]))
//
// Use this family for individual-level binary outcomes.  The latent theta
// is on the log-odds scale: theta = 0 maps to Pr(y=1) = 0.5.
//
// There is no sigma parameter — the Bernoulli distribution has fixed
// variance conditional on the probability.
//
// For aggregated binary data (e.g., "15 of 40 articles"), use the
// binomial family instead, which is more efficient and informative.

functions {
  real partial_sum_lpmf(array[] int y_slice,
                        int start, int end,
                        array[] int y,
                        vector mu) {
    return bernoulli_logit_lpmf(y[start:end] | mu[start:end]);
  }
}

data {
  int<lower=1> n_states;
  int<lower=1> n_obs;
  int<lower=1> n_units;

  // --- Binary observations ---
  array[n_obs] int<lower=0,upper=1> y;
  array[n_obs] int<lower=1,upper=n_states> obs_to_state;

  // --- State structure ---
  array[n_states] int<lower=0,upper=n_states> state_prev;
  vector<lower=0>[n_states] delta_t;

  // --- Distribution parameters ---
  real<lower=1> nu_state;
  real<lower=0> scale_state;
  int<lower=1> grainsize;
  int<lower=0,upper=1> compute_gq;
}

transformed data {
  vector[n_states] state_nu;
  array[n_obs] int seq_obs;

  for (s in 1:n_states) {
    state_nu[s] = (state_prev[s] == 0) ? 4.0 : nu_state;
  }

  for (i in 1:n_obs) {
    seq_obs[i] = i;
  }

  int n_gq = compute_gq ? n_obs : 0;
}

parameters {
  vector[n_states] theta_raw;
  real<lower=0> innov;
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

  // --- Observation model (parallelised) ---
  vector[n_obs] mu;
  for (i in 1:n_obs) {
    mu[i] = theta[obs_to_state[i]];
  }

  target += reduce_sum(partial_sum_lpmf,
                       seq_obs,
                       grainsize,
                       y, mu);
}

generated quantities {
  vector[n_gq] log_lik;
  array[n_gq] int y_rep;

  if (compute_gq) {
    for (i in 1:n_obs) {
      real mu_i = theta[obs_to_state[i]];
      log_lik[i] = bernoulli_logit_lpmf(y[i] | mu_i);
      y_rep[i] = bernoulli_logit_rng(mu_i);
    }
  }
}
