// dlvm_cumulative.stan — Dynamic Latent Variable Model (Cumulative/Ordinal observation family)
//
// Observation model: Ordered logistic (cumulative logit) for ordinal data.
//   Pr(y[i] = k) = inv_logit(c[k] - theta[s(i)]) - inv_logit(c[k-1] - theta[s(i)])
//
// Use this family for ordinal/categorical response data with K ordered levels
// coded as integers 1..K.  The model estimates K-1 cutpoints that partition
// the latent scale into K probability bins.
//
// There is no sigma parameter — the observation noise is absorbed into the
// logistic link's inherent variance (pi^2/3).
//
// NOTE: reduce_sum is not used here because ordered_logistic requires
// passing an `ordered` vector of cutpoints, which cannot be used as a
// function argument in Stan's type system.  The loop is fast in practice
// since ordered_logistic_lpmf is O(1) per observation.

data {
  int<lower=1> n_states;
  int<lower=1> n_obs;
  int<lower=1> n_units;

  // --- Ordinal observations ---
  int<lower=2> K;                                    // number of ordinal categories
  array[n_obs] int<lower=1,upper=K> y;               // observed ordinal responses (1..K)
  array[n_obs] int<lower=1,upper=n_states> obs_to_state;

  // --- State structure ---
  array[n_states] int<lower=0,upper=n_states> state_prev;
  vector<lower=0>[n_states] delta_t;

  // --- Distribution parameters ---
  real<lower=1> nu_state;
  real<lower=0> scale_state;
  int<lower=1> grainsize;       // unused but kept for API consistency
  int<lower=0,upper=1> compute_gq;
}

transformed data {
  vector[n_states] state_nu;

  for (s in 1:n_states) {
    state_nu[s] = (state_prev[s] == 0) ? 4.0 : nu_state;
  }

  int n_gq = compute_gq ? n_obs : 0;
}

parameters {
  vector[n_states] theta_raw;
  real<lower=0> innov;
  ordered[K-1] cutpoints;
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
  cutpoints ~ normal(0, 3);

  // --- Observation model (direct loop) ---
  for (i in 1:n_obs) {
    y[i] ~ ordered_logistic(theta[obs_to_state[i]], cutpoints);
  }
}

generated quantities {
  vector[n_gq] log_lik;
  array[n_gq] int y_rep;

  if (compute_gq) {
    for (i in 1:n_obs) {
      real mu_i = theta[obs_to_state[i]];
      log_lik[i] = ordered_logistic_lpmf(y[i] | mu_i, cutpoints);
      y_rep[i] = ordered_logistic_rng(mu_i, cutpoints);
    }
  }
}
