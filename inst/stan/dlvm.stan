// dlvm.stan — Dynamic Robust Latent Variable Model
//
// A generalised state-space model with:
//   - Non-centered student-t random walk (robust innovations)
//   - Continuous-time support via delta_t scaling
//   - Optional known measurement error (SEM)
//   - Observation-count scaling: sigma is per-document noise,
//     means of n documents have noise sigma/sqrt(n)
//   - Dual observation modes: pre-aggregated means OR individual observations
//   - Within-chain parallelisation via reduce_sum
//   - Generated quantities for LOO-CV and posterior predictive checks
//
// Mathematical specification:
//   State:   theta[t] = theta[t-1] + innov * sqrt(delta_t) * z[t]
//            z[t] ~ student_t(nu_state, 0, scale_state)
//   Obs:     y[i] ~ student_t(nu_obs, theta[cy_id[i]], sigma_eff[i])
//            sigma_eff[i] = sqrt(sigma^2/obs_n[i] + se[i]^2)  if SEM known
//            sigma_eff[i] = sigma / sqrt(obs_n[i])             if SEM unknown

functions {
  /**
   * Partial log-likelihood for reduce_sum parallelisation.
   * Evaluates the student-t log-density for a slice of observations.
   *
   * @param slice_obs  Dummy integer array (required by reduce_sum interface)
   * @param start      Start index of the slice (1-indexed)
   * @param end        End index of the slice (1-indexed)
   * @param y          Full vector of observed values
   * @param mu         Full vector of location parameters (thetas mapped to obs)
   * @param sigma_eff  Full vector of effective scale parameters
   * @param nu_obs     Degrees of freedom for observation distribution
   * @return           Log-likelihood contribution from this slice
   */
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
  // --- Dimensions ---
  int<lower=1> n_states;                   // total latent state positions (unit × time-slot)
  int<lower=1> n_obs;                      // number of observations (individual docs or means)
  int<lower=1> n_units;                    // number of units (countries, actors, etc.)

  // --- Observations ---
  vector[n_obs] y;                         // observed values (means or individual scores)
  array[n_obs] int<lower=1,upper=n_states> obs_to_state;  // maps each obs to its latent state

  // --- Measurement error (optional) ---
  int<lower=0,upper=1> has_se;             // 1 if standard errors are provided
  vector<lower=0>[n_obs] se;               // known SEs (zeros if not provided)

  // --- Observation counts per cell (for means-mode sigma/sqrt(n) scaling) ---
  // In individual mode, all values are 1 and sigma_eff == sigma.
  // In means mode, obs_n[i] is the number of raw observations aggregated
  // into datum i, so that sigma_eff = sigma / sqrt(obs_n[i]).
  array[n_obs] int<lower=1> obs_n;

  // --- State structure ---
  array[n_states] int<lower=0,upper=n_states> state_prev;  // predecessor index (0 = no predecessor)
  vector<lower=0>[n_states] delta_t;       // time gap to predecessor (0 for initial states)

  // --- Distribution parameters (passed as data to preserve user's choices) ---
  real<lower=1> nu_obs;                    // df for observation distribution (default: 4)
  real<lower=1> nu_state;                  // df for state innovation distribution (default: 4)
  real<lower=0> scale_state;               // scale for state innovation prior (default: 4)

  // --- Parallelisation tuning ---
  int<lower=1> grainsize;                  // reduce_sum grainsize (1 = auto)

  // --- Optional generated quantities ---
  int<lower=0,upper=1> compute_gq;         // 1 = compute log_lik/y_rep; 0 = skip for speed
}

transformed data {
  // Pre-compute per-state innovation df
  // Initial states (no predecessor) get nu=4 regardless of nu_state
  // (matches the original model's logic for initialising new chains)
  vector[n_states] state_nu;
  array[n_obs] int seq_obs;                // dummy index array for reduce_sum
  vector[n_obs] sqrt_obs_n;                // pre-computed sqrt(obs_n) for efficiency

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
  vector[n_states] theta_raw;             // non-centered innovation draws
  real<lower=0> innov;                    // global innovation scale
  real<lower=0> sigma;                    // per-document observation noise scale
}

transformed parameters {
  // --- Build latent state trajectories ---
  vector[n_states] theta;

  for (s in 1:n_states) {
    if (state_prev[s] == 0) {
      // Initial state for a unit: draw from prior directly
      theta[s] = theta_raw[s];
    } else {
      // Random walk step, scaled by sqrt(delta_t) for continuous time
      theta[s] = theta[state_prev[s]] + theta_raw[s] * innov * sqrt(delta_t[s]);
    }
  }
}

model {
  // --- State process priors ---
  theta_raw ~ student_t(state_nu, 0, scale_state);
  innov ~ normal(0, 1);
  sigma ~ normal(0, 1);

  // --- Observation model (parallelised) ---
  // Map latent states to observation locations
  vector[n_obs] mu;
  vector[n_obs] sigma_eff;

  for (i in 1:n_obs) {
    mu[i] = theta[obs_to_state[i]];
  }

  // Compute effective observation scale: sigma/sqrt(n), with optional SE
  // sigma is per-document noise; the mean of n docs has noise sigma/sqrt(n)
  if (has_se) {
    for (i in 1:n_obs) {
      sigma_eff[i] = sqrt(square(sigma) / obs_n[i] + square(se[i]));
    }
  } else {
    for (i in 1:n_obs) {
      sigma_eff[i] = sigma / sqrt_obs_n[i];
    }
  }

  // Parallelised log-likelihood via reduce_sum
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
