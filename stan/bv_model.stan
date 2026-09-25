// Statistical model from JDG v7.26. Modern Stan array syntax.
data {
  int<lower=0> N;
  int<lower=0> J;
  array[N] int<lower=1, upper=J> subj;
  vector[N] y;
  array[N] int<lower=0, upper=1> censored;
  real log_loq;

  real<lower=0> prior_cv_between;
  real<lower=0> prior_cv_method;
  real<lower=0> prior_cv_method_sd;
  real<lower=0> weight_cv_inter;
  real<lower=0> prior_cv_within;
  real<lower=0> weight_cv_within;

  real<lower=0> mu_mean;

  real<lower=2.1> prior_nu;
  real<lower=0> prior_nu_sd;
}

parameters {
  real mu;
  vector[J] mu_subject;
  real<lower=0> cv_between;

  real log_mu_cv_within;
  real<lower=0> sigma_log_cv_within;
  vector[J] z_cv_within;

  real<lower=0> sigma_method;
  real<lower=2.1> nu_residual;
}

transformed parameters {
  real sigma_between = cv_between;
  vector<lower=0>[J] cv_within_subject;
  for (j in 1:J) {
    cv_within_subject[j] = exp(log_mu_cv_within +
                               sigma_log_cv_within * z_cv_within[j]);
  }
}

model {
  cv_between           ~ normal(prior_cv_between, weight_cv_inter);
  mu                   ~ normal(log(mu_mean), 2);

  log_mu_cv_within     ~ normal(log(prior_cv_within), weight_cv_within);
  sigma_log_cv_within  ~ normal(0, weight_cv_within);
  z_cv_within          ~ normal(0, 1);

  sigma_method         ~ normal(prior_cv_method, prior_cv_method_sd);

  // Røraas hierarchy: normally distributed homeostatic setpoints; adaptive
  // Student-t residual at the sample/observation level.
  nu_residual          ~ normal(prior_nu, prior_nu_sd);
  mu_subject           ~ normal(mu, sigma_between);

  for (i in 1:N) {
    // cv_within_subject and sigma_method are normal-equivalent log-SDs.
    // Stan's Student-t parameter is a scale, not its marginal SD.
    real sigma_obs_ne = sqrt(cv_within_subject[subj[i]]^2 + sigma_method^2);
    real sigma_obs_t = sigma_obs_ne * sqrt((nu_residual - 2) / nu_residual);
    if (censored[i] == 1) {
      target += student_t_lcdf(log_loq | nu_residual,
                               mu_subject[subj[i]], sigma_obs_t);
    } else {
      y[i] ~ student_t(nu_residual, mu_subject[subj[i]], sigma_obs_t);
    }
  }
}

generated quantities {
  vector[J] cv_within_original_scale;
  real dCVP_20; real dCVP_50; real dCVP_80;
  vector[1000] dCVP_sim;

  real mu_cv_within     = exp(log_mu_cv_within + sigma_log_cv_within^2 / 2);
  real sigma_cv_within  = sqrt((exp(sigma_log_cv_within^2) - 1) *
                               exp(2 * log_mu_cv_within +
                                   sigma_log_cv_within^2));

  for (j in 1:J) {
    cv_within_original_scale[j] = sqrt(exp(cv_within_subject[j]^2) - 1);
  }

  for (n in 1:1000) {
    real sim_log_sigma = normal_rng(log_mu_cv_within, sigma_log_cv_within);
    real sim_sigma     = exp(sim_log_sigma);
    dCVP_sim[n]        = sqrt(exp(sim_sigma^2) - 1);
  }

  dCVP_20 = quantile(dCVP_sim, 0.2);
  dCVP_50 = quantile(dCVP_sim, 0.5);
  dCVP_80 = quantile(dCVP_sim, 0.8);
}
