data{
  vector[78] factor_v;
  array[78] int time;
  array[78] int sample_id;
}
parameters{
  vector[26] z1;
  vector[26] z2;
  real<lower=0> sigma_b1;
  real<lower=0> sigma_b2;
  real beta1_bar;
  real beta2_bar;
  real alpha;
  real<lower=0> sigma;
}
model{
  vector[78] mu;
  sigma ~ exponential(1);
  alpha ~ normal(60, 20);
  sigma_b1 ~ exponential(1);
  sigma_b2 ~ exponential(1);
  beta1_bar ~ normal(0, 2);
  beta2_bar ~ normal(0, 2);
  z1 ~ normal(0, 1);
  z2 ~ normal(0, 1);
  for (i in 1:78) {
    mu[i] = alpha + 
      (beta1_bar + z1[sample_id[i]] * sigma_b1) * (time[i] - 11) * step(11 - time[i]) +
      (beta2_bar + z2[sample_id[i]] * sigma_b2) * (time[i] - 11) * step(time[i] - 11);
  }
  factor_v ~ student_t(2, mu, sigma);
}
generated quantities{
  vector[78] log_lik;
  vector[78] mu;
  for (i in 1:78) {
    mu[i] = alpha + 
      (beta1_bar + z1[sample_id[i]] * sigma_b1) * (time[i] - 11) * step(11 - time[i]) +
      (beta2_bar + z2[sample_id[i]] * sigma_b2) * (time[i] - 11) * step(time[i] - 11);
  }
  for (i in 1:78) log_lik[i] = student_t_lpdf(factor_v | 2, mu[i], sigma);
}
