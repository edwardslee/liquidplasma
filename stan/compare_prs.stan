data{
  int N;
  vector[N] prs;
  array[N] int product;
}
parameters{
  vector[2] mu;
  vector<lower=0>[2] sigma;
}
model{
  sigma ~ exponential(0.2);
  mu ~ normal(90, 20);
  for (i in 1:N) {
    prs[i] ~ student_t(2, mu[product[i]], sigma[product[i]]);  
  }
  
}
generated quantities{
  vector[N] log_lik;
  real mu_diff;
  
  for (i in 1:N) log_lik[i] = student_t_lpdf(prs[i] | 2 , mu[product[i]], sigma[product[i]]);
  mu_diff = mu[1] - mu[2];
}
