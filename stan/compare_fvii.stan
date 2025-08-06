data{
  vector[50] fvii;
  array[50] int product;
}
parameters{
  vector[2] mu;
  vector<lower=0>[2] sigma;
}
model{
  sigma ~ exponential(0.2);
  mu ~ normal(70, 20);
  for (i in 1:50) {
    fvii[i] ~ student_t(2, mu[product[i]], sigma[product[i]]);  
  }
  
}
generated quantities{
  vector[50] log_lik;
  real mu_diff;
  
  for ( i in 1:50 ) log_lik[i] = student_t_lpdf(fvii[i] | 2 , mu[product[i]], sigma[product[i]]);
  mu_diff = mu[1] - mu[2];
}
