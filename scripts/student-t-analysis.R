library(tidyverse)
library(rethinking)
library(cmdstanr)
library(readxl)

# import data
df <- read_xlsx("data/liquid_plasma_data.xlsx", na = "NA")


#--------------------
#    fibrinogen 

# filter for fibrinogen data from liquid plasma
df_fibrinogen <- df |>
  filter(product == "liquid plasma",
         factor == "fibrinogen") |>
  select(sample_id, day, factor, value) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day)

# day 15 and day 26 data
df_fibrinogen_15 <- df_fibrinogen |> filter(time == 15)
df_fibrinogen_26 <- df_fibrinogen |> filter(time == 26)

# fit t distribution to day 15 data with fixed nu = 2
fit_fibrinogen_15 <- ulam(
  alist(
    fibrinogen ~ dstudent(2, mu, sigma),
    mu ~ dnorm(280, 50),
    sigma ~ dexp(1)
  ),
  data = df_fibrinogen_15, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# fit t distribution to day 15 data and fitting for nu
fit_fibrinogen_15_nu <- ulam(
  alist(
    fibrinogen ~ dstudent(nu, mu, sigma),
    mu ~ dnorm(280, 50),
    sigma ~ dexp(1),
    nu ~ dexp(0.1)
  ),
  data = df_fibrinogen_15, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# compare two models -- fixed nu has better predictive fit
compare(fit_fibrinogen_15, fit_fibrinogen_15_nu)


# fit t distribution to day 26 data with fixed nu = 2
fit_fibrinogen_26 <- ulam(
  alist(
    fibrinogen ~ dstudent(2, mu, sigma),
    mu ~ dnorm(280, 50),
    sigma ~ dexp(1)
  ),
  data = df_fibrinogen_26, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# fit t distribution to day 26 data and fitting for nu
fit_fibrinogen_26_nu <- ulam(
  alist(
    fibrinogen ~ dstudent(nu, mu, sigma),
    mu ~ dnorm(280, 50),
    sigma ~ dexp(1),
    nu ~ dexp(0.1)
  ),
  data = df_fibrinogen_26, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# compare two models -- fixed nu has better predictive fit
compare(fit_fibrinogen_26, fit_fibrinogen_26_nu)



#--------------------
# protein c

# filter for protein c data from liquid plasma
df_prc <- df |>
  filter(product == "liquid plasma",
         factor == "protein c") |>
  select(sample_id, day, factor, value) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         prc = `protein c`)

# day 15 and day 26 data
df_prc_15 <- df_prc |> filter(time == 15)
df_prc_26 <- df_prc |> filter(time == 26)

# fit t distribution to day 15 data with fixed nu = 2
fit_prc_15 <- ulam(
  alist(
    prc ~ dstudent(2, mu, sigma),
    mu ~ dnorm(100, 20),
    sigma ~ dexp(1)
  ),
  data = df_prc_15, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# fit t distribution to day 15 data and fitting for nu
fit_prc_15_nu <- ulam(
  alist(
    prc ~ dstudent(nu, mu, sigma),
    mu ~ dnorm(100, 20),
    sigma ~ dexp(1),
    nu ~ dexp(0.05)
  ),
  data = df_prc_15, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# compare two models -- fixed nu has better predictive fit
compare(fit_prc_15, fit_prc_15_nu)


# fit t distribution to day 26 data with fixed nu = 2
fit_prc_26 <- ulam(
  alist(
    prc ~ dstudent(2, mu, sigma),
    mu ~ dnorm(100, 20),
    sigma ~ dexp(1)
  ),
  data = df_prc_26, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# fit t distribution to day 26 data and fitting for nu
fit_prc_26_nu <- ulam(
  alist(
    prc ~ dstudent(nu, mu, sigma),
    mu ~ dnorm(100, 20),
    sigma ~ dexp(1),
    nu ~ dexp(0.05)
  ),
  data = df_prc_26, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# compare two models -- fixed nu has better predictive fit
compare(fit_prc_26, fit_prc_26_nu)


#--------------------
# protein s

# filter for protein s data from liquid plasma
df_prs <- df |>
  filter(product == "liquid plasma",
         factor == "protein s") |>
  select(sample_id, day, factor, value) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         prs = `protein s`)

# day 15 and day 26 data
df_prs_15 <- df_prs |> filter(time == 15)
df_prs_26 <- df_prs |> filter(time == 26)

# fit t distribution to day 15 data with fixed nu = 2
fit_prs_15 <- ulam(
  alist(
    prs ~ dstudent(2, mu, sigma),
    mu ~ dnorm(100, 20),
    sigma ~ dexp(1)
  ),
  data = df_prs_15, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# fit t distribution to day 15 data and fitting for nu
fit_prs_15_nu <- ulam(
  alist(
    prs ~ dstudent(nu, mu, sigma),
    mu ~ dnorm(100, 20),
    sigma ~ dexp(1),
    nu ~ dexp(0.5)
  ),
  data = df_prs_15, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# compare two models -- fixed nu has better predictive fit
compare(fit_prs_15, fit_prs_15_nu)


# fit t distribution to day 26 data with fixed nu = 2
fit_prs_26 <- ulam(
  alist(
    prs ~ dstudent(2, mu, sigma),
    mu ~ dnorm(100, 20),
    sigma ~ dexp(1)
  ),
  data = df_prs_26, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# fit t distribution to day 26 data and fitting for nu
fit_prs_26_nu <- ulam(
  alist(
    prs ~ dstudent(nu, mu, sigma),
    mu ~ dnorm(100, 20),
    sigma ~ dexp(1),
    nu ~ dexp(0.5)
  ),
  data = df_prs_26, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# compare two models -- fixed nu has better predictive fit
compare(fit_prs_26, fit_prs_26_nu)


#--------------------
# factor v

# filter for factor v data from liquid plasma
df_fv <- df |>
  filter(product == "liquid plasma",
         factor == "factor v") |>
  select(sample_id, day, factor, value) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         fv = `factor v`)

# day 15 and day 26 data
df_fv_15 <- df_fv |> filter(time == 15)
df_fv_26 <- df_fv |> filter(time == 26)

# fit t distribution to day 15 data with fixed nu = 2
fit_fv_15 <- ulam(
  alist(
    fv ~ dstudent(2, mu, sigma),
    mu ~ dnorm(85, 15),
    sigma ~ dexp(1)
  ),
  data = df_fv_15, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# fit t distribution to day 15 data and fitting for nu
fit_fv_15_nu <- ulam(
  alist(
    fv ~ dstudent(nu, mu, sigma),
    mu ~ dnorm(85, 15),
    sigma ~ dexp(1),
    nu ~ dexp(0.5)
  ),
  data = df_fv_15, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# compare two models -- fixed nu has better predictive fit
compare(fit_fv_15, fit_fv_15_nu)


# fit t distribution to day 26 data with fixed nu = 2
fit_fv_26 <- ulam(
  alist(
    fv ~ dstudent(2, mu, sigma),
    mu ~ dnorm(85, 15),
    sigma ~ dexp(1)
  ),
  data = df_fv_26, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# fit t distribution to day 26 data and fitting for nu
fit_fv_26_nu <- ulam(
  alist(
    fv ~ dstudent(nu, mu, sigma),
    mu ~ dnorm(85, 15),
    sigma ~ dexp(1),
    nu ~ dexp(0.5)
  ),
  data = df_fv_26, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# compare two models -- fixed nu has better predictive fit
compare(fit_fv_26, fit_fv_26_nu)


#--------------------
# factor vii

# filter for factor vii data from liquid plasma
df_fvii <- df |>
  filter(product == "liquid plasma",
         factor == "factor vii") |>
  select(sample_id, day, factor, value) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         fvii = `factor vii`)

# day 15 and day 26 data
df_fvii_15 <- df_fvii |> filter(time == 15)
df_fvii_26 <- df_fvii |> filter(time == 26)

# fit t distribution to day 15 data with fixed nu = 2
fit_fvii_15 <- ulam(
  alist(
    fvii ~ dstudent(2, mu, sigma),
    mu ~ dnorm(100, 25),
    sigma ~ dexp(1)
  ),
  data = df_fvii_15, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# fit t distribution to day 15 data and fitting for nu
fit_fvii_15_nu <- ulam(
  alist(
    fvii ~ dstudent(nu, mu, sigma),
    mu ~ dnorm(100, 25),
    sigma ~ dexp(1),
    nu ~ dexp(0.5)
  ),
  data = df_fvii_15, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# compare two models -- fixed nu has better predictive fit
compare(fit_fvii_15, fit_fvii_15_nu)


# fit t distribution to day 26 data with fixed nu = 2
fit_fvii_26 <- ulam(
  alist(
    fvii ~ dstudent(2, mu, sigma),
    mu ~ dnorm(100, 25),
    sigma ~ dexp(1)
  ),
  data = df_fvii_26, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# fit t distribution to day 26 data and fitting for nu
fit_fvii_26_nu <- ulam(
  alist(
    fvii ~ dstudent(nu, mu, sigma),
    mu ~ dnorm(100, 25),
    sigma ~ dexp(1),
    nu ~ dexp(0.5)
  ),
  data = df_fvii_26, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# compare two models -- fixed nu has better predictive fit
compare(fit_fvii_26, fit_fvii_26_nu)




#--------------------
# factor viii

# filter for factor viii data from liquid plasma
df_fviii <- df |>
  filter(product == "liquid plasma",
         factor == "factor viii") |>
  select(sample_id, day, factor, value) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         fviii = `factor viii`)

# day 15 and day 26 data
df_fviii_15 <- df_fviii |> filter(time == 15)
df_fviii_26 <- df_fviii |> filter(time == 26)

# fit t distribution to day 15 data with fixed nu = 2
fit_fviii_15 <- ulam(
  alist(
    fviii ~ dstudent(2, mu, sigma),
    mu ~ dnorm(80, 20),
    sigma ~ dexp(1)
  ),
  data = df_fviii_15, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# fit t distribution to day 15 data and fitting for nu
fit_fviii_15_nu <- ulam(
  alist(
    fviii ~ dstudent(nu, mu, sigma),
    mu ~ dnorm(80, 20),
    sigma ~ dexp(1),
    nu ~ dexp(0.5)
  ),
  data = df_fviii_15, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# compare two models -- fixed nu has better predictive fit
compare(fit_fviii_15, fit_fviii_15_nu)


# fit t distribution to day 26 data with fixed nu = 2
fit_fviii_26 <- ulam(
  alist(
    fviii ~ dstudent(2, mu, sigma),
    mu ~ dnorm(80, 20),
    sigma ~ dexp(1)
  ),
  data = df_fviii_26, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# fit t distribution to day 26 data and fitting for nu
fit_fviii_26_nu <- ulam(
  alist(
    fviii ~ dstudent(nu, mu, sigma),
    mu ~ dnorm(80, 20),
    sigma ~ dexp(1),
    nu ~ dexp(0.5)
  ),
  data = df_fviii_26, chains = 2, log_lik = TRUE, iter = 10000, cores = 2
)

# compare two models -- fixed nu has better predictive fit
compare(fit_fviii_26, fit_fviii_26_nu)
