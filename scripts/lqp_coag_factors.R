library(tidyverse)
library(rethinking)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(shinystan)
library(readxl)

# import data
df <- read_xlsx("data/liquid_plasma_data_updated.xlsx", na = "NA")

df_results <- tibble(factor = rep(c('fibrinogen',
                                    'protein c',
                                    'protein s',
                                    'factor v',
                                    'factor vii',
                                    'factor viii'), each = 3),
                     product = "liquid plasma",
                     time = rep(c(15, 26, 27), 6),
                     mean = 0,
                     ci_lo = 0, ci_hi = 0)

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
df_fibrinogen_27 <- df_fibrinogen |> filter(time == 27)

# fit t distribution to day 15 data with fixed nu = 2
fit_fibrinogen_15 <- ulam(
  alist(
    fibrinogen ~ dstudent(2, mu, sigma),
    mu ~ dnorm(280, 50),
    sigma ~ dexp(1)
  ),
  data = df_fibrinogen_15, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


# fit t distribution to day 26 data with fixed nu = 2
fit_fibrinogen_26 <- ulam(
  alist(
    fibrinogen ~ dstudent(2, mu, sigma),
    mu ~ dnorm(280, 50),
    sigma ~ dexp(1)
  ),
  data = df_fibrinogen_26, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


# fit t distribution to day 27 data with fixed nu = 2
fit_fibrinogen_27 <- ulam(
  alist(
    fibrinogen ~ dstudent(2, mu, sigma),
    mu ~ dnorm(280, 50),
    sigma ~ dexp(1)
  ),
  data = df_fibrinogen_27, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)

df_post <- precis(fit_fibrinogen_15, prob = 0.95) |> as.data.frame()
df_results[1, "mean"] <- df_post[rownames(df_post) == "mu", "mean"]
df_results[1, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[1, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

df_post <- precis(fit_fibrinogen_26, prob = 0.95) |> as.data.frame()
df_results[2, "mean"] <- df_post[rownames(df_post) == "mu", "mean"]
df_results[2, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[2, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

df_post <- precis(fit_fibrinogen_27, prob = 0.95) |> as.data.frame()
df_results[3, "mean"] <- df_post[rownames(df_post) == "mu", "mean"]
df_results[3, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[3, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]


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
df_prc_27 <- df_prc |> filter(time == 27)

# fit t distribution to day 15 data with fixed nu = 2
fit_prc_15 <- ulam(
  alist(
    prc ~ dstudent(2, mu, sigma),
    mu ~ dnorm(90, 20),
    sigma ~ dexp(1)
  ),
  data = df_prc_15, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


# fit t distribution to day 26 data with fixed nu = 2
fit_prc_26 <- ulam(
  alist(
    prc ~ dstudent(2, mu, sigma),
    mu ~ dnorm(90, 20),
    sigma ~ dexp(1)
  ),
  data = df_prc_26, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


# fit t distribution to day 26 data with fixed nu = 2
fit_prc_27 <- ulam(
  alist(
    prc ~ dstudent(2, mu, sigma),
    mu ~ dnorm(90, 20),
    sigma ~ dexp(1)
  ),
  data = df_prc_27, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


df_post <- precis(fit_prc_15, prob = 0.95) |> as.data.frame()
df_results[4, "mean"] <- df_post[rownames(df_post) == "mu", "mean"]
df_results[4, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[4, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

df_post <- precis(fit_prc_26, prob = 0.95) |> as.data.frame()
df_results[5, "mean"] <- df_post[rownames(df_post) == "mu", "mean"]
df_results[5, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[5, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

df_post <- precis(fit_prc_27, prob = 0.95) |> as.data.frame()
df_results[6, "mean"] <- df_post[rownames(df_post) == "mu", "mean"]
df_results[6, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[6, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

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
df_prs_27 <- df_prs |> filter(time == 27)

# fit t distribution to day 15 data with fixed nu = 2
fit_prs_15 <- ulam(
  alist(
    prs ~ dstudent(2, mu, sigma),
    mu ~ dnorm(80, 20),
    sigma ~ dexp(1)
  ),
  data = df_prs_15, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


# fit t distribution to day 26 data with fixed nu = 2
fit_prs_26 <- ulam(
  alist(
    prs ~ dstudent(2, mu, sigma),
    mu ~ dnorm(80, 20),
    sigma ~ dexp(1)
  ),
  data = df_prs_26, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


# fit t distribution to day 26 data with fixed nu = 2
fit_prs_27 <- ulam(
  alist(
    prs ~ dstudent(2, mu, sigma),
    mu ~ dnorm(80, 20),
    sigma ~ dexp(1)
  ),
  data = df_prs_27, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


df_post <- precis(fit_prs_15, prob = 0.95) |> as.data.frame()
df_results[7, "mean"] <- df_post[rownames(df_post) == "mu", "mean"]
df_results[7, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[7, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

df_post <- precis(fit_prs_26, prob = 0.95) |> as.data.frame()
df_results[8, "mean"] <- df_post[rownames(df_post) == "mu", "mean"]
df_results[8, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[8, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

df_post <- precis(fit_prs_27, prob = 0.95) |> as.data.frame()
df_results[9, "mean"] <- df_post[rownames(df_post) == "mu", "mean"]
df_results[9, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[9, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

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
df_fv_27 <- df_fv |> filter(time == 27)

# fit t distribution to day 15 data with fixed nu = 2
fit_fv_15 <- ulam(
  alist(
    fv ~ dstudent(2, mu, sigma),
    mu ~ dnorm(60, 20),
    sigma ~ dexp(1)
  ),
  data = df_fv_15, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


# fit t distribution to day 26 data with fixed nu = 2
fit_fv_26 <- ulam(
  alist(
    fv ~ dstudent(2, mu, sigma),
    mu ~ dnorm(60, 20),
    sigma ~ dexp(1)
  ),
  data = df_fv_26, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


# fit t distribution to day 26 data with fixed nu = 2
fit_fv_27 <- ulam(
  alist(
    fv ~ dstudent(2, mu, sigma),
    mu ~ dnorm(60, 20),
    sigma ~ dexp(1)
  ),
  data = df_fv_27, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


df_post <- precis(fit_fv_15, prob = 0.95) |> as.data.frame()
df_results[10, "mean"] <- df_post[rownames(df_post) == "mu", "mean"]
df_results[10, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[10, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

df_post <- precis(fit_fv_26, prob = 0.95) |> as.data.frame()
df_results[11, "mean"] <- df_post[rownames(df_post) == "mu", "mean"]
df_results[11, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[11, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

df_post <- precis(fit_fv_27, prob = 0.95) |> as.data.frame()
df_results[12, "mean"] <- df_post[rownames(df_post) == "mu", "mean"]
df_results[12, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[12, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

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
df_fvii_27 <- df_fvii |> filter(time == 27)

# fit t distribution to day 15 data with fixed nu = 2
fit_fvii_15 <- ulam(
  alist(
    fvii ~ dstudent(2, mu, sigma),
    mu ~ dnorm(70, 20),
    sigma ~ dexp(1)
  ),
  data = df_fvii_15, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


# fit t distribution to day 26 data with fixed nu = 2
fit_fvii_26 <- ulam(
  alist(
    fvii ~ dstudent(2, mu, sigma),
    mu ~ dnorm(70, 20),
    sigma ~ dexp(1)
  ),
  data = df_fvii_26, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


# fit t distribution to day 26 data with fixed nu = 2
fit_fvii_27 <- ulam(
  alist(
    fvii ~ dstudent(2, mu, sigma),
    mu ~ dnorm(70, 20),
    sigma ~ dexp(1)
  ),
  data = df_fvii_27, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)



df_post <- precis(fit_fvii_15, prob = 0.95) |> as.data.frame()
df_results[13, "mean"]  <- df_post[rownames(df_post) == "mu", "mean"]
df_results[13, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[13, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

df_post <- precis(fit_fvii_26, prob = 0.95) |> as.data.frame()
df_results[14, "mean"]  <- df_post[rownames(df_post) == "mu", "mean"]
df_results[14, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[14, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

df_post <- precis(fit_fvii_27, prob = 0.95) |> as.data.frame()
df_results[15, "mean"]  <- df_post[rownames(df_post) == "mu", "mean"]
df_results[15, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[15, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

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
df_fviii_27 <- df_fviii |> filter(time == 27)

# fit t distribution to day 15 data with fixed nu = 2
fit_fviii_15 <- ulam(
  alist(
    fviii ~ dstudent(2, mu, sigma),
    mu ~ dnorm(50, 20),
    sigma ~ dexp(1)
  ),
  data = df_fviii_15, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


# fit t distribution to day 26 data with fixed nu = 2
fit_fviii_26 <- ulam(
  alist(
    fviii ~ dstudent(2, mu, sigma),
    mu ~ dnorm(50, 20),
    sigma ~ dexp(1)
  ),
  data = df_fviii_26, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)


# fit t distribution to day 26 data with fixed nu = 2
fit_fviii_27 <- ulam(
  alist(
    fviii ~ dstudent(2, mu, sigma),
    mu ~ dnorm(50, 20),
    sigma ~ dexp(1)
  ),
  data = df_fviii_27, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)

df_post <- precis(fit_fviii_15, prob = 0.95) |> as.data.frame()
df_results[16, "mean"]  <- df_post[rownames(df_post) == "mu", "mean"]
df_results[16, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[16, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

df_post <- precis(fit_fviii_26, prob = 0.95) |> as.data.frame()
df_results[17, "mean"]  <- df_post[rownames(df_post) == "mu", "mean"]
df_results[17, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[17, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]

df_post <- precis(fit_fviii_27, prob = 0.95) |> as.data.frame()
df_results[18, "mean"]  <- df_post[rownames(df_post) == "mu", "mean"]
df_results[18, "ci_lo"] <- df_post[rownames(df_post) == "mu", "2.5%"]
df_results[18, "ci_hi"] <- df_post[rownames(df_post) == "mu", "97.5%"]


#--------------------
#    thawed plasma
df_fibrinogen_tp <- df |>
  filter(product == "thawed plasma",
         factor == "fibrinogen") |>
  select(sample_id, day, factor, value) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day)

df_prc_tp <- df |>
  filter(product == "thawed plasma",
         factor == "protein c") |>
  select(sample_id, day, factor, value) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         prc = `protein c`)

df_prs_tp <- df |>
  filter(product == "thawed plasma",
         factor == "protein s") |>
  select(sample_id, day, factor, value) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         prs = `protein s`)

df_fv_tp <- df |>
  filter(product == "thawed plasma",
         factor == "factor v",
         sample_source == "bag") |>
  select(sample_id, day, factor, value) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         fv = `factor v`)

df_fvii_tp <- df |>
  filter(product == "thawed plasma",
         factor == "factor vii") |>
  select(sample_id, day, factor, value) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         fvii = `factor vii`)

df_fviii_tp <- df |>
  filter(product == "thawed plasma",
         factor == "factor viii") |>
  select(sample_id, day, factor, value) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         fviii = `factor viii`)

# fit t distribution to day 15 data with fixed nu = 2
fit_fibrinogen_tp <- ulam(
  alist(
    fibrinogen ~ dstudent(2, mu, sigma),
    mu ~ dnorm(280, 50),
    sigma ~ dexp(1)
  ),
  data = df_fibrinogen_tp, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)

fit_prc_tp <- ulam(
  alist(
    prc ~ dstudent(2, mu, sigma),
    mu ~ dnorm(90, 20),
    sigma ~ dexp(1)
  ),
  data = df_prc_tp, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)

fit_prs_tp <- ulam(
  alist(
    prs ~ dstudent(2, mu, sigma),
    mu ~ dnorm(80, 20),
    sigma ~ dexp(1)
  ),
  data = df_prs_tp, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)

fit_fv_tp <- ulam(
  alist(
    fv ~ dstudent(2, mu, sigma),
    mu ~ dnorm(60, 20),
    sigma ~ dexp(1)
  ),
  data = df_fv_tp, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)

fit_fvii_tp <- ulam(
  alist(
    fvii ~ dstudent(2, mu, sigma),
    mu ~ dnorm(70, 20),
    sigma ~ dexp(1)
  ),
  data = df_fvii_tp, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)

fit_fviii_tp <- ulam(
  alist(
    fviii ~ dstudent(2, mu, sigma),
    mu ~ dnorm(50, 20),
    sigma ~ dexp(1)
  ),
  data = df_fviii_tp, chains = 4, log_lik = TRUE, iter = 10000, cores = 2
)

df_post <- precis(fit_fibrinogen_tp, prob = 0.95) |> as.data.frame()
df_results <- bind_rows(df_results,
                        tibble(
                          factor = "fibrinogen",
                          product = "thawed_plasma",
                          time = 5,
                          mean = df_post[rownames(df_post) == "mu", "mean"],
                          ci_lo = df_post[rownames(df_post) == "mu", "2.5%"],
                          ci_hi = df_post[rownames(df_post) == "mu", "97.5%"]
                        )
)

df_post <- precis(fit_prc_tp, prob = 0.95) |> as.data.frame()
df_results <- bind_rows(df_results,
                        tibble(
                          factor = "protein c",
                          product = "thawed_plasma",
                          time = 5,
                          mean = df_post[rownames(df_post) == "mu", "mean"],
                          ci_lo = df_post[rownames(df_post) == "mu", "2.5%"],
                          ci_hi = df_post[rownames(df_post) == "mu", "97.5%"]
                        )
)

df_post <- precis(fit_prs_tp, prob = 0.95) |> as.data.frame()
df_results <- bind_rows(df_results,
                        tibble(
                          factor = "protein s",
                          product = "thawed_plasma",
                          time = 5,
                          mean = df_post[rownames(df_post) == "mu", "mean"],
                          ci_lo = df_post[rownames(df_post) == "mu", "2.5%"],
                          ci_hi = df_post[rownames(df_post) == "mu", "97.5%"]
                        )
)

df_post <- precis(fit_fv_tp, prob = 0.95) |> as.data.frame()
df_results <- bind_rows(df_results,
                        tibble(
                          factor = "factor v",
                          product = "thawed_plasma",
                          time = 5,
                          mean = df_post[rownames(df_post) == "mu", "mean"],
                          ci_lo = df_post[rownames(df_post) == "mu", "2.5%"],
                          ci_hi = df_post[rownames(df_post) == "mu", "97.5%"]
                        )
)

df_post <- precis(fit_fvii_tp, prob = 0.95) |> as.data.frame()
df_results <- bind_rows(df_results,
                        tibble(
                          factor = "factor vii",
                          product = "thawed_plasma",
                          time = 5,
                          mean = df_post[rownames(df_post) == "mu", "mean"],
                          ci_lo = df_post[rownames(df_post) == "mu", "2.5%"],
                          ci_hi = df_post[rownames(df_post) == "mu", "97.5%"]
                        )
)

df_post <- precis(fit_fviii_tp, prob = 0.95) |> as.data.frame()
df_results <- bind_rows(df_results,
                        tibble(
                          factor = "factor viii",
                          product = "thawed_plasma",
                          time = 5,
                          mean = df_post[rownames(df_post) == "mu", "mean"],
                          ci_lo = df_post[rownames(df_post) == "mu", "2.5%"],
                          ci_hi = df_post[rownames(df_post) == "mu", "97.5%"]
                        )
)




#------------------
#      plotting
df_lp <- df |>
  filter(product == "liquid plasma",
         day != 27) |>
  rename(time = day)

ggplot(df_lp) +
  geom_jitter(aes(as.factor(time), value), width = 0.2, height = 0) +
  geom_point(aes(as.factor(time), mean), data = df_results, color = "red") + 
  geom_linerange(aes(x = as.factor(time), ymin = ci_lo, ymax = ci_hi), data = df_results) +
  facet_wrap(~factor, scales = "free") +
  coord_cartesian(ylim = c(0, NA)) +
  theme_light() +
  xlab("Day of Storage") +
  ylab("Activity Level (%)")
