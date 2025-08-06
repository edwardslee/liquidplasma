library(tidyverse)
library(rethinking)
library(cmdstanr)
library(readxl)

# import data
df <- read_xlsx("data/liquid_plasma_data.xlsx", na = "NA")

df_results <- tibble(factor = c("fibrinogen", "protein c", "protein s", "factor v", "factor vii", "factor viii"),
                     mean = NULL,
                     CI_lo = NULL,
                     CI_hi = NULL)

#--------------------
#    fibrinogen 

# filter for fibrinogen data
df_fibrinogen_compare <- df |>
  filter(factor == "fibrinogen") |>
  select(sample_id, product, day, factor, value) |>
  filter(day == 26 | day == 5) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day) |>
  mutate(product = ifelse(product == "liquid plasma", 1, 2))


mod_compare_fibrinogen <- cmdstan_model("stan/compare_fibrinogen.stan")

data_list_compare_fibrinogen <- list(
  fibrinogen = df_fibrinogen_compare$fibrinogen,
  product = df_fibrinogen_compare$product
)

fit_fibrinogen_compare <- mod_compare_fibrinogen$sample(
  data = data_list_compare_fibrinogen,
  # seed = 1024,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.95
)


df_post <- precis(fit_fibrinogen_compare, prob = 0.95) |> as.data.frame()

df_results[1, "mean"]  <- df_post[rownames(df_post) == "mu_diff", "mean"]
df_results[1, "CI_lo"] <- df_post[rownames(df_post) == "mu_diff", "2.5%"]
df_results[1, "CI_hi"] <- df_post[rownames(df_post) == "mu_diff", "97.5%"]


#--------------------
#    protein c

# filter for protein c data
df_prc_compare <- df |>
  filter(factor == "protein c") |>
  select(sample_id, product, day, factor, value) |>
  filter(day == 26 | day == 5) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         prc = `protein c`) |>
  mutate(product = ifelse(product == "liquid plasma", 1, 2))


mod_compare_prc <- cmdstan_model("stan/compare_prc.stan")

data_list_compare_prc <- list(
  prc = df_prc_compare$prc,
  product = df_prc_compare$product
)

fit_prc_compare <- mod_compare_prc$sample(
  data = data_list_compare_prc,
  # seed = 1024,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.95
)

df_post <- precis(fit_prc_compare, prob = 0.95) |> as.data.frame()

df_results[2, "mean"]  <- df_post[rownames(df_post) == "mu_diff", "mean"]
df_results[2, "CI_lo"] <- df_post[rownames(df_post) == "mu_diff", "2.5%"]
df_results[2, "CI_hi"] <- df_post[rownames(df_post) == "mu_diff", "97.5%"]



#--------------------
#    protein s

# filter for protein c data
df_prs_compare <- df |>
  filter(factor == "protein s") |>
  select(sample_id, product, day, factor, value) |>
  filter(day == 26 | day == 5) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         prs = `protein s`) |>
  mutate(product = ifelse(product == "liquid plasma", 1, 2))

mod_compare_prs <- cmdstan_model("stan/compare_prs.stan")

data_list_compare_prs <- list(
  prs = df_prs_compare$prs,
  product = df_prs_compare$product
)

fit_prs_compare <- mod_compare_prs$sample(
  data = data_list_compare_prs,
  # seed = 1024,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.95
)

df_post <- precis(fit_prs_compare, prob = 0.95) |> as.data.frame()

df_results[3, "mean"]  <- df_post[rownames(df_post) == "mu_diff", "mean"]
df_results[3, "CI_lo"] <- df_post[rownames(df_post) == "mu_diff", "2.5%"]
df_results[3, "CI_hi"] <- df_post[rownames(df_post) == "mu_diff", "97.5%"]



#--------------------
#    factor v

# filter for factor v data
df_fv_compare <- df |>
  filter(factor == "factor v") |>
  select(sample_id, product, day, factor, value) |>
  filter(day == 26 | day == 5) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         fv = `factor v`) |>
  mutate(product = ifelse(product == "liquid plasma", 1, 2))

mod_compare_fv <- cmdstan_model("stan/compare_fv.stan")

data_list_compare_fv <- list(
  fv = df_fv_compare$fv,
  product = df_fv_compare$product
)

fit_fv_compare <- mod_compare_fv$sample(
  data = data_list_compare_fv,
  # seed = 1024,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.95
)

df_post <- precis(fit_fv_compare, prob = 0.95) |> as.data.frame()

df_results[4, "mean"]  <- df_post[rownames(df_post) == "mu_diff", "mean"]
df_results[4, "CI_lo"] <- df_post[rownames(df_post) == "mu_diff", "2.5%"]
df_results[4, "CI_hi"] <- df_post[rownames(df_post) == "mu_diff", "97.5%"]



#--------------------
#    factor vii

# filter for factor v data
df_fvii_compare <- df |>
  filter(factor == "factor vii") |>
  select(sample_id, product, day, factor, value) |>
  filter(day == 26 | day == 5) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         fvii = `factor vii`) |>
  mutate(product = ifelse(product == "liquid plasma", 1, 2))

mod_compare_fvii <- cmdstan_model("stan/compare_fvii.stan")

data_list_compare_fvii <- list(
  fvii = df_fvii_compare$fvii,
  product = df_fvii_compare$product
)

fit_fvii_compare <- mod_compare_fvii$sample(
  data = data_list_compare_fvii,
  # seed = 1024,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.95
)

df_post <- precis(fit_fvii_compare, prob = 0.95) |> as.data.frame()

df_results[5, "mean"]  <- df_post[rownames(df_post) == "mu_diff", "mean"]
df_results[5, "CI_lo"] <- df_post[rownames(df_post) == "mu_diff", "2.5%"]
df_results[5, "CI_hi"] <- df_post[rownames(df_post) == "mu_diff", "97.5%"]

# compare to frequentist method
t.test(fvii ~ product, data = df_fvii_compare)





#--------------------
#    factor viii

# filter for factor viii data
df_fviii_compare <- df |>
  filter(factor == "factor viii") |>
  select(sample_id, product, day, factor, value) |>
  filter(day == 26 | day == 5) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         fviii = `factor viii`) |>
  mutate(product = ifelse(product == "liquid plasma", 1, 2))

mod_compare_fviii <- cmdstan_model("stan/compare_fviii.stan")

data_list_compare_fviii <- list(
  fviii = df_fviii_compare$fviii,
  product = df_fviii_compare$product
)

fit_fviii_compare <- mod_compare_fviii$sample(
  data = data_list_compare_fviii,
  # seed = 1024,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.95
)

df_post <- precis(fit_fviii_compare, prob = 0.95) |> as.data.frame()

df_results[6, "mean"]  <- df_post[rownames(df_post) == "mu_diff", "mean"]
df_results[6, "CI_lo"] <- df_post[rownames(df_post) == "mu_diff", "2.5%"]
df_results[6, "CI_hi"] <- df_post[rownames(df_post) == "mu_diff", "97.5%"]

# compare to frequentist method
df_frequentist <- bind_rows(
  t.test(fibrinogen ~ product, data = df_fibrinogen_compare) |> broom::tidy(),
  t.test(prc ~ product, data = df_prc_compare) |> broom::tidy(),
  t.test(prs ~ product, data = df_prs_compare) |> broom::tidy(),
  t.test(fv ~ product, data = df_fv_compare) |> broom::tidy(),
  t.test(fvii ~ product, data = df_fvii_compare) |> broom::tidy(),
  t.test(fviii ~ product, data = df_fviii_compare)  |> broom::tidy()
)


print(df_results)
print(df_frequentist)