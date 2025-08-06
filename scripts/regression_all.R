library(tidyverse)
library(rethinking)
library(cmdstanr)
library(readxl)

df <- read_xlsx("data/liquid_plasma_data.xlsx", na = "NA")


df_results <- tibble(factor = rep(c('fibrinogen',
                                    'protein c',
                                    'protein s',
                                    'factor v',
                                    'factor vii',
                                    'factor viii'), each = 2),
                     product = "liquid plasma",
                     parameter = rep(c("beta1", "beta2"), 6),
                     mean = 0,
                     ci_lo = 0, ci_hi = 0)

# ------------------
#      fibrinogen
#

# liquid plasma data, convert days to start from day 0
df_lp <- df |>
  filter(product == "liquid plasma",
         factor == "fibrinogen") |>
  select(sample_id, day, factor, value) |>
  mutate(day = day - 15) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day)

data_list <- list(
  fibrinogen = df_lp$fibrinogen,
  time = df_lp$time,
  sample_id = as.factor(df_lp$sample_id)
)

# piecewise regression with breakpoint at day of expiration
mod_fib <- cmdstan_model("stan/fibrinogen_piecewise_regression.stan")

fit_fibrinogen_pw <- mod_fib$sample(
  data = data_list,
  # seed = 1024,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.95
)

df_post <- precis(fit_fibrinogen_pw, prob = 0.95) |> as.data.frame()
df_results[1, "mean"]  <- df_post[rownames(df_post) == "beta1_bar", "mean"]
df_results[1, "ci_lo"] <- df_post[rownames(df_post) == "beta1_bar", "2.5%"]
df_results[1, "ci_hi"] <- df_post[rownames(df_post) == "beta1_bar", "97.5%"]

df_results[2, "mean"]  <- df_post[rownames(df_post) == "beta2_bar", "mean"]
df_results[2, "ci_lo"] <- df_post[rownames(df_post) == "beta2_bar", "2.5%"]
df_results[2, "ci_hi"] <- df_post[rownames(df_post) == "beta2_bar", "97.5%"]

df_frequentist <- lm(fibrinogen ~ time, df_lp) |>
  broom::tidy(conf.int = TRUE) |>
  (\(x) `[`(x, 2, ))()

# ------------------
#      protein c
#

# liquid plasma data, convert days to start from day 0
df_lp <- df |>
  filter(product == "liquid plasma",
         factor == "protein c") |>
  select(sample_id, day, factor, value) |>
  mutate(day = day - 15) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         protein_c = `protein c`)

data_list <- list(
  protein_c = df_lp$protein_c,
  time = df_lp$time,
  sample_id = as.factor(df_lp$sample_id)
)


# piecewise regression with breakpoint at day of expiration
mod_prc <- cmdstan_model("stan/proteinc_piecewise_regression.stan")

fit_prc_pw <- mod_prc$sample(
  data = data_list,
  # seed = 1024,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.95
)

df_post_prc <- precis(fit_prc_pw, prob = 0.95) |> as.data.frame()
df_results[3, "mean"]  <- df_post_prc[rownames(df_post) == "beta1_bar", "mean"]
df_results[3, "ci_lo"] <- df_post_prc[rownames(df_post) == "beta1_bar", "2.5%"]
df_results[3, "ci_hi"] <- df_post_prc[rownames(df_post) == "beta1_bar", "97.5%"]

df_results[4, "mean"]  <- df_post_prc[rownames(df_post) == "beta2_bar", "mean"]
df_results[4, "ci_lo"] <- df_post_prc[rownames(df_post) == "beta2_bar", "2.5%"]
df_results[4, "ci_hi"] <- df_post_prc[rownames(df_post) == "beta2_bar", "97.5%"]

df_frequentist <- bind_rows(df_frequentist,
                            lm(protein_c ~ time, df_lp) |>
                              broom::tidy(conf.int = TRUE) |>
                              (\(x) `[`(x, 2, ))()
)                            
# ------------------
#      protein s
#

# liquid plasma data, convert days to start from day 0
df_lp <- df |>
  filter(product == "liquid plasma",
         factor == "protein s") |>
  select(sample_id, day, factor, value) |>
  mutate(day = day - 15) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         protein_s = `protein s`)

data_list <- list(
  protein_s = df_lp$protein_s,
  time = df_lp$time,
  sample_id = as.factor(df_lp$sample_id)
)


# piecewise regression with breakpoint at day of expiration
mod_prs <- cmdstan_model("stan/proteins_piecewise_regression.stan")

fit_prs_pw <- mod_prs$sample(
  data = data_list,
  # seed = 1024,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.95
)

df_post_prc <- precis(fit_prs_pw, prob = 0.95) |> as.data.frame()
df_results[5, "mean"]  <- df_post_prc[rownames(df_post_prc) == "beta1_bar", "mean"]
df_results[5, "ci_lo"] <- df_post_prc[rownames(df_post_prc) == "beta1_bar", "2.5%"]
df_results[5, "ci_hi"] <- df_post_prc[rownames(df_post_prc) == "beta1_bar", "97.5%"]

df_results[6, "mean"]  <- df_post_prc[rownames(df_post_prc) == "beta2_bar", "mean"]
df_results[6, "ci_lo"] <- df_post_prc[rownames(df_post_prc) == "beta2_bar", "2.5%"]
df_results[6, "ci_hi"] <- df_post_prc[rownames(df_post_prc) == "beta2_bar", "97.5%"]


df_frequentist <- bind_rows(df_frequentist,
                            lm(protein_s ~ time, df_lp) |>
                              broom::tidy(conf.int = TRUE) |>
                              (\(x) `[`(x, 2, ))()
)   

# ------------------
#      factor v
#

# liquid plasma data, convert days to start from day 0
df_lp <- df |>
  filter(product == "liquid plasma",
         factor == "factor v") |>
  select(sample_id, day, factor, value) |>
  mutate(day = day - 15) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         factor_v = `factor v`)

data_list <- list(
  factor_v = df_lp$factor_v,
  time = df_lp$time,
  sample_id = as.factor(df_lp$sample_id)
)


# piecewise regression with breakpoint at day of expiration
mod_fv <- cmdstan_model("stan/factorv_piecewise_regression.stan")

fit_fv_pw <- mod_fv$sample(
  data = data_list,
  # seed = 1024,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.95
)

df_post_fv <- precis(fit_fv_pw, prob = 0.95) |> as.data.frame()
df_results[7, "mean"]  <- df_post_fv[rownames(df_post_fv) == "beta1_bar", "mean"]
df_results[7, "ci_lo"] <- df_post_fv[rownames(df_post_fv) == "beta1_bar", "2.5%"]
df_results[7, "ci_hi"] <- df_post_fv[rownames(df_post_fv) == "beta1_bar", "97.5%"]

df_results[8, "mean"]  <- df_post_fv[rownames(df_post_fv) == "beta2_bar", "mean"]
df_results[8, "ci_lo"] <- df_post_fv[rownames(df_post_fv) == "beta2_bar", "2.5%"]
df_results[8, "ci_hi"] <- df_post_fv[rownames(df_post_fv) == "beta2_bar", "97.5%"]


df_frequentist <- bind_rows(df_frequentist,
                            lm(factor_v ~ time, df_lp) |>
                              broom::tidy(conf.int = TRUE) |>
                              (\(x) `[`(x, 2, ))()
)   

# ------------------
#      factor vii
#

# liquid plasma data, convert days to start from day 0
df_lp <- df |>
  filter(product == "liquid plasma",
         factor == "factor vii") |>
  select(sample_id, day, factor, value) |>
  mutate(day = day - 15) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         factor_vii = `factor vii`)

data_list <- list(
  factor_vii = df_lp$factor_vii,
  time = df_lp$time,
  sample_id = as.factor(df_lp$sample_id)
)


# piecewise regression with breakpoint at day of expiration
mod_fvii <- cmdstan_model("stan/factorvii_piecewise_regression.stan")

fit_fvii_pw <- mod_fvii$sample(
  data = data_list,
  # seed = 1024,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.95
)

df_post_fvii <- precis(fit_fvii_pw, prob = 0.95) |> as.data.frame()
df_results[9, "mean"]  <- df_post_fvii[rownames(df_post_fvii) == "beta1_bar", "mean"]
df_results[9, "ci_lo"] <- df_post_fvii[rownames(df_post_fvii) == "beta1_bar", "2.5%"]
df_results[9, "ci_hi"] <- df_post_fvii[rownames(df_post_fvii) == "beta1_bar", "97.5%"]

df_results[10, "mean"]  <- df_post_fvii[rownames(df_post_fvii) == "beta2_bar", "mean"]
df_results[10, "ci_lo"] <- df_post_fvii[rownames(df_post_fvii) == "beta2_bar", "2.5%"]
df_results[10, "ci_hi"] <- df_post_fvii[rownames(df_post_fvii) == "beta2_bar", "97.5%"]

df_frequentist <- bind_rows(df_frequentist,
                            lm(factor_vii ~ time, df_lp) |>
                              broom::tidy(conf.int = TRUE) |>
                              (\(x) `[`(x, 2, ))()
)   


# ------------------
#      factor viii
#

# liquid plasma data, convert days to start from day 0
df_lp <- df |>
  filter(product == "liquid plasma",
         factor == "factor viii") |>
  select(sample_id, day, factor, value) |>
  mutate(day = day - 15) |>
  pivot_wider(names_from = factor, values_from = value) |>
  rename(time = day,
         factor_viii = `factor viii`)

data_list <- list(
  factor_viii = df_lp$factor_viii,
  time = df_lp$time,
  sample_id = as.factor(df_lp$sample_id)
)


# piecewise regression with breakpoint at day of expiration
mod_fviii <- cmdstan_model("stan/factorviii_piecewise_regression.stan")

fit_fviii_pw <- mod_fviii$sample(
  data = data_list,
  # seed = 1024,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  parallel_chains = 4,
  refresh = 500,
  adapt_delta = 0.95
)

df_post_fviii <- precis(fit_fviii_pw, prob = 0.95) |> as.data.frame()
df_results[11, "mean"]  <- df_post_fviii[rownames(df_post_fvii) == "beta1_bar", "mean"]
df_results[11, "ci_lo"] <- df_post_fviii[rownames(df_post_fvii) == "beta1_bar", "2.5%"]
df_results[11, "ci_hi"] <- df_post_fviii[rownames(df_post_fvii) == "beta1_bar", "97.5%"]

df_results[12, "mean"]  <- df_post_fviii[rownames(df_post_fvii) == "beta2_bar", "mean"]
df_results[12, "ci_lo"] <- df_post_fviii[rownames(df_post_fvii) == "beta2_bar", "2.5%"]
df_results[12, "ci_hi"] <- df_post_fviii[rownames(df_post_fvii) == "beta2_bar", "97.5%"]


df_frequentist <- bind_rows(df_frequentist,
                            lm(factor_viii ~ time, df_lp) |>
                              broom::tidy(conf.int = TRUE) |>
                              (\(x) `[`(x, 2, ))()
)   

# frequentist method
df_frequentist$factor <- c("fibrinogen", "protein c", "protein s", "factor v", "factor vii", "factor viii")


# --------------------------
# posterior predictive plots
df <- read_xlsx("data/liquid_plasma_data.xlsx", na = "NA")

draws_array <- fit_fibrinogen_pw$draws(variables = c("beta1_bar", "beta2_bar", "sigma_b1", "sigma_b2", "alpha", "sigma"))
time <- seq(0, 16.5, 0.5)
alpha_level <- 0.6
size_level <- 1.7
width_level <- .3

posterior_means <- function(draws_array, time) {
  df_factor <- tibble(time = time)
  beta1_bar <- draws_array[, , "beta1_bar"] |> as.vector()
  beta2_bar <- draws_array[, , "beta2_bar"] |> as.vector()
  sigma_b1 <- draws_array[, , "sigma_b1"] |> as.vector()
  sigma_b2 <- draws_array[, , "sigma_b2"] |> as.vector()
  alpha <- draws_array[, , "alpha"] |> as.vector()
  sigma <- draws_array[, , "sigma"] |> as.vector()
  
  df_factor <- df_factor |> mutate(
    factor_avg = purrr::map(time, function(x) if(x <= 11) {alpha + beta1_bar * (x - 11)} else {alpha + beta2_bar * (x - 11)}),
    mean_factor = map_dbl(factor_avg, function(x) mean(x)),
    ci_lo = map_dbl(factor_avg, function(x) quantile(x, probs = c(0.025))),
    ci_hi = map_dbl(factor_avg, function(x) quantile(x, probs = c(0.975)))
  )
  
  df_factor |> select(-factor_avg) |> mutate(time = time + 15)
}

# fibrinogen
draws_array_fibrinogen <- fit_fibrinogen_pw$draws(variables = c("beta1_bar", "beta2_bar", "sigma_b1", "sigma_b2", "alpha", "sigma"))

df_factor_fibrinogen  <- posterior_means(draws_array_fibrinogen, time) 
p_fibrinogen <- df |>
  filter(factor == "fibrinogen", product == "liquid plasma") |>
  ggplot() +
  geom_ribbon(aes(time, ymin = ci_lo, ymax = ci_hi), data = df_factor_fibrinogen, alpha = 0.3) +
  geom_jitter(aes(day, value), width = width_level, height = 0, alpha = alpha_level, stroke = NA, size = size_level) +
  geom_line(aes(time, mean_factor), data = df_factor_fibrinogen) +
  theme_light() +
  xlab("Day of Storage") +
  ylab("Activity (mg/dL)") +
  coord_cartesian(xlim = c(15, 28), ylim = c(0, NA))

# protein c
draws_array_prc <- fit_prc_pw$draws(variables = c("beta1_bar", "beta2_bar", "sigma_b1", "sigma_b2", "alpha", "sigma"))

df_factor_prc <- posterior_means(draws_array_prc, time) 
p_prc <- df |>
  filter(factor == "protein c", product == "liquid plasma") |>
  ggplot() +
  geom_ribbon(aes(time, ymin = ci_lo, ymax = ci_hi), data = df_factor_prc, alpha = 0.3) +
  geom_jitter(aes(day, value), width = width_level, height = 0, alpha = alpha_level, stroke = NA, size = size_level) +
  geom_line(aes(time, mean_factor), data = df_factor_prc) +
  theme_light() +
  xlab("Day of Storage") +
  ylab("Activity (%)") +
  coord_cartesian(xlim = c(15, 28), ylim = c(0, NA))
p_prc


# protein s
draws_array_prs <- fit_prs_pw$draws(variables = c("beta1_bar", "beta2_bar", "sigma_b1", "sigma_b2", "alpha", "sigma"))

df_factor_prs <- posterior_means(draws_array_prs, time) 
p_prs <- df |>
  filter(factor == "protein s", product == "liquid plasma") |>
  ggplot() +
  geom_ribbon(aes(time, ymin = ci_lo, ymax = ci_hi), data = df_factor_prs, alpha = 0.3) +
  geom_jitter(aes(day, value), width = width_level, height = 0, alpha = alpha_level, stroke = NA, size = size_level) +
  geom_line(aes(time, mean_factor), data = df_factor_prs) +
  theme_light() +
  xlab("Day of Storage") +
  ylab("Activity (%)") +
  coord_cartesian(xlim = c(15, 28), ylim = c(0, NA))
p_prs



# factor v
draws_array_fv <- fit_fv_pw$draws(variables = c("beta1_bar", "beta2_bar", "sigma_b1", "sigma_b2", "alpha", "sigma"))

df_factor_fv <- posterior_means(draws_array_fv, time) 
p_fv <- df |>
  filter(factor == "factor v", product == "liquid plasma") |>
  ggplot() +
  geom_ribbon(aes(time, ymin = ci_lo, ymax = ci_hi), data = df_factor_fv, alpha = 0.3) +
  geom_jitter(aes(day, value), width = width_level, height = 0, alpha = alpha_level, stroke = NA, size = size_level) +
  geom_line(aes(time, mean_factor), data = df_factor_fv) +
  theme_light() +
  xlab("Day of Storage") +
  ylab("Activity (%)") +
  coord_cartesian(xlim = c(15, 28), ylim = c(0, NA))
p_fv



# factor vii
draws_array_fvii <- fit_fvii_pw$draws(variables = c("beta1_bar", "beta2_bar", "sigma_b1", "sigma_b2", "alpha", "sigma"))

df_factor_fvii <- posterior_means(draws_array_fvii, time) 
p_fvii <- df |>
  filter(factor == "factor vii", product == "liquid plasma") |>
  ggplot() +
  geom_ribbon(aes(time, ymin = ci_lo, ymax = ci_hi), data = df_factor_fvii, alpha = 0.3) +
  geom_jitter(aes(day, value), width = width_level, height = 0, alpha = alpha_level, stroke = NA, size = size_level) +
  geom_line(aes(time, mean_factor), data = df_factor_fvii) +
  theme_light() +
  xlab("Day of Storage") +
  ylab("Activity (%)") +
  coord_cartesian(xlim = c(15, 28), ylim = c(0, NA))
p_fvii



# factor viii
draws_array_fviii <- fit_fviii_pw$draws(variables = c("beta1_bar", "beta2_bar", "sigma_b1", "sigma_b2", "alpha", "sigma"))

df_factor_fviii <- posterior_means(draws_array_fviii, time) 
p_fviii <- df |>
  filter(factor == "factor viii", product == "liquid plasma") |>
  ggplot() +
  geom_ribbon(aes(time, ymin = ci_lo, ymax = ci_hi), data = df_factor_fviii, alpha = 0.3) +
  geom_jitter(aes(day, value), width = width_level, height = 0, alpha = alpha_level, stroke = NA, size = size_level) +
  geom_line(aes(time, mean_factor), data = df_factor_fviii) +
  theme_light() +
  xlab("Day of Storage") +
  ylab("Activity (%)") +
  coord_cartesian(xlim = c(15, 28), ylim = c(0, NA))
p_fviii

p_fibrinogen <- p_fibrinogen +
  facet_wrap(~factor) + xlab(NULL) + theme(axis.text = element_text(color="black")) +
  scale_x_continuous(
    breaks = seq(0, 40, by = 4),           # Major ticks every 4 units
    minor_breaks = seq(0, 40, by = 2)      # Minor ticks every 2 unit
  )
p_prc <- p_prc + facet_wrap(~factor) + xlab(NULL) + ylab(NULL) + theme(axis.text = element_text(color="black")) +
  scale_x_continuous(
    breaks = seq(0, 40, by = 4),           # Major ticks every 4 units
    minor_breaks = seq(0, 40, by = 2)      # Minor ticks every 2 unit
  )
p_prs <- p_prs + facet_wrap(~factor) + xlab(NULL) + ylab(NULL) + theme(axis.text = element_text(color="black")) +
  scale_x_continuous(
    breaks = seq(0, 40, by = 4),           # Major ticks every 4 units
    minor_breaks = seq(0, 40, by = 2)      # Minor ticks every 2 unit
  )
p_fv <- p_fv + facet_wrap(~factor) + theme(axis.text = element_text(color="black")) +
  scale_x_continuous(
    breaks = seq(0, 40, by = 4),           # Major ticks every 4 units
    minor_breaks = seq(0, 40, by = 2)      # Minor ticks every 2 unit
  )
p_fvii <- p_fvii + facet_wrap(~factor) + ylab(NULL) + theme(axis.text = element_text(color="black")) +
  scale_x_continuous(
    breaks = seq(0, 40, by = 4),           # Major ticks every 4 units
    minor_breaks = seq(0, 40, by = 2)      # Minor ticks every 2 unit
  )
p_fviii <- p_fviii + facet_wrap(~factor) + ylab(NULL) + theme(axis.text = element_text(color="black")) +
  scale_x_continuous(
    breaks = seq(0, 40, by = 4),           # Major ticks every 4 units
    minor_breaks = seq(0, 40, by = 2)      # Minor ticks every 2 unit
  )

fig_b <- egg::ggarrange(p_fibrinogen, p_prc, p_prs, p_fv, p_fvii, p_fviii, nrow = 2)

ggsave(filename = "figures/fig1b1.pdf", plot = fig_b, device = "pdf", colormodel = "cmyk",
       width = 6.3, height = 3.8, units = "in")


# figure a
width_size <- 0.2

fig_a <- df |> 
  mutate(factor = factor(factor, levels = c("fibrinogen", "protein c", "protein s", "factor v", "factor vii", "factor viii"))) |>
  mutate(product_label = ifelse(product == "liquid plasma", "LQP", "TP"),
         product_label = paste0(product_label, day),
         product_label = as.factor(product_label)) |>
  ggplot() +
  geom_boxplot(aes(product_label, value), outlier.shape = NA, width = width_size + .25) +
  geom_jitter(aes(product_label, value), width = width_size, height = 0, stroke = NA, size = 1.4, alpha = 0.45) + 
  facet_wrap(~factor, scales = "free") +
  theme_light() +
  coord_cartesian(ylim = c(0, NA)) + theme(axis.text = element_text(color="black")) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

ggsave(filename = "figures/fig1a1.pdf", plot = fig_a, device = "pdf", colormodel = "cmyk",
       width = 6.3, height = 4, units = "in")
