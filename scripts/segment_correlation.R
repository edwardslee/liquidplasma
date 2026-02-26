library(tidyverse)
library(ggpmisc)
library(rethinking)

# lm_formula <- y ~ x

df <- read_csv("data/lqp_correlation.csv")

df_seg <- df |> filter(source == "seg")
df_bag <- df |> filter(source == "bag")


df_join <- left_join(df_seg, df_bag,
                     by = c("factor", "bag"),
                     suffix = c(".seg", ".bag"))

# EDA plot
# df_join |>
#   filter(bag != "D", bag != "H", bag != "G", bag != "I") |>
# ggplot(aes(x = value.seg, y = value.bag)) +
#   # points
#   geom_point() +
#   # y = x diagonal
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
#   # per-facet linear regression line
#   geom_smooth(method = "lm", se = TRUE, color = "blue") +
#   # equation + R^2 per-facet; parse = TRUE makes it render as math
#   stat_poly_eq(
#     aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#     formula = lm_formula,
#     parse = TRUE,
#     size = 3,
#     label.x.npc = "left",  # position in npc units per facet
#     label.y.npc = 0.95
#   ) +
#   facet_wrap(~ factor, scales = "free") +
#   coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
#   labs(x = "value.bag", y = "value.seg", 
#        title = "value.seg vs value.bag: diagonal (y=x) + per-facet LM and equation")


# pearson correlation
# cor(df_bag |> filter(factor == "Factor V", bag != "D") |> select(value),
#     df_seg |> filter(factor == "Factor V", bag != "D") |> select(value))
# 
# cor(df_bag |> filter(factor == "Factor VII", bag != "D") |> select(value),
#     df_seg |> filter(factor == "Factor VII", bag != "D") |> select(value))
# 
# cor(df_bag |> filter(factor == "Factor VIII", bag != "D") |> select(value),
#     df_seg |> filter(factor == "Factor VIII", bag != "D") |> select(value))



df_join_rethink <- rename(df_join, value_bag = value.bag, value_seg = value.seg) |>
  filter(bag != "D", bag != "H", bag != "G", bag != "I") |>
  select(bag, factor, value_bag, value_seg)

# frequentist approach
# lm(value_bag ~ value_seg, df_join_rethink |> filter(factor == "Factor V")) |>
#   broom::tidy(conf.int = TRUE)
# lm(value_bag ~ value_seg, df_join_rethink |> filter(factor == "Factor VII")) |>
#   broom::tidy(conf.int = TRUE)
# lm(value_bag ~ value_seg, df_join_rethink |> filter(factor == "Factor VIII")) |>
#   broom::tidy(conf.int = TRUE)
# lm(value_bag ~ value_seg, df_join_rethink |> filter(factor == "Fibrinogen")) |>
#   broom::tidy(conf.int = TRUE)
# lm(value_bag ~ value_seg, df_join_rethink |> filter(factor == "Protein C")) |>
#   broom::tidy(conf.int = TRUE)
# lm(value_bag ~ value_seg, df_join_rethink |> filter(factor == "Protein S")) |>
#   broom::tidy(conf.int = TRUE)


df_vii <- df_join_rethink |> filter(bag != "I", factor == "Factor VII") |> select(value_bag, value_seg)

m_vii <- ulam(
  alist(
    value_bag ~ dnorm(mu, sd),
    mu <- beta * value_seg + alpha,
    beta ~ dnorm(1, 0.1),
    alpha ~ dnorm(0, 20),
    sd ~ dexp(1)
  ),
  data = df_vii, iter = 5000, chains = 4
)

df_viii <- df_join_rethink |> filter(bag != "D", factor == "Factor VIII") |> select(value_bag, value_seg)

m_viii <- ulam(
  alist(
    value_bag ~ dnorm(mu, sd),
    mu <- beta * value_seg + alpha,
    beta ~ dnorm(1, 0.1),
    alpha ~ dnorm(0, 20),
    sd ~ dexp(1)
  ),
  data = df_viii, iter = 5000, chains = 4
)

precis(m_viii)

df_ps <- df_join_rethink |> filter(bag != "D", factor == "Protein S") |> select(value_bag, value_seg)

m_ps <- ulam(
  alist(
    value_bag ~ dnorm(mu, sd),
    mu <- beta * value_seg + alpha,
    beta ~ dnorm(1, 0.08),
    alpha ~ dnorm(0, 3),
    sd ~ dexp(1)
  ),
  data = df_ps, iter = 5000, chains = 4
)


df_pc <- df_join_rethink |> filter(bag != "D", factor == "Protein C") |> select(value_bag, value_seg)

m_pc <- ulam(
  alist(
    value_bag ~ dnorm(mu, sd),
    mu <- beta * value_seg + alpha,
    beta ~ dnorm(1, 0.1),
    alpha ~ dnorm(0, 5),
    sd ~ dexp(1)
  ),
  data = df_pc, iter = 5000, chains = 4
)

precis(m_pc)



df_fib <- df_join_rethink |> filter(bag != "D", factor == "Fibrinogen") |> select(value_bag, value_seg)

m_fib <- ulam(
  alist(
    value_bag ~ dnorm(mu, sd),
    mu <- beta * value_seg + alpha,
    beta ~ dnorm(1, 0.2),
    alpha ~ dnorm(0, 20),
    sd ~ dexp(1)
  ),
  data = df_fib, iter = 5000, chains = 4
)




df_v <- df_join_rethink |> filter(bag != "D", bag != "I", factor == "Factor V") |> select(value_bag, value_seg)

m_v <- ulam(
  alist(
    value_bag ~ dnorm(mu, sd),
    mu <- beta * value_seg + alpha,
    beta ~ dnorm(1, 0.25),
    alpha ~ dnorm(0, 20),
    sd ~ dexp(1)
  ),
  data = df_v, iter = 5000, chains = 4
)

precis(m_v, prob = 0.95)
precis(m_vii, prob = 0.95)
precis(m_viii, prob = 0.95)
precis(m_fib, prob = 0.95)
precis(m_pc, prob = 0.95)
precis(m_ps, prob = 0.95)

post_v <- extract.samples(m_v, n = 5000)
post_fv <- tibble(fv = seq(min(df_v$value_seg)-5, max(df_v$value_seg)+5, by = 2))
post_fv <- post_fv |>
  mutate(post = purrr::map(fv, function(x) post_v$beta * x + post_v$alpha)) |>
  mutate(mean = map_dbl(post, function(x) x |> mean()),
         ci_lo = map_dbl(post, function(x) quantile(x, probs = 0.025)),
         ci_hi = map_dbl(post, function(x) quantile(x, probs = 0.975)),
         factor = "Factor V") 

pv <- ggplot(df_v) +
  geom_point(aes(value_seg, value_bag)) +
  geom_line(aes(fv, mean), data = post_fv) +
  geom_ribbon(aes(x = fv, ymin = ci_lo, ymax = ci_hi), data = post_fv, alpha = 0.3) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  facet_wrap(~factor) +
  # y = x diagonal
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  theme_light() +
  xlab("Segment Activity (%)") + ylab("Bag Activity (%)") +
  theme(strip.text = element_text(colour = 'black'))+ theme(axis.text = element_text(size = 7),
                                                            axis.title = element_text(size = 9))


post_vii <- extract.samples(m_vii, n = 5000)
post_fvii <- tibble(fvii = seq(min(df_vii$value_seg)-5, max(df_vii$value_seg)+5, by = 2))
post_fvii <- post_fvii |>
  mutate(post = purrr::map(fvii, function(x) post_vii$beta * x + post_vii$alpha)) |>
  mutate(mean = map_dbl(post, function(x) x |> mean()),
         ci_lo = map_dbl(post, function(x) quantile(x, probs = 0.025)),
         ci_hi = map_dbl(post, function(x) quantile(x, probs = 0.975)),
         factor = "Factor VII")

pvii <- ggplot(df_vii) +
  geom_point(aes(value_seg, value_bag)) +
  geom_line(aes(fvii, mean), data = post_fvii) +
  geom_ribbon(aes(x = fvii, ymin = ci_lo, ymax = ci_hi), data = post_fvii, alpha = 0.3) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  facet_wrap(~factor) +
  # y = x diagonal
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  theme_light() +
  xlab("Segment Activity (%)") + ylab("Bag Activity (%)") +
  theme(strip.text = element_text(colour = 'black'))+ theme(axis.text = element_text(size = 7),
                                                            axis.title = element_text(size = 9))



post_viii <- extract.samples(m_viii, n = 5000)
post_fviii <- tibble(fviii = seq(min(df_viii$value_seg)-5, max(df_viii$value_seg)+5, by = 2))
post_fviii <- post_fviii |>
  mutate(post = purrr::map(fviii, function(x) post_viii$beta * x + post_viii$alpha)) |>
  mutate(mean = map_dbl(post, function(x) x |> mean()),
         ci_lo = map_dbl(post, function(x) quantile(x, probs = 0.025)),
         ci_hi = map_dbl(post, function(x) quantile(x, probs = 0.975)),
         factor = "Factor VIII")

pviii <- ggplot(df_viii) +
  geom_point(aes(value_seg, value_bag)) +
  geom_line(aes(fviii, mean), data = post_fviii) +
  geom_ribbon(aes(x = fviii, ymin = ci_lo, ymax = ci_hi), data = post_fviii, alpha = 0.3) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  facet_wrap(~factor) +
  # y = x diagonal
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  theme_light() +
  xlab("Segment Activity (%)") + ylab("Bag Activity (%)") +
  theme(strip.text = element_text(colour = 'black'))+ theme(axis.text = element_text(size = 7),
                                                            axis.title = element_text(size = 9))




post_s <- extract.samples(m_ps, n = 5000)
post_ps <- tibble(ps = seq(min(df_ps$value_seg), max(df_ps$value_seg)+5, by = 2))
post_ps <- post_ps |>
  mutate(post = purrr::map(ps, function(x) post_s$beta * x + post_s$alpha)) |>
  mutate(mean = map_dbl(post, function(x) x |> mean()),
         ci_lo = map_dbl(post, function(x) quantile(x, probs = 0.025)),
         ci_hi = map_dbl(post, function(x) quantile(x, probs = 0.975)),
         factor = "Protein S")

pps <- ggplot(df_ps) +
  geom_point(aes(value_seg, value_bag)) +
  geom_line(aes(ps, mean), data = post_ps) +
  geom_ribbon(aes(x = ps, ymin = ci_lo, ymax = ci_hi), data = post_ps, alpha = 0.3) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  facet_wrap(~factor) +
  # y = x diagonal
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  theme_light() +
  xlab("Segment Activity (%)") + ylab("Bag Activity (%)") +
  theme(strip.text = element_text(colour = 'black'))+ theme(axis.text = element_text(size = 7),
                                                            axis.title = element_text(size = 9))
pps


post_c <- extract.samples(m_pc, n = 5000)
post_pc <- tibble(pc = seq(min(df_pc$value_seg), max(df_pc$value_seg)+5, by = 2))
post_pc <- post_pc |>
  mutate(post = purrr::map(pc, function(x) post_c$beta * x + post_c$alpha)) |>
  mutate(mean = map_dbl(post, function(x) x |> mean()),
         ci_lo = map_dbl(post, function(x) quantile(x, probs = 0.025)),
         ci_hi = map_dbl(post, function(x) quantile(x, probs = 0.975)),
         factor = "Protein C")

ppc <- ggplot(df_pc) +
  geom_point(aes(value_seg, value_bag)) +
  geom_line(aes(pc, mean), data = post_pc) +
  geom_ribbon(aes(x = pc, ymin = ci_lo, ymax = ci_hi), data = post_pc, alpha = 0.3) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  facet_wrap(~factor) +
  # y = x diagonal
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  theme_light() +
  xlab("Segment Activity (%)") + ylab("Bag Activity (%)") +
  theme(strip.text = element_text(colour = 'black'))+ theme(axis.text = element_text(size = 7),
                                                            axis.title = element_text(size = 9))
ppc

post_f <- extract.samples(m_fib, n = 5000)
post_fib <- tibble(fib = seq(min(df_fib$value_seg), max(df_fib$value_seg)+5, by = 2))
post_fib <- post_fib |>
  mutate(post = purrr::map(fib, function(x) post_f$beta * x + post_f$alpha)) |>
  mutate(mean = map_dbl(post, function(x) x |> mean()),
         ci_lo = map_dbl(post, function(x) quantile(x, probs = 0.025)),
         ci_hi = map_dbl(post, function(x) quantile(x, probs = 0.975)),
         factor = "Fibrinogen")

pfib <- ggplot(df_fib) +
  geom_point(aes(value_seg, value_bag)) +
  geom_line(aes(fib, mean), data = post_fib) +
  geom_ribbon(aes(x = fib, ymin = ci_lo, ymax = ci_hi), data = post_fib, alpha = 0.3) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  facet_wrap(~factor) +
  # y = x diagonal
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  theme_light() +
  xlab("Segment Activity (mg/dL)") + ylab("Bag Activity (mg/dL)") +
  theme(strip.text = element_text(colour = 'black')) + theme(axis.text = element_text(size = 7),
                                                             axis.title = element_text(size = 9))
pfib



ggsave(filename = "figures/correlation.pdf", egg::ggarrange(pfib, pv, pvii, pviii, pps, ppc),
       width = 3.8,
       height = 5.5,
       units = "in",
       device = "pdf",
       colormode = "cmyk")
