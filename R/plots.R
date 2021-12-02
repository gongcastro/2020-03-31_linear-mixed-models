library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(gganimate)
library(purrr)
library(stringr)
library(truncnorm)
library(rlang)
library(lme4)
library(wesanderson)
library(here)

# load functions
source(here("R", "utils.R"))

theme_set(theme_custom())

# generate naming reaction times (RT) ----
dat <- simulate_rts(
  intercepts = c(0.1, 0.3, 0.55, 0.9, 1.5),
  sd = 0.01,
  corr = 2,
  n_conditions = 2,
  n_trials = 5
)

# aggregatted data
dat_sum <- dat %>%
  group_by(participant, frequency) %>%
  summarise(rt = mean(rt), .groups = "drop") %>%
  ungroup()


# fit linear models ----
fit_all <- lm(rt ~ frequency, data = dat) # only fixed effects
fit_sum <- lm(rt ~ frequency, data = dat_sum) # only fixed effects (aggregated data)
fit_int <- lmer(rt ~ frequency + (1 | participant), data = dat) # fixed effects + random intercepts
fit_slo <- lmer(rt ~ frequency + (1 + frequency|participant), data = dat) # fixed effects + random intercepts and slopes


# visualise predictions of fixed effects-only model ----
ggplot(dat, aes(frequency, rt)) +
  geom_jitter(alpha = 0.8, size = 3, width = 0.1) +
  stat_summary(
    aes(y = predict(fits$all), group = 1), fun.data = "mean_se", 
    geom = "ribbon", alpha = 0.5, colour = NA
  ) +
  stat_summary(
    aes(y = predict(fits$all), group = 1), fun = "mean", geom = "line", 
    size = 1, colour = "black"
  ) +
  stat_summary(aes(group = 1), fun = "mean", geom = "line", size = 1) +
  geom_text(
    x = 1.5, y = 2, 
    label = paste0("Intercept: ", round(fit_all$coefficients[1], 2), "\nSlope: ", round(fit_all$coefficients[2], 2))
  ) +
  labs(x = "Frequency", y = "Reaction time (s)", colour = "Participant", subtitle = "No random effects") +
  scale_colour_manual(values = wes_palette("Darjeeling1", 5, type = "discrete"))

ggsave(here("img", "plot_model.png"))


# visualise predictions of fixed effects-only model fitted on aggregated data ---
dat %>%
  group_by(frequency, participant) %>%
  summarise(rt = mean(rt), .groups = "drop") %>%
  ggplot(aes(frequency, rt, colour = participant)) +
  geom_jitter(alpha = 0.8, size = 3, width = 0.1) +
  stat_summary(aes(y = predict(fit_sum), group = 1), fun.data = "mean_se", geom = "ribbon", alpha = 0.5, colour = NA) +
  stat_summary(aes(y = predict(fit_sum), group = 1), fun = "mean", geom = "line", size = 1) +
  stat_summary(aes(y = fitted(fit_sum)), fun = "mean", geom = "point", size = 4, colour = "black") +
  annotate(
    x = 1.5, y = 2, geom = "text",
    label = paste0("Intercept: ", round(model$coefficients[1], 2), "\nSlope: ", round(model$coefficients[2], 2))
  ) +
  labs(x = "Frequency", y = "Reaction time (s)", colour = "Participant", subtitle = "No random effects, aggregated by participant") +
  scale_colour_manual(values = wes_palette("Darjeeling1", 5, type = "discrete"))

ggsave(here("img", "plot_model_agg.png"))


# visualise predicitons of model with random intercepts ----

# get model coefficients
fit_int_coefs <- data.frame(
  participant = paste0("participant", 1:5),
  x1 = 2.3,
  y1 = seq(1, 3, length.out = 5),
  intercept = round(coef(fit_int)$participant[, 1], 2),
  slope = round(coef(fit_int)$participant[, 2], 2)
)

ggplot(dat, aes(frequency, rt, colour = participant)) +
  geom_jitter(alpha = 0.8, size = 3, width = 0.1) +
  stat_summary(aes(y = predict(fit_int), group = participant), fun = "mean", geom = "point", size = 4) +
  stat_summary(aes(y = predict(fit_int), group = participant), fun = "mean", geom = "line", size = 1) +
  stat_summary(aes(y = predict(fit_int), group = 1), fun.data = "mean_se", geom = "ribbon", colour = NA, alpha = 0.5) +
  stat_summary(aes(y = predict(fit_int), group = 1), fun = "mean", geom = "point", size = 4) +
  stat_summary(aes(y = predict(fit_int), group = 1), fun = "mean", geom = "line", size = 1) +
  annotate(
    x = 1.5, y = 3, colour = "black", geom = "label",
    label = paste0("Intercept: ", round(fixef(fit_int)[1], 2), "\nCoefficient: ", round(fixef(fit_int)[2], 2)),
  ) +
  geom_text(data = fit_int_coefs, aes(x = x1, y = y1, label = paste0("Intercept: ", intercept, "\nSlope: ", slope)),
            colour = "black", show.legend = FALSE, size = 3) +
  labs(x = "Frequency", y = "Reaction time (s)", colour = "Participant",
       subtitle = "Random intercepts by participant") +
  scale_colour_manual(values = wes_palette("Darjeeling1", 5, type = "discrete"))

ggsave(here("img", "plot_model_int.png"))


# visualise predictions of model with random intercepts and slopes ----

# get model coefficients
fit_slo_coefs <- data.frame(
  participant = paste0("participant", 1:5),
  x1 = 2.3,
  y1 = seq(1, 3.5, length.out = 5),
  intercept = round(coef(fit_slo)$participant[, 1], 2),
  slope = round(coef(fit_slo)$participant[, 2], 2)
)

ggplot(dat, aes(frequency, rt, colour = participant)) +
  geom_jitter(alpha = 0.8, size = 3, width = 0.1) +
  stat_summary(aes(y = fitted(fit_slo), group = participant), fun = "mean", geom = "point", size = 4) +
  stat_summary(aes(y = fitted(fit_slo), group = participant), fun = "mean", geom = "line", size = 1) +
  stat_summary(aes(y = predict(fit_slo), group = 1), fun.data = "mean_se", geom = "ribbon", colour = NA, alpha = 0.5) +
  
  stat_summary(aes(y = fitted(fit_slo), group = 1), fun = "mean", geom = "point", size = 4) +
  stat_summary(aes(y = fitted(fit_slo), group = 1), fun = "mean", geom = "line", size = 1) +
  annotate(
    x = 1.5, y = 3, colour = "black", geom = "label",
    label = paste0("Intercept: ", round(fixef(fit_slo)[1], 2), "\nCoefficient: ", round(fixef(fit_slo)[1], 2))) +
  geom_text(
    data = model_slo_coefs,
    aes(x = x1, y = y1, label = paste0("Intercept: ", intercept, "\nSlope: ", slope)), 
    colour = "black", show.legend = FALSE, size = 3
  ) +
  labs(x = "Frequency", y = "Reaction time (s)", colour = "Participant", subtitle = "Random intercepts by participant") +
  scale_colour_manual(values = wes_palette("Darjeeling1", 5, type = "discrete"))

ggsave(here("img", "plot_model_slo.png"))

# create animation ----

# model coefficients
coefs <- map(list(fit_all, fit_int, fit_slo), summary) %>%
  map(get_coefs) %>%
  set_names(c("Only fixed effects", "Adding random intercepts", "Adding random intercepts and slopes")) %>%
  bind_rows(.id = "model")

# get model predictions
predictions <- dat %>%
  mutate(
    fit_all = predict(fit_all, newdata = .),
    fit_int = predict(fit_int, newdata = .),
    fit_slo = predict(fit_slo, newdata = .)
  ) %>%
  pivot_longer(fit_all:fit_slo, names_to = "model", values_to = "fit") %>%
  mutate(
    model = case_when(
      model=="fit_all" ~ "Only fixed effects",
      model=="fit_int" ~ "Adding random intercepts",
      model=="fit_slo" ~ "Adding random intercepts and slopes"
    )
  )

# animate plots
predictions_animate <- ggplot(predictions, aes(x = frequency, y = rt)) +
  geom_line(aes(y = fit, group = participant, colour = participant)) +
  geom_point(aes(colour = participant), alpha = 0.8, size = 4, position = position_jitter(width = 0.1, seed = 888)) +
  geom_text(data = coefs, aes(x = 0.5, y = 4.2, label = model), hjust = 0, show.legend = FALSE, size = 6) + 
  stat_summary(aes(y = fit, group = 1), fun.data = "mean_se", geom = "ribbon", colour = NA, alpha = 0.5) +
  stat_summary(aes(y = fit, group = 1), fun = "mean", geom = "line", size = 2) +
  stat_summary(aes(y = fit, group = 1), fun = "mean", geom = "point", size = 4) +
  labs(x = "Frequency", y = "Reaction time (s)", colour = "Participant") +
  scale_colour_manual(values = wes_palette("Darjeeling1", 5, type = "discrete")) +
  transition_states(states = model, transition_length = 3, state_length = 2, wrap = TRUE)

animate(
  predictions_animate, height = 500, width = 500, 
  fps = 30, rewind = TRUE, duration = 10, device = "png"
)

anim_save(filename = here("img", "plot_model.gif"))
