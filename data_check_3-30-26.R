library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)

path <- "height_data.xlsx"

# D1: side stems recorded as counts (num_main + num_small), no individual heights
d1 <- read_excel(path, sheet = "D1 3 2 26", skip = 2) |>
  transmute(
    date,
    plantid,
    main_height,
    num_side_stems = num_main + num_small,
    mean_side_height = NA_real_,
    bush_score = NA_real_
  )

# D2–D5: individual side stem heights recorded as s_1...s_n
read_stem_sheet <- function(sheet) {
  df <- read_excel(path, sheet = sheet)

  # Standardize column names (handles "bush score" -> "bush_score", etc.)
  names(df) <- str_replace_all(str_to_lower(names(df)), " ", "_")

  # Standardize plant id column name
  if ("plant" %in% names(df)) df <- rename(df, plantid = plant)

  # Add bush_score if absent
  if (!"bush_score" %in% names(df)) df$bush_score <- NA_real_

  # Fill date downward (some sheets only fill it on first row)
  df <- fill(df, date)

  s_cols <- names(df)[str_starts(names(df), "s_")]

  df |>
    mutate(
      num_side_stems = rowSums(!is.na(across(all_of(s_cols)))),
      mean_side_height = if_else(
        num_side_stems > 0,
        rowMeans(across(all_of(s_cols)), na.rm = TRUE),
        NA_real_
      )
    ) |>
    select(date, plantid, main_height, num_side_stems, mean_side_height, bush_score)
}

data_sheets <- setdiff(excel_sheets(path), c("check", "meta_data"))
stem_sheets <- setdiff(data_sheets, "D1 3 2 26")

height_data <- bind_rows(d1, map(stem_sheets, read_stem_sheet)) |>
  mutate(
    treatment = case_when(
      str_starts(plantid, "th") ~ "top_half_cut",
      str_starts(plantid, "lc") ~ "leaf_cut",
      str_starts(plantid, "c")  ~ "control"
    ),
    treatment = factor(treatment, levels = c("control", "leaf_cut", "top_half_cut"))
  )

# ── Exploratory plots ──────────────────────────────────────────────────────────

# Main stem height by collection week and treatment group.
# Includes all 5 weeks (D1 = pre-treatment baseline through D5).
# Boxes show the distribution within each group; points are individual plants.
height_data |>
  mutate(week = factor(format(as.Date(date), "%b %d"),
                       levels = format(sort(unique(as.Date(date))), "%b %d"))) |>
  ggplot(aes(x = week, y = main_height, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(aes(),
             position = position_jitterdodge(jitter.width = 0.1),
             alpha = 0.4, size = 1) +
  labs(x = "Collection week", y = "Main stem height (cm)",
       fill = "Treatment",
       title = "Main stem height by week and treatment group") +
  theme_classic() +
  theme(panel.grid = element_blank())

# Main stem height over time by treatment.
# top_half_cut plants stay low (main stem was removed); control and leaf_cut
# grow upward over time, with leaf_cut trending slightly higher.
ggplot(height_data, aes(x = date, y = main_height, color = treatment)) +
  geom_point(alpha = 0.4, position = position_jitter(width = 0.5)) +
  geom_smooth(se = FALSE) +
  labs(x = "Date", y = "Main stem height (cm)", color = "Treatment",
       title = "Main stem height over time by treatment")

# Side stem count over time by treatment.
# top_half_cut plants maintain more side stems throughout — consistent with
# apical dominance being released after the top of the stem was removed.
ggplot(height_data, aes(x = date, y = num_side_stems, color = treatment)) +
  geom_point(alpha = 0.4, position = position_jitter(width = 0.5)) +
  geom_smooth(se = FALSE) +
  labs(x = "Date", y = "Number of side stems", color = "Treatment",
       title = "Side stem count over time by treatment")

# Mean side stem height over time by treatment.
# D1 (Mar 02) has no individual side stem measurements, so this plot
# starts from D2 (Mar 03). control plants have the tallest side stems by late
# March, with leaf_cut and top_half_cut lower and more similar to each other.
ggplot(height_data |> filter(!is.na(mean_side_height)),
       aes(x = date, y = mean_side_height, color = treatment)) +
  geom_point(alpha = 0.4, position = position_jitter(width = 0.5)) +
  geom_smooth(se = FALSE) +
  labs(x = "Date", y = "Mean side stem height (cm)", color = "Treatment",
       title = "Mean side stem height over time by treatment")

# ── Late March outlier check (D5, 2026-03-23) ─────────────────────────────────

# IQR-based flagging: no statistical outliers detected. The high values
# (lc4 = 89.6 cm main height, c11 = 64 cm mean side height) are the tallest
# plants but within the spread of their groups. All tall plants are
# leaf_cut or control — no data entry concerns flagged.

flag_outliers <- function(x) {
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  x > q[2] + 1.5 * iqr
}

height_data |>
  filter(date == max(date)) |>
  mutate(
    outlier_main   = flag_outliers(main_height),
    outlier_side_h = flag_outliers(mean_side_height),
    outlier_side_n = flag_outliers(num_side_stems)
  ) |>
  filter(outlier_main | outlier_side_h | outlier_side_n) |>
  select(date, plantid, treatment, main_height, mean_side_height, num_side_stems,
         starts_with("outlier"))

# Top 10 tallest plants in late March for reference
height_data |>
  filter(date == max(date)) |>
  arrange(desc(main_height)) |>
  select(plantid, treatment, main_height, mean_side_height, num_side_stems) |>
  head(10)

# ── Model data preparation ─────────────────────────────────────────────────────

# D1 (March 2) is pre-treatment; D2 (March 3) is the first post-treatment
# measurement. D1 is therefore used as a baseline covariate to control for
# pre-existing size differences between plants, not as a time point in the
# growth model. Controlling for baseline reduces noise and isolates the
# treatment effect on growth.

# Baseline heights from D1 (pre-treatment)
baseline <- height_data |>
  filter(date == min(date)) |>
  select(plantid, baseline_height = main_height)

# Post-treatment data (D2–D5): day 0 = March 3 (first post-treatment collection)
# Days: D2 = 0, D3 = 6, D4 = 13, D5 = 20
model_data <- height_data |>
  filter(date > min(date)) |>
  mutate(day = as.numeric(as.Date(date) - as.Date("2026-03-03"))) |>
  left_join(baseline, by = "plantid")

# ── Mixed models ───────────────────────────────────────────────────────────────

library(lme4)
library(lmerTest) # Satterthwaite p-values for lmer models

# -- Model 1: Main stem height (control vs. leaf_cut only) ---------------------
# top_half_cut excluded: its main stem was physically removed, so near-zero
# variance in that group violates homoscedasticity and the comparison is
# biologically trivial. The meaningful question is whether leaf damage
# affects main stem elongation in otherwise intact plants.
#
# Fixed effects: baseline height (size covariate) + treatment * day
# Random effect: plant-level intercept (repeated measures)
fit_main_cl <- lmer(
  main_height ~ baseline_height + treatment * day + (1 | plantid),
  data = model_data |> filter(treatment %in% c("control", "leaf_cut"))
)
summary(fit_main_cl)
# Key result: leaf_cut:day = +0.45 cm/day (p = 0.011)
# Leaf-damaged plants grow ~0.45 cm/day faster than controls after treatment.
# Baseline height is a strong predictor (p < 0.001), confirming size differences
# at treatment time warranted inclusion as a covariate.

# -- Model 2: Mean side stem height (all three groups) -------------------------
# D1 has no individual side stem measurements, so no direct baseline for this
# outcome. baseline_height (main stem at D1) included as a plant size proxy.
fit_side_height <- lmer(
  mean_side_height ~ baseline_height + treatment * day + (1 | plantid),
  data = model_data
)
summary(fit_side_height)
# Key result: top_half_cut:day = -0.57 cm/day (p < 0.001)
# Top-half-cut side stems grow more slowly than control side stems.
# leaf_cut:day is not significant (p = 0.11) — leaf damage does not detectably
# alter side stem elongation rate.

# -- Model 3: Side stem count (all three groups) --------------------------------
# Count outcome — overdispersion present in raw data (variance/mean > 1).
# Negative binomial GLMM used; the large dispersion parameter at convergence
# (effectively Poisson) suggests random effects absorb most between-plant
# variance. baseline_side_count (num_main + num_small from D1) included directly.
#
# Note: MASS::glmer.nb masks dplyr::select — use dplyr::select explicitly.
baseline_stems <- d1 |>
  dplyr::select(plantid, baseline_side_count = num_side_stems)

model_data_stems <- model_data |>
  left_join(baseline_stems, by = "plantid")

fit_side_count <- glmer.nb(
  num_side_stems ~ baseline_side_count + treatment * day + (1 | plantid),
  data = model_data_stems
)
summary(fit_side_count)
# Key result: treatment × day interactions not significant for either leaf_cut
# (p = 0.48) or top_half_cut (p = 0.12) after controlling for baseline stem
# count and plant random effects. A small overall decline in side stem count
# across all groups over time (day: p = 0.03).

# -- Model 4: Total side stem height (all three groups) ------------------------
# total_side_height = mean_side_height * num_side_stems (sum of all side stem
# heights per plant per week). Captures combined count + size signal.
# Log-transformed due to right skew; no zeros present so log is clean.
model_data <- model_data |>
  mutate(total_side_height = mean_side_height * num_side_stems)

fit_total_side <- lmer(
  log(total_side_height) ~ baseline_height + treatment * day + (1 | plantid),
  data = model_data
)
summary(fit_total_side)
# Key result: total side stem height increases over time across all groups
# (day: p < 0.001; ~4.4% increase per day on log scale).
# treatment × day interactions are NOT significant: leaf_cut (p = 0.93),
# top_half_cut (p = 0.71). Unlike mean_side_height (where top_half_cut was
# detectably slower), the compound outcome of count × height does not differ
# in growth trajectory between groups. High within-group variance likely
# limits power for this outcome.

# ── Model prediction plots ─────────────────────────────────────────────────────
# All predictions are population-level (re.form = NA) at mean baseline values.

mean_baseline       <- mean(model_data$baseline_height, na.rm = TRUE)
mean_baseline_stems <- mean(model_data_stems$baseline_side_count, na.rm = TRUE)

# Main stem height — control vs. leaf_cut
grid_main <- expand.grid(
  treatment       = factor(c("control", "leaf_cut"),
                           levels = c("control", "leaf_cut", "top_half_cut")),
  day             = seq(0, 20, by = 0.5),
  baseline_height = mean_baseline
)
grid_main$pred <- predict(fit_main_cl, newdata = grid_main, re.form = NA)

ggplot(model_data |> filter(treatment %in% c("control", "leaf_cut")),
       aes(x = day, y = main_height, color = treatment)) +
  geom_point(alpha = 0.3) +
  geom_line(data = grid_main, aes(y = pred), linewidth = 1) +
  labs(x = "Days post-treatment", y = "Main stem height (cm)",
       color = "Treatment",
       title = "Main stem height — model predictions vs. raw data",
       subtitle = "Predictions at mean baseline height; points are individual plants")

# Mean side stem height — all three groups
grid_side_h <- expand.grid(
  treatment       = factor(c("control", "leaf_cut", "top_half_cut"),
                           levels = c("control", "leaf_cut", "top_half_cut")),
  day             = seq(0, 20, by = 0.5),
  baseline_height = mean_baseline
)
grid_side_h$pred <- predict(fit_side_height, newdata = grid_side_h, re.form = NA)

ggplot(model_data |> filter(!is.na(mean_side_height)),
       aes(x = day, y = mean_side_height, color = treatment)) +
  geom_point(alpha = 0.3) +
  geom_line(data = grid_side_h, aes(y = pred), linewidth = 1) +
  labs(x = "Days post-treatment", y = "Mean side stem height (cm)",
       color = "Treatment",
       title = "Mean side stem height — model predictions vs. raw data",
       subtitle = "Predictions at mean baseline height; points are individual plants")

# Side stem count — all three groups (predictions back-transformed from log scale)
grid_side_n <- expand.grid(
  treatment           = factor(c("control", "leaf_cut", "top_half_cut"),
                               levels = c("control", "leaf_cut", "top_half_cut")),
  day                 = seq(0, 20, by = 0.5),
  baseline_side_count = mean_baseline_stems
)
grid_side_n$pred <- predict(fit_side_count, newdata = grid_side_n,
                            re.form = NA, type = "response")

ggplot(model_data_stems, aes(x = day, y = num_side_stems, color = treatment)) +
  geom_point(alpha = 0.3) +
  geom_line(data = grid_side_n, aes(y = pred), linewidth = 1) +
  labs(x = "Days post-treatment", y = "Number of side stems",
       color = "Treatment",
       title = "Side stem count — model predictions vs. raw data",
       subtitle = "Predictions at mean baseline count; points are individual plants")

# Total side stem height — all three groups (back-transformed from log scale)
grid_total <- expand.grid(
  treatment       = factor(c("control", "leaf_cut", "top_half_cut"),
                           levels = c("control", "leaf_cut", "top_half_cut")),
  day             = seq(0, 20, by = 0.5),
  baseline_height = mean_baseline
)
grid_total$pred <- exp(predict(fit_total_side, newdata = grid_total, re.form = NA))

ggplot(model_data |> filter(!is.na(total_side_height)),
       aes(x = day, y = total_side_height, color = treatment)) +
  geom_point(alpha = 0.3) +
  geom_line(data = grid_total, aes(y = pred), linewidth = 1) +
  labs(x = "Days post-treatment", y = "Total side stem height (cm)",
       color = "Treatment",
       title = "Total side stem height — model predictions vs. raw data",
       subtitle = "Predictions at mean baseline height; points are individual plants") +
  theme_classic()

# ── Bush score confounding check ───────────────────────────────────────────────
# bush_score is only recorded in D4 and D5. If bushy plants cluster in one
# treatment group, it could confound the later time point results.

# Counts of bush_score by treatment
height_data |>
  filter(!is.na(bush_score)) |>
  count(treatment, bush_score) |>
  tidyr::pivot_wider(names_from = bush_score, values_from = n,
                     values_fill = 0, names_prefix = "score_")

# Chi-squared test of independence between bush_score and treatment
# p = 0.20 — no evidence that bush phenotype is unevenly distributed across
# treatments. top_half_cut has somewhat fewer bushy plants (5 vs 9–11),
# which is biologically plausible given stem removal, but not significant.
bush_tab <- height_data |>
  filter(!is.na(bush_score)) |>
  mutate(bush_score = factor(bush_score)) |>
  with(table(treatment, bush_score))

chisq.test(bush_tab)

# ── Biological interpretation note ────────────────────────────────────────────
# The leaf_cut acceleration in main stem growth (+0.45 cm/day, ~40% faster
# than control) is consistent with a compensatory growth response. Partial
# defoliation is known to trigger meristematic growth in some species via
# jasmonate signaling and altered source-sink dynamics — the plant may
# prioritize rapid stem elongation to restore photosynthetic capacity.
# This is worth noting as a hypothesis, but 3 weeks of data and n ~13/group
# warrants caution before strong claims. Follow-up collections will be
# informative.
