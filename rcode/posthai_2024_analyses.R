library(tictoc)
library(tidyverse)
library(dplyr)
library(gridExtra)
library(patchwork)
library(broom)
library(performance)
library(easystats)
library(gtsummary)
library(emmeans)
library(dlookr)
library(boot)
library(lme4)
library(lmerTest)
library(mgcv)
library(report)
library(gratia)
library(ggplot2)
library(sjPlot)
library(see)
library(viridis)
library(patchwork)
library(effectsize)
library(parameters)
library(boot)

# Define default parameters.
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
save_dir <- "/Users/dschonhaut/box/projects/leads_tau_spread/analysis/r_models/gam/idv_rois"
today <- format(Sys.Date(), "%Y-%m-%d")
overwrite = FALSE
verbose = TRUE
tau_subjs <- c(
  'LDS0070166', 'LDS0070174', 'LDS0070199', 'LDS0070263',
  'LDS0070302', 'LDS0070313', 'LDS0070327', 'LDS0070328',
  'LDS0070344', 'LDS0070410', 'LDS0070417', 'LDS0070426',
  'LDS0070428', 'LDS0070472', 'LDS0100297', 'LDS0100330',
  'LDS0110021', 'LDS0110022', 'LDS0110040', 'LDS0110041',
  'LDS0110078', 'LDS0110079', 'LDS0110219', 'LDS0110257',
  'LDS0110295', 'LDS0110458', 'LDS0110481', 'LDS0180096',
  'LDS0180153', 'LDS0180160', 'LDS0180177', 'LDS0180214',
  'LDS0220031', 'LDS0220071', 'LDS0220159', 'LDS0220202',
  'LDS0220252', 'LDS0220273', 'LDS0220290', 'LDS0220334',
  'LDS0220429', 'LDS0220432', 'LDS0220502', 'LDS0350171',
  'LDS0350210', 'LDS0350390', 'LDS0360121', 'LDS0360283',
  'LDS0360291', 'LDS0360320', 'LDS0360324', 'LDS0360331',
  'LDS0360354', 'LDS0370006', 'LDS0370007', 'LDS0370008',
  'LDS0370010', 'LDS0370019', 'LDS0370038', 'LDS0370042',
  'LDS0370061', 'LDS0370065', 'LDS0370089', 'LDS0370094',
  'LDS0370097', 'LDS0370100', 'LDS0370108', 'LDS0370110',
  'LDS0370147', 'LDS0370163', 'LDS0370167', 'LDS0370172',
  'LDS0370180', 'LDS0370221', 'LDS0370222', 'LDS0370239',
  'LDS0370285', 'LDS0370286', 'LDS0370304', 'LDS0370321',
  'LDS0370325', 'LDS0370329', 'LDS0370335', 'LDS0370338',
  'LDS0370343', 'LDS0370392', 'LDS0370397', 'LDS0370399',
  'LDS0370400', 'LDS0370427', 'LDS0370435', 'LDS0370449',
  'LDS0370454', 'LDS0370465', 'LDS0370516', 'LDS0370533',
  'LDS0670111', 'LDS0730044', 'LDS0730055', 'LDS0730103',
  'LDS0730141', 'LDS0730146', 'LDS0730150', 'LDS0730161',
  'LDS0730228', 'LDS0730280', 'LDS0730284', 'LDS0730353',
  'LDS0730372', 'LDS0730402', 'LDS0730422', 'LDS0730474',
  'LDS0730478', 'LDS0730521', 'LDS0980440', 'LDS0990406',
  'LDS0990407', 'LDS0990408', 'LDS0990409', 'LDS0990430',
  'LDS0990444', 'LDS0990493', 'LDS1770107', 'LDS1770112',
  'LDS1770129', 'LDS1770145', 'LDS1770181', 'LDS1770209',
  'LDS1770213', 'LDS1770223', 'LDS1770241', 'LDS1770266',
  'LDS1770292', 'LDS1770293', 'LDS1770300', 'LDS1770305',
  'LDS1770310', 'LDS1770362', 'LDS1770375', 'LDS1770394',
  'LDS1770413', 'LDS1770501', 'LDS3600030', 'LDS3600032',
  'LDS3600056', 'LDS3600114', 'LDS3600115', 'LDS3600116',
  'LDS3600133', 'LDS3600140', 'LDS3600187', 'LDS3600262',
  'LDS3600265', 'LDS3600281', 'LDS3600326', 'LDS3600350',
  'LDS3600361', 'LDS3600401', 'LDS3600415', 'LDS3600419',
  'LDS3600421', 'LDS3600455', 'LDS3600484', 'LDS9410149',
  'LDS9410298'
)

# Load dataframes.
tau_suvrsf <- file.path(data_dir, "tau-rois-agg_eoad-long_formatted_2024-01-24.csv")
tau_suvrs <- read_csv(tau_suvrsf)
tau_suvrs_widef <- file.path(data_dir, "tau-rois-agg_eoad-wide_formatted_2023-12-28.csv")
tau_suvrs_wide <- read_csv(tau_suvrs_widef)
tau_suvrs_wide$parc <- "metarois"

# Clean dataframes.
transform_df <- function(df) {
  df %>%
    # filter(parc %in% c("metarois")) %>%
    mutate(across(where(is.character), as.factor)) %>%
    mutate(apoe4_allelesf = factor(apoe4_alleles, ordered = FALSE)) %>%
    mutate(apoe4_allelesof = factor(apoe4_alleles, ordered = TRUE)) %>%
    { contrasts(.$sex) <- 'contr.sum'; . }
}
tau_suvrs <- transform_df(tau_suvrs)
tau_suvrs_wide <- transform_df(tau_suvrs_wide)

# Create binned variable columns.
process_cats <- function(df) {
  # Function to bin age
  bin_age <- function(age) {
    age <- as.integer(age)
    if (age >= 45 & age < 55) {
      return("45-54")
    } else if (age >= 55 & age < 60) {
      return("55-59")
    } else if (age >= 60 & age <= 65) {
      return("60-65")
    } else {
      return("?")
    }
  }

  # Function to bin fbb_cl
  bin_fbb_cl <- function(cl) {
    if (cl < 80) {
      return("<80")
    } else if (cl >= 80 & cl <= 120) {
      return("80-120")
    } else {
      return(">120")
    }
  }

  # Function to bin mmse
  bin_mmse <- function(mmse) {
    if (mmse >= 25) {
      return("25-30")
    } else if (mmse < 25 & mmse >= 20) {
      return("20-24")
    } else {
      return("<20")
    }
  }

  # Applying the functions and converting to categorical variables
  df %>%
    mutate(age_at_ftp_bl_bin = sapply(age_at_ftp_bl, bin_age),
           fbb_cl_bl_bin = sapply(fbb_cl_bl, bin_fbb_cl),
           mmse_bl_bin = sapply(mmse_bl, bin_mmse)) %>%
    mutate(fbb_cl_bl_bin = fct_relevel(fbb_cl_bl_bin, c("<80", "80-120", ">120")),
           mmse_bl_bin = fct_relevel(mmse_bl_bin, c("<20", "20-24", "25-30")))
}

tau_suvrs <- process_cats(tau_suvrs)
tau_suvrs_wide <- process_cats(tau_suvrs_wide)

# Merge the tau_suvrs (1 row per subj/visit/roi) dataframe with the mixed model fits.
right_cols <- c("subj", "parc", "roi", "fixed_icpt", "fixed_time",
                "rdm_icpt", "rdm_time", "icpt", "time")
if (!all(c("icpt", "time") %in% colnames(tau_suvrs))) {
  tau_suvrs <- tau_suvrs %>%
    left_join(tau_mms %>% select(all_of(right_cols)),
              by = c("subj", "parc", "roi"))
}
if (!all(c("icpt", "time") %in% colnames(tau_suvrs))) {
  tau_suvrs_wide <- tau_suvrs_wide %>%
    left_join(tau_suvrs_wide %>% select(all_of(right_cols)),
              by = c("subj", "parc", "roi"))
}

tau_suvrs2 <- tau_suvrs %>%
  filter(visit > 1) %>%
  mutate(suvr_last = suvr - suvr_chg_from_last) %>%
  group_by(subj, roi) %>%
  summarize(
    ftp_visits = unique(ftp_visits),
    suvr_bl = unique(suvr_bl),
    suvr_last = mean(suvr_last),
    suvr_annchg_from_last = mean(suvr_annchg_from_last),
    icpt = unique(icpt),
    time = unique(time),
    age_at_ftp_bl_mc = unique(age_at_ftp_bl_mc),
    fbb_cl_bl_mc = unique(fbb_cl_bl_mc),
    cdr_sb_bl_mc = unique(cdr_sb_bl_mc),
    apoe4_alleles_mc = unique(apoe4_alleles_mc),
    sex = unique(sex)
  ) %>%
  ungroup()

# LMEMs
df <- tau_suvrs %>% filter(roi %in% rois) %>% tibble()
xmod <- lmer(
  suvr ~ roi * ftp_yrs_from_bl + (roi + ftp_yrs_from_bl|subj),
  data=df, REML=T
)

joint_tests(xmod)
pairs(emmeans(xmod, ~ roi * ftp_yrs_from_bl, at = list(ftp_yrs_from_bl = c(0))))
emm <- emmeans(xmod, pairwise ~ ftp_yrs_from_bl | roi, at = list(ftp_yrs_from_bl = c(0, 1)))
emt <- emtrends(xmod, pairwise ~ roi, var = "ftp_yrs_from_bl", adjust="none")
plot(emt, comparisons=T, adjust="none") + theme_classic()
pwpp(emt, adjust="none") + theme_classic()
emmip(xmod, roi ~ ftp_yrs_from_bl, cov.reduce = range) + theme_classic()


# Visualize data
esquisser(tau_suvrs)
esquisser(tau_suvrs2)

model <- lmer(suvr ~ ftp_yrs_from_bl + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
              data = df_roi)

#!! Start analysis for AD/PD 2024
df <- tau_suvrs
antiamy_subjs = c(
  "LDS0360378",
  "LDS9410139",
  "LDS9410287",
  "LDS9410376",
  "LDS9410396",
  "LDS9410450",
  "LDS9410459"
)
df %>% filter(!subj %in% antiamy_subjs) %>% dim()
contrasts(tau_suvrs$sex) <- "contr.sum"

roi <- "temporal"
df_roi <- df %>% filter(
  dx == "EOAD",
  ftp_visits > 1,
  roi == "temporal"
)

model <- lmer(suvr ~ ftp_yrs_from_bl + (1 + ftp_yrs_from_bl | subj), data = df_roi, REML = TRUE)
summary(model)
confint(model, method = "Wald")
report(model)
#!! End analysis for AD/PD 2024

# Age, sex, etc.
# Define inputs
df <- as.tibble(tau_suvrs_wide)
df$sex <- factor(as.character(df$sex), levels=c("m", "f"))

bs <- "cs"
rois <- c("frontal", "parietal", "temporal")
covs <- c("age_at_ftp_bl_mc", "fbb_cl_bl_mc", "severity_bl", "converted", "cdr_sb_bl_mc", "mmse_bl_mc", "apoe4_alleles_mc", "sex")

# Fit models for each ROI and covariate combination using nested map functions
models_bl <- expand.grid(roi=rois, cov=covs) %>%
  mutate(
    dv = paste0("suvr_bl_", roi),
    formula = paste0(dv, "~", cov),
    model = pmap(list(formula, list(df)), ~ lm(as.formula(..1), data = ..2))
  )
models_gam <- expand.grid(roi=rois, cov=covs) %>%
  mutate(
    dv = paste0("suvr_annchg_from_last_", roi),
    iv = paste0("s(suvr_bl_", roi, ", bs=bs)"),
    formula = paste0(dv, "~", iv, "+", cov),
    model = pmap(list(formula, list(df)), ~ gam(as.formula(..1), data = ..2, method="REML"))
  )
models_bl$model %>%
  map(summary)

# *****
# Fit models for each ROI and covariate combination using nested map functions
bs2 <- "re"
covs2 <- c("age_at_ftp_bl_bin", "apoe4_allelesf", "fbb_cl_bl_bin", "mmse_bl_bin", "sex")
models_gam2 <- expand.grid(roi=rois, cov=covs2) %>%
  mutate(
    dv = paste0("suvr_annchg_from_last_", roi),
    iv1 = paste0("s(suvr_bl_", roi, ", bs=bs)"),
    iv2 = paste0("s(suvr_bl_", roi, ", by=", cov, ", m=1, id=1", ", bs=bs)"),
    iv3 = paste0("s(", cov, ", bs=bs2)"),
    formula = paste0(dv, "~", iv1, " + ", iv2, " + ", iv3),
    model = pmap(list(formula, list(df)), ~ gam(as.formula(..1), data = ..2, method="REML"))
  )

models_gam2$model %>%
  map(summary)


draw(models_gam2$model[[3]])

# *****

# -----
models_bl_full <- expand.grid(roi=rois) %>%
  mutate(
    dv = paste0("suvr_bl_", roi),
    cov = paste("age_at_ftp_bl_mc", "fbb_cl_bl_mc", "cdr_sb_bl_mc", "apoe4_alleles_mc", "sex", sep=" + "),
    formula = paste0(dv, "~", cov),
    model = pmap(list(formula, list(df)), ~ lm(as.formula(..1), data = ..2))
  )
models_gam_full <- expand.grid(roi=rois) %>%
  mutate(
    dv = paste0("suvr_annchg_from_last_", roi),
    iv = paste0("s(suvr_bl_", roi, ", bs=bs)"),
    cov = paste("age_at_ftp_bl_mc", "fbb_cl_bl_mc", "cdr_sb_bl_mc", "apoe4_alleles_mc", "sex", sep=" + "),
    formula = paste0(dv, "~", iv, "+", cov),
    model = pmap(list(formula, list(df)), ~ gam(as.formula(..1), data = ..2, method="REML"))
  )
models_cog_full <- expand.grid(roi=rois) %>%
  mutate(
    dv = "cdr_sb_annchg",
    iv1 = "s(cdr_sb_bl_mc, bs=bs)",
    iv2 = paste0("s(suvr_bl_", roi, ", bs=bs)"),
    iv3 = paste0("s(suvr_annchg_from_last_", roi, ", bs=bs)"),
    cov = paste("age_at_ftp_bl_mc", "fbb_cl_bl_mc", "apoe4_alleles_mc", "sex", sep=" + "),
    formula = paste0(dv, "~", iv1, "+", cov),
    model = pmap(list(formula, list(df)), ~ gam(as.formula(..1), data = ..2, method="REML"))
  )
# Print the summary info for each fitted model
ii = 1; format_table(model_parameters(models_bl_full$model[[ii]], ci_method="wald"), digits=4, ci_digits=4, p_digits=3, stars=T) %>% mutate(roi=models_bl_full$roi[[ii]])

models_gam_full$model %>%
  map(summary)
ii = 3; format_table(model_parameters(models_gam_full$model[[ii]], ci_method="wald", drop="Intercept|suvr"), digits=4, ci_digits=4, p_digits=3, stars=T) %>% mutate(roi=models_gam_full$roi[[ii]]) %>% select(-Component)

models_cog_full$model %>%
  map(summary)
ii = 1; format_table(model_parameters(models_cog_full$model[[ii]], ci_method="wald", drop="Intercept|suvr"), digits=4, ci_digits=4, p_digits=3, stars=T) %>% mutate(roi=models_cog_full$roi[[ii]]) %>% select(-Component)

summary(models_gam_full$model[[ii]])

map(models_bl_full$model, ~ {
  params <- model_parameters(.x, ci_method = "wald")
  params <- standardize_names(params)
  formatted_table <- format_table(params, digits = 4, ci_digits=4, stars = T)
  as.data.frame(formatted_table)
}) %>% bind_rows() %>% select(-Component) %>% View()

map(models_gam$model, ~ {
  params <- model_parameters(.x, ci_method = "wald", drop="Intercept|suvr")
  params <- standardize_names(params)
  formatted_table <- format_table(params, digits = 4, ci_digits=4, stars = T)
  as.data.frame(formatted_table)
}) %>% bind_rows() %>% select(-Component) %>% View()

map(models_gam_full$model, ~ {
  params <- model_parameters(.x, ci_method = "wald", drop="Intercept|suvr")
  params <- standardize_names(params)
  formatted_table <- format_table(params, digits = 4, ci_digits=4, stars = T)
  as.data.frame(formatted_table)
}) %>% bind_rows() %>% select(-Component) %>% View()

map(models_cog_full$model, ~ {
  params <- model_parameters(.x, ci_method = "wald", drop="Intercept")
  params <- standardize_names(params)
  formatted_table <- format_table(params, digits = 4, ci_digits=4, stars = T)
  as.data.frame(formatted_table)
}) %>% bind_rows() %>% select(-Component) %>% View()

m <- models_gam_full$model[[3]]

# Graph model results
colors <- c(temporal = "#2E45B8", parietal = "#3EBCD2", frontal = "#FF4983")

plot_model_results <- function(models, var_name, df) {
  # Preparing data for plotting'
  data_list <- lapply(models, function(model) {
    df$predicted <- predict(model, newdata = df, type = "terms")[, var_name]
    df$residuals <- residuals(model)
    df$partial_residuals <- df$residuals + df$predicted
    df$region <- as.character(model$call$data)[2]  # Assuming the region is the second argument in the model formula
    return(df)
  })
  return(df)
}

plot_data <- do.call(rbind, data_list)

# Bootstrapping for confidence intervals
get_ci <- function(model, var_name, data) {
  boot_fun <- function(data, indices) {
    d <- data[indices, ]
    fit <- gam(as.formula(model$call$formula), data = d, method = "REML")
    predict(fit, newdata = data, type = "terms")[, var_name]
  }
  results <- boot(data, boot_fun, R = 1000)  # Adjust R for more or fewer bootstrap iterations
  boot.ci(results, type = "bca")
}

ci_list <- lapply(models, function(model) get_ci(model, var_name, df))
ci_data <- do.call(rbind, ci_list)
ci_data <- mutate(ci_data, region = names(models))

# Plotting
gg <- ggplot(plot_data, aes_string(x = var_name, y = "partial_residuals", color = "region")) +
  geom_point(alpha = 0.5) +
  geom_line(data = ci_data, aes_string(y = "fit", ymin = "lower", ymax = "upper"), stat = "smooth", method = "loess", se = TRUE) +
  scale_color_manual(values = colors) +
  labs(title = paste("Partial Residuals and Model Slopes for", var_name), x = var_name, y = "Partial Contribution") +
  theme_minimal()

return(gg)
}

var_name = "age_at_ftp_bl_mc"
ggplot_obj <- plot_model_results(models_gam_full$model, var_name, df[complete.cases(df$apoe4_alleles),])
print(ggplot_obj)

## Model longitudinal CDR-SB
cdrf <- file.path(data_dir, "cdr-sb_eoad-long_formatted_165subjs_2024-01-24.csv")
cdr <- read_csv(cdrf)
cdr$subj <- factor(cdr$subj)
cdr$visit <- factor(cdr$visit, levels=c(1, 2, 3, 4, 5), ordered=TRUE)
cdr$sex <- factor(cdr$sex)
contrasts(cdr$sex) <- "contr.sum"
cdr$cdr_sb <- as.numeric(cdr$cdr_sb)
cdr$cdr_sb_last <- as.numeric(cdr$cdr_sb_last)
cdr$cdr_sb_bl <- as.numeric(cdr$cdr_sb_bl)
cdr$cdr_yrs_from_bl <- as.numeric(cdr$cdr_yrs_from_bl)
cdr$cdr_sb_annchg_from_last <- as.numeric(cdr$cdr_sb_annchg_from_last)
cdr$apoe4_alleles <- as.numeric(cdr$apoe4_alleles)
cdr$apoe4_allelesf <- factor(cdr$apoe4_alleles, levels=c(0, 1, 2), ordered=FALSE)
cdr$apoe4_allelesof <- factor(cdr$apoe4_alleles, levels=c(0, 1, 2), ordered=TRUE)
cdr$age_at_cdr <- as.numeric(cdr$age_at_cdr)
cdr$age_at_cdr_bl <- as.numeric(cdr$age_at_cdr_bl)
cdr$age_at_cdr_last <- as.numeric(cdr$age_at_cdr_last)
mean_age_cdr <- (cdr %>% filter(visit>1) %>% summarize(age_at_cdr_last = mean(age_at_cdr_last)))$age_at_cdr_last
cdr$age_at_cdr_last_mc <- cdr$age_at_cdr_last - mean_age_cdr

cdr_mod <- gam(
  cdr_sb_annchg_from_last ~ s(cdr_sb_last, bs="cs") + s(subj, bs="re"),
  data=cdr %>% filter(visit>1),
  method="REML"
)

plot(cdr_mod, select=1, residuals=T, shift=fixef(cdr_mod)["(Intercept)"], xlim=c(0, 18), ylim=c(-2, 10), seWithMean=T, rug=F)
plot(cdr_mod_full, select=1, residuals=T, shift=fixef(cdr_mod_full)["(Intercept)"], xlim=c(0, 18), ylim=c(-2, 10), seWithMean=T, rug=F)
plot(cdr_mod_full2, select=1, residuals=T, shift=fixef(cdr_mod_full)["(Intercept)"], xlim=c(0, 18), ylim=c(-2, 10), seWithMean=T, rug=F)
plot(cdr_mod_full, select=2, residuals=F, shift=fixef(cdr_mod_full)["(Intercept)"], xlim=c(0, 18), ylim=c(1, 3.5), seWithMean=T, rug=F)
plot(cdr_mod_full, select=3)
plot(cdr_mod_full, select=4)
plot(cdr_mod_full, select=5)

cdr_mod2 <- gamm(
  cdr_sb_annchg_from_last ~ s(cdr_sb_last, bs="cs"),
  random=list(subj=~1, visit=~1),
  data=cdr %>% filter(visit>1),
  method="REML"
)

cdr_mod_full <- gam(
  cdr_sb_annchg_from_last ~ age_at_cdr_last_mc + apoe4_alleles + sex + s(cdr_sb_last, bs="cs") + s(subj, bs="re"),
  data=cdr %>% filter(visit>1),
  method="REML"
)

cdr_mod_full2 <- gam(
  cdr_sb_annchg_from_last ~ s(cdr_sb_last, bs="cs") + s(age_at_cdr_last, bs="cs") + s(apoe4_alleles, k=2, bs="cs") + s(sex, bs="re") + s(subj, bs="re"),
  data=cdr %>% filter(visit>1),
  method="REML"
)

cdr_mod_full3 <- gamm(
  cdr_sb_annchg_from_last ~ age_at_cdr_last + apoe4_alleles + sex + s(cdr_sb_last, bs="cs"),
  random=list(subj=~1, visit=~1),
  data=cdr %>% filter(visit>1),
  method="REML"
)

## dlaa;d
metarois <- c("temporal", "parietal", "frontal")
df <- tau_suvrs %>%
  filter(
    subj %in% tau_subjs,
    roi %in% metarois,
    severity_bl %in% c("MCI", "Dementia"),
    dx == "EOAD",
    ftp_visits > 1
  )
df$subj <- factor(as.character(df$subj))
df$visit <- factor(df$visit, levels=c(1, 2, 3, 4), ordered=TRUE)
mean_age <- (df %>% filter(roi=="frontal", visit>1) %>% summarize(age_at_ftp_last = mean(age_at_ftp_last)))$age_at_ftp_last
df$age_at_ftp_last_mc <- df$age_at_ftp_last - mean_age
# df$sex <- factor(df$sex)

## Model ROIs like CDR...
roi_mods_bl_full <- list(
  lm(suvr ~ age_at_ftp_bl_mc + apoe4_alleles_mc + cdr_sb_bl_mc + fbb_cl_bl_mc + sex,
     data=tau_suvrs %>% filter(roi=="frontal", visit==1)),
  lm(suvr ~ age_at_ftp_bl_mc + apoe4_alleles_mc + cdr_sb_bl_mc + fbb_cl_bl_mc + sex,
     data=tau_suvrs %>% filter(roi=="parietal", visit==1)),
  lm(suvr ~ age_at_ftp_bl_mc + apoe4_alleles_mc + cdr_sb_bl_mc + fbb_cl_bl_mc + sex,
     data=tau_suvrs %>% filter(roi=="temporal", visit==1))
)

rois <- c("frontal", "parietal", "temporal")
roi_mods_gam <- list(
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(subj, bs="re"),
      data=df %>% filter(roi=="frontal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(subj, bs="re"),
      data=df %>% filter(roi=="parietal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(subj, bs="re"),
      data=df %>% filter(roi=="temporal", visit>1),
      method="REML")
)
roi_mods_gam_full <- list(
  gam(suvr_annchg_from_last ~ age_at_ftp_last_mc + apoe4_alleles + sex + s(suvr_last, bs="cs") + s(subj, bs="re"),
      data=df %>% filter(roi=="frontal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ age_at_ftp_last_mc + apoe4_alleles + sex + s(suvr_last, bs="cs") + s(subj, bs="re"),
      data=df %>% filter(roi=="parietal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ age_at_ftp_last_mc + apoe4_alleles + sex + s(suvr_last, bs="cs") + s(subj, bs="re"),
      data=df %>% filter(roi=="temporal", visit>1),
      method="REML")
)
xmod <- roi_mods_gam[[3]]
plot(xmod, select=1, residuals=T, shift=fixef(xmod)["(Intercept)"], xlim=c(1.22, 3.5), ylim=c(-0.2, 0.2), seWithMean=T, rug=F)

roi_mods_gam_full <- list(
  gam(suvr_annchg_from_last ~ apoe4_alleles + sex + s(suvr_last, bs="cs") + s(age_at_ftp_last, bs="cs") + s(subj, bs="re"),
      data=df %>% filter(roi=="frontal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ apoe4_alleles + sex + s(suvr_last, bs="cs") + s(age_at_ftp_last, bs="cs") + s(subj, bs="re"),
      data=df %>% filter(roi=="parietal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ apoe4_alleles + sex + s(suvr_last, bs="cs") + s(age_at_ftp_last, bs="cs") + s(subj, bs="re"),
      data=df %>% filter(roi=="temporal", visit>1),
      method="REML")
)
roi_mods_gam_full <- list(
  gam(suvr_annchg_from_last ~ apoe4_allelesf + sex + s(suvr_last, bs="cs") + s(age_at_ftp_last, bs="cs") + s(subj, bs="re"),
      data=df %>% filter(roi=="frontal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ apoe4_allelesf + sex + s(suvr_last, bs="cs") + s(age_at_ftp_last, bs="cs") + s(subj, bs="re"),
      data=df %>% filter(roi=="parietal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ apoe4_allelesf + sex + s(suvr_last, bs="cs") + s(age_at_ftp_last, bs="cs") + s(subj, bs="re"),
      data=df %>% filter(roi=="temporal", visit>1),
      method="REML")
)
roi_mods_gam_full2 <- list(
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(age_at_ftp_last, bs="cs") + s(apoe4_alleles, bs="cs", k=2) + s(sex, bs="re") + s(subj, bs="re"),
      data=df %>% filter(roi=="frontal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(age_at_ftp_last, bs="cs") + s(apoe4_alleles, bs="cs", k=2) + s(sex, bs="re") + s(subj, bs="re"),
      data=df %>% filter(roi=="parietal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(age_at_ftp_last, bs="cs") + s(apoe4_alleles, bs="cs", k=2) + s(sex, bs="re") + s(subj, bs="re"),
      data=df %>% filter(roi=="temporal", visit>1),
      method="REML")
)
roi_mods_gam_full <- list(
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(age_at_ftp_last, bs="cs") + s(apoe4_allelesf, bs="re") + s(sex, bs="re") + s(subj, bs="re"),
      data=df %>% filter(roi=="frontal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(age_at_ftp_last, bs="cs") + s(apoe4_allelesf, bs="re") + s(sex, bs="re") + s(subj, bs="re"),
      data=df %>% filter(roi=="parietal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(age_at_ftp_last, bs="cs") + s(apoe4_allelesf, bs="re") + s(sex, bs="re") + s(subj, bs="re"),
      data=df %>% filter(roi=="temporal", visit>1),
      method="REML")
)

roi_ <- "pericalcarine"
xmod <- gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(subj, bs="re"),
            data=tau_suvrs %>% filter(roi==roi_, visit>1),
            method="REML")
summary(xmod)
plot(xmod, select=1, shift=fixef(xmod)["(Intercept)"], seWithMean=T,
     residuals=T, shade=T, xlim=c(1, 3), ylim=c(-0.1, 0.2),
     xlab=str_c(str_to_upper(roi_), " SUVR last"),
     ylab=str_c(str_to_upper(roi_), " ΔSUVR/year"))

roi_mods_gam_full <- list(
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(subj, bs="re") + age_at_ftp_bl_mc + apoe4_alleles_mc + cdr_sb_bl_mc + fbb_cl_bl_mc + sex,
      data=tau_suvrs %>% filter(roi=="frontal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(subj, bs="re") + age_at_ftp_bl_mc + apoe4_alleles_mc + cdr_sb_bl_mc + fbb_cl_bl_mc + sex,
      data=tau_suvrs %>% filter(roi=="parietal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(subj, bs="re") + age_at_ftp_bl_mc + apoe4_alleles_mc + cdr_sb_bl_mc + fbb_cl_bl_mc + sex,
      data=tau_suvrs %>% filter(roi=="temporal", visit>1),
      method="REML")
)

map(roi_mods_gam_full, ~ {
  params <- model_parameters(.x, ci_method = "wald")
  params <- standardize_names(params)
  formatted_table <- format_table(params, digits = 4, ci_digits=4, stars = T)
  as.data.frame(formatted_table)
}) %>% bind_rows() %>% select(-Component) %>% View()

ii <- 3
plot.gam(roi_mods_gam[[ii]], shade=T, residuals=T, shift=fixef(roi_mods_gam[[ii]])["(Intercept)"], select=1)
plot.gam(roi_mods_gam_full[[ii]], shade=T, residuals=T, shift=fixef(roi_mods_gam_full[[ii]])["(Intercept)"], select=1)

## Predict!!
xstep <- 0.001 # predict in increments of this SUVR
cutoffs <- tibble()
pred <- tibble()

# CDR-SB
# Predict ΔCDR-SB ~ last CDR-SB
df <- cdr %>% filter(visit>1) %>% tibble()
dv <- "cdr_sb_annchg_from_last"
iv <- "cdr_sb_last"
smooth_iv <- "s(cdr_sb_last)"
xmod <- cdr_mod

# Calculate cutoff values
# cutoffs_ <- df %>%
#   summarize(
#     iv = iv,
#     n_obs = n(),
#     n_subj = n_distinct(subj),
#     xmin = round(max(min(!!sym(iv)), quantile(!!sym(iv), 0.25) - (1.5 * IQR(!!sym(iv)))), 1),
#     xmax = round(min(max(!!sym(iv)), quantile(!!sym(iv), 0.75) + (1.5 * IQR(!!sym(iv)))), 1),
#     ymin = round(max(min(!!sym(dv)), quantile(!!sym(dv), 0.25) - (1.5 * IQR(!!sym(dv)))), 2),
#     ymax = round(min(max(!!sym(dv)), quantile(!!sym(dv), 0.75) + (1.5 * IQR(!!sym(dv)))), 2)
#   )
cutoffs_ <- df %>%
  summarize(
    iv = iv,
    n_obs = n(),
    n_subj = n_distinct(subj),
    xmin = round(min(!!sym(iv)), 3),
    xmax = round(max(!!sym(iv)), 3),
    ymin = round(min(!!sym(dv)), 3),
    ymax = round(max(!!sym(dv)), 3)
  )
xvals <- seq(cutoffs_$xmin, cutoffs_$xmax, by = xstep)
cutoffs <- rbind(cutoffs, cutoffs_)

# Get model estimates
pred_ <- with(df, expand.grid(xvals))
names(pred_) <- iv
pred_ <- cbind(
  pred_,
  data.frame(
    predict(
      xmod,
      newdata = pred_,
      type = "response",
      se.fit = TRUE,
      terms = c("(Intercept)", smooth_iv),
      newdata.guaranteed = TRUE
    )
  )
)
pred_ <- pred_ %>%
  rename(iv_value=iv) %>%
  mutate(iv = iv, .before = "iv_value") %>%
  mutate(upper = fit + (2 * se.fit), lower = fit - (2 * se.fit)) %>%
  tibble()
pred <- rbind(pred, pred_)

# ROIs
# Predict ΔSUVR ~ last SUVR, separately in each region
df <- tau_suvrs %>% filter(visit>1, roi %in% rois) %>% tibble()
df$roi <- factor(as.character(df$roi), levels=rois)
iroi <- 1
for (roi_ in rois) {
  # Create variable names based on the current ROI
  dv <- "suvr_annchg_from_last"
  iv <- "suvr_last"
  iv_roi <- str_c(iv, roi_, sep="_")
  smooth_iv <- "s(suvr_last)"
  xmod <- roi_mods_gam[[iroi]]
  iroi = iroi + 1

  # Calculate cutoffs for the current ROI
  # cutoffs_ <- df %>%
  #   filter(roi == roi_) %>%
  #   summarise(
  #     iv = iv_roi,
  #     n_obs = n(),
  #     n_subj = n_distinct(subj),
  #     xmin = round(max(min(!!sym(iv)), quantile(!!sym(iv), 0.25) - (1.5 * IQR(!!sym(iv)))), 1),
  #     xmax = round(min(max(!!sym(iv)), quantile(!!sym(iv), 0.75) + (1.5 * IQR(!!sym(iv)))), 1),
  #     ymin = round(max(min(!!sym(dv)), quantile(!!sym(dv), 0.25) - (1.5 * IQR(!!sym(dv)))), 2),
  #     ymax = round(min(max(!!sym(dv)), quantile(!!sym(dv), 0.75) + (1.5 * IQR(!!sym(dv)))), 2)
  #   )
  cutoffs_ <- df %>%
    filter(roi == roi_) %>%
    summarise(
      iv = iv_roi,
      n_obs = n(),
      n_subj = n_distinct(subj),
      xmin = round(min(!!sym(iv)), 3),
      xmax = round(max(!!sym(iv)), 3),
      ymin = round(min(!!sym(dv)), 3),
      ymax = round(max(!!sym(dv)), 3)
    )
  xvals <- seq(cutoffs_$xmin, cutoffs_$xmax, by = xstep)
  cutoffs <- rbind(cutoffs, cutoffs_)

  # Get model estimates
  pred_ <- with(df, expand.grid(xvals))
  names(pred_) <- iv
  pred_ <- cbind(
    pred_,
    data.frame(
      predict(
        xmod,
        newdata = pred_,
        type = "response",
        se.fit = TRUE,
        terms = c("(Intercept)", smooth_iv),
        newdata.guaranteed = TRUE
      )
    )
  )
  pred_ <- pred_ %>%
    rename(iv_value=iv) %>%
    mutate(iv = iv_roi, .before = "iv_value") %>%
    mutate(upper = fit + (2 * se.fit), lower = fit - (2 * se.fit)) %>%
    tibble()
  pred <- rbind(pred, pred_)
}

# Save the output.
pred_outputf <- file.path(data_dir, "tau-metarois-gam_preds_all-eoad-gt1ftp_2024-01-05.csv")
write_csv(pred, pred_outputf)
print(str_c("Saved", pred_outputf))

## Freesurfer ROIs
df <- tau_suvrs %>% filter(parc=="fsroi_bilat", severity_bl %in% c("MCI", "Dementia"))
df$subj <- factor(as.character(df$subj), levels=unique(as.character(df$subj)))
df$roi <- factor(as.character(df$roi), levels=unique(as.character(df$roi)))
df$severity_bl <- factor(as.character(df$severity_bl), levels=c("MCI", "Dementia"))

rois <- unique(df$roi)
fit_model_for_roi <- function(roi) {
  df_ <- df[df$roi == roi, ]
  mod <- lmer(
    "suvr ~ severity_bl * ftp_yrs_from_bl + (1 | subj) + (0 + ftp_yrs_from_bl | subj)",
    data = df_
  )
  return(mod)
}
fsroi_mods <- lapply(rois, fit_model_for_roi)
fsroi_mods <- tibble(roi = rois, model = fsroi_mods)

roi_ <- "caudalmiddlefrontal"
df_ <- df %>% filter(roi==roi_)
xmod <- (fsroi_mods %>% filter(roi==roi_))$model[[1]]
parameters(fsroi_mods$model[[1]], effects="fixed")
random <- estimate_grouplevel(xmod)
reshaped <- reshape_grouplevel(random, indices = "Coefficient")

output <- tibble()
roi_ <- "superiorfrontal"
df_ <- df[df$roi==roi_,]
xmod <- lmer(
  "suvr ~ severity_bl * ftp_yrs_from_bl + (1 + ftp_yrs_from_bl|subj)",
  data=df_
)
emm <- emmeans(xmod, pairwise ~ severity_bl, infer = c(T, T), adjust = "none")
pairs(emmeans(xmod, ~ severity_bl * ftp_yrs_from_bl, at = list(ftp_yrs_from_bl = c(0))))
emt <- emtrends(xmod, pairwise ~ severity_bl, var = "ftp_yrs_from_bl", infer = c(T, T), adjust = "none")
summary(emt, infer=c(T, T))
emmip(xmod, ftp_yrs_from_bl ~ severity_bl, cov.reduce = range)
summary(xmod)
report(xmod)
get_emtrends(xmod, at = c("severity_bl", "ftp_yrs_from_bl = c(0, 1)"), infer=T)
ggplot(df_, aes(x=ftp_yrs_from_bl, y=suvr)) +
  geom_point(aes(color=severity_bl)) +
  geom_smooth(aes(color=severity_bl), method="lm", se=F) +
  geom_smooth(color="black", method="lm", se=F) +
  # guides(color=F) +
  theme_classic()

# Baseline: MCI - Dementia
get_emcontrasts(xmod, contrast="severity_bl", at="ftp_yrs_from_bl=0", adjust="none", infer=T)

# Rate of change: MCI - Dementia


# ----------------
# Initialize an empty dataframe to store the results..
mm_fits <- data.frame(
  subj = character(),
  roi = character(),
  fixed_icpt = numeric(),
  fixed_time = numeric(),
  rdm_icpt = numeric(),
  rdm_time = numeric(),
  stringsAsFactors = FALSE
)

# For each roi in df$roi, fit a linear mixed-effects model like this:
# suvr ~ ftp_yrs_from_bl + (1|subj) + (0 + ftp_yrs_from_bl|subj)
# then use the fitted model to append a row to the mm_fits dataframe
# for each subject and roi.
df <- tau_suvrs %>% filter(parc %in% c("fsroi_bilat", "metarois"), severity_bl %in% c("MCI", "Dementia"))
df$subj <- factor(as.character(df$subj), levels=unique(as.character(df$subj)))
df$parc <- factor(as.character(df$parc), levels=unique(as.character(df$parc)))
df$roi <- factor(as.character(df$roi), levels=unique(as.character(df$roi)))
df$severity_bl <- factor(as.character(df$severity_bl), levels=c("MCI", "Dementia"))
rois <- sort(unique(df$roi))
for (roi in rois) {
  # Select data for the current ROI
  df_roi <- df[df$roi == roi, ]

  # Fit the linear mixed-effects model
  # mm <- lmer(
  #   suvr ~ ftp_yrs_from_bl + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
  #   data = df_roi
  # )
  mm <- lmer(
    suvr_roi_mc ~ ftp_yrs_from_bl + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
    data = df_roi
  )

  # Extract fixed and random effects
  mm_fixefs <- fixef(mm)
  mm_ranefs <- ranef(mm)$subj

  # Create a dataframe for the current roi's results
  mm_fits_roi <- data.frame(
    subj = rownames(mm_ranefs),
    roi = roi,
    fixed_icpt = mm_fixefs[["(Intercept)"]],
    fixed_time = mm_fixefs[["ftp_yrs_from_bl"]],
    rdm_icpt = mm_ranefs[, "(Intercept)"],
    rdm_time = mm_ranefs[, "ftp_yrs_from_bl"]
  )

  # Append the results to the 'mm_fits' dataframe
  mm_fits <- bind_rows(mm_fits, mm_fits_roi)
}

mm_fits <- mm_fits %>%
  mutate(
    icpt = fixed_icpt + rdm_icpt,
    time = fixed_time + rdm_time
  )
mm_fits <- mm_fits %>%
  inner_join(
    df %>%
      distinct(subj, roi, .keep_all=TRUE),
    by = c("subj", "roi")
  ) %>%
  select(
    subj, parc, roi, dx, ftp_visits,
    age_at_ftp_bl, age_at_ftp_bl_mc, sex,
    apoe4_alleles, apoe4_alleles_mc,
    fbb_cl_bl, fbb_cl_bl_mc,
    mmse_bl, mmse_bl_mc,
    cdr_sb_bl, cdr_sb_bl_mc,
    fixed_icpt, fixed_time,
    rdm_icpt, rdm_time,
    icpt, time
  )

# Preview summary stats
mm_fits %>%
  group_by(roi) %>%
  summarize(
    n_subj = n_distinct(subj),
    min_icpt = min(icpt),
    max_icpt = max(icpt),
    mean_icpt = mean(icpt),
    sd_icpt = sd(icpt),
    min_time = min(time),
    max_time = max(time),
    mean_time = mean(time),
    sd_time = sd(time)
  ) %>%
  arrange(-mean_time) %>%
  view()

# Save the dataframe with subject- and ROI-specific SUVR intercepts and slopes
mm_outputf <- file.path(data_dir, "tau-rois-mms_all-eoad-gt1ftp_2023-10-19.csv")
write_csv(
  mm_fits,
  mm_outputf,
)
