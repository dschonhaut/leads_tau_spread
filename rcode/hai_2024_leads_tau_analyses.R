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

# Load dataframes.
tau_suvrsf <- file.path(data_dir, "tau-rois-agg_eoad-long_formatted_2023-12-28.csv")
tau_suvrs <- read_csv(tau_suvrsf)
tau_suvrs_widef <- file.path(data_dir, "tau-rois-agg_eoad-wide_formatted_2023-12-28.csv")
tau_suvrs_wide <- read_csv(tau_suvrs_widef)
tau_suvrs_wide$parc <- "metarois"
# tau_mmsf <- file.path(data_dir, "tau-rois-mms_all-eoad-gt1ftp_2023-10-21.csv")
# tau_mms <- read_csv(tau_mmsf)

# Clean dataframes.
transform_df <- function(df) {
  df %>%
    # filter(parc %in% c("metarois")) %>%
    mutate(across(where(is.character), as.factor)) %>%
    mutate(apoe4_allelesf = factor(apoe4_alleles, ordered = FALSE)) %>%
    mutate(apoe4_alleles = factor(apoe4_alleles, ordered = TRUE)) %>%
    { contrasts(.$sex) <- 'contr.sum'; . }
}
tau_suvrs <- transform_df(tau_suvrs)
tau_suvrs_wide <- transform_df(tau_suvrs_wide)
# tau_mms <- transform_df(tau_mms)

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
tau_mms <- process_cats(tau_mms)
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

df <- as.tibble(tau_suvrs_wide)
mod_name <- paste("mod_long_full_reml", roi, sep="_")
gmod <- gam(
  suvr_annchg_from_last_temporal ~
    # Main effects
    s(suvr_bl_temporal),
    # s(age_at_ftp_bl_mc) +
    # s(fbb_cl_bl_mc) +
    # s(cdr_sb_bl_mc) +
    # apoe4_alleles +
    # sex,
  data=df,
  method="REML"
)
bs = "cs"
gmod_t_by_t <- gam(suvr_annchg_from_last_temporal ~ s(suvr_bl_temporal, bs = bs), data=df, method="REML")
gmod_p_by_p <- gam(suvr_annchg_from_last_parietal ~ s(suvr_bl_parietal, bs = bs), data=df, method="REML")
gmod_f_by_f <- gam(suvr_annchg_from_last_frontal ~ s(suvr_bl_frontal, bs = bs), data=df, method="REML")

gmod_t_by_p <- gam(suvr_annchg_from_last_temporal ~ s(suvr_bl_parietal, bs = bs), data=df, method="REML")
gmod_t_by_f <- gam(suvr_annchg_from_last_temporal ~ s(suvr_bl_frontal, bs = bs), data=df, method="REML")
gmod_p_by_t <- gam(suvr_annchg_from_last_parietal ~ s(suvr_bl_temporal, bs = bs), data=df, method="REML")
gmod_p_by_f <- gam(suvr_annchg_from_last_parietal ~ s(suvr_bl_frontal, bs = bs), data=df, method="REML")
gmod_f_by_t <- gam(suvr_annchg_from_last_frontal ~ s(suvr_bl_temporal, bs = bs), data=df, method="REML")
gmod_f_by_p <- gam(suvr_annchg_from_last_frontal ~ s(suvr_bl_parietal, bs = bs), data=df, method="REML")

gmod_t_by_tp <- gam(suvr_annchg_from_last_temporal ~ s(suvr_bl_temporal, bs = bs) + s(suvr_bl_parietal, bs = bs), data=df, method="REML")
gmod_t_by_tf <- gam(suvr_annchg_from_last_temporal ~ s(suvr_bl_temporal, bs = bs) + s(suvr_bl_frontal, bs = bs), data=df, method="REML")
gmod_p_by_pt <- gam(suvr_annchg_from_last_parietal ~ s(suvr_bl_parietal, bs = bs) + s(suvr_bl_temporal, bs = bs), data=df, method="REML")
gmod_p_by_pf <- gam(suvr_annchg_from_last_parietal ~ s(suvr_bl_parietal, bs = bs) + s(suvr_bl_frontal, bs = bs), data=df, method="REML")
gmod_f_by_ft <- gam(suvr_annchg_from_last_frontal ~ s(suvr_bl_frontal, bs = bs) + s(suvr_bl_temporal, bs = bs), data=df, method="REML")
gmod_f_by_fp <- gam(suvr_annchg_from_last_frontal ~ s(suvr_bl_frontal, bs = bs) + s(suvr_bl_parietal, bs = bs), data=df, method="REML")

AIC(
  gmod_t_by_t, gmod_t_by_p, gmod_t_by_f, gmod_t_by_tp, gmod_t_by_tf,
  gmod_p_by_p, gmod_p_by_t, gmod_p_by_f, gmod_p_by_pt, gmod_p_by_pf,
  gmod_f_by_f, gmod_f_by_t, gmod_f_by_p, gmod_f_by_ft, gmod_f_by_fp
)

# Don't run these model comparison tests
# if using select=T or a shrinkage smooth
anova.gam(gmod_t_by_t, gmod_t_by_tp, test = "F") #*
anova.gam(gmod_t_by_t, gmod_t_by_tf, test = "F")
anova.gam(gmod_p_by_p, gmod_p_by_pt, test = "F") #*
anova.gam(gmod_p_by_p, gmod_p_by_pf, test = "F")
anova.gam(gmod_f_by_f, gmod_f_by_ft, test = "F")
anova.gam(gmod_f_by_f, gmod_f_by_fp, test = "F") #*

summary(gmod_t_by_t)
summary(gmod_t_by_tp)
summary(gmod_t_by_tf)
summary(gmod_p_by_p)
summary(gmod_p_by_pt)
summary(gmod_p_by_pf)
summary(gmod_f_by_f)
summary(gmod_f_by_ft)
summary(gmod_f_by_fp)

p1 <- (
  draw(
  gmod_t_by_t,
  constant=fixef(gmod_t_by_t),
  rug=F,
  smooth_col="#2E45B8",
  ci_col="#2E45B8",
  nrow=1,
  ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p2 <- (
  draw(
    gmod_t_by_p,
    constant=fixef(gmod_t_by_p),
    rug=F,
    smooth_col="#3EBCD2",
    ci_col="#3EBCD2",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p3 <- (
  draw(
    gmod_t_by_f,
    constant=fixef(gmod_t_by_f),
    rug=F,
    smooth_col="#FF4983",
    ci_col="#FF4983",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p1 + p2 + p3 + plot_layout(ncol = 3)

p4 <- (
  draw(
    gmod_p_by_t,
    constant=fixef(gmod_p_by_t),
    rug=F,
    smooth_col="#2E45B8",
    ci_col="#2E45B8",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p5 <- (
  draw(
    gmod_p_by_p,
    constant=fixef(gmod_p_by_p),
    rug=F,
    smooth_col="#3EBCD2",
    ci_col="#3EBCD2",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p6 <- (
  draw(
    gmod_p_by_f,
    constant=fixef(gmod_p_by_f),
    rug=F,
    smooth_col="#FF4983",
    ci_col="#FF4983",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p4 + p5 + p6 + plot_layout(ncol = 3)

p7 <- (
  draw(
    gmod_f_by_t,
    constant=fixef(gmod_f_by_t),
    rug=F,
    smooth_col="#2E45B8",
    ci_col="#2E45B8",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p8 <- (
  draw(
    gmod_f_by_p,
    constant=fixef(gmod_f_by_p),
    rug=F,
    smooth_col="#3EBCD2",
    ci_col="#3EBCD2",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p9 <- (
  draw(
    gmod_f_by_f,
    constant=fixef(gmod_f_by_f),
    rug=F,
    smooth_col="#FF4983",
    ci_col="#FF4983",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p7 + p8 + p9 + plot_layout(ncol = 3)

p1 + p5 + p9 + plot_layout(ncol = 3)

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
esquisser(tau_mms)
esquisser(tau_suvrs2)

tau_mms %>%
 filter(parc %in% "metarois") %>%
 ggplot() +
 aes(x = icpt, y = time, colour = ftp_visits) +
 geom_point(shape = "circle", size = 2L) +
 scale_color_viridis_c(option = "viridis", direction = 1) +
 theme_linedraw() +
 facet_wrap(vars(roi))

View(tau_suvrs)

df <- tau_suvrs
contrasts(tau_suvrs$sex) <- "contr.sum"
roi <- "frontal"
df_roi <- df[df$roi == roi, ]
model <- lmer(suvr ~ ftp_yrs_from_bl + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
              data = df_roi)


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

# Model longitudinal CDR-SB
cdrf <- file.path(data_dir, "cdr-sb_eoad-long_formatted_193subjs_2023-12-30.csv")
cdr <- read_csv(cdrf)
cols <- c("cdr_sb", "cdr_sb_bl", "cdr_yrs_from_bl")
cdr$subj <- as.factor(cdr$subj)
cdr$cdr_sb <- as.numeric(cdr$cdr_sb)
cdr$cdr_sb_last <- as.numeric(cdr$cdr_sb_last)
cdr$cdr_sb_bl <- as.numeric(cdr$cdr_sb_bl)
cdr$cdr_yrs_from_bl <- as.numeric(cdr$cdr_yrs_from_bl)
cdr$cdr_sb_annchg_from_last <- as.numeric(cdr$cdr_sb_annchg_from_last)

cdr_mod <- gam(
  cdr_sb_annchg_from_last ~ s(cdr_sb_last, bs="cs") + s(subj, bs="re"),
  data=cdr %>% filter(visit>1),
  method="REML"
)


# Model ROIs like CDR...
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
      data=tau_suvrs %>% filter(roi=="frontal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(subj, bs="re"),
      data=tau_suvrs %>% filter(roi=="parietal", visit>1),
      method="REML"),
  gam(suvr_annchg_from_last ~ s(suvr_last, bs="cs") + s(subj, bs="re"),
      data=tau_suvrs %>% filter(roi=="temporal", visit>1),
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

### Predict!!
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
cutoffs_ <- df %>%
  summarize(
    iv = iv,
    n_obs = n(),
    n_subj = n_distinct(subj),
    xmin = round(max(min(!!sym(iv)), quantile(!!sym(iv), 0.25) - (1.5 * IQR(!!sym(iv)))), 1),
    xmax = round(min(max(!!sym(iv)), quantile(!!sym(iv), 0.75) + (1.5 * IQR(!!sym(iv)))), 1),
    ymin = round(max(min(!!sym(dv)), quantile(!!sym(dv), 0.25) - (1.5 * IQR(!!sym(dv)))), 2),
    ymax = round(min(max(!!sym(dv)), quantile(!!sym(dv), 0.75) + (1.5 * IQR(!!sym(dv)))), 2)
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
  cutoffs_ <- df %>%
    filter(roi == roi_) %>%
    summarise(
      iv = iv_roi,
      n_obs = n(),
      n_subj = n_distinct(subj),
      xmin = round(max(min(!!sym(iv)), quantile(!!sym(iv), 0.25) - (1.5 * IQR(!!sym(iv)))), 1),
      xmax = round(min(max(!!sym(iv)), quantile(!!sym(iv), 0.75) + (1.5 * IQR(!!sym(iv)))), 1),
      ymin = round(max(min(!!sym(dv)), quantile(!!sym(dv), 0.25) - (1.5 * IQR(!!sym(dv)))), 2),
      ymax = round(min(max(!!sym(dv)), quantile(!!sym(dv), 0.75) + (1.5 * IQR(!!sym(dv)))), 2)
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
pred_outputf <- file.path(data_dir, "tau-metarois-gam_preds_all-eoad-gt1ftp_2024-01-01.csv")
write_csv(pred, pred_outputf)
print(str_c("Saved", pred_outputf))

## Freesurfer ROIs
df <- tau_suvrs %>% filter(parc=="fsroi_bilat", severity_bl %in% c("MCI", "Dementia"))
df$subj <- factor(as.character(df$subj), levels=unique(as.character(df$subj)))
df$roi <- factor(as.character(df$roi), levels=unique(as.character(df$roi)))
df$severity_bl <- factor(as.character(df$severity_bl), levels=c("MCI", "Dementia"), ordered=T)

xmod <- lmer()
rois <- unique(df$roi)
fsroi_mods <- expand.grid(roi_=rois) %>%
  mutate(
    formula = "suvr ~ severity_bl + ftp_yrs_from_bl + (1 + ftp_yrs_from_bl|subj)"
    model = pmap(list(formula, list(df[df$roi==roi_,])), ~ lmer(as.formula(..1), data = ..2))
  )
