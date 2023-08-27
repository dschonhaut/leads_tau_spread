library(tidyverse)
library(easystats)
library(tictoc)
library(mgcv)
library(gratia)
library(gtsummary)
library(emmeans)
# library(sjPlot)
# library(sjstats)
library(ggeffects)

# Setup

## Define plotting variables
blue1 <- "#D6DBF5"
blue2 <- "#2E45B8"
blue3 <- "#141F52"
cyan1 <- "#6FE4FB"
cyan2 <- "#3EBCD2"
cyan3 <- "#005F73"
green1 <- "#D2F9F0"
green2 <- "#1DC9A4"
green3 <- "#0E6452"
yellow1 <- "#FCDE83"
yellow2 <- "#F9C31F"
yellow3 <- "#CFA300"
beige1 <- "#EDE3D5"
beige2 <- "#D4C6B6"
beige3 <- "#BDAF97"
orange1 <- "#FCB583"
orange2 <- "#F97A1F"
orange3 <- "#C75400"
red1 <- "#FFA39F"
red2 <- "#E3120B"
red3 <- "#861320"
pink1 <- "#FFA4C2"
pink2 <- "#FF4983"
pink3 <- "#A80055"
violet1 <- "#E6D4FA"
violet2 <- "#B38FE7"
violet3 <- "#7D55C7"
white <- "#FFFFFF"
gray1 <- "#F2F2F2"
gray2 <- "#D9D9D9"
gray3 <- "#B3B3B3"
gray4 <- "#595959"
gray5 <- "#333333"
gray6 <- "#1A1A1A"
black <- "#000000"
pal24 <- c(
  blue2, orange2, cyan2, pink2, green2, violet2, yellow2, red2,
  blue3, orange3, cyan3, pink3, green3, violet3, yellow3, red3,
  blue1, orange1, cyan1, pink1, green1, violet1, yellow1, red1
)

## Define parameters
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
data_file <- file.path(data_dir, "tau-eoad-re_2023-05-16.csv")
save_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/analysis/r_models/gam/bward_sel/dv_chg_suvr_re"
today <- format(Sys.Date(), "%Y-%m-%d")
overwrite = FALSE
verbose = TRUE

## Load and format the data
check_cols <- c("subj", "roi", "suvr_bl_re", "chg_yr_re",
                "age_at_ftp_bl", "sex", "fbb_cl_bl",
                "apoe4_alleles", "cdr_sb_bl")
df <- read.csv(data_file)
df <- as.data.frame(unclass(df), stringsAsFactors = TRUE)
df <- df[complete.cases(df[,check_cols]), check_cols]
df <- data.frame(lapply(df, function(x) if (is.factor(x)) droplevels(x) else x))
df$age_at_ftp_bl_mc <- df$age_at_ftp_bl - mean(df$age_at_ftp_bl)
df$apoe4_alleles <- factor(df$apoe4_alleles, ordered = TRUE)
df$apoe4_alleles_int <- as.numeric(as.character(df$apoe4_alleles))
contrasts(df$roi) <- 'contr.sum'
contrasts(df$sex) <- 'contr.sum'

## Load the fitted model
mod_long_full <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r9_lin_age_reml_2023-06-24.rds"))
)
sum_mod_long_full <- readRDS(
  file.path(save_dir, paste0("sum_mod_long_full_r9_lin_age_reml_2023-06-24.rds"))
)

## Define groups of regions to plot
roi_map <- list(
  "all" = list(
    "title" = "All regions",
    "rois" = c(
      "all", "acc", "amy", "cun", "erc", "fsf", "ifg", "ins", "ipc",
      "itg", "locc", "ling", "mfg", "mtg", "ofc", "pcc", "phg",
      "prcu", "v1", "ssm", "sfg", "spc", "stg", "smg", "tp"
    ),
    "labels" = c(
      "Global", "ACC", "Amygdala", "Cuneus", "Entorhinal", "Fusiform",
      "IFG", "Insula", "IPC", "ITG", "Lat. Occipital", "Lingual",
      "MFG", "MTG", "OFC", "PCC", "PHG", "Precuneus", "Prim. Visual",
      "Sensorimotor", "SFG", "SPC", "STG", "Supramarginal", "Temp. Pole"
    )
  ),
  "global" = list(
    "title" = "Global",
    "rois" = c("all"),
    "labels" = c("Global")
  ),
  "mtl" = list(
    "title" = "MTL",
    "rois" = c("all", "amy", "erc", "phg"),
    "labels" = c("Global", "Amygdala", "Entorhinal", "PHG")
  ),
  "temporal" = list(
    "title" = "Temporal",
    "rois" = c("fsf", "itg", "mtg", "stg", "tp"),
    "labels" = c("Fusiform", "ITG", "MTG", "STG", "Temp. Pole")
  ),
  "parietal" = list(
    "title" = "Parietal",
    "rois" = c("pcc", "prcu", "smg", "ipc", "spc", "ssm"),
    "labels" = c("PCC", "Precuneus", "Supramarginal", "IPC", "SPC", "Sensorimotor")
  ),
  "occipital" = list(
    "title" = "Occipital",
    "rois" = c("ling", "locc", "cun", "v1"),
    "labels" = c("Lingual", "Lat. Occipital", "Cuneus", "Prim. Visual")
  ),
  "frontal" = list(
    "title" = "Frontal",
    "rois" = c("sfg", "mfg", "ifg", "acc", "ofc", "ins"),
    "labels" = c("SFG", "MFG", "IFG", "ACC", "OFC", "Insula")
  ),
  "meta_temporal" = list(
    "title" = "Meta-temporal ROI",
    "rois" = c("amy", "erc", "phg", "fsf", "itg", "mtg"),
    "labels" = c("Amygdala", "Entorhinal", "PHG", "Fusiform", "ITG", "MTG")
  ),
  "slow" = list(
    "title" = "Region",
    "rois" = c(
      "all", "acc", "amy", "erc", "ins", "pcc", "phg", "stg",
      "tp"
    ),
    "labels" = c(
      "Global", "ACC", "Amygdala", "Entorhinal", "Insula", "PCC", "PHG", "STG",
      "Temp. Pole"
    )
  ),
  "average" = list(
    "title" = "Region",
    "rois" = c(
      "all", "fsf", "ling", "ofc", "v1", "ssm",
      "spc", "smg"
    ),
    "labels" = c(
      "Global", "Fusiform", "Lingual", "OFC", "Prim. Visual", "Sensorimotor",
      "SPC", "Supramarginal"
    )
  ),
  "fast" = list(
    "title" = "Region",
    "rois" = c(
      "all", "cun", "ifg", "ipc", "itg", "locc", "mfg", "mtg",
      "prcu", "sfg"
    ),
    "labels" = c(
      "Global", "Cuneus", "IFG", "IPC", "ITG", "Lat. Occipital", "MFG", "MTG",
      "Precuneus", "SFG"
    )
  ),
  "late" = list(
    "title" = "Region",
    "rois" = c(
      "all", "ofc", "v1", "ssm"
    ),
    "labels" = c(
      "Global", "OFC", "Prim. Visual", "Sensorimotor"
    )
  )
)

## Summarize the subject demographics.
(df %>%
    filter(roi == "itg") %>%
    select(age_at_ftp_bl, sex, apoe4_alleles, fbb_cl_bl, cdr_sb_bl) %>%
    tbl_summary(
      digits = all_continuous() ~ 1,
      label = list(
        age_at_ftp_bl ~ "Age",
        sex ~ "Sex",
        apoe4_alleles ~ "APOE-ε4 alleles",
        fbb_cl_bl ~ "Aβ-PET Centiloids",
        cdr_sb_bl ~ "CDR-SB"
      )
    ) %>%
    bold_labels()
)

keep_cols <- c("suvr_bl_re", "age_at_ftp_bl", "fbb_cl_bl", "cdr_sb_bl")
chart.Correlation(df[keep_cols])

## Check the model
par(mfrow = c(2, 2))
gam.check(mod_long_full_r9_lin_age_reml)
AIC(mod_long_full_r9_lin_age_reml)
sum_mod_long_full_r9_lin_age_reml

model_parameters(mod_long_full)
plot(ggpredict(mod_long_full, terms = "age_at_ftp_bl_mc"))
draw(mod_long_full, select = 1)

est_rel <- estimate_relation(mod_long_full, include_random = FALSE)

deriv <- estimate_slopes(mod_long_full,
                         trend = "suvr_bl_re",
                         at = "suvr_bl_re",
                         length = 100)
plot(deriv)
summary(deriv)


## Plot model results
par(mfrow = c(3, 3))
plot.gam(mod_long_full, pages = 3)

mod_long_full %>%
  tbl_regression(
    digits = all_continuous() ~ 4,
    tidy_fun = tidy_gam
    # pvalue_fun = ~style_pvalue(.x, digits = 4)
  ) %>%
  add_q() %>%
  bold_p(q = TRUE)
  tab_options(
    table.font.size = "small"
  )

ref_grid(mod_long_full)

## Look at contrasting means
emm_long_full_age <- emmeans(mod_long_full, ~ age_at_ftp_bl_mc,
                             at = list(age_at_ftp_bl_mc = c(
                               50 - mean(df$age_at_ftp_bl),
                               55 - mean(df$age_at_ftp_bl),
                               60 - mean(df$age_at_ftp_bl),
                               65 - mean(df$age_at_ftp_bl)
                             )),
                             adjust = "none", infer = TRUE)
emt_age <- emtrends(mod_long_full, "roi", var = "age_at_ftp_bl_mc",
                    adjust = "BH", infer = TRUE)
data.frame(emt_age) %>%
  arrange(age_at_ftp_bl_mc.trend) %>%
  select(-SE, -df) %>%
  rename(
    chg_yr_re = age_at_ftp_bl_mc.trend,
    lower = lower.CL,
    upper = upper.CL,
    t = t.ratio,
    p_fdr = p.value
  ) # %>%
  # save_csv(
  #   file.path(data_dir, paste0("emt_age", "_", today, ".csv")),
  #   overwrite = overwrite
  # )

emt_apoe <- emtrends(mod_long_full, "roi", var = "apoe4_alleles_int",
                     adjust = "BH", infer = TRUE)
data.frame(emt_apoe) %>%
  arrange(apoe4_alleles_int.trend) %>%
  select(-SE, -df) %>%
  rename(
    chg_yr_re = apoe4_alleles_int.trend,
    lower = lower.CL,
    upper = upper.CL,
    t = t.ratio,
    p_fdr = p.value
  ) # %>%
  # save_csv(
  #   file.path(data_dir, paste0("emt_apoe", "_", today, ".csv")),
  #   overwrite = overwrite
  # )

emm_long_full_apoe <- emmeans(mod_long_full, ~ apoe4_alleles | roi,
                              adjust = "none", infer = TRUE)
plot(emm_long_full_apoe)


pairs(emm, simple = "roi")
plot(emm)

emm_roi <- emmeans(mod_long_full, ~ roi,
                   adjust = "none", infer = TRUE, weights = "proportional")
# Calculate mean values of suvr_bl_re within each roi level
mean_suvr_roi <- emmeans(mod_long_full, ~ roi, at = list(suvr_bl_re = "range"))

# Create a custom contrast matrix
custom_contrast <- cbind(1, -mean_suvr_roi$emmean$suvr_bl_re)

# Calculate marginal means using the custom contrast
emm_roi <- emmeans(mod_long_full, ~ roi, adjust = "none", infer = TRUE)

emm_apoe <- emmeans(mod_long_full, pairwise ~ apoe4_alleles,
                    adjust = "none")
emm_apoe_roi <- emmeans(mod_long_full, pairwise ~ apoe4_alleles, by = "roi",
                        adjust = "none")
emm_age <- emmeans(mod_long_full, ~ age_at_ftp_bl,
                   at = list(age_at_ftp_bl = c(45, 50, 55, 60, 65)),
                   adjust = "none", infer = TRUE)
emm_age_roi <- emmeans(mod_long_full, ~ age_at_ftp_bl, by = "roi",
                       at = list(age_at_ftp_bl = c(45, 50, 55, 60, 65)),
                       adjust = "none", infer = TRUE)

plot(emm_age) +
  xlab("FTP SUVR Δ/year") +
  ylab("Age at baseline") +
  coord_flip() +
  theme_modern() +
  theme(
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length =  unit(5, "pt"),
    axis.text.x = element_text(margin = margin(t = 5)),
    axis.text.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(vjust = 6),
    axis.title.y = element_text(vjust = -2)
  )

plot(emm_roi)

plot(
  xmod, select = 1, shift=icpt,
  xlim = c(xmin, xmax), ylim = c(-0.06, 0.12),
  residuals = TRUE, rug = FALSE, shade = TRUE,
  xlab = "Baseline FTP SUVR", ylab = "FTP SUVR Δ/year"
)

emm_roi <- emmeans(ref_grid(mod_long_full), "roi", null = mean(df$chg_yr_re), adjust = "BH", infer = TRUE)
emm_roi <- test(emm_roi, null = mean(df$chg_yr_re), adjust = "BH")

aov_mod_long_full <- anova(mod_long_full)
ovv <- data.frame(aov_mod_long_full$s.table)[2:25,]
alpha = 0.05
ovv$p_fdr <- p.adjust(ovv$p.value, method = "fdr")
ovv$sig_fdr <- ovv$p_fdr < alpha
View(ovv)

### Plot baseline SUVR ~ change in SUVR: Global
# Get model residuals
xmod <- mod_long_full
mod_terms <- predict.gam(xmod, type="terms")
icpt <- fixef(xmod)["(Intercept)"]
dat <- data.frame(df)
dat$mod <- "mod_long_full"
dat$resid <- residuals.gam(xmod, type = "response")
dat$partial <- (dat$resid
               + icpt
               + mod_terms[,c("s(suvr_bl_re)")]
               + mod_terms[,c("age_at_ftp_bl_mc")]
               + mod_terms[,c("apoe4_alleles_int")]
               )

# Get model estimates
xmod_names <- paste(names(fixef(xmod)), smooths(xmod))
xmin <- round(
  max(
    min(df$suvr_bl_re),
    quantile(df$suvr_bl_re, 0.25) - (1.5 * IQR(df$suvr_bl_re))
    ),
  1)
xmax <- round(
  min(
    max(df$suvr_bl_re),
    quantile(df$suvr_bl_re, 0.75) + (1.5 * IQR(df$suvr_bl_re))
    ),
  1)
xvals <- seq(xmin, xmax, by = 0.01)
pred <- with(df,
              expand.grid(
                suvr_bl_re = xvals,
                age_at_ftp_bl_mc = mean(df$age_at_ftp_bl_mc),
                apoe4_alleles_int = mean(df$apoe4_alleles_int),
                roi = "acc",
                subj = "LDS0070166")
)
pred <- cbind(
  pred,
  data.frame(
    predict(
      xmod,
      newdata = pred,
      type = "response",
      se.fit = TRUE,
      terms = c("(Intercept)", "s(suvr_bl_re)", "age_at_ftp_bl_mc", "apoe4_alleles_int"),
      newdata.guaranteed = TRUE
    )
  )
)
pred$roi <- "all"
pred <- transform(pred, upper = fit + (2 * se.fit), lower = fit - (2 * se.fit))
pred$deriv_fit <- c(diff(pred$fit), diff(pred$fit)[length(diff(pred$fit))])

roi <- "global"
ymin <- round(
  max(
    min(df$chg_yr_re),
    quantile(df$chg_yr_re, 0.25) - (1.5 * IQR(df$chg_yr_re))
  ),
  2)
ymax <- round(
  min(
    max(df$chg_yr_re),
    quantile(df$chg_yr_re, 0.75) + (1.5 * IQR(df$chg_yr_re))
  ),
  2)

# Make plot
plot.gam(mod_long_full,
         select = 1,
         shade = TRUE,
         rug = FALSE,
         seWithMean = TRUE,
         residuals = TRUE,
         shift = fixef(mod_long_full)["(Intercept)"],
         xlim = c(xmin, xmax),
         ylim = c(ymin, ymax),
         xlab = "Baseline FTP SUVR",
         ylab = "FTP SUVR Δ/year"
         )

ggplot(pred) +
  geom_hline(yintercept = 0, color = black) +
  geom_point(
    data = dat, aes(x = suvr_bl_re, y = partial),
    # data = dat, aes(x = suvr_bl_re, y = chg_yr_re),
    size = 0.05, alpha = 0.33, color = blue2,
  ) +
  geom_line(
    aes(x = suvr_bl_re, y = fit, color = roi),
    size = 1.2
  ) +
  geom_ribbon(
    aes(x = suvr_bl_re, y = fit, ymin = lower, ymax = upper, fill = roi),
    alpha = 0.2, color=NA
  ) +
  scale_x_continuous(
    limits = c(xmin, xmax),
    breaks = seq(1, 3.5, by = 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(ymin, ymax),
    breaks = seq(ymin, ymax, by = 0.04),
    expand = c(0, 0)
  ) +
  scale_color_discrete(
    type = c(black), # append(black, pal24),
    name = "Region",
    limits = roi_map[[roi]][["rois"]],
    labels = roi_map[[roi]][["labels"]]
  ) +
  scale_fill_discrete(
    type = c(black), # append(black, pal24),
    name = "Region",
    limits = roi_map[[roi]][["rois"]],
    labels = roi_map[[roi]][["labels"]]
  ) +
  labs(
    x = "Baseline FTP SUVR",
    y = "FTP SUVR Δ/year",
  ) +
  theme_modern() +
  theme(
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length =  unit(5, "pt"),
    axis.text.x = element_text(margin = margin(t = 5)),
    axis.text.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(vjust = 6),
    axis.title.y = element_text(vjust = -4),
    plot.title = element_text(size = rel(1.75), hjust = 0.5)
  )

### Plot baseline SUVR ~ change in SUVR: Regional
# Get model residuals
xmod <- mod_long_full
mod_terms <- predict.gam(xmod, type="terms")
icpt <- fixef(xmod)["(Intercept)"]
dat <- data.frame(df)
dat$mod <- "mod_long_full"
dat$resid <- residuals.gam(xmod, type = "response")
dat$partial <- (dat$resid
                + icpt
                + mod_terms[,c("s(suvr_bl_re)")]
)

# Get model estimates for each region
xmod_names <- paste(names(fixef(xmod)), smooths(xmod))
cutoffs <- df %>%
  group_by(roi) %>%
  summarize(
    n = n(),
    xmin = round(max(min(suvr_bl_re), quantile(suvr_bl_re, 0.25) - (1.5 * IQR(suvr_bl_re))), 1),
    xmax = round(min(max(suvr_bl_re), quantile(suvr_bl_re, 0.75) + (1.5 * IQR(suvr_bl_re))), 1),
    ymin = round(max(min(chg_yr_re), quantile(chg_yr_re, 0.25) - (1.5 * IQR(chg_yr_re))), 2),
    ymax = round(min(max(chg_yr_re), quantile(chg_yr_re, 0.75) + (1.5 * IQR(chg_yr_re))), 2)
  )
cutoffs <- rbind(
  cutoffs,
  df %>%
    summarize(
      roi = "all",
      n = n(),
      xmin = round(max(min(suvr_bl_re), quantile(suvr_bl_re, 0.25) - (1.5 * IQR(suvr_bl_re))), 1),
      xmax = round(min(max(suvr_bl_re), quantile(suvr_bl_re, 0.75) + (1.5 * IQR(suvr_bl_re))), 1),
      ymin = round(max(min(chg_yr_re), quantile(chg_yr_re, 0.25) - (1.5 * IQR(chg_yr_re))), 2),
      ymax = round(min(max(chg_yr_re), quantile(chg_yr_re, 0.75) + (1.5 * IQR(chg_yr_re))), 2)
    )
)
xmin <- min(cutoffs$xmin)
xmax <- max(cutoffs$xmax)
xvals <- seq(xmin, xmax, by = 0.01)
pred <- with(df,
             expand.grid(
               suvr_bl_re = xvals,
               age_at_ftp_bl_mc = mean(df$age_at_ftp_bl_mc),
               apoe4_alleles_int = mean(df$apoe4_alleles_int),
               roi = unique(df$roi),
               subj = "LDS0070166")
)
pred <- pred %>%
  left_join(
    cutoffs[, c("roi", "xmin", "xmax")],
    by = join_by(roi),
    relationship = "many-to-one"
  ) %>%
  filter(suvr_bl_re >= xmin, suvr_bl_re <= xmax) %>%
  select(-xmin, -xmax)
pred <- cbind(
  pred,
  data.frame(
    predict(
      xmod,
      newdata = pred,
      type = "response",
      se.fit = TRUE,
      terms = c("(Intercept)", colnames(mod_terms)[1:30]),
      newdata.guaranteed = TRUE
    )
  )
)
pred <- transform(pred, upper = fit + (2 * se.fit), lower = fit - (2 * se.fit))

# Add model estimates for the global mean across regions
pred_ <- with(df,
             expand.grid(
               suvr_bl_re = xvals,
               age_at_ftp_bl_mc = mean(df$age_at_ftp_bl_mc),
               apoe4_alleles_int = mean(df$apoe4_alleles_int),
               roi = "acc",
               subj = "LDS0070166")
)
pred_ <- cbind(
  pred_,
  data.frame(
    predict(
      xmod,
      newdata = pred_,
      type = "response",
      se.fit = TRUE,
      terms = c("(Intercept)", "s(suvr_bl_re)", "age_at_ftp_bl_mc", "apoe4_alleles_int"),
      newdata.guaranteed = TRUE
    )
  )
)
pred_$roi <- "all"
pred_ <- transform(pred_, upper = fit + (2 * se.fit), lower = fit - (2 * se.fit))
pred <- rbind(pred, pred_)

# Save the output.
write.csv(pred,
          file.path(data_dir, "mod-long-full_2023-07-23.csv"),
          row.names=F)

# Plot the data
roi <- "late"

select_rois <- roi_map[[roi]][["rois"]]
pred_sub <- pred %>%
  filter(roi %in% select_rois) %>%
  left_join(
    cutoffs[, c("roi", "xmin", "xmax")],
    by = join_by(roi),
    relationship = "many-to-one"
    ) %>%
  filter(suvr_bl_re >= xmin, suvr_bl_re <= xmax) %>%
  select(-xmin, -xmax)
ymin <- -0.07 # plyr::round_any(min(pred_sub$lower), 0.01, f = floor) # -0.04
ymax <- 0.09 # plyr::round_any(max(pred_sub$upper), 0.01, f = ceiling) # 0.1

# Make plot
ggplot(pred_sub) +
  geom_hline(yintercept = 0, color = black) +
  geom_line(
    aes(x = suvr_bl_re, y = fit, color = roi),
    size = 1.2
  ) +
  geom_ribbon(
    aes(x = suvr_bl_re, y = fit, ymin = lower, ymax = upper, fill = roi),
    alpha = 0.15, color=NA
  ) +
  scale_x_continuous(
    limits = c(xmin, xmax),
    breaks = seq(1, 4.5, by = 0.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(ymin, ymax),
    # breaks = seq(ymin, ymax, by = 0.04),
    breaks = seq(-0.06, 0.09, by = 0.03),
    expand = c(0, 0)
  ) +
  scale_color_discrete(
    # type = append(black, pal24),
    type = c(black, cyan2, pink2, green2),
    name = "Region",
    limits = roi_map[[roi]][["rois"]],
    labels = roi_map[[roi]][["labels"]]
  ) +
  scale_fill_discrete(
    # type = append(black, pal24),
    type = c(black, cyan2, pink2, green2),
    name = "Region",
    limits = roi_map[[roi]][["rois"]],
    labels = roi_map[[roi]][["labels"]]
  ) +
  labs(
    x = "Baseline FTP SUVR",
    y = "FTP SUVR Δ/year",
  ) +
  theme_modern() +
  theme(
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length =  unit(5, "pt"),
    axis.text.x = element_text(margin = margin(t = 5)),
    axis.text.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(vjust = 6),
    axis.title.y = element_text(vjust = -4),
    plot.title = element_text(size = rel(1.75), hjust = 0.5)
  )

### Plot baseline SUVR ~ Age at baseline
# Get model residuals
xmod <- mod_long_full
mod_terms <- predict.gam(xmod, type="terms")
icpt <- fixef(xmod)["(Intercept)"]
dat <- data.frame(df)
dat$mod <- "mod_long_full"
dat$resid <- residuals.gam(xmod, type = "response")
dat$partial <- (dat$resid
                + icpt
                + mod_terms[,c("s(suvr_bl_re)")]
                + mod_terms[,c("age_at_ftp_bl_mc")]
                + mod_terms[,c("apoe4_alleles_int")]
)

# Get model estimates
xmod_names <- paste(names(fixef(xmod)), smooths(xmod))
xmin <- round(
  max(
    min(df$age_at_ftp_bl_mc),
    quantile(df$age_at_ftp_bl_mc, 0.25) - (1.5 * IQR(df$age_at_ftp_bl_mc))
  ),
  1)
xmax <- round(
  min(
    max(df$age_at_ftp_bl_mc),
    quantile(df$age_at_ftp_bl_mc, 0.75) + (1.5 * IQR(df$age_at_ftp_bl_mc))
  ),
  1)
xvals <- seq(xmin, xmax, by = 0.1)
pred <- with(df,
             expand.grid(
               suvr_bl_re = mean(df$suvr_bl_re),
               age_at_ftp_bl_mc = xvals,
               apoe4_alleles_int = mean(df$apoe4_alleles_int),
               roi = "acc",
               subj = "LDS0070166")
)
pred <- cbind(
  pred,
  data.frame(
    predict(
      xmod,
      newdata = pred,
      type = "response",
      se.fit = TRUE,
      # terms = c("(Intercept)", "s(suvr_bl_re)", "age_at_ftp_bl_mc", "apoe4_alleles_int"),
      terms = c("(Intercept)", "age_at_ftp_bl_mc", "apoe4_alleles_int"),
      newdata.guaranteed = TRUE
    )
  )
)
pred$roi <- "all"
pred <- transform(pred, upper = fit + (2 * se.fit), lower = fit - (2 * se.fit))

roi <- "global"
ymin <- round(
  max(
    min(df$chg_yr_re),
    quantile(df$chg_yr_re, 0.25) - (1.5 * IQR(df$chg_yr_re))
  ),
  2)
ymax <- round(
  min(
    max(df$chg_yr_re),
    quantile(df$chg_yr_re, 0.75) + (1.5 * IQR(df$chg_yr_re))
  ),
  2)

# Make plot
ggplot(pred) +
  geom_hline(yintercept = 0, color = black) +
  geom_point(
    data = dat, aes(x = age_at_ftp_bl_mc, y = partial),
    # data = dat, aes(x = age_at_ftp_bl_mc, y = chg_yr_re),
    size = 0.05, alpha = 0.33, color = blue2,
  ) +
  geom_line(
    aes(x = age_at_ftp_bl_mc, y = fit, color = roi),
    size = 1.2
  ) +
  geom_ribbon(
    aes(x = age_at_ftp_bl_mc, y = fit, ymin = lower, ymax = upper, fill = roi),
    alpha = 0.2, color=NA
  ) +
  scale_x_continuous(
    limits = c(xmin, xmax),
    breaks = seq(50, 65, by = 5) - mean(df$age_at_ftp_bl),
    labels = seq(50, 65, by = 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(ymin, ymax),
    breaks = seq(ymin, ymax, by = 0.04),
    expand = c(0, 0)
  ) +
  scale_color_discrete(
    type = c(black), # append(black, pal24),
    name = "Region",
    limits = roi_map[[roi]][["rois"]],
    labels = roi_map[[roi]][["labels"]]
  ) +
  scale_fill_discrete(
    type = c(black), # append(black, pal24),
    name = "Region",
    limits = roi_map[[roi]][["rois"]],
    labels = roi_map[[roi]][["labels"]]
  ) +
  labs(
    x = "Age at baseline",
    y = "FTP SUVR Δ/year",
  ) +
  theme_modern() +
  theme(
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length =  unit(5, "pt"),
    axis.text.x = element_text(margin = margin(t = 5)),
    axis.text.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(vjust = 6),
    axis.title.y = element_text(vjust = -4),
    plot.title = element_text(size = rel(1.75), hjust = 0.5)
  )

### Plot baseline SUVR ~ APOE4 alleles
# Get model residuals
xmod <- mod_long_full
mod_terms <- predict.gam(xmod, type="terms")
icpt <- fixef(xmod)["(Intercept)"]
dat <- data.frame(df)
dat$mod <- "mod_long_full"
dat$resid <- residuals.gam(xmod, type = "response")
dat$partial <- (dat$resid
                + icpt
                + mod_terms[,c("s(suvr_bl_re)")]
                + mod_terms[,c("age_at_ftp_bl_mc")]
                + mod_terms[,c("apoe4_alleles_int")]
)

# Get model estimates
xmod_names <- paste(names(fixef(xmod)), smooths(xmod))
pred <- with(df,
             expand.grid(
               suvr_bl_re = mean(df$suvr_bl_re),
               age_at_ftp_bl_mc = 0,
               apoe4_alleles_int = unique(df$apoe4_alleles_int),
               roi = "acc",
               subj = "LDS0070166")
)
pred <- cbind(
  pred,
  data.frame(
    predict(
      xmod,
      newdata = pred,
      type = "response",
      se.fit = TRUE,
      terms = c("(Intercept)", "age_at_ftp_bl_mc", "apoe4_alleles_int"),
      newdata.guaranteed = TRUE
    )
  )
)
pred$roi <- "all"
pred <- transform(pred, upper = fit + (2 * se.fit), lower = fit - (2 * se.fit))

roi <- "global"
ymin <- round(
  max(
    min(df$chg_yr_re),
    quantile(df$chg_yr_re, 0.25) - (1.5 * IQR(df$chg_yr_re))
  ),
  2)
ymax <- round(
  min(
    max(df$chg_yr_re),
    quantile(df$chg_yr_re, 0.75) + (1.5 * IQR(df$chg_yr_re))
  ),
  2)

# Make plot
ggplot(pred) +
  geom_hline(yintercept = 0, color = black) +
  geom_point(
    data = dat, aes(x = apoe4_alleles_int, y = partial),
    # data = dat, aes(x = apoe4_alleles_int, y = chg_yr_re),
    size = 0.05, alpha = 0.33, color = blue2,
    position = position_jitter(w = 0.25, h = 0)
  ) +
  geom_point(
    aes(x = apoe4_alleles_int, y = fit),
    size = 3, color = black
  ) +
  geom_errorbar(
    aes(x = apoe4_alleles_int, ymin = lower, ymax = upper),
    color = black, width = 0.1, size = 0.75
  ) +
  scale_x_continuous(
    limits = c(-0.5, 2.5),
    breaks = c(0, 1, 2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(ymin, ymax),
    breaks = seq(ymin, ymax, by = 0.04),
    expand = c(0, 0)
  ) +
  # scale_color_discrete(
  #   type = reds, # append(black, pal24),
  #   name = "Region",
  #   limits = roi_map[[roi]][["rois"]],
  #   labels = roi_map[[roi]][["labels"]]
  # ) +
  # scale_fill_discrete(
  #   type = pal24, # append(black, pal24),
  #   name = "Region",
  #   limits = roi_map[[roi]][["rois"]],
  #   labels = roi_map[[roi]][["labels"]]
  # ) +
  labs(
    x = "APOE-ε4 alleles",
    y = "FTP SUVR Δ/year",
  ) +
  theme_modern() +
  theme(
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length =  unit(5, "pt"),
    axis.text.x = element_text(margin = margin(t = 5)),
    axis.text.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(vjust = 6),
    axis.title.y = element_text(vjust = -4),
    plot.title = element_text(size = rel(1.75), hjust = 0.5)
  )
