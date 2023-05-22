library(tidyverse)
library(lme4)
library(lmerTest)
library(mgcv)
library(gratia)
library(viridis)
library(broom)

# ------------------------------------------------------------------------------
# Import and format the data
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
tau <- read.csv(file.path(data_dir, "tau-rois-agg_2023-05-12.csv"))
tau_bl <- tau[(tau["visit"]==1) & (tau["parc"]=="roi24"),]
tau_eoad_bl <- tau_bl[(tau_bl["dx"]=="EOAD"),]
remove(tau)
tau_eoad_long <- read.csv(file.path(data_dir, "tau-rois-agg_eoad-long_2023-05-12.csv"))
tau_eoad_long <- tau_eoad_long[(tau_eoad_long["parc"]=="roi24"),]

# Convert strings to factors
tau_bl <- as.data.frame(unclass(tau_bl), stringsAsFactors = TRUE)
tau_eoad_bl <- as.data.frame(unclass(tau_eoad_bl), stringsAsFactors = TRUE)
tau_eoad_long <- as.data.frame(unclass(tau_eoad_long), stringsAsFactors = TRUE)
tau_bl$apoe4_allelesf <- factor(tau_bl$apoe4_alleles, ordered = TRUE)
tau_eoad_bl$apoe4_allelesf <- factor(tau_eoad_bl$apoe4_alleles, ordered = TRUE)
tau_eoad_long$apoe4_allelesf <- factor(tau_eoad_long$apoe4_alleles, ordered = TRUE)

eoad_subjs <- levels(as.factor(tau_eoad_long$subj))
rois <- levels(as.factor(tau_eoad_long$roi))

# Fit random effects model
# ------------------------------------------------------------------------------
# Create an empty dataframe to store the results
df <- tau_eoad_long
df <- df[(df$ftp_visits>1),]
df$subjroi <- factor(paste(df$subj, df$roi, sep = '_'))
results <- data.frame(subjroi = character(),
                      parc = character(),
                      icpt = numeric(),
                      slope = numeric(),
                      stringsAsFactors = FALSE)

# Loop over unique values of df$parc
unique_parc <- unique(df$parc)
for (atlas in unique_parc) {
  # Subset the dataframe for the current parc
  df_ <- df[df$parc == atlas, ]

  # Fit the random effects model
  re_mod <- lmer(suvr ~ 0 + (1|subjroi) + (0 + ftp_yrs_from_bl|subjroi), df_, REML = TRUE)

  # Extract the random effects of the fitted model
  random_effects <- ranef(re_mod)$subjroi

  # Extract the intercept and slope from the random effects
  icpt <- random_effects[[1]]
  slope <- random_effects[[2]]

  # Append the results to the output dataframe
  results <- rbind(results, data.frame(subjroi = rownames(random_effects),
                                       parc = atlas,
                                       icpt = icpt,
                                       slope = slope,
                                       stringsAsFactors = F))
}
# Save the output.
write.csv(results,
          file.path(data_dir, "tau-rois-agg_eoad-subj-random-slopes-icpts_gt1visit_2023-05-14.csv"),
          row.names=F)
# ------------------------------------------------------------------------------













# ------------------------------------------------------------------------------
df <- tau_eoad_long
df$subjroi <- factor(paste(df$subj, df$roi, sep = '_'))
mod0 <- lmer(suvr ~ 0 +
               (1|subjroi) +
               (0+ftp_yrs_from_bl|subjroi),
             data=df,
             REML=T)
summary(mod0)

df <- tau_eoad_long[(tau_eoad_long$ftp_visits>1),]
df$subjroi <- factor(paste(df$subj, df$roi, sep = '_'))
mod0_ <- lmer(suvr ~ 0 +
               (1|subjroi) +
               (0+ftp_yrs_from_bl|subjroi),
             data=df,
             REML=T)
summary(mod0_)
write.csv(ranef(mod0_)$subjroi,
          file.path(data_dir, "tau-rois-agg_eoad-subj-random-slopes-icpts_gt1visit_2023-05-13.csv"))


# Compare EOAD vs. Abeta- CON at baseline:
# Linear regression controlling age + sex
check_cols <- c("subj", "roi", "chg_suvr", "suvr_bl", "age_at_ftp_bl", "ftp_yrs_from_bl")
df <- tau_eoad_long[(tau_eoad_long$roi=="erc"),]
df <- df[complete.cases(df[,check_cols]),]
mod <- bam(
  chg_suvr ~
    s(suvr_bl) +
    s(ftp_yrs_from_bl) +
    s(subj, bs = "re"),
  data = df,
  REML = TRUE,
  family = gaussian,
  select = TRUE
)
layout(matrix(1:4,ncol=2,byrow=TRUE))
gam.check(mod) + abline(0, 1, col="red")
layout(1)
summary(mod)
layout(matrix(1:6,ncol=3,byrow=TRUE))
plot(mod, residuals=TRUE, seWithMean=TRUE)
layout(1)






mod1 <- bam(
  chg_suvr ~
    # Smooths
    s(suvr_bl) +
    s(age_at_ftp_bl) +
    # ti(suvr_bl, age_at_ftp_bl) +
    # Covariates
    s(ftp_yrs_from_bl) +
    s(subj, bs = "re"),
  data = df,
  # REML = TRUE,
  family = gaussian,
  select = TRUE
)
layout(matrix(1:4,ncol=2,byrow=TRUE))
gam.check(mod1) + abline(0, 1, col="red")
layout(1)
summary(mod1)
AIC(mod, mod1)
anova(mod1)
layout(matrix(1:6,ncol=3,byrow=TRUE))
plot(mod1, residuals=TRUE, seWithMean=TRUE)
layout(1)

# asdfjkl;
check_cols <- c("subj", "roi", "chg_suvr", "suvr_bl", "age_at_ftp_bl", "ftp_yrs_from_bl")
df <- tau_eoad_long[(tau_eoad_long$roi=="erc"),]
df <- tau_eoad_long
df <- df[complete.cases(df[,check_cols]),]
mod1 <- bam(
  suvr ~
    s(ftp_yrs_from_bl) +
    s(roi, bs="re") +
    s(ftp_yrs_from_bl, by=roi) +
    # s(age_at_ftp_bl) +
    # ti(ftp_yrs_from_bl, age_at_ftp_bl) +
    s(subj, bs="re"),
  data = df,
  family = gaussian,
  select = TRUE,
  nthreads=c(4, 1)
)
layout(matrix(1:4,ncol=2,byrow=TRUE))
gam.check(mod1) + abline(0, 1, col="red")
layout(1)
summary(mod1)
AIC(mod1)
layout(matrix(1:6,ncol=3,byrow=T))
plot.gam(
  mod1,
  residuals=F,
  seWithMean=T,
  # shift=fixef(mod1)["(Intercept)"],
  shade=T,
  scheme=2,
  pages=4
  )
layout(1)
plot.gam(mod1, residuals=T, seWithMean=T, shift=fixef(mod1)["(Intercept)"], shade=T, scheme=2, select=2)


mod2 <- bam(
  suvr ~
    s(ftp_yrs_from_bl) +
    s(age_at_ftp_bl) +
    ti(ftp_yrs_from_bl, age_at_ftp_bl) +
    s(roi, bs="re") +
    s(subj, bs="re"),
  data = df,
  family = gaussian,
  select = TRUE,
  nthreads=c(4, 1)
)
layout(matrix(1:4,ncol=2,byrow=TRUE))
gam.check(mod2) + abline(0, 1, col="red")
layout(1)
summary(mod2)
AIC(mod1, mod2)
anova(mod1, mod2, test="F")
layout(matrix(1:6,ncol=3,byrow=T))
plot.gam(
  mod2,
  residuals=F,
  seWithMean=T,
  # shift=fixef(mod2)["(Intercept)"],
  shade=T,
  scheme=2,
  pages=4
)
layout(1)
plot.gam(
  mod2,
  residuals=T,
  seWithMean=T,
  shift=fixef(mod2)["(Intercept)"],
  shade=T,
  scheme=2,
  select=2
  )


# asdfjkl;!!!
check_cols <- c("subj", "roi", "chg_suvr", "suvr_bl", "age_at_ftp_bl", "ftp_yrs_from_bl")
df <- tau_eoad_long[(tau_eoad_long$visit>1),]
# df <- tau_eoad_long
df <- df[complete.cases(df[,check_cols]),]
mod1 <- bam(
  chg_suvr ~
    # LINEAR
    #ftp_yrs_from_bl +
    # SMOOTHS
    s(roi, bs="re") +
    # te(suvr_bl, ftp_yrs_from_bl, bs=c("cr", "re"), by=roi) +
    s(suvr_bl, by=ftp_yrs_from_bl, bs="cr", k=24) +
    ti(suvr_bl, ftp_yrs_from_bl, bs=c("cr", "re"), by=roi, k=12) +
    # s(suvr_bl, by=roi) +

    # COVARIATES
    # s(ftp_yrs_from_bl, bs="re") +
    # ti(ftp_yrs_from_bl, suvr_bl) +
    # s(ftp_yrs_from_bl, by=roi) +
    s(subj, bs="re") +
    s(subj, bs="re", by=ftp_yrs_from_bl),
  data = df,
  family = gaussian,
  select = TRUE,
  na.action="na.omit",
  nthreads=c(16, 1)
)
layout(matrix(1:4,ncol=2,byrow=TRUE))
gam.check(mod1) + abline(0, 1, col="red")
layout(1)
summary(mod1)
AIC(mod1)
# AIC(mod1, mod1)
# anova(mod1, mod1, test="F")
# layout(matrix(1:6,ncol=3,byrow=T))
plot.gam(
  mod1,
  residuals=F,
  seWithMean=T,
  # shift=fixef(mod1)["(Intercept)"],
  shade=T,
  scheme=2,
  pages=4
)
layout(1)
plot.gam(
  mod1,
  residuals=F,
  seWithMean=T,
  shift=fixef(mod1)["(Intercept)"],
  shade=T,
  scheme=2,
  select=2
)


mod_age <- bam(
  chg_suvr ~
    # SMOOTHS
    s(suvr_bl) +
    s(age_at_ftp_bl) +
    s(roi, bs="re") +
    s(suvr_bl, by=roi) +
    s(age_at_ftp_bl, by=roi) +
    ti(suvr_bl, age_at_ftp_bl) +
    # COVARIATES
    s(ftp_yrs_from_bl) +
    ti(ftp_yrs_from_bl, suvr_bl) +
    s(ftp_yrs_from_bl, by=roi) +
    ti(ftp_yrs_from_bl, age_at_ftp_bl) +
    s(subj, bs="re"),
  data = df,
  family = gaussian,
  select = TRUE,
  na.action="na.omit",
  nthreads=c(16, 1)
)
layout(matrix(1:4,ncol=2,byrow=TRUE))
gam.check(mod_age) + abline(0, 1, col="red")
layout(1)
summary(mod_age)
AIC(mod1, mod_age)
anova(mod1, mod_age, test="F")
# layout(matrix(1:6,ncol=3,byrow=T))
plot.gam(
  mod_age,
  residuals=F,
  seWithMean=T,
  # shift=fixef(mod_age)["(Intercept)"],
  shade=T,
  scheme=2,
  pages=6
)
layout(1)
plot.gam(
  mod_age,
  residuals=F,
  seWithMean=T,
  shift=fixef(mod_age)["(Intercept)"],
  shade=T,
  scheme=2,
  select=2
)
