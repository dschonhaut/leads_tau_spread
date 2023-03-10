library(gamm4)
library(ggplot2)
library(tidyverse)

# ------------------------------------------------------------------------------
# Import and format the data
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
filepath <- file.path(data_dir, "tau_rois_spec.csv")
tau_rois <- read.csv(filepath)
tau_rois <- subset(tau_rois, select = -X)

# Replace missing values with NA
cols_to_replace <- c("minority", "apoe_genotype", "apoe4_alleles", "cdrsob_baseline", "mmse_baseline")
tau_rois[cols_to_replace] <- lapply(tau_rois[cols_to_replace], function(x) {
  x[x == ""] <- NA
  x
})

# Subset the dataframe
tau_rois_eoad <- tau_rois[(tau_rois$dx == "EOAD") & (tau_rois$pet_visits > 1), ]

# Convert strings to factors
tau_rois_eoad <- as.data.frame(unclass(tau_rois_eoad), stringsAsFactors = TRUE)
eoad_subjs <- levels(as.factor(tau_rois_eoad$subj))
rois <- levels(as.factor(tau_rois_eoad$roi))

# ------------------------------------------------------------------------------
