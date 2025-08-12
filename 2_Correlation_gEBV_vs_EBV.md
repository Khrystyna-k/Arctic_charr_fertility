```{r}
# ------------------------------------------------------------------------------
# EBV vs GEBV Correlation for Motility Traits in Rainbow Trout
# Author: Khrystyna
# Description: This script loads EBVs and GEBVs, matches IDs, filters, merges,
#              calculates correlations, and plots results for five traits.
# ------------------------------------------------------------------------------

# Required libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(psych)
library(tidyr)

# Load unique animal identifiers (change path as needed)
animal_ids <- read.table("mot_conc_unique.txt", header = TRUE)

# Load pedigree file
ped_file <- read.table("renadd03.ped", header = FALSE)
colnames(ped_file) <- c("level", "2", "3", "4", "5", "6", "7", "8", "9", "pittag.male")

# -------------------------- Load EBVs --------------------------

# Helper function to load, join, and filter EBV data
load_ebv <- function(file_path, ped_data, animal_ids, trait_col, effect_val = "3") {
  df <- read.table(file_path, header = TRUE) %>%
    inner_join(ped_data, by = "level") %>%
    filter(pittag.male %in% animal_ids$pittag.male, effect == effect_val)
  return(df)
}

VCL <- load_ebv("blup_results/solutions_VCL", ped_file, animal_ids, "VCL")
VAP <- load_ebv("blup_results/solutions_VAP", ped_file, animal_ids, "VAP")
VSL <- load_ebv("blup_results/solutions_VSL", ped_file, animal_ids, "VSL")
TM  <- load_ebv("blup_results/solutions_tot_mot", ped_file, animal_ids, "TM")

# SD uses a different PED file and effect type
ped_fileSD <- read.table("blup_results/renadd04SD.ped", header = FALSE)
colnames(ped_fileSD) <- c("level", "2", "3", "4", "5", "6", "7", "8", "9", "pittag.male")
SD <- read.table("blup_results/solutions_SD", header = TRUE) %>%
  inner_join(ped_fileSD, by = "level") %>%
  filter(pittag.male %in% animal_ids$X1pittag.male, effect == "4")

# -------------------------- Load GEBVs --------------------------

# Helper function for GEBV loading
load_gebv <- function(file_path, ped_data, animal_ids, effect_val = "3") {
  df <- read.table(file_path, header = TRUE) %>%
    inner_join(ped_data, by = "level") %>%
    filter(pittag.male %in% animal_ids$X1pittag.male, effect == effect_val)
  return(df)
}

gTM  <- load_gebv("gblup_results/solutions_snp_TM", ped_file, animal_ids)
gVCL <- load_gebv("gblup_results/solutions_snp_VCL", ped_file, animal_ids)
gVAP <- load_gebv("gblup_results/solutions_snp_VAP", ped_file, animal_ids)
gVSL <- load_gebv("gblup_results/solutions_snp_VSL", ped_file, animal_ids)

# SD uses a different PED file
gped_fileSD <- read.table("gblup_results/renadd04SD.ped", header = FALSE)
colnames(gped_fileSD) <- c("level", "2", "3", "4", "5", "6", "7", "8", "9", "pittag.male")
gSD <- read.table("gblup_results/solutions_SD", header = TRUE) %>%
  inner_join(gped_fileSD, by = "level") %>%
  filter(pittag.male %in% animal_ids$X1pittag.male, effect == "4")

# -------------------------- Merge EBVs and GEBVs --------------------------

# Format EBV and GEBV data frames
all_ebvs <- data.frame(
  pittag.male = SD$pittag.male,
  SD = SD$solution,
  TM = TM$solution,
  VCL = VCL$solution,
  VAP = VAP$solution,
  VSL = VSL$solution
)

all_gebvs <- data.frame(
  pittag.male = gSD$pittag.male,
  SD = gSD$solution,
  TM = gTM$solution,
  VCL = gVCL$solution,
  VAP = gVAP$solution,
  VSL = gVSL$solution
)

# Merge into one data frame
all <- inner_join(all_ebvs, all_gebvs, by = "pittag.male", suffix = c(".EBV", ".GEBV"))

# -------------------------- Correlation Analysis --------------------------

# Calculate Pearson correlations between EBV and GEBV
cor_results <- data.frame(
  Trait = c("SD", "TM", "VCL", "VAP", "VSL"),
  Correlation = sapply(c("SD", "TM", "VCL", "VAP", "VSL"), function(trait) {
    cor(all[[paste0(trait, ".EBV")]], all[[paste0(trait, ".GEBV")]])
  })
)

print(cor_results)

# -------------------------- Plotting --------------------------

# Shared theme for all plots
My_Theme <- theme(
  legend.position = "none",
  text = element_text(size = 12, color = 'black'),
  plot.title = element_text(hjust = 0.5),
  axis.text = element_text(size = 14)
)

# Plot each trait
plot_trait_corr <- function(trait, label_x = NULL, label_y = NULL) {
  ggplot(all, aes_string(paste0(trait, ".EBV"), paste0(trait, ".GEBV"))) +
    geom_point(shape = 21, fill = "steelblue3", color = "black", size = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    stat_cor(method = "pearson", label.x = label_x, label.y = label_y, size = 5) +
    theme_classic2() +
    My_Theme +
    labs(title = trait, x = " ", y = " ")
}

e1 <- plot_trait_corr("SD", -0.6, 0.7)
e3 <- plot_trait_corr("TM", NULL, 10)
e7 <- plot_trait_corr("VCL")
e8 <- plot_trait_corr("VAP")
e9 <- plot_trait_corr("VSL")

# Save arranged figure
png("Ebv_cor_GEBV.png", width = 9, height = 11, units = "in", res = 300)
fig <- ggarrange(e1, e3, e7, e8, e9, nrow = 5, ncol = 1)
annotate_figure(fig, 
  left = textGrob("GEBVs", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
  bottom = textGrob("EBVs", gp = gpar(cex = 1.3))
)
dev.off()


```

