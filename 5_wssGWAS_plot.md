
```{r}
# Load required library
library(CMplot)

# -------------------------------
# Load and Prepare SNP Mapping Data
# -------------------------------

# Read SNP map file with chromosome and base pair position
map2 <- read.table("/home/rabu0002/kalamri-k/GWAS_blupf90_mot2/map_chr_codes_5191.txt", header = TRUE)
colnames(map2) <- c("SNP", "CHR", "BP")

# Combine map with trait-specific variance estimates
# Traits: SC = Sperm Concentration, TM = Total Motility, VCL/VAP/VSL = sperm velocity metrics
var_5 <- data.frame(map2, SD_v$SD, TM_v$TM, VCL_v$VCL, VAP_v$VAP, VSL_v$VSL)
colnames(var_5) <- c("SNP", "Chromosome", "Position", "SC", "TM", "VCL", "VAP", "VSL")

# -------------------------------
# Manhattan Plots of SNP Variance Contributions
# -------------------------------

# Save to PDF for high-resolution printing
pdf(file = "wssGBLUP_var_LG.pdf", width = 9, height = 10)
par(mfrow = c(5, 1), mar = c(6, 10, 0, 0))  # 5 stacked plots

# Function to reduce repetition in CMplot calls
plot_trait_variance <- function(data, trait_col, label) {
  CMplot(data[, c("SNP", "Chromosome", "Position", trait_col)],
         plot.type = "m", LOG10 = FALSE, cex = 0.5, multracks = FALSE,
         threshold = 1, threshold.lty = 2, threshold.lwd = 2, threshold.col = "red",
         signal.cex = 0.8, signal.col = "red",
         col = c("navyblue", "orange1"), ylim = c(0, 4),
         cex.axis = 1.5, cex.lab = 1, main = "", file.output = FALSE, verbose = TRUE)
  mtext(label, side = 2, line = 4)
}

plot_trait_variance(var_5, "SC", "SC")
plot_trait_variance(var_5, "TM", "TM")
plot_trait_variance(var_5, "VCL", "VCL")
plot_trait_variance(var_5, "VAP", "VAP")
plot_trait_variance(var_5, "VSL", "VSL")

dev.off()

# -------------------------------
# QQ Plots of Trait P-values
# -------------------------------

# Save all QQ plots to one PNG
png("wssGBLUP_QQ.png", width = 900, height = 700, pointsize = 16)

# Function for QQ plots
plot_qq <- function(data, trait_col, trait_name) {
  CMplot(data[, c("SNP", "CHR", "BP", trait_col)],
         plot.type = "q", LOG10 = FALSE, cex = 0.5,
         threshold = bonf, threshold.lty = 2, threshold.lwd = 2,
         threshold.col = "red", signal.cex = 2, signal.col = "red",
         ylim = c(0, 6), ylab = "", main = trait_name,
         main.cex = 1.5, main.font = 1, file.output = FALSE, verbose = TRUE)
}

plot_qq(pval_5, "SC", "SC")
plot_qq(pval_5, "TM", "TM")
plot_qq(pval_5, "VCL", "VCL")
plot_qq(pval_5, "VAP", "VAP")
plot_qq(pval_5, "VSL", "VSL")

dev.off()

```

