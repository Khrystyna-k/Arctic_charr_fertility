
```{r}
# Load required library
library(dplyr)

# Step 1: Load the base SNP map (BLUPf90 format)
map <- read.table("/home/rabu0002/kalamri-k/GWAS_blupf90_dens2/ac_map_allP.txt_clean", header = FALSE)

# Step 2: Reorder and rename columns for compatibility
# Original columns assumed to be: V1 = CHR, V2 = SNP, V3 = ?, V4 = BP
map <- map[, -3]                         # Remove the unnecessary 3rd column
map <- map[, c(2, 1, 3)]                 # Reorder to: SNP, CHR, BP
colnames(map) <- c("SNP", "CHR", "BP")  # Assign column names

# Step 3: Replace chromosome 40 with label 'Pseudo' (e.g., unplaced scaffolds)
map$CHR[map$CHR == 40] <- "Pseudo"

# Step 4: Load SNP annotation file (with SNP names for BP positions)
all_codes <- read.table("/home/rabu0002/kalamri-k/Motility_geno_2020/map_codes_allSNP.txt", header = TRUE)

# Step 5: Merge annotation with map based on physical position
map_new <- map %>%
  left_join(all_codes, by = c("BP" = "Position"))

# Step 6: Keep only rows with matching SNPs from original map
map_new <- map_new[map_new$SNP %in% map$SNP, ]

# Step 7: Remove duplicate SNPs if any
map_new <- map_new[!duplicated(map_new$SNP), ]

# Step 8: Replace missing SNP names with 'Pseudo'
map_new$Name[is.na(map_new$Name)] <- "Pseudo"

# Step 9: Create final map object with renamed SNPs
map_new2 <- map_new[, c("SNP", "Name", "BP")]
colnames(map_new2) <- c("SNP", "CHR", "BP")

# Step 10: Plot SNP density per chromosome using CMplot
# The color gradient represents bin density (low = green, high = red)
CMplot(
  map_new2,
  type = "p",              # plot points
  plot.type = "d",         # density plot
  bin.size = 1e6,          # 1 Mb bin size
  chr.den.col = c("darkgreen", "yellow", "red"),
  file = 'tiff',           # output format
  memo = "",               # optional label in filename
  dpi = 300,               # high-res output
  main = " ",              # plot title (left blank)
  file.output = TRUE,      # save output
  verbose = TRUE,          # print progress
  width = 9, height = 6    # figure size in inches
)

```

