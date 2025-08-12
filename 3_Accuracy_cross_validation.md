```{r}
# ------------------------------------------------------------------
# Cross-validation and accuracy estimation for BLUP and GBLUP models
# Trait: Total sperm motility (TM)
# Author: Khrystyna
# 
# ------------------------------------------------------------------

# Load required libraries
library(ggplot2)
library(GGally)
library(ggpubr)
library(egg)
library(plotrix)
library(dplyr)

# -------------------------------
# 1. Data Preparation
# -------------------------------

# Read phenotype file and select necessary columns
pheno <- read.table("../mot_conc_unique.txt", header = TRUE)
pheno_TM <- pheno[, c("fam_id", "pittag.male", "perc_mot")]

# Read genotype file to retain matching individuals
geno <- read.table("../SNP_file39.txt_clean", header = FALSE)
id_geno <- data.frame(geno$V1)
pheno_id <- pheno_TM %>% inner_join(id_geno, by = c("pittag.male" = "geno.V1"))

# -------------------------------
# 2. Create CV Groups (3-fold)
# -------------------------------

set.seed(4500)
fold <- 3
N <- nrow(pheno_id)
group_size <- floor(N / fold)
rem <- 1:N
id <- list()

for (i in 1:fold) {
  tst <- sample(rem, size = group_size, replace = FALSE)
  rem <- setdiff(rem, tst)
  id[[i]] <- tst
}

# Generate phenotypic replicates with missing TM values for test sets
for (rep in 1:3) {
  pheno_i <- pheno_id
  pheno_i[id[[rep]], "perc_mot"] <- 0
  write.table(pheno_i, file = paste0("pheno_TM_rep1_", rep, ".txt"), row.names = FALSE,
              quote = FALSE, sep = " ", col.names = TRUE)
}

# -------------------------------
# 3. Function to Validate Accuracy
# -------------------------------

validate_accuracy <- function(rep_tag = "1", y_file = "y_star_residual.txt", solution_prefix = "solutions_BLUP_site") {

  y_star <- read.table(y_file, header = FALSE, sep = "	", stringsAsFactors = FALSE)
  blup_id <- read.table("renadd02.ped", header = FALSE)
  id_ystar <- blup_id %>% inner_join(y_star, by = "V1") %>%
    transmute(Id = V10, Id_num = V1, y_star = V11, residual = V12)

  acc_mat <- matrix(NA, nrow = 3, ncol = 1)
  blup_cv <- data.frame()

  for (i in 1:3) {
    # Read phenotype
    val_file <- paste0("pheno_TM_rep", rep_tag, "_", i, ".txt")
    val <- read.table(val_file, header = TRUE, stringsAsFactors = FALSE)
    val$pittag.male <- as.character(val$pittag.male)

    val <- val %>%
      inner_join(id_ystar, by = c("pittag.male" = "Id")) %>%
      filter(perc_mot == 0) %>%
      select(Id = pittag.male, Id_num, TM_star = y_star)

    # Read solution
    sol_file <- paste0(solution_prefix, i, ".txt")
    solution <- read.table(sol_file, header = TRUE, sep = "	", stringsAsFactors = FALSE)
    solution <- solution[solution$trait == 1 & solution$effect == 2, ]

    # Merge and correlate
    blup <- val %>% left_join(solution, by = c("Id_num" = "level"))
    acc_mat[i, 1] <- cor(blup$TM_star, blup$solution)
    blup_cv <- rbind(blup_cv, blup)
  }

  print(paste("Mean accuracy (rep", rep_tag, "):", round(mean(acc_mat[, 1], na.rm = TRUE), 2)))
  print(paste("Standard error:", std.error(acc_mat[, 1])))
}

# Run validation for all 5 replicates
validate_accuracy("1", "y_star_residual.txt", "solutions_BLUP_site")
validate_accuracy("2", "y_star_residual.txt", "solutions2_BLUP_site")
validate_accuracy("3", "y_star_residual.txt", "solutions3_BLUP_site")
validate_accuracy("4", "y_star_residual.txt", "solutions4_BLUP_site")
validate_accuracy("5", "y_star_residual.txt", "solutions5_BLUP_site")

# -------------------------------
# 4. Repeat for GBLUP
# -------------------------------

validate_accuracy_gblup <- function(rep_tag = "1", y_file = "y_star_Gresidual.txt", solution_prefix = "solutions1_GBLUP_site") {

  y_star <- read.table(y_file, header = FALSE, sep = "	", stringsAsFactors = FALSE)
  blup_id <- read.table("renadd02.ped", header = FALSE)
  id_ystar <- blup_id %>% inner_join(y_star, by = "V1") %>%
    transmute(Id = V10, Id_num = V1, y_star = V11, residual = V12)

  acc_mat <- matrix(NA, nrow = 3, ncol = 1)
  blup_cv <- data.frame()

  for (i in 1:3) {
    val_file <- paste0("pheno_TM_rep", rep_tag, "_", i, ".txt")
    val <- read.table(val_file, header = TRUE, stringsAsFactors = FALSE)
    val$pittag.male <- as.character(val$pittag.male)

    val <- val %>%
      inner_join(id_ystar, by = c("pittag.male" = "Id")) %>%
      filter(perc_mot == 0) %>%
      select(Id = pittag.male, Id_num, TM_star = y_star)

    sol_file <- paste0(solution_prefix, i, ".txt")
    solution <- read.table(sol_file, header = TRUE, sep = "	", stringsAsFactors = FALSE)
    solution <- solution[solution$trait == 1 & solution$effect == 2, ]

    blup <- val %>% left_join(solution, by = c("Id_num" = "level"))
    acc_mat[i, 1] <- cor(blup$TM_star, blup$solution)
    blup_cv <- rbind(blup_cv, blup)
  }

  print(paste("Mean GBLUP accuracy (rep", rep_tag, "):", round(mean(acc_mat[, 1], na.rm = TRUE), 2)))
  print(paste("Standard error:", std.error(acc_mat[, 1])))
}

# Run GBLUP validations
validate_accuracy_gblup("1", "y_star_Gresidual.txt", "solutions1_GBLUP_site")
validate_accuracy_gblup("2", "y_star_Gresidual.txt", "solutions2_GBLUP_site")
validate_accuracy_gblup("3", "y_star_Gresidual.txt", "solutions3_GBLUP_site")
validate_accuracy_gblup("4", "y_star_Gresidual.txt", "solutions4_GBLUP_site")
validate_accuracy_gblup("5", "y_star_Gresidual.txt", "solutions5_GBLUP_site")

```

