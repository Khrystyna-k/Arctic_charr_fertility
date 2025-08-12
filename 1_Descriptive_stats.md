```{r}

# Load required libraries
library(readxl)
library(dplyr)
library(psych)


# Load and clean the motility data
motility <- read_excel("Motility_results_summary_2020.xlsx") %>%
  select(-c(3, 7))  # Remove unneeded columns by position

# Automatically calculate sperm counts per motility category
motility <- motility %>%
  mutate(across(ends_with("%"), ~ round(.x * `No of sperm` / 100, 2), 
                .names = "{str_remove(.col, ' %')}_count"))

# Summarize total sperm and average velocity traits per fish
sum_no_sperm <- motility %>%
  group_by(id, type, Control) %>%
  summarize(
    total_sperm   = sum(`No of sperm`, na.rm = TRUE),
    motility_count = sum(Motility_count, na.rm = TRUE),
    rapid_count    = sum(Rapid progressive_count, na.rm = TRUE),
    medium_count   = sum(`Medium progressive_count`, na.rm = TRUE),
    nonprog_count  = sum(`Non progressive_count`, na.rm = TRUE),
    VCL            = mean(`VCL um/s`, na.rm = TRUE),
    VAP            = mean(`VAP um/s`, na.rm = TRUE),
    VSL            = mean(`VSL um/s`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    perc_mot     = round(motility_count / total_sperm * 100, 0),
    perc_rapid   = round(rapid_count    / total_sperm * 100, 0),
    perc_medium  = round(medium_count   / total_sperm * 100, 0),
    perc_nonprog = round(nonprog_count  / total_sperm * 100, 0),
    across(c(VCL, VAP, VSL), round, 2)
  )
  
# Subset and describe data for motility
motility_perc <- sum_no_sperm[c(1:4, 12:15)]
motility_perc_breed <- motility_perc %>% filter(type == "breed")
motility_per_cool <- motility_perc %>% filter(type == "cool", Control == "no")
motility_per_control <- motility_perc %>% filter(Control == "control")

mot_means_breed <- describeBy(motility_perc_breed)
mot_means_cool <- describeBy(motility_per_cool)
mot_means_control <- describeBy(motility_per_control)
motility_perc_all <- describeBy(motility_perc)

# Statistical tests for cooling treatment
cool <- sum_no_sperm %>% filter(type == "cool")
describeBy(cool, cool$Control)



# Analysis for breed group by exercise status
breed <- sum_no_sperm %>% filter(type == "breed")
breed$id <- as.numeric(breed$id)
coonc_all6$id <- as.numeric(coonc_all6$id)
breed_add_exercised <- breed %>% inner_join(coonc_all6, "id")
describeBy(breed_add_exercised, "execised")


# Export descriptive statistics tables
tab_df(mot_means_breed[5:8,], show.rownames = TRUE,
       title = "Descriptive statistics for motility (breed group)",
       file = "Stat_motility_breed.doc")

tab_df(mot_means_cool[5:8,], show.rownames = TRUE,
       title = "Descriptive statistics for motility (cool group)",
       file = "Stat_motility_cool.doc")

tab_df(mot_means_control[5:8,], show.rownames = TRUE,
       title = "Descriptive statistics for motility (control group)",
       file = "Stat_motility_cotrol.doc")

tab_df(motility_perc_all[5:8,], show.rownames = TRUE,
       title = "Descriptive statistics for motility (all)",
       file = "Stat_motility_all.doc")

# Velocity summary and export
velocity_means <- sum_no_sperm[c(1:4, 9, 10, 11)]
vel_breed <- velocity_means %>% filter(type == "breed")
vel_cool <- velocity_means %>% filter(type == "cool", Control == "no")
vel_control <- velocity_means %>% filter(Control == "control")

vel_means_breed <- describeBy(vel_breed)
vel_means_cool <- describeBy(vel_cool)
vel_means_control <- describeBy(vel_control)
velocity_all_mean <- describeBy(velocity_means)

tab_df(vel_means_breed[5:7,], show.rownames = TRUE,
       title = "Descriptive statistics for velocity (breed group)",
       file = "Stat_velocity_breed.doc")

tab_df(vel_means_cool[5:7,], show.rownames = TRUE,
       title = "Descriptive statistics for velocity (cool group)",
       file = "Stat_velocity_cool.doc")

tab_df(vel_means_control[5:7,], show.rownames = TRUE,
       title = "Descriptive statistics for velocity (control group)",
       file = "Stat_velocity_control.doc")

tab_df(velocity_all_mean[5:7,], show.rownames = TRUE,
       title = "Descriptive statistics for velocity (all)",
       file = "Stat_velocity_all.doc")


```

