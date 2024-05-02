#### adjusted from Henk van den Toorn #####

library(lme4)
library(tidyverse)
library(magrittr)
library(readxl)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Read Genes output table can be found in PRIDE submission 
df <- read.table("report.gg_matrix.tsv",header = TRUE, sep="\t", na.strings=c("","NA"))

# Read Important Genes file
important_genes <- read_excel("Important proteins.xlsx") %>% as.data.frame()

df <- df %>% pivot_longer(cols = -Genes,
                          names_to = "File.Name",
                          values_to = "intensity")

# Remove rows containing "FT" or "T10" as are not samples of interest
df<- df[!grepl("_FT_|_T10_", df$File.Name),]

#Clean File Names
df<- df %>% mutate(File.Name = case_when(grepl("CTRL_P3_V1", File.Name) ~ "C1_T0",
                                       grepl("CTRL_P3_V2", File.Name) ~"C1_T1",
                                       grepl("CTRL_P3_V3", File.Name) ~"C1_T2",
                                       grepl("CTRL_P4_V1", File.Name) ~"C2_T0",
                                       grepl("CTRL_P4_V2", File.Name) ~"C2_T1",
                                       grepl("CTRL_P4_V3", File.Name) ~"C2_T2",
                                       grepl("P1_T0", File.Name) ~"P1_T0",
                                       grepl("P1_T1", File.Name) ~"P1_T1",
                                       grepl("P1_T2", File.Name) ~"P1_T2",
                                       grepl("P1_T3", File.Name) ~"P1_T3",
                                       grepl("P1_T4", File.Name) ~"P1_T4",
                                       grepl("P1_T5", File.Name) ~"P1_T5",
                                       grepl("P1_T6", File.Name) ~"P1_T6",
                                       grepl("P1_T7", File.Name) ~"P1_T7",
                                       grepl("P1_T8", File.Name) ~"P1_T8",
                                       grepl("P1_T9", File.Name) ~"P1_T9",
                                       grepl("P2_T0", File.Name) ~"P2_T0",
                                       grepl("P2_T1", File.Name) ~"P2_T1",
                                       grepl("P2_T2", File.Name) ~"P2_T2",
                                       grepl("P2_T3", File.Name) ~"P2_T3",
                                       grepl("P2_T4", File.Name) ~"P2_T4"))

#Calculate median of injections
df <- df %>%
  group_by(File.Name, Genes) %>%
  summarise(Median_Intensity = median(intensity, na.rm = TRUE))

#log2 transform intensities
df <- df %>%
  filter(!is.na(Median_Intensity)) %>%
  mutate(Log2_Median_Intensity = log2(Median_Intensity))

# Standards are in g/ml, the table contains ranges
standards <- list( "ALB"      = c(3.5, 5.2) * 1e-2,
                   "FBA"      = c(2.0, 4.0) * 1e-3,
                   "A2M"      = c(0.9, 4.0) * 1e-3,
                   "SERPINA1" = c(7.8, 20)  * 1e-4,
                   "HP"       = c(3.0, 22)  * 1e-4,
                   "TTR"      = c(2.8, 3.5) * 1e-4,
                   "CP"       = c(1.5, 6.0) * 1e-4,
                   "F2"       = c(1.0, 1.0) * 1e-4,
                   "KLKB1"    = c(5.0, 5.0) * 1e-5,
                   "C6"       = c(4.8, 6.4) * 1e-5,
                   "C9"       = c(4.7, 6.9) * 1e-5,
                   "F12"      = c(2.9, 2.9) * 1e-5,
                   "C1R"      = c(2.5, 3.8) * 1e-5,
                   "CFP"      = c(2.4, 3.2) * 1e-5,
                   "C2"       = c(2.2, 3.4) * 1e-5,
                   "VWF"      = c(7, 7)     * 1e-6,
                   "F10"      = c(5.0, 5.0) * 1e-6,
                   "F9"       = c(4.0, 4.0) * 1e-6,
                   "TFRC"     = c(0.8, 1.8) * 1e-6,
                   "F7"       = c(1.0, 1.0) * 1e-6,
                   "MBL2"     = c(0.3, 4.1) * 1e-6,
                   "B2M"      = c(8.0, 24)  * 1e-7,
                   "F8"       = c(1.0, 1.0) * 1e-7,
                   "IGF2"     = c(9.9, 50)  * 1e-8,
                   "MB"       = c(6.0, 85)  * 1e-9,
                   "PRL"      = c(1.0, 7.0) * 1e-9,
                   "INS"      = c(2.0, 8.4) * 1e-10)
  
# Convert data into tibble, and units into g/ml
standards <- as_tibble(t(as_tibble(standards)), rownames="Genes") %>%
  rename("lower_ref_conc" = "V1", "upper_ref_conc" = "V2") %>%
  mutate(refconc = (lower_ref_conc + upper_ref_conc) / 2)
print(standards)

df %<>% mutate(File.Name = as.factor(File.Name))

calibration <- standards %>%
  right_join(df) %>%
  filter(!is.na(File.Name)) %>%
  mutate(log_refconc = log(refconc))

fitted <- calibration %>%
  group_by(File.Name) %>%
  nest() %>%
  mutate(fit = map(data, ~lm(.x$log_refconc ~ .x$Log2_Median_Intensity))) %>%
  mutate(log_concentration = map2(fit, data, ~ predict(.x, .y))) %>%
  select(-fit) %>% unnest(c(data, log_concentration)) %>%
  mutate(concentration = exp(log_concentration))

fitted<- fitted %>% 
  mutate(concentration = concentration * 100000) %>%  #convert concentration to mg/dl
  select(Genes,concentration,File.Name)%>%
  group_by(Genes) %>%
  pivot_wider(names_from = "File.Name",
              values_from = "concentration")

#keep Genes that are present in Important Genes file
merged_df <- subset(fitted, (Genes %in% important_genes$Genes))

# Find the lowest value in a dataframe
lowest_value <- min(unlist(merged_df), na.rm = TRUE)

# Impute the lowest value to all NA 
merged_df[is.na(merged_df)] <- as.numeric(lowest_value)
merged_df <- merged_df %>% mutate_at(vars(-1), as.numeric)

#save new files
openxlsx::write.xlsx(merged_df,"Supplementary Table 1.xlsx")

