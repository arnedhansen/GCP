
#### Code for Eye Event Comparison Table ####

# Packages
library(readxl)
library(openxlsx)
library(tidyverse)
library(dplyr)
library(writexl)


#### 1. Define electrodes
electrodes <- c('P1', 'P2', 'P4', 'P6', 'PO3', 'PO4', 'POz', 'Pz', 'PPO1', 'PPO2', 'POO4h', 'PPO5h', 'PPO6h')


#### 2. Define data type: Singletaper, Multitaper, NoOPTICAT_multi, NoOPTICAT_single, blk.xlsx, blk_ET.xlsx, deviation.xlsx
type = "deviation.xlsx"

if (type == "blk.xlsx") {
  type = "blk"
} else if (type == "blk_ET.xlsx") {
  type = "blk_ET"
} else if (type == "deviation.xlsx") {
  type = "deviation"
}


#### 3. Read data

peaks_P1 <- read.xlsx(paste("Input/PeakMeans/SubMeans_", type, "_P1.xlsx", sep = ""))
peaks_P2 <- read.xlsx(paste("Input/PeakMeans/SubMeans_", type, "_P2.xlsx", sep = ""))
peaks_P4 <- read.xlsx(paste("Input/PeakMeans/SubMeans_", type, "_P4.xlsx", sep = ""))
peaks_P6 <- read.xlsx(paste("Input/PeakMeans/SubMeans_", type, "_P6.xlsx", sep = ""))
peaks_PO3 <- read.xlsx(paste("Input/PeakMeans/SubMeans_", type, "_PO3.xlsx", sep = ""))
peaks_PO4 <- read.xlsx(paste("Input/PeakMeans/SubMeans_", type, "_PO4.xlsx", sep = ""))
peaks_POz <- read.xlsx(paste("Input/PeakMeans/SubMeans_", type, "_POz.xlsx", sep = ""))
peaks_Pz <- read.xlsx(paste("Input/PeakMeans/SubMeans_", type, "_Pz.xlsx", sep = ""))
peaks_PPO1 <- read.xlsx(paste("Input/PeakMeans/SubMeans_", type, "_PPO1.xlsx", sep = ""))
peaks_PPO2 <- read.xlsx(paste("Input/PeakMeans/SubMeans_", type, "_PPO2.xlsx", sep = ""))
peaks_POO4h <- read.xlsx(paste("Input/PeakMeans/SubMeans_", type, "_POO4h.xlsx", sep = ""))
peaks_PPO5h <- read.xlsx(paste("Input/PeakMeans/SubMeans_", type, "_PPO5h.xlsx", sep = ""))
peaks_PPO6h <- read.xlsx(paste("Input/PeakMeans/SubMeans_", type, "_PPO6h.xlsx", sep = ""))


#### 4. Create table with means

# Initialize an empty data frame
Peaks <- NA
Peaks <- data.frame(
  Electrodes = character(),  
  FreqLow_M = numeric(),         
  PowLow_M = numeric(),       
  FreqHigh_M = numeric(),
  PowHigh_M = numeric(),
  FreqLow_SD = numeric(),         
  PowLow_SD = numeric(),       
  FreqHigh_SD = numeric(),
  PowHigh_SD = numeric())


# Iterate over the electrodes
for (electrode in electrodes) {
 
   # Form the variable name dynamically (assuming you have objects like t_PPO1, t_P2, etc.)
  electrode_var <- get(paste0("peaks_", electrode))
  
  # Select the second row and bind it to the Power data frame
  Peaks[electrode, 2:5] <- c(round(apply(electrode_var, 1, mean), 2))
  Peaks[electrode, 6:9] <- c(round(apply(electrode_var, 1, sd), 2))
}

Peaks$Electrodes <- electrodes

# Create to separate data frames
Peaks_Freq <- Peaks[, c("Electrodes", "FreqLow_M", "FreqLow_SD", "FreqHigh_M", "FreqHigh_SD")]
Peaks_Pow <- Peaks[, c("Electrodes", "PowLow_M", "PowLow_SD", "PowHigh_M", "PowHigh_SD")]


#### 5. Include t-values

t_P1 <- read.xlsx(paste("Input/TTestsEEG/TTestresults_", type, "_P1.xlsx", sep = ""))
t_P2 <- read.xlsx(paste("Input/TTestsEEG/TTestresults_", type, "_P2.xlsx", sep = ""))
t_P4 <- read.xlsx(paste("Input/TTestsEEG/TTestresults_", type, "_P4.xlsx", sep = ""))
t_P6 <- read.xlsx(paste("Input/TTestsEEG/TTestresults_", type, "_P6.xlsx", sep = ""))
t_PO3 <- read.xlsx(paste("Input/TTestsEEG/TTestresults_", type, "_PO3.xlsx", sep = ""))
t_PO4 <- read.xlsx(paste("Input/TTestsEEG/TTestresults_", type, "_PO4.xlsx", sep = ""))
t_POz <- read.xlsx(paste("Input/TTestsEEG/TTestresults_", type, "_POz.xlsx", sep = ""))
t_Pz <- read.xlsx(paste("Input/TTestsEEG/TTestresults_", type, "_Pz.xlsx", sep = ""))
t_PPO1 <- read.xlsx(paste("Input/TTestsEEG/TTestresults_", type, "_PPO1.xlsx", sep = ""))
t_PPO2 <- read.xlsx(paste("Input/TTestsEEG/TTestresults_", type, "_PPO2.xlsx", sep = ""))
t_POO4h <- read.xlsx(paste("Input/TTestsEEG/TTestresults_", type, "_POO4h.xlsx", sep = ""))
t_PPO5h <- read.xlsx(paste("Input/TTestsEEG/TTestresults_", type, "_PPO5h.xlsx", sep = ""))
t_PPO6h <- read.xlsx(paste("Input/TTestsEEG/TTestresults_", type, "_PPO6h.xlsx", sep = ""))


# Initialize an empty data frame
Frequeny <- NA
Frequency <- data.frame()
Power <- NA
Power <- data.frame()

# Iterate over the electrodes
for (electrode in electrodes) {
  # Form the variable name dynamically (assuming you have objects like t_P1, t_P2, etc.)
  electrode_var <- get(paste0("t_", electrode))
  
  # Bind rows to dataframes
  Frequency <- rbind(Frequency, round(electrode_var[1, ], 3))
  Power <- rbind(Power, round(electrode_var[2, ], 3))
} 

# Frequency data frame
Frequency$Electrodes <- electrodes
colnames(Frequency) <- c("H", "p_value", "CI_low", "CI_high", "t_value", "SD", "Electrodes")
Frequency <- Frequency[, c("Electrodes", "t_value", "CI_low", "CI_high", "p_value", "SD")]

# Power data frame
Power$Electrodes <- electrodes
colnames(Power) <- c("H", "p_value", "CI_low", "CI_high", "t_value", "SD", "Electrodes")
Power <- Power[, c("Electrodes", "t_value", "CI_low", "CI_high", "p_value", "SD")]


#### 6. Combine the data framese
Frequency_Complete <- inner_join(Peaks_Freq, Frequency, by = "Electrodes")
Power_Complete <- inner_join(Peaks_Pow, Power, by = "Electrodes")


#### 7. Save
write_xlsx(Frequency_Complete, paste("Output/Table_Frequency_", type, ".xlsx", sep = ""), col_names = TRUE)
write_xlsx(Power_Complete, paste("Output/Table_Power_", type, ".xlsx", sep = ""), col_names = TRUE)



