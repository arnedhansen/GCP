
### Code for Comparing the EEG and Saccade Results Directly ###


# Packages
library(readxl)
library(openxlsx)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(forcats)
library(grid)

# Read EEG data for all electrodes
t_single_P1 <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_P1.xlsx")
t_single_P2 <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_P2.xlsx")
t_single_P4 <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_P4.xlsx")
t_single_P6 <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_P6.xlsx")
t_single_PO3 <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_PO3.xlsx")
t_single_PO4 <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_PO4.xlsx")
t_single_POz <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_POz.xlsx")
t_single_Pz <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_Pz.xlsx")
t_single_PPO1 <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_PPO1.xlsx")
t_single_PPO2 <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_PPO2.xlsx")
t_single_POO4h <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_POO4h.xlsx")
t_single_PPO5h <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_PPO5h.xlsx")
t_single_PPO6h <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_PPO6h.xlsx")

# Read ET data
t_eyes <- read.xlsx("Input/Numbers_Eyeevents_Means.xlsx")



#### 1. Calculate differences between contrast conditions for all electrodes

# Power
diff_powSO_P1 <- t(t_single_P1[4,] - t_single_P1[2,])
diff_powSO_P2 <- t(t_single_P2[4,] - t_single_P2[2,])
diff_powSO_P4 <- t(t_single_P4[4,] - t_single_P4[2,])
diff_powSO_P6 <- t(t_single_P6[4,] - t_single_P6[2,])
diff_powSO_PO3 <- t(t_single_PO3[4,] - t_single_PO3[2,])
diff_powSO_PO4 <- t(t_single_PO4[4,] - t_single_PO4[2,])
diff_powSO_POz <- t(t_single_POz[4,] - t_single_POz[2,])
diff_powSO_Pz <- t(t_single_Pz[4,] - t_single_Pz[2,])
diff_powSO_PPO1 <- t(t_single_PPO1[4,] - t_single_PPO1[2,])
diff_powSO_PPO2 <- t(t_single_PPO2[4,] - t_single_PPO2[2,])
diff_powSO_POO4h <- t(t_single_POO4h[4,] - t_single_POO4h[2,])
diff_powSO_PPO5h <- t(t_single_PPO5h[4,] - t_single_PPO5h[2,])
diff_powSO_PPO6h <- t(t_single_PPO6h[4,] - t_single_PPO6h[2,])

# Frequency
diff_freqSO_P1 <- t(t_single_P1[3,] - t_single_P1[1,])
diff_freqSO_P2 <- t(t_single_P2[3,] - t_single_P2[1,])
diff_freqSO_P4 <- t(t_single_P4[3,] - t_single_P4[1,])
diff_freqSO_P6 <- t(t_single_P6[3,] - t_single_P6[1,])
diff_freqSO_PO3 <- t(t_single_PO3[3,] - t_single_PO3[1,])
diff_freqSO_PO4 <- t(t_single_PO4[3,] - t_single_PO4[1,])
diff_freqSO_POz <- t(t_single_POz[3,] - t_single_POz[1,])
diff_freqSO_Pz <- t(t_single_Pz[3,] - t_single_Pz[1,])
diff_freqSO_PPO1 <- t(t_single_PPO1[3,] - t_single_PPO1[1,])
diff_freqSO_PPO2 <- t(t_single_PPO2[3,] - t_single_PPO2[1,])
diff_freqSO_POO4h <- t(t_single_POO4h[3,] - t_single_POO4h[1,])
diff_freqSO_PPO5h <- t(t_single_PPO5h[3,] - t_single_PPO5h[1,])
diff_freqSO_PPO6h <- t(t_single_PPO6h[3,] - t_single_PPO6h[1,])

# Saccadic frequency (on average across trials and subjects)
diff_sacc <- t_eyes$Saccade_high - t_eyes$Saccade_low



#### 2. Save differences between contrast conditions in data frames for all electrodes

# Power
diff_electrodes_pow <- data.frame(Saccade_Diff = diff_sacc,
                                  P1 = diff_powSO_P1,
                                  P2 = diff_powSO_P2,
                                  P4 = diff_powSO_P4,
                                  P6 = diff_powSO_P6,
                                  PO3 = diff_powSO_PO3,
                                  PO4 = diff_powSO_PO4,
                                  POz = diff_powSO_POz,
                                  Pz = diff_powSO_Pz,
                                  PPO1 = diff_powSO_PPO1,
                                  PPO2 = diff_powSO_PPO2,
                                  POO4h = diff_powSO_POO4h,
                                  PPO5h = diff_powSO_PPO5h,
                                  PPO6h = diff_powSO_PPO6h)

colnames(diff_electrodes_pow) <- c('Saccade_Diff', 'P1', 'P2', 'P4', 'P6', 'PO3', 
                                   'PO4', 'POz', 'Pz', 'PPO1', 'PPO2', 
                                   'POO4h', 'PPO5h', 'PPO6h')

# Frequency
diff_electrodes_freq <- data.frame(Saccade_Diff = diff_sacc,
                                  P1 = diff_freqSO_P1,
                                  P2 = diff_freqSO_P2,
                                  P4 = diff_freqSO_P4,
                                  P6 = diff_freqSO_P6,
                                  PO3 = diff_freqSO_PO3,
                                  PO4 = diff_freqSO_PO4,
                                  POz = diff_freqSO_POz,
                                  Pz = diff_freqSO_Pz,
                                  PPO1 = diff_freqSO_PPO1,
                                  PPO2 = diff_freqSO_PPO2,
                                  POO4h = diff_freqSO_POO4h,
                                  PPO5h = diff_freqSO_PPO5h,
                                  PPO6h = diff_freqSO_PPO6h)

colnames(diff_electrodes_freq) <- c('Saccade_Diff', 'P1', 'P2', 'P4', 'P6', 'PO3', 
                                    'PO4', 'POz', 'Pz', 'PPO1', 'PPO2', 
                                    'POO4h', 'PPO5h', 'PPO6h')




#### 3. Calculate sum of participants with effects according to average for all electrodes

# Power
P1_sump <- length(which(diff_electrodes_pow$Saccade_Diff < 0 & diff_electrodes_pow$P1 > 0))
P2_sump <- length(which(diff_electrodes_pow$Saccade_Diff < 0 & diff_electrodes_pow$P2 > 0))
P4_sump <- length(which(diff_electrodes_pow$Saccade_Diff < 0 & diff_electrodes_pow$P4 > 0))
P6_sump <- length(which(diff_electrodes_pow$Saccade_Diff < 0 & diff_electrodes_pow$P6 > 0))
PO3_sump <- length(which(diff_electrodes_pow$Saccade_Diff < 0 & diff_electrodes_pow$PO3 > 0))
PO4_sump <- length(which(diff_electrodes_pow$Saccade_Diff < 0 & diff_electrodes_pow$PO4 > 0))
POz_sump <- length(which(diff_electrodes_pow$Saccade_Diff < 0 & diff_electrodes_pow$POz > 0))
Pz_sump <- length(which(diff_electrodes_pow$Saccade_Diff < 0 & diff_electrodes_pow$Pz > 0))
PPO1_sump <- length(which(diff_electrodes_pow$Saccade_Diff < 0 & diff_electrodes_pow$PPO1 > 0))
PPO2_sump <- length(which(diff_electrodes_pow$Saccade_Diff < 0 & diff_electrodes_pow$PPO2 > 0))
POO4h_sump <- length(which(diff_electrodes_pow$Saccade_Diff < 0 & diff_electrodes_pow$POO4h > 0))
PPO5h_sump <- length(which(diff_electrodes_pow$Saccade_Diff < 0 & diff_electrodes_pow$PPO5h > 0))
PPO6h_sump <- length(which(diff_electrodes_pow$Saccade_Diff < 0 & diff_electrodes_pow$PPO6h > 0))

# Frequency
P1_sumf <- length(which(diff_electrodes_freq$Saccade_Diff < 0 & diff_electrodes_freq$P1 > 0))
P2_sumf <- length(which(diff_electrodes_freq$Saccade_Diff < 0 & diff_electrodes_freq$P2 > 0))
P4_sumf <- length(which(diff_electrodes_freq$Saccade_Diff < 0 & diff_electrodes_freq$P4 > 0))
P6_sumf <- length(which(diff_electrodes_freq$Saccade_Diff < 0 & diff_electrodes_freq$P6 > 0))
PO3_sumf <- length(which(diff_electrodes_freq$Saccade_Diff < 0 & diff_electrodes_freq$PO3 > 0))
PO4_sumf <- length(which(diff_electrodes_freq$Saccade_Diff < 0 & diff_electrodes_freq$PO4 > 0))
POz_sumf <- length(which(diff_electrodes_freq$Saccade_Diff < 0 & diff_electrodes_freq$POz > 0))
Pz_sumf <- length(which(diff_electrodes_freq$Saccade_Diff < 0 & diff_electrodes_freq$Pz > 0))
PPO1_sumf <- length(which(diff_electrodes_freq$Saccade_Diff < 0 & diff_electrodes_freq$PPO1 > 0))
PPO2_sumf <- length(which(diff_electrodes_freq$Saccade_Diff < 0 & diff_electrodes_freq$PPO2 > 0))
POO4h_sumf <- length(which(diff_electrodes_freq$Saccade_Diff < 0 & diff_electrodes_freq$POO4h > 0))
PPO5h_sumf <- length(which(diff_electrodes_freq$Saccade_Diff < 0 & diff_electrodes_freq$PPO5h > 0))
PPO6h_sumf <- length(which(diff_electrodes_freq$Saccade_Diff < 0 & diff_electrodes_freq$PPO6h > 0))



#### 4. Plot number of subjects with coherent effect patterns for all electrodes

# Power: Summary of sums
smmry_pow <- data.frame(P1 = c(P1_sump, P1_sump-10),
                        P2 = c(P2_sump, P2_sump-10),
                        P4 = c(P4_sump, P4_sump-10),
                        P6 = c(P6_sump, P6_sump-10),
                        PO3 = c(PO3_sump, PO3_sump-10),
                        PO4 = c(PO4_sump, PO4_sump-10),
                        POz = c(POz_sump, POz_sump-10),
                        Pz = c(Pz_sump, Pz_sump-10),
                        PPO1 = c(PPO1_sump, PPO1_sump-10),
                        PPO2 = c(PPO2_sump, PPO2_sump-10),
                        POO4h = c(POO4h_sump, POO4h_sump-10),
                        PPO5h = c(PPO5h_sump, PPO5h_sump-10),
                        PPO6h = c(PPO6h_sump, PPO6h_sump-10))

# Power: Preparation
smmry_pow <- gather(smmry_pow, key = "Electrode", value = "Number")
smmry_pow$Group <- c(rep(c("Coherent", "Incoherent"), nrow(smmry_pow)/2))
smmry_pow$Electrode <- factor(smmry_pow$Electrode, levels = c('P1', 'P2', 'P4', 'P6', 'PO3', 'PO4', 'POz', 'Pz', 'PPO1', 'PPO2', 'POO4h', 'PPO5h', 'PPO6h'))

# Power: Figure
Power_Numbers <- ggplot(smmry_pow, aes(x = factor(Electrode), y = Number, fill = Group)) +
  geom_bar(stat = "identity", position = "identity", color = "black", width = 0.7) +
  geom_hline(yintercept = 0, color = "black", size = 0.5)  +
  scale_fill_manual(values = c("Coherent" = "plum4", "Incoherent" = "bisque3"), name = "Value") +
  ylab('Number of Subjects') +
  scale_y_continuous(breaks = seq(-8, 8, by = 2)) +
  theme_classic() +
  theme(text = element_text(size = 55),
        axis.title.y = element_text(margin = margin(0, 17, 0, 0)),
        axis.text.y = element_text(margin = margin(0, 7, 0, 0)),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(size = 0.15),
        legend.position = "none")

# Power: Save
ggsave("Boxplots/Numberplots_Power.pdf", Power_Numbers, height = 15, width = 30)


# Frequency: Summary of sums
smmry_freq <- data.frame(P1 = c(P1_sumf, P1_sumf-10),
                         P2 = c(P2_sumf, P2_sumf-10),
                         P4 = c(P4_sumf, P4_sumf-10),
                         P6 = c(P6_sumf, P6_sumf-10),
                         PO3 = c(PO3_sumf, PO3_sumf-10),
                         PO4 = c(PO4_sumf, PO4_sumf-10),
                         POz = c(POz_sumf, POz_sumf-10),
                         Pz = c(Pz_sumf, Pz_sumf-10),
                         PPO1 = c(PPO1_sumf, PPO1_sumf-10),
                         PPO2 = c(PPO2_sumf, PPO2_sumf-10),
                         POO4h = c(POO4h_sumf, POO4h_sumf-10),
                         PPO5h = c(PPO5h_sumf, PPO5h_sumf-10),
                         PPO6h = c(PPO6h_sumf, PPO6h_sumf-10))

# Frequency: Preparation
smmry_freq <- gather(smmry_freq, key = "Electrode", value = "Number")
smmry_freq$Group <- c(rep(c("Coherent", "Incoherent"), nrow(smmry_freq)/2))
smmry_freq$Electrode <- factor(smmry_freq$Electrode, levels = c('P1', 'P2', 'P4', 'P6', 'PO3', 'PO4', 'POz', 'Pz', 'PPO1', 'PPO2', 'POO4h', 'PPO5h', 'PPO6h'))

# Frequency: Figure
Frequency_Numbers <- ggplot(smmry_freq, aes(x = factor(Electrode), y = Number, fill = Group)) +
  geom_bar(stat = "identity", position = "identity", color = "black", width = 0.7) +
  geom_hline(yintercept = 0, color = "black", size = 0.5)  +
  scale_fill_manual(values = c("Coherent" = "plum4", "Incoherent" = "bisque3"), name = "Value") +
  ylab('Number of Subjects') +
  scale_y_continuous(breaks = seq(-8, 8, by = 2)) +
  theme_classic() +
  theme(text = element_text(size = 55),
        axis.title.y = element_text(margin = margin(0, 17, 0, 0)),
        axis.text.y = element_text(margin = margin(0, 7, 0, 0)),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(size = 0.15),
        legend.position = "none")

# Frequency: Save 
save("Boxplots/Numberplots_Frequency.pdf", Frequency_Numbers, height = 15, width = 30)



#### 5. Scatter plots for exemplary electrodes (PPO2 and POO4h)

diff_sacc <- t_eyes$Saccade_high - t_eyes$Saccade_low

# GPP: Prepare
diff_powSO_PPO2 <- t(t_single_PPO2[4,] - t_single_PPO2[2,])
ComparisonPow_PPO2 <- data.frame(SaccDiff = diff_sacc,
                                 PowDiff = diff_powSO_PPO2)

# GPP: Figure
Power <-
  ggplot(ComparisonPow_PPO2, aes(x = SaccDiff, y = X4)) + 
  geom_point(size = 10, col = 'plum4') +
  theme_classic() +
  ylab('GPP Difference [dB]') +
  xlab('Average Saccade Amount Difference') +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, margin = margin(0, 0, 15, 0)), 
        text = element_text(size = 40),
        axis.title.y = element_text(margin = margin(0, 17, 0, 0) , size = 40),
        axis.text.y = element_text(margin = margin(0, 7, 0, 0)),
        axis.title.x = element_text(margin = margin(17, 0, 0, 0), size = 40),
        axis.line = element_line(size = 0.15)) 

# GPP: Save
ggsave("Boxplots/PowerEEGET.pdf", Power, height = 12, width = 15)



# GPF: Prepare
ComparisonFreq_POO4h <- data.frame(SaccDiff = diff_sacc,
                         POO4h = diff_freqSO_POO4h)

# GPF: Figure
Frequency <-
ggplot(ComparisonFreq_POO4h, aes(x = SaccDiff, y = X3)) + 
  geom_point(size = 10, col = 'plum4') +
  theme_classic() +
  ylab('GPF Difference [Hz]') +
  xlab('Average Saccade Amount Difference') +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, margin = margin(0, 0, 15, 0)), 
        text = element_text(size = 40),
        axis.title.y = element_text(margin = margin(0, 17, 0, 0), size = 40),
        axis.text.y = element_text(margin = margin(0, 7, 0, 0)),
        axis.title.x = element_text(margin = margin(17, 0, 0, 0), size = 40) ,
        axis.line = element_line(size = 0.15))

# GPF: Save
ggsave("/Boxplots/FrequencyEEGET.pdf", Frequency, height = 12, width = 15)



#### 6. Plot both figures together

# Create layout function (something I found in the internet)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

# New page with specified layout
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2))) # grid with 1 rosw, 2 columns

# Print the plots
print(Frequency, vp = vplayout(1, 1))
print(Power, vp = vplayout(1, 2))




