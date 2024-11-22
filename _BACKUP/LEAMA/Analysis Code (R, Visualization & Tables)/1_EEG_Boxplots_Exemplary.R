
#### Code for Exemplary GPP/GPF Result Boxplots (PPO2 and POO4h) ####

# Packages
library(readxl)
library(openxlsx)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(forcats)
library(grid)

# Read data
POO4h <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_POO4h.xlsx")
PPO2 <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_PPO2.xlsx")


#### 1. Create tables
Frequency <- data.frame(Subject = rep(c(1:10),2),
                        Contrast = c(rep(c("Low"), 10), rep(c("High"), 10)),
                        Frequency = c(t(POO4h[1,]), t(POO4h[3,])))

Power <- data.frame(Subject = rep(c(1:10),2),
                    Contrast = c(rep(c("Low"), 10), rep(c("High"), 10)),
                    Power = c(t(PPO2[2,]), t(PPO2[4,])))

Power <- Power %>% mutate(Contrast = fct_relevel(Contrast, "High", after = Inf))
Frequency <- Frequency %>% mutate(Contrast = fct_relevel(Contrast, "High", after = Inf))


#### 2. Boxplots

# Power (PPO2)
Pow <-
  ggplot(Power, aes(x = Contrast, y = Power)) + 
  geom_boxplot(aes(fill = factor(Contrast)), outlier.shape = NA) +
  geom_point(alpha = 0.7, shape = 1, size = 4) + 
  geom_line(aes(group = factor(Subject)), color = "darkgrey", size = 0.3) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 7, color = "black") +
  theme_minimal() +
  ylab('Power [dB]') +
  ylim(0.5, 4.5) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, margin = margin(0, 0, 15, 0)), 
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(margin = margin(0, 17, 0, 0)),
        axis.text.y = element_text(margin = margin(0, 7, 0, 0)),
        text = element_text(size = 25)) +
  scale_fill_manual('Contrast', values = c('antiquewhite2', 'plum4')) +
  ggtitle("Power (PPO2)")

# Frequency (POO4h)
Freq <-
  ggplot(Frequency, aes(x = Contrast, y = Frequency)) + 
  geom_boxplot(aes(fill = factor(Contrast)), outlier.shape = NA) +
  geom_point(alpha = 0.7, shape = 1, size = 4) + 
  geom_line(aes(group = factor(Subject)), color = "darkgrey", size = 0.3) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 7, color = "black") +
  theme_minimal() +
  ylab('Frequency [Hz]') +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, margin = margin(0, 0, 15, 0)), 
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(margin = margin(0, 17, 0, 0)),
        axis.text.y = element_text(margin = margin(0, 7, 0, 0)),
        text = element_text(size = 25)) +
  scale_fill_manual('Contrast', values = c('antiquewhite2', 'plum4')) +
  ggtitle("Frequency (POO4h)")



#### 3. Plot figures together

# Create layout function (something I found in the internet)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

# New page with specified layout
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2))) # grid with 1 row, 2 columns

# Print the plots
print(Freq, vp = vplayout(1, 1))
print(Pow, vp = vplayout(1, 2))


#### 4. Save
ggsave('Boxplots/PowerPPO2.pdf', Pow, height = 10, width = 8)
ggsave('Boxplots/FrequencyPOO4h.pdf', Freq, height = 10, width = 8)


