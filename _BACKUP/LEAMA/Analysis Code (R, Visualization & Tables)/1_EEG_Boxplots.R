
#### Code for GPP/GPF Result Boxplots ####


# Packages
library(readxl)
library(openxlsx)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(forcats)
library(grid)

# Read data for POO4h and PPO2
t_singlePOO4h <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_POO4h.xlsx")
t_multiPOO4h <- read.xlsx("Input/PeakMeans/SubMeans_Multitaper_POO4h.xlsx")
t_noopti_sPOO4h <- read.xlsx("Input/PeakMeans/SubMeans_NoOPTICAT_single_POO4h.xlsx")
t_noopti_mPOO4h <- read.xlsx("Input/PeakMeans/SubMeans_NoOPTICAT_multi_POO4h.xlsx")
t_singlePPO2 <- read.xlsx("Input/PeakMeans/SubMeans_Singletaper_PPO2.xlsx")
t_multiPPO2 <- read.xlsx("Input/PeakMeans/SubMeans_Multitaper_PPO2.xlsx")
t_noopti_sPPO2 <- read.xlsx("Input/PeakMeans/SubMeans_NoOPTICAT_single_PPO2.xlsx")
t_noopti_mPPO2 <- read.xlsx("Input/PeakMeans/SubMeans_NoOPTICAT_multi_PPO2.xlsx")

names(t_singlePOO4h) <- NULL
names(t_multiPOO4h) <- NULL
names(t_noopti_sPOO4h) <- NULL
names(t_noopti_mPOO4h) <- NULL
names(t_singlePPO2) <- NULL
names(t_multiPPO2) <- NULL
names(t_noopti_sPPO2) <- NULL
names(t_noopti_mPPO2) <- NULL


#### 1. Create table
OPTICAT = c(rep("With OPTICAT", 40), rep("Without OPTICAT", 40))
Tapers = c(rep("Singletaper", 20), rep("Multitaper", 20), rep("Singletaper", 20), rep("Multitaper", 20))
Contr = c(rep(c(rep("Low", 10), rep("High", 10)),4))
subject = c(rep(1:10, 8))

table <- NA
table <- data.frame(Subject = subject,
                    Opticat = OPTICAT,
                    TaperType = Tapers,
                    Contrast = Contr)

table$Frequency <- c(t_singlePOO4h[1,], t_singlePOO4h[3,], t_multiPOO4h[1,], t_multiPOO4h[3,], 
                     t_noopti_sPOO4h[1,], t_noopti_sPOO4h[3,], t_noopti_mPOO4h[1,], t_noopti_mPOO4h[3,])
table$Power <- c(t_singlePPO2[2,],t_singlePPO2[4,], t_multiPPO2[2,], t_multiPPO2[4,], 
                 t_noopti_sPPO2[2,], t_noopti_sPPO2[4,], t_noopti_mPPO2[2,], t_noopti_mPPO2[4,])

write.xlsx(table, "Output/AllConditions.xlsx")


#### 2. Read the table (to skip previous steps)
table <- read_excel("Output/AllConditions.xlsx")


#### 3. Prepare table
table$Frequency <- as.numeric(table$Frequency)
table$Power <- as.numeric(table$Power)
table$Subject <- factor(table$Subject)
table$Opticat <- factor(table$Opticat)
table$Contrast <- factor(table$Contrast)
table$TaperType <- factor(table$TaperType)

table <- table %>% mutate(TaperType = fct_relevel(TaperType, "Multitaper", after = Inf))
table <- table %>% mutate(Contrast = fct_relevel(Contrast, "High", after = Inf))

tableSingle <- table[!table$TaperType == "Multitaper",]
tableMulti <- table[!table$TaperType == "Singletaper",]


#### 4. Power Boxplots

# Singletaper, for OPTICAT / no OPTICAT
PowSingle <-
  ggplot(tableSingle, aes(x = Contrast, y = Power)) + 
  geom_boxplot(aes(fill = factor(Contrast)), outlier.shape = NA) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "black") +
  facet_grid(~Opticat) +
  geom_point(alpha = 0.7, shape = 1, size = 1.5) + 
  geom_line(aes(group = factor(Subject)), color = "darkgrey", size = 0.3) +
  theme_minimal() +
  ylim(-0.2, 4) +
  ylab('Power [dB]') +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, margin = margin(0, 0, 15, 0)), 
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(margin = margin(0, 17, 0, 0)),
        axis.text.y = element_text(margin = margin(0, 7, 0, 0)),
        text = element_text(size = 15)) +
  scale_fill_manual('Contrast', values = c('antiquewhite2', 'plum4')) +
  ggtitle("Singletaper")

# Multitaper, for OPTICAT / no OPTICAT
PowMulti <-
  ggplot(tableMulti, aes(x = Contrast, y = Power)) + 
  geom_boxplot(aes(fill = factor(Contrast)), outlier.shape = NA) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "black") +
  facet_grid(~Opticat) +
  geom_point(alpha = 0.7, shape = 1, size = 1.5) + 
  geom_line(aes(group = factor(Subject)), color = "darkgrey", size = 0.3) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, margin = margin(0, 0, 15, 0)), 
        axis.ticks.y = element_line(color = 'white'),
        axis.line.y = element_line(color = 'white'),
        axis.title.y = element_text(color = 'white', margin = margin(0, 15, 0, 0)),
        axis.text.y = element_text(color = 'white', margin = margin(0, 9, 0, 0)),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        text = element_text(size = 15),
        # legend.title = element_text(size = 10),
        # legend.position = c(0.87, 0.9),
        # legend.background = element_rect(fill = "white", colour = "white")
        legend.position = "none") +
  ylim(-0.2, 4) +
  ylab('Power [dB]') +
  scale_fill_manual('Contrast', values = c('antiquewhite2', 'plum4')) +
  ggtitle("Multitaper")



#### 5. Frequency Boxplots

# Singletaper, for OPTICAT / no OPTICAT
FreqSingle <-
  ggplot(tableSingle, aes(x = Contrast, y = Frequency)) + 
  geom_boxplot(aes(fill = factor(Contrast)), outlier.shape = NA) +
  facet_grid(~Opticat) +
  geom_point(alpha = 0.7, shape = 1, size = 1.5) + 
  geom_line(aes(group = factor(Subject)), color = "darkgrey", size = 0.3) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "black") +
  theme_minimal() +
  ylim(30, 70) +
  ylab('Frequency [Hz]') +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, margin = margin(0, 0, 15, 0)), 
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(margin = margin(0, 15, 0, 0)),
        text = element_text(size = 15)) +
  scale_fill_manual('Contrast', values = c('antiquewhite2', 'plum4')) +
  ggtitle("Singletaper")

# Multitaper, for OPTICAT / no OPTICAT
FreqMulti <-
  ggplot(tableMulti, aes(x = Contrast, y = Frequency)) + 
  geom_boxplot(aes(fill = factor(Contrast)), outlier.shape = NA) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "black") +
  facet_grid(~Opticat) +
  geom_point(alpha = 0.7, shape = 1, size = 1.5) + 
  geom_line(aes(group = factor(Subject)), color = "darkgrey", size = 0.3) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, margin = margin(0, 0, 15, 0)), 
        axis.ticks.y = element_line(color = 'white'),
        axis.line.y = element_line(color = 'white'),
        axis.title.y = element_text(color = 'white', margin = margin(0, 15, 0, 0)),
        axis.text.y = element_text(color = 'white'),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        text = element_text(size = 15),
        # legend.title = element_text(size = 10),
        # legend.position = c(0.87, 0.9),
        # legend.background = element_rect(fill = "white", colour = "white")
        legend.position = "none") +
  ylim(30, 70) +
  scale_fill_manual('Contrast', values = c('antiquewhite2', 'plum4')) +
  ggtitle("Multitaper")



#### 6. Plot all figures together 

#### Power

# Create layout function (something I found in the internet)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

# New page with specified layout
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2))) # grid with 1 row, 2 columns

# Print the plots
print(PowSingle, vp = vplayout(1, 1))
print(PowMulti, vp = vplayout(1, 2))


#### Frequency

# New page with specified layout
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2))) # grid with 1 row, 2 columns

# Print the plots
print(FreqSingle, vp = vplayout(1, 1))
print(FreqMulti, vp = vplayout(1, 2))





