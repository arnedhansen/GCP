
#### Code for Eye Event Comparison Table ####


# Packages
library(readxl)
library(openxlsx)
library(tidyverse)
library(dplyr)
library(writexl)

# Read data
t <- read.xlsx("Input/Numbers_Eyeevents_Means.xlsx")


#### 1. Create table with all means

Blink_t <- t.test(t$Blinks_high, t$Blinks_low, paired = TRUE)
Fixation_t <- t.test(t$Fixation_high, t$Fixation_low, paired = TRUE)
Saccade_t <- t.test(t$Saccade_high, t$Saccade_low, paired = TRUE)


#### 2. Create table with all means
Eye_Summary <- data_frame(Event = c("Blink", "Fixation", "Saccade"),
                          MeanLow = round(c(mean(t$Blinks_low), mean(t$Fixation_low), mean(t$Saccade_low)), 2),
                          SDLow = round(c(sd(t$Blinks_low), sd(t$Fixation_low), sd(t$Saccade_low)), 2),
                          MeanHigh = round(c(mean(t$Blinks_high), mean(t$Fixation_high), mean(t$Saccade_high)), 2),
                          SDHigh = round(c(sd(t$Blinks_high), sd(t$Fixation_high), sd(t$Saccade_high)), 2),
                          t_value = round(c(Blink_t$statistic, Fixation_t$statistic, Saccade_t$statistic), 3),
                          CIlow = round(c(Blink_t$conf.int[1], Fixation_t$conf.int[1], Saccade_t$conf.int[1]), 3),
                          CIhigh = round(c(Blink_t$conf.int[2], Fixation_t$conf.int[2], Saccade_t$conf.int[2]), 3),
                          p_value = round(c(Blink_t$p.value, Fixation_t$p.value, Saccade_t$p.value), 3),
                          SD = round(c(Blink_t$stderr*sqrt(10), Fixation_t$stderr*sqrt(10), Saccade_t$stderr*sqrt(10)), 3))

#### 3. Save
write_xlsx(Eye_Summary, 'Output/Table_ETEvents.xlsx', col_names = TRUE)
