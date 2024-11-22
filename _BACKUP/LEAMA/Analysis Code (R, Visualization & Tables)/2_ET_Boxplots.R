
### Code for Eye Tracker Result Boxplots ###


# Packages
library(readxl)
library(openxlsx)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(forcats)
library(grid)

# Read data
t <- read.xlsx("Input/Numbers_Eyeevents_Means.xlsx")


#### 1. Create data frame with eye event numbers

# Create columns
sub <- rep(1:nrow(t), 2)
contrast <- c(rep('low', nrow(t)), rep('high', nrow(t)))
sacc <- c(t$Saccade_low, t$Saccade_high)
fix <- c(t$Fixation_low, t$Fixation_high)
blk <- c(t$Blinks_low, t$Blinks_high)

eye <- data.frame(Subject = sub,
                  Contrast = contrast,
                  Saccade = sacc,
                  Fixation = fix,
                  Blink = blk)

eye <- eye %>% mutate(Contrast = fct_relevel(Contrast, "high", after = Inf))


#### 2. Create data frame difference between eye event numbers

subdiff <- rep(1:nrow(t), 3)
events <- c(rep('Saccadic', nrow(t)), rep('Fixational', nrow(t)), rep('Blinks', nrow(t)))
diffs <- c(t$Saccade_high-t$Saccade_low, t$Fixation_high-t$Fixation_low, t$Blinks_high-t$Blinks_low)

Difference <- data.frame(Subject = subdiff,
                         Event = events,
                         Difference = diffs)
Difference$Event <- factor(Difference$Event, levels = c("Saccadic", "Fixational", "Blinks"))


#### 3. Boxplots

# Saccades
Sacc <-
  ggplot(eye, aes(x = Contrast, y = Saccade)) + 
  geom_boxplot(aes(fill = factor(Contrast)), outlier.shape = NA) +
  geom_point(alpha = 0.7, shape = 1, size = 1.5) + 
  geom_line(aes(group = factor(Subject)), color = "darkgrey", size = 0.3) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "black") +
  theme_minimal() +
  ylab("Average Amount per Trial") +
  ylim(0, 2.5) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, margin = margin(0, 0, 15, 0)), 
        legend.position = "none",
        axis.ticks.x = element_line(color = 'white'),
        axis.title.x = element_text(color = 'white', margin = margin(15,0,0,0)),
        axis.text.x = element_text(color = 'white'),
        axis.line.x = element_blank(),
        axis.title.y = element_text(margin = margin(0, 17, 0, 0)),
        text = element_text(size = 15)) +
  scale_fill_manual('Contrast', values = c('antiquewhite2', 'plum4')) +
  ggtitle("Saccadic Events")

# Fixations
Fix <-
  ggplot(eye, aes(x = Contrast, y = Fixation)) + 
  geom_boxplot(aes(fill = factor(Contrast)), outlier.shape = NA) +
  geom_point(alpha = 0.7, shape = 1, size = 1.5) + 
  geom_line(aes(group = factor(Subject)), color = "darkgrey", size = 0.3) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "black") +
  theme_minimal() +
  ylim(0, 2.5) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, margin = margin(0, 0, 15, 0)), 
        axis.ticks = element_line(color = 'white'),
        axis.line = element_line(color = 'white'),
        axis.title.y = element_text(margin = margin(0, 15, 0, 0)),
        axis.title = element_text(color = 'white'),
        axis.title.x = element_text(margin = margin(15,0,0,0)),
        axis.text = element_text(color = 'white'),
        axis.line.x = element_blank(),
        text = element_text(size = 15),
        # legend.title = element_text(size = 10),
        # legend.position = c(0.87, 0.9),
        # legend.background = element_rect(fill = "white", colour = "white")
        legend.position = "none") +
  scale_fill_manual('Contrast', values = c('antiquewhite2', 'plum4')) +
  ggtitle("Fixational Events")

# Blinks
Blk <-
  ggplot(eye, aes(x = Contrast, y = Blink)) + 
  geom_boxplot(aes(fill = factor(Contrast)), outlier.shape = NA) +
  geom_point(alpha = 0.7, shape = 1, size = 1.5) + 
  geom_line(aes(group = factor(Subject)), color = "darkgrey", size = 0.3) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "black") +
  theme_minimal() +
  ylim(0, 2.5) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, margin = margin(0, 0, 15, 0)), 
        axis.ticks = element_line(color = 'white'),
        axis.line = element_line(color = 'white'),
        axis.title.y = element_text(margin = margin(0, 15, 0, 0)),
        axis.title = element_text(color = 'white'),
        axis.title.x = element_text(margin = margin(15,0,0,0)),
        axis.text = element_text(color = 'white'),
        axis.line.x = element_blank(),
        text = element_text(size = 15),
        # legend.title = element_text(size = 10),
        # legend.position = c(0.87, 0.9),
        # legend.background = element_rect(fill = "white", colour = "white")) +
        legend.position = "bottom",
        legend.direction = "horizontal") +
  scale_fill_manual('Contrast', values = c('antiquewhite2', 'plum4')) +
  ggtitle("Blinks")

# Differences
Diff <-
  ggplot(Difference, aes(x = factor(Event), y = Difference)) + 
  geom_boxplot(aes(fill = factor(Event)), outlier.shape = NA) +
  geom_point(alpha = 0.7, shape = 1, size = 1.5) + 
  geom_line(aes(group = factor(Subject)), color = "darkgrey", size = 0.3) +
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 4, color = "black") +
  theme_minimal() +
  xlab("Event Type") +
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, margin = margin(0, 0, 15, 0)), 
        axis.ticks.y = element_line(color = 'white'),
        axis.line.y = element_line(color = 'white'),
        axis.title.y = element_text(color = 'white', margin = margin(0, 15, 0, 0)),
        # axis.text.y = element_text(color = 'black'),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(margin = margin(15,0,0,0)),
        # axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 15),
        legend.position = "none") +
  scale_fill_manual('Contrast', values = c('antiquewhite3', 'antiquewhite3', 'antiquewhite3')) +
  ggtitle("Differences [high - low]")


#### 4. Plot all figures together

# Create layout function (something I found in the internet)
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

# New page with specified layout
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 4))) # grid with 1 row, 2 columns

# Print the plots
print(Sacc, vp = vplayout(1, 1))
print(Fix, vp = vplayout(1, 2))
print(Blk, vp = vplayout(1, 3))
print(Diff, vp = vplayout(1, 4))
