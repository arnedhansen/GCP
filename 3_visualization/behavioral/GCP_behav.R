# GCP behavioral data visualizations

# Setup ----------------------------------------------------------------------------------------------------------------------------------
pacman::p_load(tidyverse, 
               janitor,  # tidy data easily
               patchwork # assemble graphs
)
#install.packages("colorspace")
library(colorspace)
library(ggdist)
library(ggplot2)
library(dplyr)

# Load GCP data
dat <- read.csv("/Volumes/methlab/Students/Arne/GCP/data/features/merged_data.csv")

# Change variable types
dat$ID <- as.integer(dat$ID)                       # ID as integer
dat$Condition <- as.factor(dat$Condition)          # Condition as factor
dat$Accuracy <- as.numeric(dat$Accuracy)           # Accuracy as numeric
dat$ReactionTime <- as.numeric(dat$ReactionTime)   # ReactionTime as numeric
dat$GazeDeviation <- as.numeric(dat$GazeDeviation) # GazeDeviation as numeric
dat$MSRate <- as.numeric(dat$MSRate)               # MSRate as numeric
dat$GammaPower <- as.numeric(dat$GammaPower)       # Probe as character
dat$GammaFreq <- as.numeric(dat$GammaFreq)         # Match as integer

# Define color palette
pal <- c("#FF8C00", "#A034F0")
# pal <- c("#ADD8E6", "#D3D3D3", "#FFC0CB") # Pastel blue, pastel black (grey), pastel red
# pal <- c("#4682B4", "#696969", "#CD5C5C") # Darker blue, darker grey, darker red
# pal <- c("#0D4F8B", "#8B0000") # Very dark blue, very dark grey, very dark red

## ACCURACY across trials RAINPLOT  --------------------------------------------------------------------------------------------------------------
# Compute accuracy average
acc_avg <- dat %>%
  group_by(ID, Condition) %>%
  summarise(
    meanAcc = base::mean(Accuracy, na.rm = TRUE),
    .groups = "drop"
  )

# Precompute sample sizes
 sample_sizes <- acc_avg |> 
   group_by(Condition) |> 
   summarise(n = n(), .groups = "drop")

# Rainplot
acc_avg |> 
  group_by(Condition) |> 
  # Define data to plot
  ggplot(aes(x = meanAcc, y = Condition)) + 
  ggdist::stat_halfeye(
    aes(color = Condition,
        fill = after_scale(lighten(color, .5))),
    adjust = .5,
    width = .5,
    height = .6,
    .width = 0, # Keep at 0 to remove middle line along data
    justification = -.4,
    point_color = NA
  ) +
  # Boxplot
  geom_boxplot(
    aes(color = stage(Condition, after_scale = darken(color, .1, space = "HLS")),
        fill = after_scale(desaturate(lighten(color, .8), .4))),
    width = .35, 
    outlier.shape = NA
  ) +
  # Display data points with jitter on y-axis
  geom_point(
    aes(color = stage(Condition, after_scale = darken(color, .1, space = "HLS"))),
    fill = "white",
    shape = 21,
    stroke = .4,
    size = 1.75,
    position = position_jitter(seed = 1, height = 0.125),
    alpha = .5
  ) +
  geom_point(
    aes(fill = Condition),
    color = "transparent",
    shape = 21,
    stroke = .4,
    size = 1.75,
    alpha = .3,
    position = position_jitter(seed = 1, height = 0.125)
  ) +
  # Add median as text
  stat_summary(
    geom = "text",
    fun = "median",
    aes(label = round(after_stat(x), 2),
        color = stage(Condition, after_scale = darken(color, .1, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    size = 4.5,
    vjust = -3.5
  ) +
  # Add sample size as text
  geom_text(
    data = sample_sizes,
    aes(x = 101, y = Condition, label = paste("n =", n), color = Condition),
    inherit.aes = FALSE, # Prevent inheriting other aesthetics
    family = "Roboto Mono",
    size = 4,
    hjust = 0
  ) +
  # Add the labels of the y-axis
  scale_y_discrete(
    labels = c("Low Contrast", "High Contrast")
  ) +
  # Set x-axis ticks
  scale_x_continuous(
    limits = c(95, 102.5),
    breaks = seq(95, 100, by = 2.5),
    expand = c(.001, .001)
  ) +
  # Disable legend
  scale_color_manual(values = pal, guide = "none") + 
  scale_fill_manual(values = pal, guide = "none") +
  # Set x- and y-axis labels and titles
  labs(
    x = "Accuracy [%]",
    y = NULL,
    title = "Accuracy",
    subtitle = "Accuracy across Trials by Condition",
  ) +
  # Set theme
  theme_minimal(base_family = "Zilla Slab", base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(family = "Roboto Mono"),
    axis.text.y = element_text(
      color = darken(pal, .1, space = "HLS"), 
      size = 15
    ),
    axis.title.x = element_text(margin = margin(t = 10), size = 16),
    plot.subtitle = element_text(
      color = "grey40", hjust = 0,
      margin = margin(0, 0, 20, 0)
    ),
    plot.title.position = "plot",
    plot.margin = margin(15, 15, 10, 15)
  )

## REACTION TIME across trials RAINPLOT  --------------------------------------------------------------------------------------------------------------
# Compute reaction time average
rt_avg <- dat %>%
  group_by(ID, Condition) %>%
  summarise(
    meanRT = base::mean(ReactionTime, na.rm = TRUE),
    .groups = "drop"
  )
rt_avg$meanRT <- rt_avg$meanRT * 1000 # Convert to ms

# Precompute sample sizes
sample_sizes <- rt_avg |> 
  group_by(Condition) |> 
  summarise(n = n(), .groups = "drop")

# Rainplot
rt_avg |> 
  group_by(Condition) |> 
  # Define data to plot
  ggplot(aes(x = meanRT, y = Condition)) + 
  ggdist::stat_halfeye(
    aes(color = Condition,
        fill = after_scale(lighten(color, .5))),
    adjust = .5,
    width = .5,
    height = .6,
    .width = 0, # Keep at 0 to remove middle line along data
    justification = -.4,
    point_color = NA
  ) +
  # Boxplot
  geom_boxplot(
    aes(color = stage(Condition, after_scale = darken(color, .1, space = "HLS")),
        fill = after_scale(desaturate(lighten(color, .8), .4))),
    width = .35, 
    outlier.shape = NA
  ) +
  # Display data points with jitter on y-axis
  geom_point(
    aes(color = stage(Condition, after_scale = darken(color, .1, space = "HLS"))),
    fill = "white",
    shape = 21,
    stroke = .4,
    size = 1.75,
    position = position_jitter(seed = 1, height = 0.125),
    alpha = .5
  ) +
  geom_point(
    aes(fill = Condition),
    color = "transparent",
    shape = 21,
    stroke = .4,
    size = 1.75,
    alpha = .3,
    position = position_jitter(seed = 1, height = 0.125)
  ) +
  # Add median as text
  stat_summary(
    geom = "text",
    fun = "median",
    aes(label = round(after_stat(x), 2),
        color = stage(Condition, after_scale = darken(color, .1, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    size = 4.5,
    vjust = -2
  ) +
  # Add sample size as text
  geom_text(
    data = sample_sizes,
    aes(x = 710, y = Condition, label = paste("n =", n), color = Condition),
    inherit.aes = FALSE, # Prevent inheriting other aesthetics
    family = "Roboto Mono",
    size = 4,
    hjust = 0
  ) +
  # Add the labels of the y-axis
  scale_y_discrete(
    labels = c("Low Contrast", "High Contrast")
  ) +
  # Set x-axis ticks
  scale_x_continuous(
    limits = c(275, 800),
    breaks = seq(300, 700, by = 100),
    expand = c(.001, .001)
  ) +
  # Disable legend
  scale_color_manual(values = pal, guide = "none") + 
  scale_fill_manual(values = pal, guide = "none") +
  # Set x- and y-axis labels and titles
  labs(
    x = "Reaction Time [ms]",
    y = NULL,
    title = "Reaction Time",
    subtitle = "Reaction Time across Trials by Condition",
  ) +
  # Set theme
  theme_minimal(base_family = "Zilla Slab", base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(family = "Roboto Mono"),
    axis.text.y = element_text(
      color = darken(pal, .1, space = "HLS"), 
      size = 15
    ),
    axis.title.x = element_text(margin = margin(t = 10), size = 16),
    plot.subtitle = element_text(
      color = "grey40", hjust = 0,
      margin = margin(0, 0, 20, 0)
    ),
    plot.title.position = "plot",
    plot.margin = margin(15, 15, 10, 15)
  )
