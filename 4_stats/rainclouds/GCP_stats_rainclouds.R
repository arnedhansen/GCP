# Raincloud plots for AOC N-back data

# Load necessary libraries
library(ggplot2)
library(ggdist)
library(dplyr)
library(colorspace)  # provides darken(), lighten(), desaturate()

# Define colour palette
pal <- c("#FF8C00", "#A034F0", "#159090") # DataViz workshop colors
pal <- c("#ADD9E6", "#99CC99", "#FFB3BF") # Light AOC pastel colors
pal <- c("#7998A1", "#6B8F6B", "#B37D86") # Dark AOC pastel colors
pal <- c("#93B8C4", "#82AD82", "#D998A2") # Perfect AOC pastel colors

# Function to add sample size as text (offset adjusted by a fraction of the range)
add_sample <- function(x) {
  offset <- 0.025 * diff(range(x))
  return(c(y = max(x) + offset, label = length(x)))
}

# Read in the data (adjust the file path if needed)
dat <- read.csv("/Volumes/methlab/Students/Arne/AOC/data/features/merged_data_nback.csv")

# Transform reaction times to milliseconds
dat$ReactionTime <- dat$ReactionTime * 1000

# Create the GazeStd column as the average of GazeStdX and GazeStdY
dat$GazeStd <- (dat$GazeStdX + dat$GazeStdY) / 2

# Convert Condition to a factor with labels "1-back", "2-back", "3-back"
dat$Condition <- factor(dat$Condition, levels = c(1, 2, 3), labels = c("1-back", "2-back", "3-back"))

# Define the list of variables, y-axis labels and save names (to be used in file naming)
variables <- c("Accuracy", "ReactionTime", "GazeDeviation", "GazeStd", "MSRate", "Fixations", "Saccades", "AlphaPower", "IAF")
titles <- c("Accuracy", "Reaction Time", "Gaze Deviation", "Gaze Standard Deviation",
            "Microsaccade Rate", "Fixations", "Saccades", "Alpha Power", "Individual Alpha Frequency")
y_labels  <- c("Accuracy [%]", "Reaction Time [ms]", "Gaze Deviation [px]", "Gaze Std [px]",
               "Microsaccade Rate [ms/s]", "Fixations", "Saccades", "Alpha Power [\u03BCV²/Hz]", "IAF [Hz]")
save_names <- c("acc", "rt", "gazedev", "gazestd", "ms", "fix", "sacc", "pow", "iaf")

# Remove outliers using the IQR method (1.5 * IQR rule)
dat <- dat %>%
  group_by(Condition) %>%
  mutate(across(all_of(variables), ~{
    lower <- quantile(.x, 0.25, na.rm = TRUE) - 1.5 * IQR(.x, na.rm = TRUE)
    upper <- quantile(.x, 0.75, na.rm = TRUE) + 1.5 * IQR(.x, na.rm = TRUE)
    replace(.x, .x < lower | .x > upper, NA)
  })) %>%
  ungroup()

# Define the output directory and create it if it doesn't exist
output_dir <- "/Volumes/methlab/Students/Arne/AOC/figures/stats/rainclouds"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Loop over each variable to create and save raincloud plots
for(i in seq_along(variables)) {
  
  var <- variables[i]
  y_lab <- y_labels[i]  # now used as the y-axis label since data is on the y-axis
  save_name <- save_names[i]
  
  p <- dat %>% 
    group_by(Condition) %>% 
    ggplot(aes(x = Condition, y = .data[[var]])) +
    
    # Raincloud (half-eye) plot with density layer
    ggdist::stat_halfeye(
      aes(color = Condition,
          fill = after_scale(lighten(color, 0.5))),
      adjust = 0.5,
      width = 0.5,
      height = 0.6,
      .width = 0,              # removes the middle line along the data
      justification = 1.55,   # negative value places the density on the left
      side = "left",
    ) +
    
    # Boxplot layer using the same colour as the density
    geom_boxplot(
      aes(color = Condition,
          fill = after_scale(lighten(color, 0.5))),
      width = 0.35, 
      outlier.shape = NA
    ) +
    
    # Jittered data points (adjusting jitter to the x-axis because x is now categorical)
    geom_point(
      aes(color = Condition,
          fill = after_scale(lighten(color, 0.5))),
      shape = 21,
      stroke = 0.4,
      size = 3,
      position = position_jitter(seed = 1, width = 0.125),
      alpha = 0.65
    ) +
    geom_point(
      aes(fill = Condition),
      color = "transparent",
      shape = 21,
      stroke = 0.4,
      size = .25,
      alpha = 0.4,
      position = position_jitter(seed = 1, width = 0.125)
    ) +
    
    # Add median value as text, nudged to the right of the boxplot
    # stat_summary(
    #   geom = "text",
    #   fun = "median",
    #   aes(label = round(after_stat(y), 2), color = Condition),
    #   family = "Roboto Mono",
    #   fontface = "bold",
    #   size = 4.5,
    #   nudge_y = 12.3
    # ) +
    
    # Commenting out sample size addition
    # stat_summary(
    #   geom = "text",
    #   fun.data = add_sample,
    #   aes(label = paste("n =", after_stat(label), color = Condition)),
    #   family = "Roboto Mono",
    #   size = 4,
    #   nudge_x = 0.3
    # ) +
    
    # Set manual colours for consistency
    scale_color_manual(values = pal, guide = "none") + 
    scale_fill_manual(values = pal, guide = "none") +
    
    # Define labels and titles (now Condition on x-axis, data on y-axis)
    labs(
      x = "Condition",
      y = y_lab,
      title = titles[i],
      subtitle = paste("N-back", titles[i], "by Condition")
    ) +
    
    theme_minimal(base_family = "Zilla Slab", base_size = 15) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),  # remove vertical grid lines
      axis.ticks = element_blank(),
      axis.text.x = element_text(family = "Roboto Mono"),
      axis.text.y = element_text(family = "Roboto Mono", size = 15),
      axis.title.y = element_text(margin = margin(r = 10), size = 16),
      plot.subtitle = element_text(
        colour = "grey40", hjust = 0,
        margin = margin(0, 0, 20, 0)
      ),
      plot.title.position = "plot",
      plot.margin = margin(15, 15, 10, 15)
    )
  
  # Adjust the y-axis
  # "Accuracy", "ReactionTime", "GazeDeviation", "GazeStd", "MSRate", "Fixations", "Saccades", "AlphaPower", "IAF"
  if(var == "Accuracy") {
    p <- p + scale_y_continuous(
      limits = c(70, 100),
      breaks = seq(75, 100, by = 5),
      expand = c(0.001, 0.001)
    ) 
  }
  if(var == "GazeDeviation") {
    p <- p + scale_y_continuous(
      limits = c(0, 60),
      breaks = seq(0, 60, by = 20),
      expand = c(0.001, 0.001)
    ) 
  }
  if(var == "MSRate") {
    p <- p + scale_y_continuous(
      limits = c(0, 3),
      breaks = seq(0, 3, by = 1),
      expand = c(0.001, 0.001)
    ) 
  }
  if(var == "Fixations") {
    p <- p + scale_y_continuous(
      limits = c(0, 7),
      breaks = seq(0, 7, by = 1),
      expand = c(0.001, 0.001)
    ) 
  }
  if(var == "IAF") {
    p <- p + scale_y_continuous(
      limits = c(8, 13),
      breaks = seq(8, 13, by = 1),
      expand = c(0.001, 0.001)
    ) 
  }
  
  
  # Save the plot as a PNG file
  ggsave(filename = file.path(output_dir, paste0("AOC_stats_rainclouds_", save_name, "_nback.png")),
         plot = p, width = 8, height = 6, dpi = 300)
}
