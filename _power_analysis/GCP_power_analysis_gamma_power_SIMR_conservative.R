## GCP Simulation for Power Analysis on GAMMA POWER using SIMR
# Tutorial: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12504
library(lme4)
library(simr)
library(ggplot2)

# Example dataset
set.seed(123)
Subject <- factor(rep(1:30, each=10))
condition <- factor(rep(c("A", "B"), 150)) # Lowest and highest condition
gamma_power <- rnorm(300, mean=0, sd=1) + ifelse(condition == "A", 1, 0)
#sd_subj <- 0.25 # Standard deviation for subject effects
#sd_res <- 1 # Standard deviations for residual effects
#subj_eff  <- rnorm(30, 0, sd_subj)
#gamma_power <- subj_eff[Subject] + rnorm(300, 0, sd_res) + ifelse(condition=="A", 1, 0)
data <- data.frame(Subject, condition, gamma_power)

# Fit the mixed model
model <- lmer(gamma_power ~ condition + (1 | Subject), data=data)
summary(model)

# Extend the model to allow for sample size calculation
extended_model <- extend(model, along = "Subject", n=100) # Extend to a larger sample size

# Set the fixed effect size for condition (Effect size of contrast condition on GAMMA POWER)
fixef(extended_model)["conditionB"] <- 0.3

# Power analysis
sim <- powerSim(extended_model, nsim = 1000, progress = FALSE)

# Compute power curve for a range of sample sizes
power_curve <- powerCurve(extended_model, along = "Subject", nsim = 1000, breaks = seq(5, 75, by=5))

# Print, plot and save the power curve
print(power_curve)
summary(power_curve)
power_curve_summary <- summary(power_curve)
power_curve_df <- data.frame(
  Subjects = power_curve_summary$nlevels,
  Power = power_curve_summary$mean,
  Lower = power_curve_summary$lower,
  Upper = power_curve_summary$upper
)
png(file = "/Volumes/methlab/Students/Arne/GCP/figures/power_analysis/GCP_power_analysis_gamma_power_conservative.png", width = 1800, height = 1400, res = 300)
ggplot(power_curve_df, aes(x = Subjects, y = Power)) +
  geom_point(size = 2, color = "blue") +
  geom_line(color = "blue", linetype = "dotted") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 3, color = "blue", alpha = 0.5) +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "grey", size = 0.5) +  
  scale_x_continuous(
    breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75)  
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 0.90, 1)  
  ) +
  labs(
    x = "Subjects",
    y = "Power"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    axis.line = element_line(colour = "black", size = 1),
    panel.grid.major = element_blank(),  # Turn off major grid lines
    panel.grid.minor = element_blank(),  # Turn off minor grid lines
    axis.ticks.length = unit(0.2, "cm"),  # Adjust the length of the tick marks
    axis.ticks = element_line(size = 0.5),  # Size of the tick marks
    axis.ticks.x = element_line(colour = "black", size = 0.5, lineend = "square"),
    axis.ticks.y = element_line(colour = "black", size = 0.5, lineend = "square"),
    axis.ticks.margin = unit(0.3, "cm"),  # Distance between ticks and labels
    axis.text.x = element_text(margin = margin(t = 10)),  # Margins for the x-axis text
    axis.text.y = element_text(margin = margin(r = 10))   # Margins for the y-axis text
  )
dev.off()
