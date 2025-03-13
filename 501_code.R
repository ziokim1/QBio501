# QBio501 Practical data analysis
# Zio Kim

####---- Stimulus Strength vs. EMP Amplitude, Angular Displacement, Latency

# Load necessary libraries
library(ggplot2)
library(ggpubr) # For adding p-values
library(dplyr)  # For data manipulation

# Export rawdata
setwd("~/Downloads/") # Set working directory
data <- read.csv("501_rawdata.csv") # Load data

# Subset data
data1 <- data[, c(2, 3)] # stimulus strength vs. EMG amp.
data2 <- data[, c(2, 4)] # stimulus strength vs. angle
data3 <- data[, c(2, 5)] # stimulus strength vs. latency

# Define a function to perform linear regression and plot
plot_linear_regression <- function(data, x_col, y_col, title, point_color, line_color, xlab, ylab) {
  # Perform linear regression
  formula <- as.formula(paste(y_col, "~", x_col)) # Dynamically create formula
  model <- lm(formula, data = data)
  
  # Extract coefficients
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  
  # Calculate Pearson correlation, R-squared, and p-value
  cor_test <- cor.test(data[[x_col]], data[[y_col]])
  r <- cor_test$estimate
  r2 <- summary(model)$r.squared
  p_value <- cor_test$p.value
  
  # Create equation text
  equation <- paste("y =", round(slope, 2), "* x +", round(intercept, 2))
  
  # Create annotation text
  annotation <- paste(
    equation, "\n",
    "r =", round(r, 2), "\n",
    "R² =", round(r2, 2), "\n",
    "p =", format.pval(p_value, digits = 2)
  )
  
  # Plot scatter plot with regression line
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(color = point_color, size = 1.5, alpha = 0.7) + # Custom point color and size
    geom_smooth(method = "lm", col = line_color, linewidth = 0.5) + # Custom regression line color and thickness
    labs(x = xlab, y = ylab, title = title) + # Custom axis labels and title
    annotate("text", x = min(data[[x_col]]), y = max(data[[y_col]]), 
             label = annotation, hjust = 0, vjust = 1, 
             size = 2.5, col = "black", family = "Helveltica") + # Annotation in Helveltica size 7
    theme_minimal(base_family = "Helveltica") + # Use Helveltica font
    theme(
      plot.title = element_text(size = 14, face = "bold"), # Title size 14, not centered
      axis.title = element_text(size = 7, face = "bold"), # Axis titles size 7
      axis.text = element_text(size = 7), # Axis text size 7
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2), # Light grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Add border
      plot.margin = unit(c(1, 1, 1, 1), "cm") # Add margins
    )
}

# Plot the three datasets with different colors and custom axis labels
plot1 <- plot_linear_regression(data1, x_col = "stimulus_strength", y_col = "emg_amplitude", 
                                title = "a", point_color = "#1f77b4", line_color = "#ff7f0e", 
                                xlab = "Stimulus Strength (a.u.)", ylab = "EMG Amplitude (mV)")
plot2 <- plot_linear_regression(data2, x_col = "stimulus_strength", y_col = "angle", 
                                title = "b", point_color = "#2ca02c", line_color = "#d62728", 
                                xlab = "Stimulus Strength (a.u.)", ylab = "Angle (°)")
plot3 <- plot_linear_regression(data3, x_col = "stimulus_strength", y_col = "latency", 
                                title = "c", point_color = "#9467bd", line_color = "#8c564b", 
                                xlab = "Stimulus Strength (a.u.)", ylab = "Latency (ms)")

# Arrange plots in a grid with a total width of 180 mm
combined_plot_1 <- ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1, widths = c(1, 1, 1))

# Save the combined plot with a width of 180 mm and high resolution
ggsave(filename = "fig2.png", plot = combined_plot_1, 
       width = 180, height = 60, units = "mm", dpi = 300)


####---- Jendrassik manoeuvre

# Subset data
data4 <- data[1:6,6:9]

# Define a function to create boxplots with individual data points and annotate p-values
create_boxplot <- function(data, x_var, y_var, title, y_lab) {
  # Perform Wilcoxon rank-sum test (Mann-Whitney U test)
  wilcox_test <- wilcox.test(data[[y_var]] ~ data[[x_var]])
  wilcox_p_value <- wilcox_test$p.value
  
  # Perform t-test
  t_test <- t.test(data[[y_var]] ~ data[[x_var]])
  t_p_value <- t_test$p.value
  
  # Calculate y-axis limits
  y_max <- max(data[[y_var]])
  y_min <- min(data[[y_var]])
  
  # Create the plot
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) +
    geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) + # Boxplot
    geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.7, color = "black") + # Align dots horizontally
    annotate("text", x = 1.5, y = y_max + 0.1 * (y_max - y_min), 
             label = paste("Wilcoxon p =", format.pval(wilcox_p_value, digits = 2)), 
             size = 2, family = "Helvetica", hjust = 0, vjust = 1) + # Add Wilcoxon p-value
    annotate("text", x = 1.5, y = y_max + 0.05 * (y_max - y_min), 
             label = paste("t-test p =", format.pval(t_p_value, digits = 2)), 
             size = 2, family = "Helvetica", hjust = 0, vjust = 1) + # Add t-test p-value
    labs(title = title, y = y_lab, x = "") + # Axis labels and title
    scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) + # Nature-style colors
    theme_minimal(base_family = "Helvetica") + # Use Helvetica font
    theme(
      plot.title = element_text(size = 14, face = "bold"), # Title size 14
      axis.text.x = element_text(size = 7, face = "bold"), # X-axis labels size 7
      axis.text.y = element_text(size = 7), # Y-axis labels size 7
      axis.title.y = element_text(size = 7, face = "bold"), # Y-axis title size 7
      panel.grid.major = element_blank(), # Remove major grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Add border
      legend.position = "none" # Remove legend
    )
}

# Create the three boxplots
plot4 <- create_boxplot(data4, x_var = "condition", y_var = "emg", 
                        title = "a", y_lab = "EMG Amplitude (mV)")
plot5 <- create_boxplot(data4, x_var = "condition", y_var = "angle.1", 
                        title = "b", y_lab = "Angle (°)")
plot6 <- create_boxplot(data4, x_var = "condition", y_var = "latency.1", 
                        title = "c", y_lab = "Latency (ms)")

# Arrange plots in a grid
combined_plot_2 <- ggarrange(plot4, plot5, plot6, ncol = 3, nrow = 1)

# Save the combined plot with high resolution
ggsave(filename = "fig3.png", plot = combined_plot_2, 
       width = 180, height = 80, units = "mm", dpi = 300)

####---- Passive vs. Voluntary stretching

# Subset the data
data5 <- data[, 10:11]

# Define the order of conditions explicitly
data5$condition.1 <- factor(data5$condition.1, levels = c("Normal", "Voluntary", "Passive"))

# Filter out NA values for statistical tests
data_filtered <- data5 %>%
  filter(!is.na(latency.2))

# Perform Wilcoxon rank-sum test (Mann-Whitney U test) for 2 conditions
wilcox_test <- wilcox.test(latency.2 ~ condition.1, data = data_filtered)
wilcox_p_value <- wilcox_test$p.value

# Perform t-test for 2 conditions
t_test <- t.test(latency.2 ~ condition.1, data = data_filtered)
t_p_value <- t_test$p.value

# Calculate y-axis limits for annotations
y_max <- max(data5$latency.2, na.rm = TRUE)
y_min <- min(data5$latency.2, na.rm = TRUE)

# Create the plot
ggplot(data5, aes(x = condition.1, y = latency.2, fill = condition.1)) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) + # Boxplot
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.7, color = "black") + # Individual data points
  annotate("text", x = 1.5, y = y_max + 0.05 * (y_max - y_min), 
           label = paste("Wilcoxon p =", format.pval(wilcox_p_value, digits = 2)), 
           size = 3, family = "Helvetica", hjust = 0, vjust = 1) + # Add Wilcoxon p-value
  annotate("text", x = 1.5, y = y_max + 0.0 * (y_max - y_min), 
           label = paste("t-test p =", format.pval(t_p_value, digits = 2)), 
           size = 3, family = "Helvetica", hjust = 0, vjust = 1) + # Add t-test p-value
  labs(title = "Achilles Tendon Reflex Latency", y = "Latency (ms)", x = "") + # Axis labels and title
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) + # Nature-style colors
  theme_minimal(base_family = "Helvetica") + # Use Helvetica font
  theme(
    plot.title = element_text(size = 14, face = "bold"), # Title size 14
    axis.text.x = element_text(size = 7, face = "bold"), # X-axis labels size 7
    axis.text.y = element_text(size = 7), # Y-axis labels size 7
    axis.title.y = element_text(size = 7, face = "bold"), # Y-axis title size 7
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), # Add border
    legend.position = "none" # Remove legend
  )


# Save the plot with high resolution
ggsave(filename = "fig4.png", width = 180, height = 80, units = "mm", dpi = 300)
