
# Load necessary libraries
library(ggplot2)
library(randomForest)

# Load your dataset
# Replace 'your_data.csv' with the actual path to your dataset
data <- read.csv('your_data.csv')

# Create interaction term
data$Age_Glucose_Interaction <- data$age * data$Glucose

# Train a Random Forest model
set.seed(42)
rf_model <- randomForest(charlson_comorbidity_index ~ age + Glucose + Age_Glucose_Interaction,
                         data = data, ntree = 100, mtry = 2, importance = TRUE)

# Feature Importance
importance <- importance(rf_model)
print(importance)

# Visualization: Interaction Effect
interaction_plot <- ggplot(data, aes(x = Age_Glucose_Interaction, y = charlson_comorbidity_index)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Effect of Age and Glucose Interaction on CCI",
       x = "Age * Glucose Interaction",
       y = "Charlson Comorbidity Index (CCI)") +
  theme_minimal()

# Save Interaction Effect Plot
ggsave("interaction_effect_plot_rf.png", interaction_plot, dpi = 300)

# Visualization: Feature Importance
importance_df <- data.frame(Feature = rownames(importance), Importance = importance[, 1])
importance_plot <- ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance in Random Forest Model",
       x = "Features",
       y = "Importance Score") +
  theme_minimal()

# Save Feature Importance Plot
ggsave("feature_importance_rf.png", importance_plot, dpi = 300)
