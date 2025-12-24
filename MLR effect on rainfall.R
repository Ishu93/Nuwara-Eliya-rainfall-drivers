cat("\f")
rm(list = ls())

library(dplyr)
library(ggplot2)
library(trend)      # For Mann-Kendall test
library(car)       # For VIF and ANOVA
library(lmtest)    # For Durbin-Watson
library(effects)   # For effect size (eta-squared)
library(urca)
library(tseries)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# --- Load CSV ---
file_path <- choose.files()
df <- read.csv(file_path)

# --- Define seasons ---
seasons <- list(
  FIM = c(3, 4),          # First Inter-monsoon
  SWM = 5:9,              # Southwest Monsoon
  SIM = c(10, 11),        # Second Inter-monsoon
  NEM = c(12, 1, 2)       # Northeast Monsoon
)

# Step 1: Identify Trends with Mann-Kendall Test
vars <- c("rainfall", "t", "RH", "windspeed")
mk_results <- lapply(vars, function(var) {
  mk_test <- mk.test(df[[var]])
  data.frame(Variable = var, P_Value = mk_test$p.value, Trend = ifelse(mk_test$p.value < 0.05, "Significant", "Non-significant"))
})
mk_results <- do.call(rbind, mk_results)
cat("\nMann-Kendall Trend Test Results:\n")
print(mk_results)

# Step 2: Detrend Data (using linear regression against time)
df$year_month <- as.numeric(df$year + (df$month - 1) / 12)  # Continuous time variable
df_detrended <- df
for (var in vars) {
  fit <- lm(df[[var]] ~ year_month, data = df)
  df_detrended[[paste0(var, "_detrended")]] <- resid(fit)
}

vars <- c("rainfall", "t", "RH", "windspeed")

test_trend <- function(x, name) {
  cat("\n-----------------------------------------\n")
  cat("Variable:", name, "\n")
  cat("-----------------------------------------\n")
  
  # Deterministic trend test
  model <- lm(x ~ df$year_month)
  cat("\nDeterministic Trend Test (linear regression slope):\n")
  print(summary(model)$coefficients[2,])  # slope significance
  
  # ADF Test (stochastic trend)
  cat("\nADF Test (H0 = Unit Root / Stochastic Trend):\n")
  adf <- ur.df(x, type="trend", selectlags="AIC")
  print(summary(adf))
  
  # KPSS Test (H0 = Trend-stationary / deterministic trend)
  cat("\nKPSS Test:\n")
  print(kpss.test(x, null="Trend"))
  
  cat("\nInterpretation Guide:\n")
  cat(" • ADF fails & KPSS rejects → Stochastic trend\n")
  cat(" • ADF rejects & KPSS fails → Deterministic trend\n")
  cat(" • Both reject → Mixed trend\n")
  cat(" • Both fail → No trend (stationary)\n\n")
}

for (v in vars) {
  test_trend(df[[v]], v)
}


# Step 3: Assign Seasons
df_detrended$season <- case_when(
  df_detrended$month %in% seasons$FIM ~ "FIM",
  df_detrended$month %in% seasons$SWM ~ "SWM",
  df_detrended$month %in% seasons$SIM ~ "SIM",
  df_detrended$month %in% seasons$NEM ~ "NEM",
  TRUE ~ NA_character_
)

# Step 4: MLR for FIM Season
cat("\n\n--- Analysis for Season: FIM ---\n")

# Subset data for the season
season_data_FIM <- df_detrended %>% filter(season == "FIM")

# MLR model for the season
mlr_model_FIM <- lm(rainfall_detrended ~ t_detrended + RH_detrended + windspeed_detrended, 
                    data = season_data_FIM)

# MLR Summary
cat("\nMLR Summary for FIM:\n")
print(summary(mlr_model_FIM))

# Store MLR results as a data frame
mlr_results_FIM <- as.data.frame(summary(mlr_model_FIM)$coefficients)

# Residual plots (normality and homoscedasticity)
par(mfrow = c(2, 2))
plot(mlr_model_FIM, which = 1:2, main = "Residuals vs Fitted and Q-Q Plot (FIM)")
plot(mlr_model_FIM, which = 3, main = "Scale-Location (FIM)")
hist(resid(mlr_model_FIM), main = "Histogram of Residuals (FIM)", xlab = "Residuals")
par(mfrow = c(1, 1))

# Multicollinearity (VIF should be <5)
cat("\nVariance Inflation Factors (VIF) for FIM:\n")
vif_vals_FIM <- vif(mlr_model_FIM)
print(vif_vals_FIM)
vif_results_FIM <- data.frame(Variable = names(vif_vals_FIM), VIF = vif_vals_FIM)

# Autocorrelation (Durbin-Watson, p<0.05 indicates issue)
cat("\nDurbin-Watson Test for Autocorrelation (FIM):\n")
print(dwtest(mlr_model_FIM))

# ANOVA on MLR Model
anova_result_FIM <- anova(mlr_model_FIM)
cat("\nANOVA Results for FIM:\n")
print(anova_result_FIM)

# Calculate Eta-Squared (η²)
ss_total_FIM <- sum((season_data_FIM$rainfall_detrended - mean(season_data_FIM$rainfall_detrended))^2)
eta_squared_FIM <- anova_result_FIM$"Sum Sq" / ss_total_FIM
anova_result_FIM$Eta_Squared <- eta_squared_FIM
cat("\nANOVA with Eta-Squared for FIM:\n")
print(anova_result_FIM)

# ANOVA Assumption Checks
# Shapiro-Wilk test for normality
cat("\nShapiro-Wilk Test for Residual Normality (FIM):\n")
print(shapiro.test(resid(mlr_model_FIM)))

# Levene’s test for homoscedasticity across predictors
#cat("\nLevene’s Test for Homoscedasticity (FIM):\n")
#print(leveneTest(resid(mlr_model_FIM) ~ season_data_FIM$t_detrended))

# Plot residuals for visual check
ggplot(data.frame(Season = "FIM", Residuals = resid(mlr_model_FIM)), 
       aes(x = Season, y = Residuals)) +
  geom_boxplot() +
  labs(title = "Residuals for FIM (Homoscedasticity Check)") +
  theme_minimal()

# Visualize Relationships (Partial Effects)
effect_plot_FIM <- plot(allEffects(mlr_model_FIM), main = "Partial Effects of Predictors on Rainfall (FIM)")
print(effect_plot_FIM)

# Save results to CSV
write.csv(mk_results, "trend_results_FIM.csv")
write.csv(anova_result_FIM, "anova_results_FIM.csv")
write.csv(vif_results_FIM, "vif_results_FIM.csv")
write.csv(mlr_results_FIM, "MLR_Results_FIM.csv")

# Step 4: MLR for SWM Season
cat("\n\n--- Analysis for Season: SWM ---\n")

# Subset data for the season
season_data_SWM <- df_detrended %>% filter(season == "SWM")

# MLR model for the season
mlr_model_SWM <- lm(rainfall_detrended ~ t_detrended + RH_detrended + windspeed_detrended, 
                    data = season_data_SWM)

# MLR Summary
cat("\nMLR Summary for SWM:\n")
print(summary(mlr_model_SWM))

# Store MLR results as a data frame
mlr_results_SWM <- as.data.frame(summary(mlr_model_SWM)$coefficients)

# Residual plots (normality and homoscedasticity)
par(mfrow = c(2, 2))
plot(mlr_model_SWM, which = 1:2, main = "Residuals vs Fitted and Q-Q Plot (SWM)")
plot(mlr_model_SWM, which = 3, main = "Scale-Location (SWM)")
hist(resid(mlr_model_SWM), main = "Histogram of Residuals (SWM)", xlab = "Residuals")
par(mfrow = c(1, 1))

# Multicollinearity (VIF should be <5)
cat("\nVariance Inflation Factors (VIF) for SWM:\n")
vif_vals_SWM <- vif(mlr_model_SWM)
print(vif_vals_SWM)
vif_results_SWM <- data.frame(Variable = names(vif_vals_SWM), VIF = vif_vals_SWM)

# Autocorrelation (Durbin-Watson, p<0.05 indicates issue)
cat("\nDurbin-Watson Test for Autocorrelation (SWM):\n")
print(dwtest(mlr_model_SWM))

# ANOVA on MLR Model
anova_result_SWM <- anova(mlr_model_SWM)
cat("\nANOVA Results for SWM:\n")
print(anova_result_SWM)

# Calculate Eta-Squared (η²)
ss_total_SWM <- sum((season_data_SWM$rainfall_detrended - mean(season_data_SWM$rainfall_detrended))^2)
eta_squared_SWM <- anova_result_SWM$"Sum Sq" / ss_total_SWM
anova_result_SWM$Eta_Squared <- eta_squared_SWM
cat("\nANOVA with Eta-Squared for SWM:\n")
print(anova_result_SWM)

# ANOVA Assumption Checks
# Shapiro-Wilk test for normality
cat("\nShapiro-Wilk Test for Residual Normality (SWM):\n")
print(shapiro.test(resid(mlr_model_SWM)))

# Levene’s test for homoscedasticity across predictors
#cat("\nLevene’s Test for Homoscedasticity (SWM):\n")
#print(leveneTest(resid(mlr_model_SWM) ~ season_data_SWM$t_detrended))

# Plot residuals for visual check
ggplot(data.frame(Season = "SWM", Residuals = resid(mlr_model_SWM)), 
       aes(x = Season, y = Residuals)) +
  geom_boxplot() +
  labs(title = "Residuals for SWM (Homoscedasticity Check)") +
  theme_minimal()

# Visualize Relationships (Partial Effects)
effect_plot_SWM <- plot(allEffects(mlr_model_SWM), main = "Partial Effects of Predictors on Rainfall (SWM)")
print(effect_plot_SWM)

# Save results to CSV
write.csv(mk_results, "trend_results_SWM.csv")
write.csv(anova_result_SWM, "anova_results_SWM.csv")
write.csv(vif_results_SWM, "vif_results_SWM.csv")
write.csv(mlr_results_SWM, "MLR_Results_SWM.csv")

# Step 4: MLR for SIM Season
cat("\n\n--- Analysis for Season: SIM ---\n")

# Subset data for the season
season_data_SIM <- df_detrended %>% filter(season == "SIM")

# MLR model for the season
mlr_model_SIM <- lm(rainfall_detrended ~ t_detrended + RH_detrended + windspeed_detrended, 
                    data = season_data_SIM)

# MLR Summary
cat("\nMLR Summary for SIM:\n")
print(summary(mlr_model_SIM))

# Store MLR results as a data frame
mlr_results_SIM <- as.data.frame(summary(mlr_model_SIM)$coefficients)

# Residual plots (normality and homoscedasticity)
par(mfrow = c(2, 2))
plot(mlr_model_SIM, which = 1:2, main = "Residuals vs Fitted and Q-Q Plot (SIM)")
plot(mlr_model_SIM, which = 3, main = "Scale-Location (SIM)")
hist(resid(mlr_model_SIM), main = "Histogram of Residuals (SIM)", xlab = "Residuals")
par(mfrow = c(1, 1))

# Multicollinearity (VIF should be <5)
cat("\nVariance Inflation Factors (VIF) for SIM:\n")
vif_vals_SIM <- vif(mlr_model_SIM)
print(vif_vals_SIM)
vif_results_SIM <- data.frame(Variable = names(vif_vals_SIM), VIF = vif_vals_SIM)

# Autocorrelation (Durbin-Watson, p<0.05 indicates issue)
cat("\nDurbin-Watson Test for Autocorrelation (SIM):\n")
print(dwtest(mlr_model_SIM))

# ANOVA on MLR Model
anova_result_SIM <- anova(mlr_model_SIM)
cat("\nANOVA Results for SIM:\n")
print(anova_result_SIM)

# Calculate Eta-Squared (η²)
ss_total_SIM <- sum((season_data_SIM$rainfall_detrended - mean(season_data_SIM$rainfall_detrended))^2)
eta_squared_SIM <- anova_result_SIM$"Sum Sq" / ss_total_SIM
anova_result_SIM$Eta_Squared <- eta_squared_SIM
cat("\nANOVA with Eta-Squared for SIM:\n")
print(anova_result_SIM)

# ANOVA Assumption Checks
# Shapiro-Wilk test for normality
cat("\nShapiro-Wilk Test for Residual Normality (SIM):\n")
print(shapiro.test(resid(mlr_model_SIM)))

# Levene’s test for homoscedasticity across predictors
cat("\nLevene’s Test for Homoscedasticity (SIM):\n")
print(leveneTest(resid(mlr_model_SIM) ~ season_data_SIM$t_detrended))

# Plot residuals for visual check
ggplot(data.frame(Season = "SIM", Residuals = resid(mlr_model_SIM)), 
       aes(x = Season, y = Residuals)) +
  geom_boxplot() +
  labs(title = "Residuals for SIM (Homoscedasticity Check)") +
  theme_minimal()

# Visualize Relationships (Partial Effects)
effect_plot_SIM <- plot(allEffects(mlr_model_SIM), main = "Partial Effects of Predictors on Rainfall (SIM)")
print(effect_plot_SIM)

# Save results to CSV
write.csv(mk_results, "trend_results_SIM.csv")
write.csv(anova_result_SIM, "anova_results_SIM.csv")
write.csv(vif_results_SIM, "vif_results_SIM.csv")
write.csv(mlr_results_SIM, "MLR_Results_SIM.csv")

# Step 4: MLR for NEM Season
cat("\n\n--- Analysis for Season: NEM ---\n")

# Subset data for the season
season_data_NEM <- df_detrended %>% filter(season == "NEM")

# MLR model for the season
mlr_model_NEM <- lm(rainfall_detrended ~ t_detrended + RH_detrended + windspeed_detrended, 
                    data = season_data_NEM)

# MLR Summary
cat("\nMLR Summary for NEM:\n")
print(summary(mlr_model_NEM))

# Store MLR results as a data frame
mlr_results_NEM <- as.data.frame(summary(mlr_model_NEM)$coefficients)

# Residual plots (normality and homoscedasticity)
par(mfrow = c(2, 2))
plot(mlr_model_NEM, which = 1:2, main = "Residuals vs Fitted and Q-Q Plot (NEM)")
plot(mlr_model_NEM, which = 3, main = "Scale-Location (NEM)")
hist(resid(mlr_model_NEM), main = "Histogram of Residuals (NEM)", xlab = "Residuals")
par(mfrow = c(1, 1))

# Multicollinearity (VIF should be <5)
cat("\nVariance Inflation Factors (VIF) for NEM:\n")
vif_vals_NEM <- vif(mlr_model_NEM)
print(vif_vals_NEM)
vif_results_NEM <- data.frame(Variable = names(vif_vals_NEM), VIF = vif_vals_NEM)

# Autocorrelation (Durbin-Watson, p<0.05 indicates issue)
cat("\nDurbin-Watson Test for Autocorrelation (NEM):\n")
print(dwtest(mlr_model_NEM))

# ANOVA on MLR Model
anova_result_NEM <- anova(mlr_model_NEM)
cat("\nANOVA Results for NEM:\n")
print(anova_result_NEM)

# Calculate Eta-Squared (η²)
ss_total_NEM <- sum((season_data_NEM$rainfall_detrended - mean(season_data_NEM$rainfall_detrended))^2)
eta_squared_NEM <- anova_result_NEM$"Sum Sq" / ss_total_NEM
anova_result_NEM$Eta_Squared <- eta_squared_NEM
cat("\nANOVA with Eta-Squared for NEM:\n")
print(anova_result_NEM)

# ANOVA Assumption Checks
# Shapiro-Wilk test for normality
cat("\nShapiro-Wilk Test for Residual Normality (NEM):\n")
print(shapiro.test(resid(mlr_model_NEM)))

# Levene’s test for homoscedasticity across predictors
cat("\nLevene’s Test for Homoscedasticity (NEM):\n")
print(leveneTest(resid(mlr_model_NEM) ~ season_data_NEM$t_detrended))

# Plot residuals for visual check
ggplot(data.frame(Season = "NEM", Residuals = resid(mlr_model_NEM)), 
       aes(x = Season, y = Residuals)) +
  geom_boxplot() +
  labs(title = "Residuals for NEM (Homoscedasticity Check)") +
  theme_minimal()

# Visualize Relationships (Partial Effects)
effect_plot_NEM <- plot(allEffects(mlr_model_NEM), main = "Partial Effects of Predictors on Rainfall (NEM)")
print(effect_plot_NEM)

# Save results to CSV
write.csv(mk_results, "trend_results_NEM.csv")
write.csv(anova_result_NEM, "anova_results_NEM.csv")
write.csv(vif_results_NEM, "vif_results_NEM.csv")
write.csv(mlr_results_NEM, "MLR_Results_NEM.csv")
