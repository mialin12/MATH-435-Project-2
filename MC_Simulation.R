library(lme4)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(MASS)
library(kableExtra)
library(knitr)
library(patchwork)
set.seed(123)

## Labor Demand
combdata <- read.csv('combdata.csv')

combdata <- combdata %>%
  mutate(
    emp_change = total_emps - lag(total_emps, n = 1),
    log_emp = log(total_emps),
    emp_change_pct = log_emp - lag(log_emp, n = 1),
    log_min_wage = log(min_wage),
    log_cpi = log(CPI),
    log_cpidelay = log(CPI_delay),
    log_labforce = log(labor_force),
    log_workpop = log(workpop),
    d_min_wage = log_min_wage - lag(log_min_wage),
    d_labforce = log_labforce - lag(log_labforce)
  ) %>%
  drop_na()

model_ldemand <- lm(emp_change_pct ~ d_min_wage + interestrate + d_labforce, data = combdata)

coefs_ldemand <- coef(model_ldemand)
vcov_ldemand <- vcov(model_ldemand)
sigma_ldemand <- sigma(model_ldemand)


# Checking for multicollinearity
# Check multicollinearity
library(car)
vif(model_ldemand) %>% kable(col.names = c('Predictor', 'VIF Score')) %>% kable_styling(bootstrap_options ='bordered')


# Check residual normality
plot(model_ldemand)

# Print summaries to verify coefficients look sensible
summary(model_ldemand)



# Simulating CUMMULATIVE employment change over a certain amount of time
# --- Simulation inputs ---
n_sims          <- 10000    # Number of iterations
interest_rate   <- 4        # Interest Rate
t_horizon       <- 120      # Number of months to simulate over (120 -> 1 year)
felten_score    <- 0.45     # Felten index for AI component

simulate_path <- function(min_wage_growth, t_horizon, ai_decline, scenario_lambda, n_sims, ...) {
  
  total_exposure <- sum(felten_score * exp(scenario_lambda *1:t_horizon))
  beta_AI_recalibrated <- ai_decline / total_exposure
  
  # Draw coefficients once per simulation run — stays fixed over the path
  B_demand     <- mvrnorm(n_sims, mu = coefs_ldemand, Sigma = vcov_ldemand)
  beta_AI_monthly <- rnorm(n_sims, mean = beta_AI_recalibrated, sd = 0.0005)
  lambda_draw  <- rnorm(n_sims, mean = scenario_lambda, sd = scenario_lambda * 0.2)

  
  # Accumulate employment changes over the horizon
  cumulative_emp_change <- rep(0, n_sims)
  
  for (t in 1:t_horizon) {
    
    # AI exposure grows each month
    ai_exposure_t <- felten_score * exp(lambda_draw * t)
    
    # Labor force change drawn from historical distribution each month
    d_labforce_t <- rnorm(n_sims,
                          mean = mean(combdata$d_labforce, na.rm = TRUE),
                          sd   = sd(combdata$d_labforce,   na.rm = TRUE))
    
    # Monthly employment change at time t
    emp_change_t <- B_demand[,"(Intercept)"] +
      B_demand[,"d_min_wage"]   * min_wage_growth +
      B_demand[,"interestrate"] * interest_rate +
      B_demand[,"d_labforce"]   * d_labforce_t +
      (beta_AI_monthly * ai_exposure_t) +
      rnorm(n_sims, 0, sigma_ldemand)
    
    # Add this month's change to the running total
    cumulative_emp_change <- cumulative_emp_change + emp_change_t
    
    if (t %in% c(1, 12, 60, 120)) {
      cat("t =", t, "| cumulative mean:", mean(cumulative_emp_change*100), "\n")
    }
  }
  
  return(cumulative_emp_change*100)
}


slow_path <- simulate_path(
  min_wage_growth = log(1.05)/120, 
  t_horizon       = 120,  
  ai_decline      = -0.03,
  scenario_lambda = 0.03/12, 
  n_sims          = 10000, 
  coefs_lf        = coefs_lf,       vcov_lf       = vcov_lf,       sigma_lf      = sigma_lf,
  coefs_ldemand   = coefs_ldemand,  vcov_ldemand  = vcov_ldemand,  sigma_ldemand = sigma_ldemand
)

cat("Mean:   ", mean(slow_path), "\n")
cat("SD:     ", sd(slow_path), "\n")
cat("Min:    ", min(slow_path), "\n")
cat("Max:    ", max(slow_path), "\n")
cat("% negative:", mean(slow_path < 0), "\n")

slow_plot <- ggplot() + geom_histogram(aes(x = slow_path), binwidth = 1, fill = '#69b3a2', colour = 'white', alpha = 0.8) + 
  labs( subtitle = 'Slow AI growth rate', x = 'Cumulative % Change', y = "Count") + 
  geom_vline(xintercept = mean(slow_path),
             color = "red",
             linewidth = 1,
             linetype = "dashed") + theme_minimal(base_size = 10) + annotate("text",
                                                                             x = mean(slow_path),
                                                                             y = Inf,
                                                                             label = paste0("Mean = ", round(mean(slow_path), 2), "%"),
                                                                             vjust = 1.5,
                                                                             hjust = -0.2,
                                                                             color = "red",
                                                                             size = 3) 




slow_plot + labs(title  = 'Employment Change over 10 years')

# No AI, current min wage assumption
path_policy <- simulate_path(
  min_wage_growth = log(1.05)/120,
  t_horizon       = 120,
  ai_decline      = 0,
  scenario_lambda = 0.03/12,
  n_sims          = 10000
)
mean(path_policy)
t.test(path_policy, conf.level = 0.95)
cat("% negative:", mean(path_policy < 0), "\n")

none_plot <- ggplot() + geom_histogram(aes(x = path_policy), binwidth = 1,  fill = '#69b3a2', colour = 'white', alpha = 0.8) + 
  labs( subtitle = 'No AI growth',  x = 'Cumulative % Change', y = "Count") + 
  geom_vline(xintercept = mean(path_policy),
             color = "red",
             linewidth = 1,
             linetype = "dashed") + theme_minimal(base_size = 10) + annotate("text",
                                                                             x = mean(path_policy),
                                                                             y = Inf,
                                                                             label = paste0("Mean = ", round(mean(path_policy), 2), "%"),
                                                                             vjust = 1.5,
                                                                             hjust = -0.2,
                                                                             color = "red",
                                                                             size = 3) 

none_plot + labs(title = 'Employment Change over 10 years')


moderate_path <- simulate_path(
  min_wage_growth = log(1.05)/120, 
  t_horizon       = 120,  
  scenario_lambda = 0.1/12,
  ai_decline     = -0.065,
  n_sims          = 10000, 
  coefs_cpi       = coefs_cpi,      vcov_mat_cpi  = vcov_mat_cpi,  sigma_cpi     = sigma_cpi,
  coefs_lf        = coefs_lf,       vcov_lf       = vcov_lf,       sigma_lf      = sigma_lf,
  coefs_ldemand   = coefs_ldemand,  vcov_ldemand  = vcov_ldemand,  sigma_ldemand = sigma_ldemand
)
cat("Mean:   ", mean(moderate_path), "\n")
cat("% negative:", mean(moderate_path < 0), "\n")
t.test(moderate_path, conf.level = 0.95)
mod_plot <- ggplot() + geom_histogram(aes(x = moderate_path), binwidth = 1,  fill = '#69b3a2', colour = 'white', alpha = 0.8) + 
  labs( subtitle = 'Moderate AI growth rate',  x = 'Cumulative % Change', y = "Count") + 
  geom_vline(xintercept = mean(moderate_path),
             color = "red",
             linewidth = 1,
             linetype = "dashed") + theme_minimal(base_size = 10) + annotate("text",
                                                                             x = mean(moderate_path),
                                                                             y = Inf,
                                                                             label = paste0("Mean = ", round(mean(moderate_path), 2), "%"),
                                                                             vjust = 1.5,
                                                                             hjust = -0.2,
                                                                             color = "red",
                                                                             size = 3) 
mod_plot + labs(title = 'Employment Change over 10 years')


fast_path <- simulate_path(
  min_wage_growth = log(1.05)/120, 
  t_horizon       = 120,  
  scenario_lambda = 0.18/12, 
  ai_decline      = -0.1,
  n_sims          = 10000, 
  coefs_cpi       = coefs_cpi,      vcov_mat_cpi  = vcov_mat_cpi,  sigma_cpi     = sigma_cpi,
  coefs_lf        = coefs_lf,       vcov_lf       = vcov_lf,       sigma_lf      = sigma_lf,
  coefs_ldemand   = coefs_ldemand,  vcov_ldemand  = vcov_ldemand,  sigma_ldemand = sigma_ldemand
)
t.test(fast_path, conf.level = 0.95)
cat("% negative:", mean(fast_path < 0), "\n")

cat("Mean:   ", mean(fast_path), "\n")
fast_plot <-  ggplot() + geom_histogram(aes(x = fast_path), binwidth = 1, fill = '#69b3a2', colour = 'white', alpha = 0.8) + 
  labs( subtitle = 'Fast AI growth rate', x = 'Cumulative % Change', y = "Count") + 
    geom_vline(xintercept = mean(fast_path),
               color = "red",
               linewidth = 1,
               linetype = "dashed") + theme_minimal(base_size = 10) + annotate("text",
                                                                               x = mean(fast_path),
                                                                               y = Inf,
                                                                               label = paste0("Mean = ", round(mean(fast_path), 2), "%"),
                                                                               vjust = 1.5,
                                                                               hjust = -0.2,
                                                                               color = "red",
                                                                               size = 3) 
fast_plot + labs(title =  'Employment Change over 10 years')


ci_vals_none <- quantile(path_policy, c(0.025, 0.975))
ci_vals_slow <- quantile(slow_path, c(0.025, 0.975))
ci_vals_mod <- quantile(moderate_path, c(0.025, 0.975))
ci_vals_fast <- quantile(fast_path, c(0.025, 0.975))

results_ai <- data.frame(AI_Scenario = c('No AI', 'Slow AI', 'Moderate AI', 'Fast AI'), 
                         Mean = c(mean(path_policy), mean(slow_path), mean(moderate_path), mean(fast_path)),
                        CI = c(sprintf("[%.2f, %.2f]", ci_vals_none[1], ci_vals_none[2]), sprintf("[%.2f, %.2f]", ci_vals_slow[1], ci_vals_slow[2]), sprintf("[%.2f, %.2f]", ci_vals_mod[1], ci_vals_mod[2]), sprintf("[%.2f, %.2f]", ci_vals_fast[1], ci_vals_fast[2])),
                          ProbNeg = c(mean(path_policy < 0), mean(slow_path < 0), mean(moderate_path < 0), mean(fast_path < 0)))                                                                                                   
results_ai %>% kable(col.names = c('AI Scenario', 'Mean Cumulative % Change', '95% Confidence Interval', 'Proportion Negative')) %>% kable_styling(bootstrap_options = 'bordered')

(none_plot | slow_plot) / 
  (mod_plot | fast_plot)

setup_df <- data.frame(Input_Variable = c('n_sims', 't_horizon', 'scenario_lambda', 'ai_decline', 'min_wage_growth', 'felten_score'), Description = c('Number of Iterations', 
                                                                                                                                      'Time Iterated (in months)',
                                                                                                                                      'AI Growth Scenario',
                                                                                                                                      '% Change in Emp from AI',
                                                                                                                                      'Minimum Wage Growth',
                                                                                                                                      'Felten Score (AI Exposure)'), 
                       Value = c(n_sims, t_horizon, 'Chosen', 'Chosen', 'Chosen', 0.45))
setup_df %>% kable(col.names = c('Input Variable', 'Description', 'Value')) %>% kable_styling(bootstrap_options = 'bordered')                                                                                                               


scenario_lambda_df <- data.frame(Lambda = c(0.03, 0.1, 0.18), Scenario = c('Slow AI Growth', 'Moderate AI Growth', 'Fast AI Growth'))
ai_decline_df <- data.frame(Decline = c(-3, -6.5, -10), Scenario = c('Low Impact', 'Moderate Impact', 'Strong Impact'))
min_wage_df <- data.frame(Min_growth = c(0, 33, 100), Scenario = c('No Increase', 'Small Increase', 'Large Increase'), Minimum_Wage = c(7.25, 10, 14.50))

scenario_lambda_df %>% kable(col.names = c('Lambda Input', 'Scenario')) %>% kable_styling(bootstrap_options = 'bordered')
ai_decline_df %>% kable(col.names = c('% Change in Employment over 10 years from AI', 
                                      'Scenario')) %>% kable_styling(bootstrap_options = 'bordered')
min_wage_df %>% kable(col.names = c('% Minimum Wage Growth', 'Scenario', 'Minimum Wage ($)')) %>% kable_styling(bootstrap_options = c('bordered'))
### Minimum Wage Scenarios

  # Hold AI constant at moderate scenario
  moderate_lambda <- 0.10/12
  moderate_ai_decline <- -0.065
  
  # Run three minimum wage scenarios
  mw_zero <- simulate_path(
    min_wage_growth = 0,
    t_horizon       = 120,
    ai_decline      = moderate_ai_decline,
    scenario_lambda = moderate_lambda,
    n_sims          = 10000
  )
  
  mw_moderate <- simulate_path(
    min_wage_growth = log(1.33)/120, # 33% increase in min. wage
    t_horizon       = 120,
    ai_decline      = moderate_ai_decline,
    scenario_lambda = moderate_lambda,
    n_sims          = 10000
  )
  
  mw_high <- simulate_path(
    min_wage_growth = log(2)/120, # 100% increase in min. wage
    t_horizon       = 120,
    ai_decline      = moderate_ai_decline,
    scenario_lambda = moderate_lambda,
    n_sims          = 10000
  )
  
  
  results_df <- data.frame(
    emp_change = c(mw_zero, mw_moderate, mw_high),
    scenario   = rep(c("No growth (0%)", 
                       "Moderate Increase (33%)", 
                       "High Increase (100%)"), each = 10000)
  )
  
  ggplot(results_df, aes(x = emp_change, fill = scenario)) +
    geom_density(alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = mean(mw_zero),     color = "#619CFF", linewidth = 0.8) +
    geom_vline(xintercept = mean(mw_moderate), color = "#00BA38", linewidth = 0.8) +
    geom_vline(xintercept = mean(mw_high),     color = "#F8766D", linewidth = 0.8) +
    labs(
      title    = "Cumulative employment change under different minimum wage policies",
      subtitle = "Holding AI growth constant at moderate scenario",
      x        = "Cumulative % change in employment over 10 years",
      y        = "Density",
      fill     = "Minimum wage scenario"
    ) +
    theme_minimal()
  
  no_lambda <- 0
  no_ai_decline <- 0
  
  mw_no_zero <- simulate_path(
    min_wage_growth = 0,
    t_horizon       = 120,
    ai_decline      = no_ai_decline,
    scenario_lambda = no_lambda,
    n_sims          = 10000
  )
  
  mw_no_moderate <- simulate_path(
    min_wage_growth = log(1.33)/120, # 33% increase in min. wage
    t_horizon       = 120,
    ai_decline      = no_ai_decline,
    scenario_lambda = no_lambda,
    n_sims          = 10000
  )
  
  mw_no_high <- simulate_path(
    min_wage_growth = log(2)/120, # 100% increase in min. wage
    t_horizon       = 120,
    ai_decline      = no_ai_decline,
    scenario_lambda = no_lambda,
    n_sims          = 10000
  )
  
  results_no_df <- data.frame(
    emp_change = c(mw_no_zero, mw_no_moderate, mw_no_high),
    scenario   = rep(c("No growth (0%)", 
                       "Moderate Increase (33%)", 
                       "High Increase (100%)"), each = 10000)
  )
  ggplot(results_no_df, aes(x = emp_change, fill = scenario)) +
    geom_density(alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = mean(mw_no_zero),     color = "#619CFF", linewidth = 0.8) +
    geom_vline(xintercept = mean(mw_no_moderate), color = "#00BA38", linewidth = 0.8) +
    geom_vline(xintercept = mean(mw_no_high),     color = "#F8766D", linewidth = 0.8) +
    labs(
      title    = "Cumulative employment change under different minimum wage policies",
      subtitle = "No AI Growth",
      x        = "Cumulative % change in employment over 10 years",
      y        = "Density",
      fill     = "Minimum wage scenario"
    ) +
    theme_minimal()
  