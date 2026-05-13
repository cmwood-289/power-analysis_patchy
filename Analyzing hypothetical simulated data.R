# Full suite of prediction testing and power analysis for hypothetical species profiles
# Mickey Pardo - map385@cornell.edu
## First coded July 2025

library(data.table)
library(ggplot2)

# Getting started ----

# You'll need to compile all your outputs into one object. This could be a single output of one big job (100+-iteration simulation - slow)
# or many smaller ones (10x 20-iteration simulations - much much faster) that are appended.

# read in the simulation output list 
simout <- readRDS("20250624_simulation_results.rds")

# if you choose to do one big job, just name that object 'simout' or change 'simout' in this script to something else

#to calculate power to detect occupancy in a given year
## for each repetition of a given scenario, extract 85% confidence interval for 
#### occupancy estimate in the year of interest
### Compare to true occupancy for that year- if true occupancy falls within 85% CI, 
#### score as 1, otherwise 0
#### calculate proportion of runs scored as 1 for a given scenario- this is the power
### repeat for each scenario
# probably best to pick a season in the middle of the study since those are more accurate if top model did not include time
#could actually skip most of the above steps and simply use the accurate_pred slot in simout

#to calculate power to detect trend in occupancy
### run linear regression on true occupancy values- slope is the true trend
### run linear regression on estimated occupancy values - slope is estimated trend
### calculate 85% CI for the estimated slope using the SE from summary() on the lm()
### if true slope is within 85% CI, score as 1, otherwise 0
### calculate proportion of runs scored as 1 for a given scenario- this is the power
### repeat for each scenario

#How do I select all the list slots that belong to the same scenario?
### there is no slot labeling the scenario
## Need to ensure that all of the following slots match:
# species.type
# lambda
# num.sites
# det.prob
# num.visits

## Storing outputs ----
# For the published analysis, there were 108 different scenarios, and 200 iterations of each scenario 
## create data frame with one column for each sim parameter

# Create a matrix with NA values
scen_matrix <- matrix(NA, nrow = 108, ncol = 5)

# Convert the matrix to a data frame
dfscen <- data.frame(scen_matrix)

#assign column names
colnames(dfscen) <- c("species.type", "lambda", "num.sites", "det.prob", "num.visits")

#set up the possible values for each column
spp.vals <- c("common generalist", "uncommon generalist", "specialist")
lambda.vals <- c(.94, 1, 1.04) #annual growth rate
nsites.vals <- c(50, 100, 200)
detprob.vals <- c(0.2, 0.4)
visits.vals <- c(2, 4)
row.num <- 1
#loop through each of the scenarios and fill out dfscen
for (s in spp.vals){
  for (l in lambda.vals){
    for (n in nsites.vals){
      for (d in detprob.vals){
        for (v in visits.vals){
          dfscen[row.num,] <- c(s,l,n,d,v)
          row.num <- row.num+1
        }
      }
    }
  }
}

#add columns to dfscen to hold RMSE and power values
dfscen$RMSE <- NA
dfscen$power_occ_s1 <- NA
dfscen$power_occ_s5 <- NA
dfscen$power_occ_s10 <- NA
dfscen$power_occ_trend <- NA
dfscen$power_lambda_trend <- NA

## Compile relevant outputs for analysis ----
#Now loop through each scenario in dfscen, subset simout to just the slots
### that match the scenario in question, and calculate power for that scenario
for (i in 1:nrow(dfscen)){
  this_row <- dfscen[i, ]
  
  #subset simout to just the slots that match
  subset_list <- Filter(function(x)
    x$species.type == this_row$species.type &
      x$lambda  == this_row$lambda  &
      x$num.sites  == this_row$num.sites &
      x$det.prob == this_row$det.prob &
      x$num.visits == this_row$num.visits,
    simout)
  #break the loop and display an error message if subset_list is not 100 elements long
  if (length(subset_list) != 100) {
    stop("Error: subset_list is not 100 elements long")
  }
  
  #calculate squared errors for each slot of subset_list and store in vector
  sqerrs <- unlist(lapply(subset_list, function(z) (z$est_occ - z$true_occ)^2))
  #calculate RMSE by first taking mean of sqerrs and then taking sqrt
  dfscen$RMSE[i] <- sqrt(mean(sqerrs))
  ## reason we first average the sqerrs of all the models and THEN take the sqrt
  #### of the whole thing rather than separately calculating the RMSE for each model
  ##### and then averaging those RMSEs is because averaging the RMSE of each model 
  ###### would underestimate variance b/c you'd be performing sqrt operation multiple times
  
  #calculate power to accurately estimate occupancy in 1st season
  ### "accurately estimate" is defined as the 85% CI of the estimate encompassing the true value
  # Extract the 1st value from each accurate_pred
  first_occ_vals <- sapply(subset_list, function(x) x$accurate_pred[1])
  #calculate proportion of values in fifth_occ_vals equal to 1
  dfscen$power_occ_s1[i] <- length(first_occ_vals[first_occ_vals==1])/length(first_occ_vals)
  
  #Now do seasons 5 and 10
  
  # Extract the 5th value from each accurate_pred
  fifth_occ_vals <- sapply(subset_list, function(x) x$accurate_pred[5])
  #calculate proportion of values in fifth_occ_vals equal to 1
  dfscen$power_occ_s5[i] <- length(fifth_occ_vals[fifth_occ_vals==1])/length(fifth_occ_vals)
  
  # Extract the 10th value from each accurate_pred
  tenth_occ_vals <- sapply(subset_list, function(x) x$accurate_pred[10])
  #calculate proportion of values in fifth_occ_vals equal to 1
  dfscen$power_occ_s10[i] <- length(tenth_occ_vals[tenth_occ_vals==1])/length(tenth_occ_vals)
  
  #calculate power to accurately estimate trend in occupancy over 10 years
  ### create vector of values to indicate whether estimated occupancy trend was accurate
  #### (defined as 85% CI encompassing the true trend)
  
  #first clean up subset_list by removing any rows where there are NA, NaN, or Inf
  ### in est_occ
  cleaned_list_occ <- Filter(function(x) all(is.finite(x$est_occ)), subset_list)
  
  accurate_occ_trend <- sapply(cleaned_list_occ, function(x) {
    
    # Predictor: time points (1 through 10)
    time_occ <- 1:10
    
    # 1. Linear regression for estimated occupancy
    est_occ_model <- lm(x$est_occ ~ time_occ)
    est_occ_slope <- coef(est_occ_model)[2]
    est_occ_se <- summary(est_occ_model)$coefficients[2, 2]
    
    # 2. Linear regression for true occupancy
    true_occ_model <- lm(x$true_occ ~ time_occ)
    true_occ_slope <- coef(true_occ_model)[2]
    
    # 3. 85% confidence interval for est_occ_slope
    t_val_occ <- qt(0.925, df = length(time_occ) - 2)  # two-sided, df = 8
    ci_lower_occ <- est_occ_slope - t_val_occ * est_occ_se
    ci_upper_occ <- est_occ_slope + t_val_occ * est_occ_se
    
    # 4. Check if true_slope is within 85% CI
    as.integer(ci_lower_occ <= true_occ_slope & true_occ_slope <= ci_upper_occ)
  })
  #calculate proportion of values in accurate_occ_trend equal to 1
  dfscen$power_occ_trend[i] <- length(accurate_occ_trend[accurate_occ_trend==1])/length(accurate_occ_trend)
  
  #first clean up subset_list by removing any rows where there are NA, NaN, or Inf
  ### in lambda_obs
  cleaned_list_lambda <- Filter(function(x) all(is.finite(x$lambda_obs)), subset_list)
  
  #calculate power to accurately estimate trend in lambda over 10 years
  accurate_lambda_trend <- sapply(cleaned_list_lambda, function(x) {
    # Predictor: time points (1 through 10)
    time_lambda <- 1:9
    
    # 1. Linear regression for observed lambda
    lambda_obs_model <- lm(x$lambda_obs ~ time_lambda)
    lambda_obs_slope <- coef(lambda_obs_model)[2]
    lambda_obs_se <- summary(lambda_obs_model)$coefficients[2, 2]
    
    # 2. Linear regression for real lambda
    lambda_real_model <- lm(x$lambda_real ~ time_lambda)
    lambda_real_slope <- coef(lambda_real_model)[2]
    
    # 3. 85% confidence interval for lambda_obs_slope
    t_val_lambda <- qt(0.925, df = length(time_lambda) - 2)  # two-sided, df = 7
    ci_lower_lambda <- lambda_obs_slope - t_val_lambda * lambda_obs_se
    ci_upper_lambda <- lambda_obs_slope + t_val_lambda * lambda_obs_se
    
    # 4. Check if lambda_real_slope is within 85% CI
    as.integer(ci_lower_lambda <= lambda_real_slope & lambda_real_slope <= ci_upper_lambda)
  })
  #calculate proportion of values in accurate_lambda_trend equal to 1
  dfscen$power_lambda_trend[i] <- length(accurate_lambda_trend[accurate_lambda_trend==1])/length(accurate_lambda_trend)
}

# Analyze your simulation results ----

#convert dfscen to data.table for easier manipulation
dtscen <- as.data.table(dfscen)

#create column in dtscen with unique ID for each survey design
dtscen[,srvdes := paste0("s",as.character(num.sites),"_v",as.character(num.visits))]
#convert relevant columns into ordered factors or numeric format
dtscen[,":=" (srvdes = factor(srvdes, levels=c("s50_v2", "s50_v4", "s100_v2", "s100_v4", "s200_v2", "s200_v4")),
              species.type = factor(species.type, levels=c("common generalist", "uncommon generalist", "specialist")),
              lambda = as.numeric(lambda),
              num.sites = as.numeric(num.sites),
              det.prob = as.numeric(det.prob),
              num.visits = as.numeric(num.visits))]

#save dtscen as .csv file
write.csv(dtscen, "20250820_power_values_each_sim_scenario.csv")

## make some plots to visualize data----

## boxplot of RMSE ~ survey design
boxplot(RMSE ~ srvdes, data=dtscen)

## boxplot of season 1 power ~ survey design
boxplot(power_occ_s1 ~ srvdes, data=dtscen)

## boxplot of season 5 power ~ survey design
boxplot(power_occ_s5 ~ srvdes, data=dtscen)

## boxplot of season 10 power ~ survey design
boxplot(power_occ_s10 ~ srvdes, data=dtscen)

## boxplot of RMSE ~ species.type
boxplot(RMSE ~ species.type, data=dtscen)

## boxplot of season 1 power ~ species.type
boxplot(power_occ_s1 ~ species.type, data=dtscen)

## boxplot of season 5 power ~ species.type
boxplot(power_occ_s5 ~ species.type, data=dtscen)

## boxplot of season 10 power ~ species.type
boxplot(power_occ_s10 ~ species.type, data=dtscen)

## histogram of RMSE to see distribution
hist(dtscen$RMSE)

## gamma regression of RMSE as function of survey design and ecology----
### using gamma because it has no upper bound (just like RMSE) and
#### is often a good fit for right-skewed data (like our RMSE values)
rmse.gamma <- glm(RMSE ~ species.type + det.prob + lambda + scale(num.sites) + scale(num.visits),
                  family = Gamma(link = "log"), data = dtscen)
summary(rmse.gamma)

## run binomial regression models of power as function of survey design and ecology ----
## first create new columns in dtscen to indicate counts of "successes" for each scenario
dtscen[, ":="(k_occ_s1 = 100*power_occ_s1,
              k_occ_s5 = 100*power_occ_s5,
              k_occ_s10 = 100*power_occ_s10)]
## response variable for each model is of the form (cbind(k, n-k))
### this is the grouped binomial format, where you provide glm() with the counts per group
#### rather than the Bernoulli outcome (0 or 1) for each individual trial

#season 1 model
## no lambda in this model b/c pop growth rate shouldn't affect first season
s1binom <- glm(cbind(k_occ_s1, 100-k_occ_s1) ~ species.type + det.prob +
                 scale(num.sites) + scale(num.visits), 
               family = binomial(link = "logit"), data = dtscen)
summary(s1binom)

#season 5 model
s5binom <- glm(cbind(k_occ_s5, 100-k_occ_s5) ~ species.type + det.prob +
                 lambda + scale(num.sites) + scale(num.visits), 
               family = binomial(link = "logit"), data = dtscen)
summary(s5binom)

#season 10 model
s10binom <- glm(cbind(k_occ_s10, 100-k_occ_s10) ~ species.type + det.prob +
                  lambda + scale(num.sites) + scale(num.visits), 
                family = binomial(link = "logit"), data = dtscen)
summary(s10binom)
#coefficient is identical for uncommon generalist and specialist. How?

#now run model for the average power across seasons 1, 5 and 10
## I DON'T KNOW IF IT IS APPROPRIATE TO AVERAGE POWER LIKE THIS - not used in final analyses
#create column for summed k
dtscen[,k_occ_sum := k_occ_s1 + k_occ_s5 + k_occ_s10]

#calculate average power in two different ways just to check they are the same
dtscen[,":=" (avg_pow_1 = k_occ_sum/3,
              avg_pow_2 = 100*((power_occ_s1 + power_occ_s5 + power_occ_s10)/3))]
#check that both methods give the same result
all.equal(dtscen$avg_pow_1, dtscen$avg_pow_2) #yes, they are the same

#binomial model for averaged power
avgpowbinom <- glm(cbind(k_occ_sum, 300-k_occ_sum) ~ species.type + det.prob +
                     lambda + scale(num.sites) + scale(num.visits),
                   family = binomial(link = "logit"), data=dtscen)
summary(avgpowbinom)

#create line graphs showing "power" (i.e., CI coverage) as a function of season, with separate
### lines for each species type and separate panel for each survey design
### Include error bars around each point on each line
# NOTE: The error bars in the plot below represent the mean 85%CI coverage
#### across different scenarios within the same species and survey design
#### (there are multiple scenarios for the same species and survey design b/c of differences in lambda)
## In other words, these error bars don't show variability WITHIN the same scenario
### due to the stochastic nature of the simulation (Monte-Carlo variability)
# Power by season (all 10) and species.type lines

## Build per-scenario, per-season power from simout ----
# (proportion of accurate_pred == 1)
scen_power_rows <- vector("list", nrow(dfscen))

for (i in 1:nrow(dfscen)) {
  #get the ith row of dfscen (ith unique scenario)
  this_row <- dfscen[i, ]
  
  #create subset of simout that just contains the 100 slots matching the ith unique scenario (ith row in dfscen)
  subset_list <- Filter(function(x)
    x$species.type == this_row$species.type &
      x$lambda      == this_row$lambda      &
      x$num.sites   == this_row$num.sites   &
      x$det.prob    == this_row$det.prob    &
      x$num.visits  == this_row$num.visits,
    simout)
  
  #break the loop if subset_list is not 100 elements long
  if (length(subset_list) != 100) stop("Error: subset_list is not 100 elements long")
  
  #create 100 x 10 matrix of 0/1 accuracy for seasons 1 to 10 across iterations
  ## each row is an iteration
  ## each column is a season
  acc_mat <- do.call(rbind, lapply(subset_list, function(x) as.numeric(x$accurate_pred)))
  # scenario-level power per season (mean accuracy over 100 iters)
  ### calculating the mean of a column in acc_mat is the same as dividing the 
  ### number of 1s by the total number of rows (i.e., power)
  ### In other words, the column mean is the power for a given season
  pow_vec <- colMeans(acc_mat, na.rm = TRUE)
  
  scen_power_rows[[i]] <- data.frame(
    species.type = as.character(this_row$species.type),
    lambda       = as.numeric(this_row$lambda),
    num.sites    = as.numeric(this_row$num.sites),
    det.prob     = as.numeric(this_row$det.prob),
    num.visits   = as.numeric(this_row$num.visits),
    season       = 1:10,
    power        = pow_vec #10 power values, one for each season
  )
}

#convert scen_power_rows from df to data.table
scen_power <- rbindlist(scen_power_rows)

#create survey-design label
scen_power[, srvdes := paste0("s", num.sites, "_v", num.visits)]

#turn srvdes and species.type into ordered factors with more reader-friendly labels
scen_power[, ":=" (srvdes = factor(srvdes,
                                   levels = c("s50_v2","s50_v4","s100_v2","s100_v4","s200_v2","s200_v4"),
                                   labels = c("50 sites, 2 visits", "50 sites, 4 visits",
                                              "100 sites, 2 visits", "100 sites, 4 visits",
                                              "200 sites, 2 visits", "200 sites, 4 visits")),
                   species.type = factor(species.type,
                                         levels = c("common generalist","uncommon generalist","specialist")
                   ))]

# Summarize by panel (srvdes), species, and season with 85% CI across scenarios
summary_by_panel <- scen_power[
  , {
    m  <- mean(power, na.rm = TRUE)
    s  <- sd(power,   na.rm = TRUE)
    n  <- sum(!is.na(power))
    tC <- qt(0.925, df = max(n - 1, 1))     # 85% CI
    se <- if (n > 0) s / sqrt(n) else NA_real_
    lo <- pmax(0, m - tC * se)
    hi <- pmin(1, m + tC * se)
    .(mean_power = m, ci_lower = lo, ci_upper = hi, n_scen = n)
  },
  by = .(srvdes, species.type, season)
]

# Plotting results ----
## Faceted plot: 6 panels (srvdes: power), 3 lines (species), with 85% CI error bars ----
ggplot(summary_by_panel,
       aes(x = season, y = mean_power, color = species.type, group = species.type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.18, alpha = 0.7) +
  facet_wrap(~ srvdes, ncol = 3) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(limits = c(0.5,1), labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Season",
       y = "Power (85% CI)",
       color = "Species",
       title = "Occupancy power by season",
       subtitle = "Mean 85% CI coverage across scenarios") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")


## Faceted plot: 6 panels (srvdes: RMSE), 3 lines (species), with 85% CI error bars ----

# create faceted line graph figure like above except instead of "power" have RMSE
#### instead of calculating the RMSE separately for each scenario at a given pt
#### and averaging multiple RMSEs, calculate a joint MSE for all iterations of
#### all scenarios that correspond to a given pt (species, season, survey design)
#### and then take the sqrt() of that aggregate MSE
#### To create error bars that represent the variability ACROSS scenarios within
#### the same species, season, and survey design, calculate the mean and SE
#### of the MSEs (not RMSE) for each scenario

# Build per-scenario, per-season MSE 
## (mean squared error over all 100 iterations of given scenario in given season)

rows <- vector("list", nrow(dfscen))
for (i in 1:nrow(dfscen)) {
  this_row <- dfscen[i, ]
  #subset simout to just the iterations of the ith scenario
  subset_list <- Filter(function(x)
    x$species.type == this_row$species.type &
      x$lambda      == this_row$lambda      &
      x$num.sites   == this_row$num.sites   &
      x$det.prob    == this_row$det.prob    &
      x$num.visits  == this_row$num.visits,
    simout)
  
  if (length(subset_list) != 100) stop("Error: subset_list is not 100 elements long")
  
  # 100 x 10 matrix of squared errors for seasons 1..10
  sqe_mat <- do.call(rbind, lapply(subset_list, function(x) {
    e <- (x$est_occ - x$true_occ)^2
    if (length(e) != 10) stop("Each est_occ/true_occ must be length 10")
    if (!all(is.finite(e))) e[!is.finite(e)] <- NA_real_
    e
  }))
  
  # Scenario-level MSE by season (mean over iterations)
  mse_by_season <- colMeans(sqe_mat, na.rm = TRUE)
  
  rows[[i]] <- data.frame(
    species.type = as.character(this_row$species.type),
    num.sites    = as.numeric(this_row$num.sites),
    num.visits   = as.numeric(this_row$num.visits),
    lambda       = as.numeric(this_row$lambda),
    det.prob     = as.numeric(this_row$det.prob),
    season       = 1:10,
    mse_scen     = mse_by_season
  )
}
#convert rows from a df into a data.table called mse_scen
mse_scen <- rbindlist(rows)

#create column in mse_scen with survey design ID
mse_scen[, srvdes := paste0("s", num.sites, "_v", num.visits)]
#convert srvdes and species.type into factors with appropriate order and labels
mse_scen[, srvdes := factor(srvdes, levels = c("s50_v2","s50_v4","s100_v2","s100_v4","s200_v2","s200_v4"),
                            labels = c("50 sites, 2 visits", "50 sites, 4 visits",
                                       "100 sites, 2 visits", "100 sites, 4 visits",
                                       "200 sites, 2 visits", "200 sites, 4 visits"))]
mse_scen[, species.type := factor(species.type,
                                  levels = c("common generalist","uncommon generalist","specialist"))]

# Aggregate across scenarios within (species, srvdes, season)
# Point estimate: RMSE = sqrt( mean_scenario(MSE_scen) )
# CI: t-interval on mean MSE across scenarios, then sqrt the endpoints

sum_rmse <- mse_scen[
  , .(
    MSE_bar = mean(mse_scen, na.rm = TRUE),
    SD_MSE  = sd(mse_scen,   na.rm = TRUE),
    S       = sum(!is.na(mse_scen))    # number of scenarios in this panel/species/season
  ),
  by = .(srvdes, species.type, season)
]

sum_rmse[, `:=`(
  tcrit  = qt(0.925, df = pmax(S - 1, 1)),            # 85% two-sided
  SE_MSE = ifelse(S > 0, SD_MSE / sqrt(S), NA_real_)
)]

# CI on MSE, then transform to RMSE
sum_rmse[, `:=`(
  MSE_lo = pmax(0, MSE_bar - tcrit * SE_MSE),
  MSE_hi = pmax(0, MSE_bar + tcrit * SE_MSE),
  RMSE   = sqrt(MSE_bar)
)]
#add upper and lower bounds to RMSE
sum_rmse[, ":="(RMSE_lo = sqrt(MSE_lo),
                RMSE_hi = sqrt(MSE_hi)
)]

## Plot: 6 panels (srvdes), lines = species, 1 pt/season on each line, 85% CI bars around pts----
ggplot(sum_rmse,
       aes(x = season, y = RMSE, color = species.type, group = species.type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = RMSE_lo, ymax = RMSE_hi), width = 0.18, alpha = 0.7) +
  facet_wrap(~ srvdes, ncol = 3) +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "Season",
       y = "RMSE",
       color = "Species",
       title = "RMSE by season",
       subtitle = "Point = sqrt(mean MSE across scenarios); bars = 85%CI across scenarios") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")



# calculate how often top model was correct ----
## for generalists, correct top model is ~time
## for specialist, correct top model is ~time + habitat

#create empty dataframe with one row for each iteration in simout
### one col for species.type
### one col for top_mod
topmod.dt <- data.table(species.type = character(length(simout)), 
                        top_mod_psi = character(length(simout)),
                        top_mod_gamma = character(length(simout)),
                        top_mod_epsilon = character(length(simout)))

#loop through simout list and fill topmod.dt
for(i in 1:length(simout)){
  topmod.dt$species.type[i] <- simout[[i]]$species.type
  #get top model for occupancy and convert from formula to character format
  topmod.dt$top_mod_psi[i] <- deparse(simout[[i]]$top_mod$psiformula)
  #get top model for colonization and convert from formula to character format
  topmod.dt$top_mod_gamma[i] <- deparse(simout[[i]]$top_mod$gammaformula)
  #get top model for extinction and convert from formula to character format
  topmod.dt$top_mod_epsilon[i] <- deparse(simout[[i]]$top_mod$epsilonformula)
}






