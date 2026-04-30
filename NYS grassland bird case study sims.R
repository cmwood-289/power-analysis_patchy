# Overview ----
# New York State grassland birds case study #

# Code: Connor Wood & Mickey Pardo
# Data and fieldwork: Timothy Boycott, Trifosa Simamora, Steve Grodsky

library(unmarked)
# function to calculate if true occupancy in each year was within 95% CI for predicted occupancy
good_pred <- function(y){ifelse(true_occ[y] >= bootres$pred_mean[y] - mar85[y] & 
                                  true_occ[y] <= bootres$pred_mean[y] + mar85[y], 1, 0)}
# returns 1 if true occ is within 85% CI
# returns 0 if true occ is not withing 85% CI


# Part 1: running the simulation ----

## Set parameters for simulation ----
{
  n.its = 100             # Number of iterations 
  yrs = 10                # Number of survey years
  territories = 800       # Number of potential bird territories in the study area
  
  sites = c(50, 100, 200) # Number of territories that are sampled
  # det = 0.165 for HOLA           # Secondary sampling period detection probabilities
  # det = 0.481 for GRSP           # Secondary sampling period detection probabilities
  visits = c(3, 5)        # Number of secondary sampling periods
  
  # occ = 0.804 for HOLA         # Starting occupancy value in the preferred habitat
  # occ = 0.161 for GRSP         # Starting occupancy value in the preferred habitat
  # occ_alt = 0.137    # occupancy rate in non-preferred habitat 
  # occ_alt = 0 for GRSP    # occupancy rate in non-preferred habitat 
  prop_habitat = (1-.33)      # Proportion of total territories that are "preferred". In this case, it likes upland better than riparian (which is 0.37)
  
  spp_type = c("horned lark")  # give qualitative names to your species types. These are NOT used in the loop
  
  lambda = c(.93, 1, 1.07)  # annual growth rate
  
  lambda.var = 0.015        # annual +/- variation in lambda
  phi.min = 0.85              # minimum site persistence
  phi.var = 0.03           # annual +/- variation in phi
  
  scenarios=length(sites)*length(det)*length(visits)*length(occ)*length(lambda)
  numboots = 75           # for derived estimates of occupancy; use ~5 to test the code but 50-100 when running for real
  
  #create empty list to store results of each run
  output.list <- vector(mode="list",
                        length=(scenarios*n.its))
}

## Main Loop ----
for(i in 1:n.its){
  if(phi.min+phi.var > min(lambda)-lambda.var ){
    print("Error: Maximum phi may be too high to impose minimum lambda")
    break
  }
  if(length(occ)!=length(spp_type)){
    print("Error: occupancy type and rate mismatch")
    print("Must have the same number of occupancy types as rates")
    break
  }
  if(i==1){
    start.time <- Sys.time()
  }
  for (o in 1:length(occ)) {
    for (l in 1:length(lambda)) {
      
      # calculates the number of territories with specialist's habitat
      hab.terr=round(territories*prop_habitat,0)
      # outside the 'specialist' loop because it's used to create the habitat covariate
      
      # Create an empty matrix for the "habitat" sites
      true.occ.hab <- matrix(NA, nrow=hab.terr, ncol=yrs)
      
      true.occ.hab[,1] <- as.numeric(runif(hab.terr) < occ[o]) # sets initial site-by-site occupancy
      occ.sim = c(mean(true.occ.hab[,1]), rep(NA,yrs-1)) # stores population-level occupancy
      
      for(y in 2:yrs){
        psi.t=runif(n=1, min=lambda[l]-lambda.var, max=lambda[l]+lambda.var)*
          occ.sim[y-1]
        if(psi.t>1){
          occ.sim[y]=1
          true.occ.hab[,y]=1
        } 
        else {
          occ.sim[y]=psi.t
          occ.t=round(psi.t*hab.terr,0) # number of territories occupied in t
          phi.t=round(sum(true.occ.hab[,y-1]==1)*
                        runif(n=1, min=phi.min-phi.var, max=phi.min+phi.var) ,0) # minimum number of territories that persist from t-1
          
          persist.t = sample(which(true.occ.hab[,y-1]==1), 
                             phi.t,
                             replace = F) # randomly selects which sites persist between years
          true.occ.hab[persist.t,y] <- 1
          remaining=c( rep(1, occ.t-phi.t), # sites occupied in addition to the minimum phi
                       rep(0, hab.terr-occ.t)) 
          true.occ.hab[-persist.t,y] <- sample(remaining,
                                               length(remaining),
                                               replace=F)
        }
      }
      
      nonhab.terr=territories-hab.terr
      if(occ_alt==0){
        true.occ.all=rbind(true.occ.hab,
                           matrix(0, nrow=nonhab.terr, ncol=yrs))
      } else {
        # Create an empty matrix for the "non-habitat" sites
        true.occ.nonhab <- matrix(NA, nrow=nonhab.terr, ncol=yrs)
        true.occ.nonhab[,1] <- as.numeric(runif(nonhab.terr) < occ[o]*occ_alt) # sets initial site-by-site occupancy
        occ.sim = c(mean(true.occ.nonhab[,1]), rep(NA,yrs-1)) # stores population-level occupancy
        
        for(y in 2:yrs){
          psi.t=runif(n=1, min=lambda[l]-lambda.var, max=lambda[l]+lambda.var)*
            occ.sim[y-1]
          if(psi.t>1){
            occ.sim[y]=1
            true.occ.nonhab[,y]=1
          } 
          else {
            occ.sim[y]=psi.t
            occ.t=round(psi.t*nonhab.terr,0) # number of territories occupied in t
            phi.t=round(sum(true.occ.nonhab[,y-1]==1)*
                          runif(n=1, min=phi.min-phi.var, max=phi.min+phi.var) ,0) # minimum number of territories that persist from t-1
            
            persist.t = sample(which(true.occ.nonhab[,y-1]==1), 
                               phi.t,
                               replace = F) # randomly selects which sites persist between years
            true.occ.nonhab[persist.t,y] <- 1
            remaining=c( rep(1, occ.t-phi.t), # sites occupied in addition to the minimum phi
                         rep(0, nonhab.terr-occ.t)) 
            true.occ.nonhab[-persist.t,y] <- sample(remaining,
                                                    length(remaining),
                                                    replace=F)
          }
        }
        true.occ.all=rbind(true.occ.hab,
                           true.occ.nonhab)
        
      }
      
      ## Subsample the true occupancy data according to different study designs ##
      # Number of survey sites
      for(s in 1:length(sites)){
        
        keepsites <- seq(1, territories, by=territories/sites[s])
        
        # Different detection probabilities
        for(d in 1:length(det)) {
          
          # Number of visits
          for(v in 1:length(visits)) {
            visit.seq = seq(1,visits[v]*yrs,1) # reference positions for encounter history construction
            
            #subset true occupancy matrix down to just the rows of interest
            true_occ_sub <- true.occ.all[keepsites,]
            
            ## create the encounter history
            encounter_history <- matrix(0, nrow = sites[s], ncol = (visits[v]*yrs))
            
            for(r in 1:nrow(encounter_history)){
              for(c in 1:yrs){
                if(true_occ_sub[r,c] == 0){
                  next
                } else {
                  encounter_history[r,seq(visit.seq[visits[v]*c-(visits[v]-1)], 
                                          visit.seq[visits[v]*c])] <- as.numeric(runif(visits[v]) < det[d])
                }
              }
            }
            
            ## create covariates
            
            #create data.frame where each row is a site, each col is a year, and all the cells in a given column are the year number
            year_site_covs <- data.frame(matrix(rep(1:yrs, each=sites[s]), nrow=sites[s], ncol=yrs, byrow=F))
            
            #create data.frame with single column equal in length to the number of sites; 
            # 1 = suitable habitat for the specialist, 0 = unsuitable habitat
            habitat <- data.frame(habitat = c(rep(1, sum(keepsites <= hab.terr)), 
                                              rep(0, nrow(encounter_history)-sum(keepsites <= hab.terr))) )
            
            ## create unmarked frame object
            um.eh <- unmarkedMultFrame(y=encounter_history, 
                                       numPrimary = yrs, 
                                       siteCovs = habitat,
                                       yearlySiteCovs = 
                                         list(time = year_site_covs)) #yrs is defined earlier in the script
            
            #fit occupancy models
            #null model
            mod1 <- colext(psiformula= ~1, 
                           gammaformula= ~1, 
                           epsilonformula= ~1, 
                           pformula= ~1, 
                           data=um.eh) 
            #time-only model
            mod2 <- colext(psiformula= ~1,
                           gammaformula= ~time,
                           epsilonformula= ~time,
                           pformula= ~1,
                           data=um.eh)
            #habitat-only model
            mod3 <- colext(psiformula= ~habitat,
                           gammaformula= ~habitat,
                           epsilonformula= ~habitat,
                           pformula= ~1,
                           data=um.eh)
            #habitat+time model
            mod4 <- colext(psiformula= ~habitat, 
                           gammaformula= ~habitat + time,
                           epsilonformula= ~habitat + time,
                           pformula= ~1, 
                           data=um.eh)
            
            # identify the top model based on AIC
            modsel_res <- modSel(fitList(m1=mod1, m2=mod2, m3=mod3, m4=mod4))
            ## modsel_res@Full$model[1] returns a character string ("m1", "m2", etc.) indicating which model
            #### is in the first row of the modsel_res@Full dataframe. The model in the first row should
            #### be the best supported model because the dataframe is sorted with the lowest AIC model in the top row
            
            # save the top model
            topmod <- if(modsel_res@Full$model[1]=="m1"){mod1} else 
            {if(modsel_res@Full$model[1]=="m2"){mod2} else 
            {if(modsel_res@Full$model[1]=="m3"){mod3} else
            {if(modsel_res@Full$model[1]=="m4"){mod4}}}}
            
            topmod_name <- if(modsel_res@Full$model[1]=="m1"){"mod1"} else 
            {if(modsel_res@Full$model[1]=="m2"){"mod2"} else
            {if(modsel_res@Full$model[1]=="m3"){"mod3"} else
            {if(modsel_res@Full$model[1]=="m4"){"mod4"}}}}
            
            #bootstrapping the top model to calculate SE around predicted mean occupancy in each year
            bootmod <- nonparboot(topmod, B = numboots) #calculates bootstrapped means and SE for predicted (estimated) occupancy in each year
            
            #create dataframe with the predicted mean, boostrapped predicted mean, and bootstrapped SE for occupancy in each year
            bootres <- data.frame(pred_mean = bootmod@projected.mean[2,],
                                  boot_mean = smoothed(bootmod)[2,],
                                  SE = bootmod@smoothed.mean.bsse[2,])
            
            #calculate vector of 85% confidence intervals around the derived occupancy estimates
            # build the inputs for the 85% CI function built at the start
            {
              true_occ <- apply(true.occ.all, MARGIN=2, FUN=sum)/nrow(true.occ.all)
              mar85 <- round(1.44*bootres$SE, 4)
            }
            
            # calculate derived annual occupancy estimates
            est_occ = round(bootres$pred_mean,3)
            
            # calculate realized lambda and observed lambda
            lambda.real=rep(NA,yrs-1)
            lambda.obs=rep(NA,yrs-1)
            for(y in 1:(yrs-1)){
              lambda.real[y] <- true_occ[y+1]/true_occ[y]
              lambda.obs[y] <- est_occ[y+1]/est_occ[y]
            }
            
            
            # Store results in the correct index of list.results
            output.list[[min(which(as.numeric(summary(output.list)[,1])==0))]] <- list(
              species.type = spp_type[o],  #  species type
              start.occ = occ[o],          # initial occupancy
              lambda = lambda[l],          # programmed annual rate of change
              lambda_real = lambda.real,   # realized lambda based on true occupancy
              lambda_obs = lambda.obs,     # observed lambda based on estimated occupancy
              num.sites = sites[s],        # surveyed sites
              det.prob = det[d],      # detection
              num.visits = visits[v], # number of surveys per year
              top_mod_nm = topmod_name,        # top model name
              top_mod_str = topmod@estimates,        # complete summary of the top model
              true_occ = true_occ,    # true occupancy
              est_occ = est_occ,  # estimated occupancy
              est_occ_SE = bootres$SE,  # estimated occupancy
              mar85 = mar85,                # 85% CI
              accurate_pred = sapply(1:yrs, FUN=good_pred) # accurate prediction each year
            )
            
            print(paste("Scenario ", 
                        min(which(as.numeric(summary(output.list)[,1])==0))+1, 
                        "/", 
                        scenarios, "for Iteration ", i)
                  , sep="")
          }
        }
      }
    }
  }
  print(paste("Iteration ", i, "/", n.its), sep="")
  if(i==1){
    print(paste("Elapsed time for one iteration: ", round(Sys.time()-start.time,2), " minutes"), sep="")
  }
  if(i==n.its){
    print(paste("Total elapsed time: ", round((Sys.time()-start.time)/60,2), " hours"), sep="")
  }
}

#save list.results as RDS file
# update name to match workspace #
saveRDS(output.list, file="Round5_HOLA_1_20260302.rds")

output.list[[1]]


### Exploring outputs ----

# Note: this is NOT the full power power analysis, just explortatory code

# if you manually parallelize the simulations using multiple R instances, you'll need to combine the reuslts into one object
outputs <- append(output.list, 
                  readRDS('D:/SoCal_Streams/Round5/Horned_Lark_2/Round5_HOLA_2_20260302.rds') ) # the other batch of outputs

### Occupancy & lambda results ----

par(mfrow=c(2,3))
{
  # occupancy #
  for(l in 1:length(lambda)){
    sample.true=matrix(NA, nrow=1, ncol=yrs)
    sample.est=matrix(NA, nrow=1, ncol=yrs)
    for(i in 1:(scenarios*n.its*2)){
      if(outputs[[i]]$lambda==lambda[l] & outputs[[i]]$num.visits==5){
        sample.true = rbind(sample.true, outputs[[i]]$true_occ)
        sample.est = rbind(sample.est, outputs[[i]]$est_occ)
      }
      if(i==scenarios*n.its*2){
        rbind(sample.true, outputs[[i]]$true_occ)
        rbind(sample.est, outputs[[i]]$est_occ)
        sample.true=sample.true[-1,]
        sample.est=sample.est[-1,]
      }
    }
    
    plot(NA, ylim=c(0,1), xlim=c(1,yrs),
         ylab='Occupancy', xlab='Time')
    for(i in 1:dim(sample.est)[1]){
      lines(sample.est[i,]~c(1:yrs), lwd=2, col=rgb(1,0,0,.05))
      lines(sample.true[i,]~c(1:yrs), lwd=2, col=rgb(0,0,1,.05))
    }
    if(l==2){
      legend('topright', legend=c("true occ", "estimated occ"),
             col=c('blue', 'red'), lty=1)
    }
    
  }
  
  # Lambda #
  for(l in 1:length(lambda)){
    lambda.real=matrix(NA, nrow=1, ncol=yrs-1)
    lambda.obs=matrix(NA, nrow=1, ncol=yrs-1)
    for(i in 1:(scenarios*n.its*2)){
      if(outputs[[i]]$lambda==lambda[l] & outputs[[i]]$num.visits==5){
        lambda.real = rbind(lambda.real, outputs[[i]]$lambda_real)
        lambda.obs = rbind(lambda.obs, outputs[[i]]$lambda_obs)
      }
      if(i==scenarios*n.its*2){
        rbind(lambda.real, outputs[[i]]$lambda_real)
        rbind(lambda.obs, outputs[[i]]$lambda_obs)
        lambda.real=lambda.real[-1,]
        lambda.obs=lambda.obs[-1,]
      }
    }
    
    plot(NA, ylim=c(0.7,1.3), xlim=c(1,yrs-1),
         ylab='Lambda', xlab='Time')
    for(i in 1:dim(lambda.obs)[1]){
      lines(lambda.obs[i,]~c(1:(yrs-1)), lwd=2, col=rgb(1,0,0,.05))
      lines(lambda.real[i,]~c(1:(yrs-1)), lwd=2, col=rgb(0,0,1,.05))
    }
    if(l==2){
      legend('topright', legend=c("true lambda", "observed lambda"),
             col=c('blue', 'red'), lty=1)
    }
    
  }
  
}

### other outputs ----
topmod.dt <- data.table::data.table(species.type = character(length(outputs)), 
                                    top_mod_psi = character(length(outputs)),
                                    top_mod_gamma = character(length(outputs)),
                                    top_mod_epsilon = character(length(outputs)))

#loop through outputs list and fill topmod.dt
for(i in 1:length(outputs)){
  topmod.dt$species.type[i] <- outputs[[i]]$species.type
  #get top model for occupancy and convert from formula to character format
  topmod.dt$top_mod_psi[i] <- deparse(outputs[[i]]$top_mod_str$psiformula)
  #get top model for colonization and convert from formula to character format
  topmod.dt$top_mod_gamma[i] <- deparse(outputs[[i]]$top_mod$gammaformula)
  #get top model for extinction and convert from formula to character format
  topmod.dt$top_mod_epsilon[i] <- deparse(outputs[[i]]$top_mod$epsilonformula)
}
str(topmod.dt)

str(outputs[[90]]$top_mod_str)

dim(topmod.dt[ topmod.dt$species.type=="common generalist"])
dim(topmod.dt[ topmod.dt$species.type=="common generalist" & topmod.dt$top_mod_gamma=="~1"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="common generalist" & topmod.dt$top_mod_gamma=="~habitat"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="common generalist" & topmod.dt$top_mod_gamma=="~time"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="common generalist" & topmod.dt$top_mod_gamma=="~habitat + time"])[1]/3.6

dim(topmod.dt[ topmod.dt$species.type=="uncommon generalist" & topmod.dt$top_mod_gamma=="~1"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="uncommon generalist" & topmod.dt$top_mod_gamma=="~habitat"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="uncommon generalist" & topmod.dt$top_mod_gamma=="~time"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="uncommon generalist" & topmod.dt$top_mod_gamma=="~habitat + time"])[1]/3.6

dim(topmod.dt[ topmod.dt$species.type=="horned lark" & topmod.dt$top_mod_gamma=="~1"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="horned lark" & topmod.dt$top_mod_gamma=="~habitat"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="horned lark" & topmod.dt$top_mod_gamma=="~time"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="horned lark" & topmod.dt$top_mod_gamma=="~habitat + time"])[1]/3.6



sample.true=matrix(NA, nrow=1, ncol=yrs)
sample.est=matrix(NA, nrow=1, ncol=yrs)
for(i in 1:(scenarios*n.its)){
  if(output.list[[i]]$start.occ==.6 & output.list[[i]]$lambda==0.97 & 
     #output.list[[i]]$det.prob==.3 & output.list[[i]]$num.visits==2 & 
     output.list[[i]]$num.sites==50){
    # sample.true[nrow(sample.true)+1,] <- output.list[[i]]$true_occ
    # sample.est[nrow(sample.est)+1,] <- output.list[[i]]$est_occ
    sample.true = rbind(sample.true, output.list[[i]]$true_occ)
    sample.est = rbind(sample.est, output.list[[i]]$est_occ)
    
  }
  if(i ==scenarios*n.its){
    # sample.true[nrow(sampl.true)+1,] <- output.list[[i]]$true_occ
    # sample.est[nrow(sample.est)+1,] <- output.list[[i]]$est_occ
    rbind(sample.true, output.list[[i]]$true_occ)
    rbind(sample.est, output.list[[i]]$est_occ)
    sample.true=sample.true[-1,]
    sample.est=sample.est[-1,]
  }
}

plot(NA, ylim=c(0,1), xlim=c(1,yrs),
     ylab='Occupancy', xlab='Time')
for(i in 1:dim(sample.est)[1]){
  lines(sample.true[i,]~c(1:yrs), lwd=2, col=rgb(0,0,0,.2))
  lines(sample.est[i,]~c(1:yrs), lwd=2, lty='dotted', col=rgb(1,0,0,.3))
}


#save list.results as RDS file
saveRDS(output.list, file="20250826_SimPar_1.rds")


# Part 2: Prediction testing and power analysis ----
library(data.table)
library(ggplot2)
library(lme4)
library(stringr)
library(grid)
library(ggh4x)
library(rlang)

n.its.total=n.its*2 # '2' because we appended 2 sets of sim outputs (see/search "outputs <- append")

# Convert the matrix to a data frame
dfscen <- data.frame("lambda"=rep(NA, scenarios),
                     "num.sites"=rep(NA, scenarios),
                     "num.visits"=rep(NA, scenarios))

row.num <- 1
#loop through each of the scenarios and fill out dfscen
for (l in lambda){
  for (n in sites){
    for (v in visits){
      dfscen[row.num,] <- c(l,n,v)
      row.num <- row.num+1
    }
  }
}

# calculate following variables at scenario level:
# RMSE
# CI coverage for occupancy estimate in each season (1-10)
# Proportion of runs in which top model matched expectation
# Proportion of runs in which pop change was "detected" (Twining and Fuller 2025 method) - not used

#add columns to dfscen to hold scenario-level assessments of survey/model adequacy
{
  dfscen$RMSE <- vector(mode="numeric", length=nrow(dfscen)) #RMSE for the occupancy model
  dfscen$cicov1 <- vector(mode="numeric", length=nrow(dfscen)) #CI coverage for occ estimate in season 1
  dfscen$cicov2 <- vector(mode="numeric", length=nrow(dfscen))
  dfscen$cicov3 <- vector(mode="numeric", length=nrow(dfscen))
  dfscen$cicov4 <- vector(mode="numeric", length=nrow(dfscen))
  dfscen$cicov5 <- vector(mode="numeric", length=nrow(dfscen))
  dfscen$cicov6 <- vector(mode="numeric", length=nrow(dfscen))
  dfscen$cicov7 <- vector(mode="numeric", length=nrow(dfscen))
  dfscen$cicov8 <- vector(mode="numeric", length=nrow(dfscen))
  dfscen$cicov9 <- vector(mode="numeric", length=nrow(dfscen))
  dfscen$cicov10 <- vector(mode="numeric", length=nrow(dfscen)) #CI coverage for occ estimate in season 10
  dfscen$topmod_good <- vector(mode="numeric", length=nrow(dfscen)) #proportion of runs where top model was correct
  dfscen$habmod_good <- vector(mode="numeric", length=nrow(dfscen)) #proportion of runs where habitat relationship was modeled correctly (ignore time)
  dfscen$cicov_trend <- vector(mode="numeric", length=nrow(dfscen)) #proportion of runs where pop trend was detected
  dfscen$y1y5ovlp <- vector(mode="numeric", length=nrow(dfscen)) #proportion of runs where CIs for year1 and year5 overlap
  dfscen$y5y10ovlp <- vector(mode="numeric", length=nrow(dfscen))
  dfscen$y1y10ovlp <- vector(mode="numeric", length=nrow(dfscen))
}


#Now loop through each scenario in dfscen, subset outputs to just the slots
### that match the scenario in question, and calculate power for that scenario
for (i in 1:nrow(dfscen)){
  this_row <- dfscen[i, ]
  
  #subset outputs to just the slots that match the ith scenario
  subset_list <- Filter(function(x)
    x$lambda  == this_row$lambda  &
      x$num.sites  == this_row$num.sites &
      x$num.visits == this_row$num.visits,
    outputs)
  #break the loop and display an error message if subset_list is not n.its.total elements long
  if (length(subset_list) != n.its.total) {
    stop(paste0("Error: subset_list is not ", n.its.total, " elements long"))
  }
  
  #calculate squared errors for each slot of subset_list and store in vector
  sqerrs <- unlist(lapply(subset_list, function(z) (z$est_occ - z$true_occ)^2))
  #calculate RMSE by first taking mean of sqerrs and then taking sqrt
  dfscen$RMSE[i] <- sqrt(mean(sqerrs))
  ## reason we first average the sqerrs of all the models and THEN take the sqrt
  #### of the whole thing rather than separately calculating the RMSE for each model
  ##### and then averaging those RMSEs is because averaging the RMSE of each model
  ###### would underestimate variance b/c you'd be performing sqrt operation multiple times
  
  #calculate CI coverage for occupancy estimate in 1st season
  ### "accurate estimate" is defined as the 85% CI of the estimate encompassing the true value
  # Extract the 1st value from "accurate_pred" slot of each run
  acc1 <- sapply(subset_list, function(x) x$accurate_pred[1])
  #calculate proportion of values in acc1 equal to 1
  dfscen$cicov1[i] <- length(acc1[acc1==1])/length(acc1)
  
  #Now do seasons 2-10
  acc2 <- sapply(subset_list, function(x) x$accurate_pred[2])
  dfscen$cicov2[i] <- length(acc2[acc2==1])/length(acc2)
  
  acc3 <- sapply(subset_list, function(x) x$accurate_pred[3])
  dfscen$cicov3[i] <- length(acc3[acc3==1])/length(acc3)
  
  acc4 <- sapply(subset_list, function(x) x$accurate_pred[4])
  dfscen$cicov4[i] <- length(acc4[acc4==1])/length(acc4)
  
  acc5 <- sapply(subset_list, function(x) x$accurate_pred[5])
  dfscen$cicov5[i] <- length(acc5[acc5==1])/length(acc5)
  
  acc6 <- sapply(subset_list, function(x) x$accurate_pred[6])
  dfscen$cicov6[i] <- length(acc6[acc6==1])/length(acc6)
  
  acc7 <- sapply(subset_list, function(x) x$accurate_pred[7])
  dfscen$cicov7[i] <- length(acc7[acc7==1])/length(acc7)
  
  acc8 <- sapply(subset_list, function(x) x$accurate_pred[8])
  dfscen$cicov8[i] <- length(acc8[acc8==1])/length(acc8)
  
  acc9 <- sapply(subset_list, function(x) x$accurate_pred[9])
  dfscen$cicov9[i] <- length(acc9[acc9==1])/length(acc9)
  
  acc10 <- sapply(subset_list, function(x) x$accurate_pred[10])
  dfscen$cicov10[i] <- length(acc10[acc10==1])/length(acc10)
  
  #create empty vector in which to store results for top model for each run (correct or incorrect)
  topmod_good <- vector(mode="integer", length=length(subset_list))
  #create empty vector in which to store results for habitat relationship for each run (modeled correctly or not)
  habmod_good <- vector(mode="integer", length=length(subset_list))
  #loop through each run in subset_list and determine if top mod matched expectation
  for(r in 1:length(subset_list)){
    topmod_good[r] <- ifelse(dfscen$lambda[i]==1 & 
                               subset_list[[r]]$top_mod_nm == "mod3", 1,
                             ifelse(dfscen$lambda[i]!=1 & 
                                      subset_list[[r]]$top_mod_nm == "mod4", 1,
                                    ifelse(dfscen$lambda[i]==1 & 
                                             subset_list[[r]]$top_mod_nm == "mod1", 1, 
                                           ifelse(dfscen$lambda[i]!=1 & 
                                                    subset_list[[r]]$top_mod_nm == "mod2", 1, 0))))
    habmod_good[r] <- ifelse(subset_list[[r]]$top_mod_nm %in% c("mod3","mod4") , 1,
                             ifelse(subset_list[[r]]$top_mod_nm %in% c("mod1","mod2"), 1, 0))
  }
  #calculate proportion of topmods that were 1 (i.e., matched expectation)
  dfscen$topmod_good[i] <- mean(topmod_good)
  #calculate proportion of top models that correctly modeled at least the habitat relationship
  dfscen$habmod_good[i] <- mean(habmod_good)
  
  #following Twining and Fuller 2025, fit linear model est_occ ~ time for each run 
  ### of each scenario with a population change over time (i.e. exclude the scenarios with no population change over time)
  ## if the CI for the time parameter estimate does not overlap 0 and is the same
  ### sign as the direction of the simulated population trend, then the trend is 
  #### "correctly" detected. The proportion of runs in which the trend is correctly
  #### detected = "power" (i.e. CI coverage) to detect population trend. 
  ### If the CI coverage is > 0.80 for a given set of ecological parameters and sampling design,
  #### then that sampling design "works" for the ecological scenario in question
  
  #CI coverage for trend is only applicable if ith scenario included population change over time
  if (dfscen$lambda[i] != 1){
    #create empty vector in which to store results for whether trend was "detected" in each run of subset_list
    trend_detected <- vector(mode="integer", length=length(subset_list))
    #loop through subset_list
    for (rr in 1:length(subset_list)){
      #create dataframe to use for lm
      lmdat <- data.frame(est_occ = subset_list[[rr]]$est_occ, time = seq(1:10))
      #fit linear model est_occ ~ time for each run
      lm_run <- lm(est_occ ~ time, data=lmdat) #warning: essentially perfect fit: summary may be unreliable
      #calculate CI around time param estimate
      CI_run <- confint(lm_run, parm = "time", level = 0.85)
      #calculate whether CI_run contained 0 (TRUE if it contains 0, FALSE if it doesn't contain 0)
      contains0 <- ifelse(CI_run[1] <= 0 & CI_run[2] >= 0, TRUE, FALSE)
      #calculate whether parameter estimate for time is the correct sign
      right_sign <- ifelse((dfscen$lambda[i] > 1 & summary(lm_run)$coefficients["time", "Estimate"] > 0) |
                             (dfscen$lambda[i] < 1 & summary(lm_run)$coefficients["time", "Estimate"] < 0), 
                           TRUE, FALSE)
      
      #if CI_run excludes zero and sign of param est is correct, 1, otherwise, 0
      trend_detected[rr] <- ifelse(!contains0 & right_sign, 1, 0)
    }
    #calculate proportion of runs where trend was detected, and save in dfscen$cicov_trend
    dfscen$cicov_trend[i] <- sum(trend_detected)/length(trend_detected)
  }else{dfscen$cicov_trend[i] <- NA} #if ith scenario had no pop change overtime, then cicov_trend is NA
  
  ## function to test if two confidence intervals for occ estimates in 2 different years overlap
  overlap <- function(list_slot, year1, year2){
    #calculate upper bound of first CI
    upper1 <- list_slot$est_occ[year1] + list_slot$mar85[year1]
    #calculate lower bound of first CI
    lower1 <- list_slot$est_occ[year1] - list_slot$mar85[year1]
    #calculate upper bound of second CI
    upper2 <- list_slot$est_occ[year2] + list_slot$mar85[year2]
    #calculate lower bound of second CI
    lower2 <- list_slot$est_occ[year2] - list_slot$mar85[year2]
    #if the two CIs overlap, return 1, otherwise return 0
    ifelse((lower1 > lower2 & lower1 < upper2) | (upper1 > lower2 & upper1 < upper2), 1, 0)
  }
  
  ## vectors to store overlap (yes or no) for each run of ith scenario
  y1y5 <- vector(mode="integer", length=length(subset_list))
  y5y10 <- vector(mode="integer", length=length(subset_list))
  y1y10 <- vector(mode="integer", length=length(subset_list))
  for (y in 1:length(subset_list)){
    y1y5[y] <- overlap(subset_list[[y]], 1, 5)
    y5y10[y] <- overlap(subset_list[[y]], 5, 10)
    y1y10[y] <- overlap(subset_list[[y]], 1, 10)
  }
  #calculate proportion of runs where CIs overlap for each pairwise comparison
  dfscen$y1y5ovlp[i] <- mean(y1y5)
  dfscen$y5y10ovlp[i] <- mean(y5y10)
  dfscen$y1y10ovlp[i] <- mean(y1y10)
  
}

# 'warning: essentially perfect fit' can be ignored


#convert dfscen to data.table for easier manipulation
dtscen <- data.table::as.data.table(dfscen)

#create column in dtscen with unique ID for each survey design
dtscen[,srvdes := paste0("s",as.character(num.sites),"_v",as.character(num.visits))]

#save dtscen as .csv file
write.csv(dtscen, "20260302_RMSE_and_CI_coverage_each_HOLA_scenario.csv")


## Inferential stats to analyze power and RMSE at the scenario level ----

## include season as an independent variable in the model
#convert relevant columns into ordered factors or numeric format
dtscen[,":=" (srvdes = factor(srvdes, levels=c("s50_v3", "s50_v5", "s100_v3", "s100_v5", "s200_v3", "s200_v5")),
              lambda = as.numeric(lambda),
              num.sites = as.numeric(num.sites),
              num.visits = as.numeric(num.visits))]

#convert dtscen to long format and add column for season
dtscen_long <- data.table::melt(
  dtscen,
  measure.vars = patterns("^cicov\\d+$"),   # columns to pivot
  variable.name = "season",                 # new column with "cicov1", "cicov2", ...
  value.name = "cicov_occ"                      # new column with the values
)

# make season numeric instead of "cicov1" etc.
dtscen_long[, season := as.integer(gsub("cicov", "", season))]

# make sure columns are in correct format
dtscen[, `:=`(
  lambda       = as.numeric(lambda),
  num.sites    = as.numeric(num.sites),
  num.visits   = as.numeric(num.visits)
)]

### gamma regression for RMSE ----

# using gamma because it has no upper bound (just like RMSE) and is often a good fit for right-skewed data (like our RMSE values)
rmse.gamma <- glm(RMSE ~ poly(scale(lambda), 2) + # deleted 'species.type' and 'det.prob'
                    scale(num.sites) + scale(num.visits),
                  family = Gamma(link = "log"), data = dtscen)
summary(rmse.gamma)

#run binomial (logistic) regression on power to detect population change over 10 years
## first create new column in dtscen to indicate count of "successes" for each scenario
## a "success" indicates that the trend WAS detected, whereas y1y10olvp is the proportion of runs where the trend was NOT detected, hence the (1-y1y10ovlp)
dtscen[,k_pow10 := n.its.total*(1-y1y10ovlp)]
## response variable for each model is of the form (cbind(k, n-k))
### this is the grouped binomial format, where you provide glm() with the counts per group
#### rather than the Bernoulli outcome (0 or 1) for each individual trial

#Model for power to detect 10 year trend
## only use scenarios with a simulated pop change
##scale num.sites and num.visits so they can be compared directly
pow10.binom <- glm(cbind(k_pow10, n.its.total-k_pow10) ~ #species.type*scale(det.prob) + # removing species type and detection
                     scale(lambda) + scale(num.sites) + scale(num.visits),
                   family = binomial(link = "logit"), data = dtscen[lambda!=1,])
summary(pow10.binom)

### calculate power for each srvdes ----
## declining pop
powtab_decline <- dtscen[lambda==0.93, .(pow10 = 1-y1y10ovlp, pow5 = 1-y1y5ovlp),
                         by = .(srvdes)]
## increasing pop
powtab_increase <- dtscen[lambda==1.07, .(pow10 = 1-y1y10ovlp, pow5 = 1-y1y5ovlp),
                          by = .(srvdes)]



