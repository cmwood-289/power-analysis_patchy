# General-purpose power analysis for various hypothetical species proviles
# Connor Wood - cmw289@cornell.edu
## First coded August 2025; last updated 13-May 2026

library(unmarked)
# function to calculate if true occupancy in each year was within 95% CI for predicted occupancy
good_pred <- function(y){ifelse(true_occ[y] >= bootres$pred_mean[y] - mar85[y] & 
                                  true_occ[y] <= bootres$pred_mean[y] + mar85[y], 1, 0)}
# returns 1 if true occ is within 85% CI
# returns 0 if true occ is not withing 85% CI

# Set parameters for simulation ----
{
  n.its = 20              # Number of iterations; total should be >=100. Testing with <5 is nice
  yrs = 10                # Number of survey years
  territories = 400       # Number of potential bird territories in the study area
  
  sites = c(50, 100, 200) # Number of territories that are sampled
  det = c(0.35, 0.65)     # Secondary sampling period detection probabilities
  visits = c(2, 4)        # Number of secondary sampling periods
  
  spp_type = c("common generalist",
               "uncommon generalist",
               "specialist")
  prop_habitat = 0.3      # Proportion of total territories that are suitable for the specialist
  
  occ = c(0.6, 0.3, 0.6)    # Starting occupancy values; occ[3] is a specialist (total rate = occ[3]*prop_habitat)
  lambda = c(.93, 1, 1.07)  # annual growth rate
  
  lambda.var = 0.015        # annual +/- variation in lambda
  phi.min = 0.85              # minimum site persistence
  phi.var = 0.03           # annual +/- variation in phi
  
  scenarios=length(sites)*length(det)*length(visits)*length(occ)*length(lambda)
  numboots = 5           # for derived estimates of occupancy; should use 50-100 when running this for real
  
  #create empty list to store results of each run
  output.list <- vector(mode="list",
                        length=(scenarios*n.its))
  test <- vector(mode="list",
                 length=(scenarios*n.its))
  
  ## Note on 'species types' ----
  # As of May 2026, 'spp_type' and 'occ' have no intrinsic connection but are linked in 'Main Loop'
  # A future update will enforce a connection, likely by requring their input as a data frame rather than separate vectors
  
  ## Note on 'specialist' profile ----
  # As of May 2026, the specialist only occupieds "habitat". 
  # The "SoCal" and "New York" case studies relax this assumption, allowing for differential occupancy across habitat types
  # The "Kenya" case study allows for species to respond to continuous habitat variables.
  # If there is interest, I will update this master code to include those possibilities here.
  
}

# Main Loop ----

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
      
      
      ## Generalists can occupy all sites ##
      if(o < 3) {
        # Create an empty matrix for this specific simulation
        true.occ.all <- matrix(NA, nrow = territories, ncol = yrs)
        
        true.occ.all[,1] <- as.numeric(runif(territories) < occ[o]) # sets initial site-by-site occupancy
        occ.sim = c(mean(true.occ.all[,1]), rep(NA,yrs-1)) # stores population-level occupancy
        
        for(y in 2:yrs){
          psi.t=runif(n=1, 
                      min=(lambda[l]-lambda.var), 
                      max=(lambda[l]+lambda.var))*
            occ.sim[y-1] # calculate psi for next year
          if(psi.t>1){
            occ.sim[y]=1
            true.occ.all[,y]=1
          } 
          else {
            occ.sim[y]=psi.t
            occ.t=round(psi.t*territories,0) # number of territories occupied in year y
            phi.t=round( sum(true.occ.all[,y-1]==1)*
                           runif(n=1, 
                                 min=(phi.min-phi.var), 
                                 max=(phi.min+phi.var)) ,0) # minimum number of territories that persist from t-1
            persist.t = sample(which(true.occ.all[,y-1]==1), 
                               phi.t,
                               replace = F) # randomly selects which sites persist between years
            true.occ.all[persist.t,y] <- 1
            remaining=c( rep(1, occ.t-phi.t), # sites occupied in addition to the minimum phi
                         rep(0, territories-occ.t)) 
            true.occ.all[-persist.t,y] <- sample(remaining,
                                                 length(remaining),
                                                 replace=F)
          }
        }
      }
      
      # calculates the number of territories with specialist's habitat
      hab.terr=round(territories*prop_habitat,0)
      # outside the 'specialist' loop because it's used to create the habitat covariate
      
      ## Specialist can only occupy its "habitat" ##
      if(o == 3) {
        # Create an empty matrix for the special habitat only
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
        true.occ.all=rbind(true.occ.hab,
                           matrix(0, nrow=(territories-hab.terr), ncol=yrs))
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
saveRDS(output.list, file="SimPar_1_20250919.rds")

output.list[[1]]

str(test[[1]]$top_mod_str@estimates)


# Exploring outputs ----
# Note: this is NOT the full power power analysis, just explortatory code

getwd()

# if you manually parallelize the simulations using multiple R instances, you'll need to combine the reuslts into one object
outputs <- append(output.list, 
                  readRDS('D:/SoCal_Streams/Round2/SimParallel2/20250827_1_SimPar_2.rds') )


topmod.dt <- data.table::data.table(species.type = character(length(outputs)), 
                                    top_mod_psi = character(length(outputs)),
                                    top_mod_gamma = character(length(outputs)),
                                    top_mod_epsilon = character(length(outputs)))

#loop through outputs list and fill topmod.dt
for(i in 1:length(outputs)){
  topmod.dt$species.type[i] <- outputs[[i]]$species.type
  #get top model for occupancy and convert from formula to character format
  topmod.dt$top_mod_psi[i] <- deparse(outputs[[i]]$top_mod$psiformula)
  #get top model for colonization and convert from formula to character format
  topmod.dt$top_mod_gamma[i] <- deparse(outputs[[i]]$top_mod$gammaformula)
  #get top model for extinction and convert from formula to character format
  topmod.dt$top_mod_epsilon[i] <- deparse(outputs[[i]]$top_mod$epsilonformula)
}
str(topmod.dt)

dim(topmod.dt[ topmod.dt$species.type=="common generalist"])
dim(topmod.dt[ topmod.dt$species.type=="common generalist" & topmod.dt$top_mod_gamma=="~1"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="common generalist" & topmod.dt$top_mod_gamma=="~habitat"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="common generalist" & topmod.dt$top_mod_gamma=="~time"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="common generalist" & topmod.dt$top_mod_gamma=="~habitat + time"])[1]/3.6

dim(topmod.dt[ topmod.dt$species.type=="uncommon generalist" & topmod.dt$top_mod_gamma=="~1"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="uncommon generalist" & topmod.dt$top_mod_gamma=="~habitat"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="uncommon generalist" & topmod.dt$top_mod_gamma=="~time"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="uncommon generalist" & topmod.dt$top_mod_gamma=="~habitat + time"])[1]/3.6

dim(topmod.dt[ topmod.dt$species.type=="specialist" & topmod.dt$top_mod_gamma=="~1"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="specialist" & topmod.dt$top_mod_gamma=="~habitat"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="specialist" & topmod.dt$top_mod_gamma=="~time"])[1]/3.6
dim(topmod.dt[ topmod.dt$species.type=="specialist" & topmod.dt$top_mod_gamma=="~habitat + time"])[1]/3.6


par(mfrow=c(3,3))

# Occupancy results
for(l in 1:length(lambda)){
  sample.true=matrix(NA, nrow=1, ncol=yrs)
  sample.est=matrix(NA, nrow=1, ncol=yrs)
  for(i in 1:(scenarios*n.its)){
    if(output.list[[i]]$species.type=="specialist" & output.list[[i]]$lambda==lambda[l] 
       & output.list[[i]]$det.prob==.4 & output.list[[i]]$num.visits==4){
      sample.true = rbind(sample.true, output.list[[i]]$true_occ)
      sample.est = rbind(sample.est, output.list[[i]]$est_occ)
    }
    if(i ==scenarios*n.its){
      rbind(sample.true, output.list[[i]]$true_occ)
      rbind(sample.est, output.list[[i]]$est_occ)
      sample.true=sample.true[-1,]
      sample.est=sample.est[-1,]
    }
  }
  
  plot(NA, ylim=c(0,1), xlim=c(1,yrs),
       ylab='Occupancy', xlab='Time')
  for(i in 1:dim(sample.est)[1]){
    lines(sample.est[i,]~c(1:yrs), lwd=2, col=rgb(1,0,0,.1))
    lines(sample.true[i,]~c(1:yrs), lwd=2, col=rgb(0,0,1,.15))
  }
  
}



par(mfrow=c(3,3))

# Lambda results
for(l in 1:length(lambda)){
  lambda.real=matrix(NA, nrow=1, ncol=yrs-1)
  lambda.obs=matrix(NA, nrow=1, ncol=yrs-1)
  for(i in 1:(scenarios*n.its)){
    if(output.list[[i]]$species.type=="specialist" & output.list[[i]]$lambda==lambda[l] 
       & output.list[[i]]$det.prob==.4 & output.list[[i]]$num.visits==4){
      lambda.real = rbind(lambda.real, output.list[[i]]$lambda_real)
      lambda.obs = rbind(lambda.obs, output.list[[i]]$lambda_obs)
    }
    if(i ==scenarios*n.its){
      rbind(lambda.real, output.list[[i]]$lambda_real)
      rbind(lambda.obs, output.list[[i]]$lambda_obs)
      lambda.real=lambda.real[-1,]
      lambda.obs=lambda.obs[-1,]
    }
  }
  
  plot(NA, ylim=c(0.7,1.3), xlim=c(1,yrs-1),
       ylab='Lambda', xlab='Time')
  for(i in 1:dim(lambda.obs)[1]){
    lines(lambda.obs[i,]~c(1:(yrs-1)), lwd=2, col=rgb(1,0,0,.1))
    lines(lambda.real[i,]~c(1:(yrs-1)), lwd=2, col=rgb(0,0,1,.15))
  }
  
}




par(mfrow=c(2,2))

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



