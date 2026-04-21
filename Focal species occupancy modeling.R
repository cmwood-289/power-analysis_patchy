# Single-species occupancy models as power analysis case studies
# Connor Wood, 3-October 2025

library(unmarked)

# SoCal riparian ####
{
  ## read in the data ####
  effort <- read.csv('Data_SoCal/20251003_effort_file_all_SoCal_streams_sites.csv')
  effort=effort[order(effort$Site_ID),]
  
  habitat <- read.csv('Data_SoCal/Habitat_covariates.csv')
  habitat=habitat[order(habitat$Site_ID),]
  habitat$HabType <- as.numeric(ifelse(grepl("-RT0", habitat$Site_ID), 1,
                                       ifelse(grepl("-UT0", habitat$Site_ID), 0, "ERR")) )
  # Riparian coded as 1, Upland coded as 0
  habitat$LowRip <- as.numeric(ifelse(habitat$Elevation_mean<600 & habitat$HabType==1, yes=1, no=0))
  sum(habitat$HabType)
  sum(habitat$HabType)/dim(habitat)[1]
  sum(habitat$LowRip)/dim(habitat)[1]
  str(habitat)


  if(sum(habitat$Site_ID==effort$Site_ID)==dim(habitat)[1]){
    print('Habitat and effort site IDs all match')
  } else {
    print('WARNING: "habitat" and "effort" site IDs do NOT match')
  }
  
  summary(habitat$Elevation_mean[habitat$HabType==0])
  summary(habitat$Elevation_mean[habitat$HabType==1])
  
  # update nomenclature to match max_score files
  habitat$Site_ID=paste(habitat$Site_ID,habitat$Swift_ID,sep='-')
  
  # Bewick's Wren
  bewr.max <- read.csv('Data_SoCal/Bewicks Wren_max_score_summary_combined.csv')
  bewr.max=bewr.max[order(bewr.max$Site_ID),]
  bewr.max=bewr.max[bewr.max$Site_ID %in% habitat$Site_ID,]
  
  if(sum(habitat$Site_ID==bewr.max$Site_ID)==dim(habitat)[1]){
    print('Habitat and BEWR site IDs all match')
  } else {
    print('WARNING: "habitat" and "BEWR" site IDs do NOT match')
  }
  
  # Least Bell's Vireo
  lbvi.max <- read.csv('Data_SoCal/Bells Vireo_max_score_summary_combined.csv')
  lbvi.max=lbvi.max[order(lbvi.max$Site_ID),]
  lbvi.max=lbvi.max[lbvi.max$Site_ID %in% habitat$Site_ID,]
  
  if(sum(habitat$Site_ID==lbvi.max$Site_ID)==dim(habitat)[1]){
    print('Habitat and lbvi site IDs all match')
  } 
  else {
    print('WARNING: "habitat" and "LBVI" site IDs do NOT match')
  }
  
  lbvi.val <- read.csv('Data_SoCal/20251021_LBVI_sites.csv') # manually confirmed locations
  lbvi.val=lbvi.val[order(lbvi.val$Site_ID),]
  lbvi.val=lbvi.val[lbvi.val$Site_ID %in% substr(habitat$Site_ID, 1, 24),]
  
  
  # remove Site_ID columns
  {
    rownames(habitat) <- habitat$Site_ID
    rownames(effort) <- effort$Site_ID
    effort = effort[2:64]
    rownames(bewr.max) <- bewr.max$Site_ID
    bewr.max = bewr.max[2:64]
    rownames(lbvi.max) <- lbvi.max$Site_ID
    lbvi.max = lbvi.max[2:64]
  }
  
  ## prep the encounter histories ####
  
  # study design parameters
  {
    all.arus=sort(unique(habitat$Site_ID))
    ssp.length=3 # matches Copper Fire analysis
    sampling.start=seq(from=1, to=63, by=ssp.length)
    sampling.stop=seq(from=ssp.length, to=63, by=ssp.length)
    sample.intervals=length(sampling.start)
    n.ssp=ceiling(63/ssp.length)
    eff=matrix(NA, nrow=length(all.arus), ncol=n.ssp)
    survey.time=matrix(NA, nrow=length(all.arus), ncol=n.ssp)
    dets.bewr=matrix(NA, nrow=length(all.arus), ncol=n.ssp)
    dets.lbvi=matrix(NA, nrow=length(all.arus), ncol=n.ssp)
  }
  
  # build effort covariate and encounter histories
  for(a in 1:length(all.arus)){
    for(s in 1:sample.intervals){
      eff[a,s]=sum(effort[a, sampling.start[s]:sampling.stop[s]], na.rm=T)
      survey.time[a,s]=s
      if(eff[a,s]>0){
        # BEWR
        if( max(bewr.max[a, sampling.start[s]:sampling.stop[s]], na.rm=T) >= 0.968756283){
          dets.bewr[a,s]=1
        } else {
          dets.bewr[a,s]=0
        }
        # LBVI
        if( max(lbvi.max[a, sampling.start[s]:sampling.stop[s]], na.rm=T) >= 0.9253987 &
            lbvi.val$valid_TP[a]==1){
          dets.lbvi[a,s]=1
        } else {
          dets.lbvi[a,s]=0
        }
      }
    }
  }
  
  ## Fit Models ####
  
  ### Bewick's Wren ####
  occ_BEWR <- unmarkedFrameOccu(y = dets.bewr, 
                                siteCovs = habitat, 
                                obsCovs = list(effort=eff, 
                                               time=survey.time) )
  sum(occ_BEWR@siteCovs$HabType)
  sum(occ_LBVI@siteCovs$LowRip)
  
  
  bewr.0 <- occu(~1 ~1, occ_BEWR)
  bewr.1 <- occu(~effort ~1, occ_BEWR)
  bewr.2 <- occu(~time ~1, occ_BEWR)
  bewr.12 <- occu(~effort + time ~1, occ_BEWR)
  
  modSel(fitList(null=bewr.0, effort=bewr.1, time=bewr.2, effort_time=bewr.12))
  backTransform(bewr.1,type='state')
  
  bewr.1.1 <- occu(~effort ~HabType, occ_BEWR)
  bewr.1.2 <- occu(~effort ~Canopy_height_mean, occ_BEWR)
  bewr.1.3 <- occu(~effort ~LowRip, occ_BEWR)
  
  modSel(fitList(null=bewr.1, habitat_type=bewr.1.1, canopy_height=bewr.1.2, low_rip=bewr.1.3))
  
  # for supplemental materials
  # write.csv(as(modSel(fitList(psi_null=bewr.1, habitat_type=bewr.1.1, canopy_height=bewr.1.2, low_rip=bewr.1.3,
  #                  null=bewr.0, effort=bewr.1, time=bewr.2, effort_time=bewr.12)), "data.frame"),
  #           file='BEWR_modsel.csv')

  summary(bewr.1.1)

  plogis(-0.2840 +0.0388*36) # p at max effort = 0.75 - updated 'hab'
  plogis(2.37) # psi in upland = 0.915
  plogis(2.37 - 2.36) # psi in upland = 0.502
  
  
  
  ### Least Bell's Vireo ####
  #dets.lbvi_og=dets.lbvi # '_og' is the unfiltered version
  #occ_LBVI_og=occ_LBVI # '_og' is the unfiltered version
  occ_LBVI <- unmarkedFrameOccu(y = dets.lbvi, 
                                siteCovs = habitat, 
                                obsCovs = list(effort=eff, 
                                               time=survey.time) )
  
  lbvi.0 <- occu(~1 ~1, occ_LBVI)
  lbvi.1 <- occu(~effort ~1, occ_LBVI)
  lbvi.2 <- occu(~time ~1, occ_LBVI)
  lbvi.12 <- occu(~effort + time ~1, occ_LBVI)
  
  modSel(fitList(null=lbvi.0, effort=lbvi.1, time=lbvi.2, effort_time=lbvi.12))
  summary(lbvi.2)
  
  #plogis(-1.3180 + 0.0574*10) # p mid-season, og model
  plogis(-1.5366 + 0.0947*10) # p mid-season = 0.357, new model
  
  
  lbvi.2.1 <- occu(~time ~HabType, occ_LBVI)
  lbvi.2.2 <- occu(~time ~Canopy_height_mean, occ_LBVI)
  lbvi.2.3 <- occu(~time ~LowRip, occ_LBVI)
  
  modSel(fitList(null=lbvi.2, habitat_type=lbvi.2.1, canopy_height=lbvi.2.2, low_rip=lbvi.2.3))
  
  # for supplemental materials
  # write.csv(as(modSel(fitList(psi_null=lbvi.2, habitat_type=lbvi.2.1, canopy_height=lbvi.2.2, low_rip=lbvi.2.3,
  #                             null=lbvi.0, effort=lbvi.1, time=lbvi.2, effort_time=lbvi.12)), "data.frame"),
  #           file="LBVI_modsel.csv")
  
  summary(lbvi.2.3)
    plogis(-1.4712  + 0.0913*1) # p early-season = 0.20
  plogis(-1.4712  + 0.0913*10) # p mid-season = 0.36
  plogis(-1.4712  + 0.0913*21) # p late-season = 0.61
  
  
  plogis(-3.77) # psi in non-LowRip = 0.023
  plogis(-3.77 + 3.39) # psi in LowRip = 0.406
  
  
  summary(habitat$HabType) # 1 = all riparian, 0 = upland - 57% riparian
  summary(habitat$LowRip) # 1 = low rip, 0 = all else - 37% low riparian
  
  plot(habitat$HabType~habitat$Elevation_mean)
  points(habitat$HabType[substr(habitat$Site_ID, 1, 24) %in% lbvi.val$Site_ID[lbvi.val$valid_TP==1] ]
                       ~habitat$Elevation_mean[substr(habitat$Site_ID, 1, 24) %in% lbvi.val$Site_ID[lbvi.val$valid_TP==1]],
                                               col='red')
  
  lbvi.val=lbvi.val[lbvi.val$Site_ID %in% substr(habitat$Site_ID, 1, 24),]
  
  summary(habitat$LowRip)

  ####

}

# Kenyan savanna ####
{
  ## Read in .csv files ####
  effort_days <- read.csv('Data_Kenya/effort_days.csv')
  rownames(effort_days) <- effort_days$ID
  effort_days <- subset(effort_days, select = -ID)
  #head(effort_days)
  
  effort_hrs <- read.csv('Data_Kenya/effort_hrs.csv')
  rownames(effort_hrs) <- effort_hrs$ID
  effort_hrs <- subset(effort_hrs, select = -ID)
  #head(effort_hrs)
  
  g2_habitat <- read.csv('Data_Kenya/RBS_Grid2_covs.csv')
  rownames(g2_habitat) <- g2_habitat$ID
  g2_habitat <- subset(g2_habitat, select = -ID)
    hist(g2_habitat$CC)
    hist(g2_habitat$CC)
    hist(g2_habitat$Sand0)
    plot(g2_habitat$CC~g2_habitat$Sand0)
    cor(g2_habitat$CC, g2_habitat$Sand0) # weak negative: -0.15
    
    summary(g2_habitat$Sand0)
    sd(g2_habitat$Sand0)
    shapiro.test(g2_habitat$Sand0)
    ks.test(g2_habitat$Sand0, "pnorm")
    ks.test(g2_habitat$Sand0, "pnorm", 
            mean=mean(g2_habitat$Sand0),
            sd=sd(g2_habitat$Sand0))
    hist(rnorm(53,mean=mean(g2_habitat$Sand0),
               sd=sd(g2_habitat$Sand0)))
    
    
  eh_ASCI <- read.csv('Data_Kenya/Ashy Cisticola.csv')
  rownames(eh_ASCI) <- eh_ASCI$ID
  eh_ASCI <- subset(eh_ASCI, select = -ID)
  #head(eh_ASCI)
  
  eh_GRWH <- read.csv('Data_Kenya/Green Woodhoopoe.csv')
  rownames(eh_GRWH) <- eh_GRWH$ID
  eh_GRWH <- subset(eh_GRWH, select = -ID)
  #head(eh_GRWH)
  

  ## Fit Models ####
  ### Ashy Cisticola #### 
  occ_ASCI <- unmarkedFrameOccu(y = eh_ASCI[,1:10], 
                                siteCovs = g2_habitat, 
                                obsCovs = list(hours=effort_hrs, 
                                               days=effort_days)) 
  
  asci.0 <- occu( ~1 ~1, occ_ASCI)
  asci.1 <- occu(~hours ~1, occ_ASCI)
  asci.2 <- occu(~days ~1, occ_ASCI)
  asci.12 <- occu(~hours+days ~1, occ_ASCI)
  modSel(fitList(null=asci.0, hours=asci.1, days=asci.2, hours_days=asci.12))
  
  asci.1.1 <- occu(~hours ~log(CC+.01), occ_ASCI)
  asci.1.2 <- occu(~hours ~Sand0, occ_ASCI)
  modSel(fitList(null=asci.1, canopy=asci.1.1, sand=asci.1.2))
  
  # for supplemental materials
  # write.csv(as( modSel(fitList(psi_null=asci.1, canopy=asci.1.1, sand=asci.1.2,
  #                              null=asci.0, hours=asci.1, days=asci.2, hours_days=asci.12)), "data.frame"),
  #           file="ASCI_modsel.csv")

  summary(asci.1.2)  
  hist(g2_habitat$Sand0)
  plogis(3.761 - 0.151*30) # psi at 30% sand = 0.31
  plogis(3.761 - 0.151*40) # psi at 40% sand = 0.09
  plogis(3.761 - 0.151*55) # psi at 55% sand = 0.011
  
  plogis(-1.0656 +0.0277*44) # p at mean effort = 0.538
  
  # low-sand specialist
 
  ### Green Woodhoopoe #### 
  occ_GRWH <- unmarkedFrameOccu(y = eh_GRWH[,1:10], 
                                siteCovs = g2_habitat, 
                                obsCovs = list(hours=effort_hrs, 
                                               days=effort_days)) 
  
  grwh.0 <- occu( ~1 ~1, occ_GRWH)
  grwh.1 <- occu(~hours ~1, occ_GRWH)
  grwh.2 <- occu(~days ~1, occ_GRWH)
  grwh.12 <- occu(~hours+days ~1, occ_GRWH)
  modSel(fitList(null=grwh.0, hours=grwh.1, days=grwh.2, hours_days=grwh.12))
  
  grwh.1.1 <- occu(~hours ~log(CC+.01), occ_GRWH)
  grwh.1.2 <- occu(~hours ~Sand0, occ_GRWH)
  modSel(fitList(null=grwh.1, canopy=grwh.1.1, sand=grwh.1.2))
  # for supplemental materials
  # write.csv(as( modSel(fitList(psi_null=grwh.1, canopy=grwh.1.1, sand=grwh.1.2,
  #                              null=grwh.0, hours=grwh.1, days=grwh.2, hours_days=grwh.12)), "data.frame"),
  #           file="GRWH_modsel.csv")
  
  
  summary(grwh.1.2)  
  hist(g2_habitat$Sand0)
  plogis(-7.696 +0.147*30) # psi at 30% sand = 0.036
  plogis(-7.696 +0.147*40) # psi at 40% sand = 0.14
  plogis(-7.696 +0.147*55) # psi at 55% sand = 0.60
  
  summary(occ_GRWH@obsCovs$hours)
  plogis(-3.1624 +0.0392*44) # p at mean effort = 0.192
  

  ## Predicting Occupancy ###
  summary(g2_habitat$Sand0)
  sand_obs=data.frame(Sand0=seq(25,60,by=.1))
  asci.pred=predict(asci.1.2, type="state", newdata=sand_obs)
  grwh.pred=predict(grwh.1.2, type="state", newdata=sand_obs)
  
  
  par(mfrow=c(1,2))
  {
    hist(g2_habitat$Sand0, main='',
         xlab='Amount of sand observed')
    
    plot(asci.pred$Predicted~sand_obs$Sand0, 
         ylim=c(0,1), type='l', lwd=8, col='bisque3',
         ylab='Estimated occupancy', xlab='Amount of sand')
    lines(asci.pred$lower~sand_obs$Sand0,
          lwd=2, lty=2, col='bisque3')
    lines(asci.pred$upper~sand_obs$Sand0,
          lwd=2, lty=2, col='bisque3')
    
    lines(grwh.pred$Predicted~sand_obs$Sand0, 
          lwd=8, col='darkgreen')
    lines(grwh.pred$lower~sand_obs$Sand0, 
          lwd=2, lty=2, col='darkgreen')
    lines(grwh.pred$upper~sand_obs$Sand0, 
          lwd=2, lty=2, col='darkgreen')
    
  }
  
}

# NYS grasslands ####
{
  ## read in the data ####
  {
    habitat <- read.csv('Data_NYS/sites_covariates.csv')
    habitat=habitat[order(habitat$site),]
    habitat$type <- ifelse(habitat$solar_or_grassland == "grassland",0, 1)
    # grassland coded as 0, solar coded as 1
    summary(habitat$type)
    str(habitat)
    sum(habitat$type)
    
    grsp_binary <- read.csv('Data_NYS/grsp_binary_encounter_history.csv')
    hola_binary <- read.csv('Data_NYS/hola_binary_encounter_history.csv')
    grsp_binary=grsp_binary[order(grsp_binary$site),]
    hola_binary=hola_binary[order(hola_binary$site),]
  }
  
  # check alignment, trim dataframes
  {
    habitat$site==bobo_binary$site
    bobo_binary$site==sosp_binary$site
    # all sites aligned
    
    rownames(habitat) <- habitat$site
    rownames(grsp_binary) <- grsp_binary$site
    rownames(hola_binary) <- hola_binary$site
    grsp_binary = grsp_binary[,3:32]
    hola_binary = hola_binary[,3:32]
  }
  
  # study design parameters
  {
    all.arus=sort(unique(habitat$site))
    ssp.length=6 
    sampling.start=seq(from=1, to=30, by=ssp.length)
    sampling.stop=seq(from=ssp.length, to=30, by=ssp.length)
    sample.intervals=length(sampling.start)
    n.ssp=ceiling(30/ssp.length)
    eff=matrix(NA, nrow=length(all.arus), ncol=n.ssp)
    survey.time=matrix(NA, nrow=length(all.arus), ncol=n.ssp)
    dets.grsp=matrix(NA, nrow=length(all.arus), ncol=n.ssp)
    dets.hola=matrix(NA, nrow=length(all.arus), ncol=n.ssp)
  }
  
  # build effort covariate and encounter histories
  for(a in 1:length(all.arus)){
    for(s in 1:sample.intervals){
      eff[a,s]=sum(!is.na(bobo_binary[a, sampling.start[s]:sampling.stop[s]]))
      survey.time[a,s]=s
      if(eff[a,s]>0){
        # BOBO
        if( max(bobo_binary[a, sampling.start[s]:sampling.stop[s]], na.rm=T) == 1){
          dets.bobo[a,s]=1
        } else {
          dets.bobo[a,s]=0
        }
        # SOSP
        if( max(sosp_binary[a, sampling.start[s]:sampling.stop[s]], na.rm=T) == 1) {
          dets.sosp[a,s]=1
        } else {
          dets.sosp[a,s]=0
        }
        # GRSP
        if( max(grsp_binary[a, sampling.start[s]:sampling.stop[s]], na.rm=T) == 1) {
          dets.grsp[a,s]=1
        } else {
          dets.grsp[a,s]=0
        }
        # HOLA
        if( max(hola_binary[a, sampling.start[s]:sampling.stop[s]], na.rm=T) == 1) {
          dets.hola[a,s]=1
        } else {
          dets.hola[a,s]=0
        }
        # RWBL
        if( max(rwbl_binary[a, sampling.start[s]:sampling.stop[s]], na.rm=T) == 1) {
          dets.rwbl[a,s]=1
        } else {
          dets.rwbl[a,s]=0
        }
        
      }
    }
  }
  
  ## Fit Models ####
  
  ### Grasshopper Sparrow ####
  occ_GRSP <- unmarkedFrameOccu(y = dets.grsp, 
                                siteCovs = habitat, 
                                obsCovs = list(effort=eff, 
                                               time=survey.time) )
  
  grsp.0 <- occu(~1 ~1, occ_GRSP)
  grsp.1 <- occu(~effort ~1, occ_GRSP)
  grsp.2 <- occu(~time ~1, occ_GRSP)
  grsp.12 <- occu(~effort + time ~1, occ_GRSP)
  
  modSel(fitList(null=grsp.0, effort=grsp.1, time=grsp.2, effort_time=grsp.12))
  backTransform(grsp.2,type='state')
  
  sort(rowSums(grsp_binary, na.rm=T))
  
  grsp.2.1 <- occu(~time ~type, occ_GRSP)
  grsp.2.2 <- occu(~time ~prop_total_grassland_500m, occ_GRSP)
  
  modSel(fitList(null=grsp.2, habitat_type=grsp.2.1, total_grassland=grsp.2.2))
  summary(grsp.2)
  
  # for supplementals
  # write.csv(as( modSel(fitList(psi_null=grsp.2, habitat_type=grsp.2.1, total_grassland=grsp.2.2,
  #                              null=grsp.0, effort=grsp.1, time=grsp.2, effort_time=grsp.12)), "data.frame"),
  #           file="GRSP_modsel.csv")
  
  plogis(-1.268  +0.397*1) # p at time 1 = 0.295 
  plogis(-1.268  +0.397*3) # p at time 3 = 0.481 
  plogis(-1.268  +0.397*5) # p at time 5 = 0.672 
  
  plogis(-1.65) # psi = 0.161
  
  ### Horned Lark ####
  occ_HOLA <- unmarkedFrameOccu(y = dets.hola, 
                                siteCovs = habitat, 
                                obsCovs = list(effort=eff, 
                                               time=survey.time) )
  
  hola.0 <- occu(~1 ~1, occ_HOLA)
  hola.1 <- occu(~effort ~1, occ_HOLA)
  hola.2 <- occu(~time ~1, occ_HOLA)
  hola.12 <- occu(~effort + time ~1, occ_HOLA)
  
  modSel(fitList(null=hola.0, effort=hola.1, time=hola.2, effort_time=hola.12))
  backTransform(hola.2,type='state')
  
  # they're preset at all but three sites
  sort(rowSums(grsp_binary, na.rm=T))
  
  hola.0.1 <- occu(~1 ~type, occ_HOLA)
  hola.0.2 <- occu(~1 ~prop_total_grassland_500m, occ_HOLA)
  
  modSel(fitList(null=hola.0, habitat_type=hola.0.1, total_grassland=hola.0.2))
  summary(hola.0.1)
  
  # for supplementals
  # write.csv(as( modSel(fitList(psi_null=hola.0, habitat_type=hola.0.1, total_grassland=hola.0.2,
  #                              null=hola.0, effort=hola.1, time=hola.2, effort_time=hola.12)), "data.frame"),
  #           file="HOLA_modsel.csv")
  
  plogis(-1.62) # p at max effort = 0.165 
  
  plogis(1.41 ) # psi in grassland = 0.804
  plogis(1.41 - 3.25 ) # psi in solar = 0.137
  summary(habitat$type)

}


