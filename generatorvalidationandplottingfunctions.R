## Functions for transient streamflow generator validation and plotting
##  D Gorelick (May 2018)
## -------------------------------------------------------------

## -------------------- plot flow exceedence curves ------------------------

exceedenceCurves = function(SynDaily, ModDaily, decade = 1, num_syn_yrs = 10, num_rlz = 10, 
                            GeneratorInputsFolder, scen = 1, setting = "daily")
{
  Decades = seq(2010,2060,10)
  
  # 1. Organize data ---------------------------
  
  M = ModDaily[,decade,]
  
  Sarr = SynDaily[,,,decade]
  S = matrix(NA, nrow = num_syn_yrs*num_rlz*365, ncol = ncol(M))
  
  for (site in 1:ncol(S)) {
    sMsort = data.frame(sort(M[,site], decreasing = TRUE))
    sMsort$P = ((1:nrow(sMsort))-0.5)/(nrow(sMsort)+1)
    colnames(sMsort)[1] = c("Flow"); sMsort$r = 0; sMsort$type = "RHESSys"
    
    Mmelt = data.frame(sMsort$P, sMsort$r, sMsort$Flow, sMsort$type)
    colnames(Mmelt) = c("P", "Var", "Val", "Type")
    
    sSsort = matrix(NA, nrow = dim(Sarr)[2], ncol = num_rlz)
    for (rlz in 1:num_rlz) {
      sSsort[,rlz] = sort(Sarr[rlz,,site], decreasing = TRUE)
    }
    sSsort = as.data.frame(sSsort)
    sSsort$SP = ((1:nrow(sSsort))-0.5)/(nrow(sSsort)+1)
    Smelt = melt(sSsort, id = "SP"); Smelt$type = "Synthetic"
    colnames(Smelt) = c("P", "Var", "Val", "Type")
    
    # 2. Plot exceedence curves ------------------
    
    temp = rbind(Smelt, Mmelt)
    
    FDC = ggplot(temp) + geom_line(aes(x = P, y = Val, group = Var, color = Type), size = 3, alpha = 0.5) + 
      theme_classic() + scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
      ggtitle(paste(Decades[decade], " flow duration curve comparison for ", SOURCES[site], sep = "")) + 
      xlab("Exceedence Probability")
    
    ggsave(paste("transientgeneratoroutput/scenario", scen, "/FDC", Decades[decade], "site", SOURCES[site], ".png", sep = ""),
           units = "in", width = 5, height = 5)
  }
}

## -------------- monthly flow duration curves ------------------------------

monthlyexceedenceCurves = function(Syn, Mod, num_syn_yrs = 10, num_rlz = 10, scen = 1, setting = "monthly")
{
  Decades = seq(2010,2060,10)
  SOURCES = c("Haw River", "Little River", "Flat River", "Cane Creek", "Jordan Evap.", "Michie Evap.", "CCR Evap.")
  #Syn = SyntheticArrayHolder
  #Mod = ModeledArray
  
  BigSet = c()
  for (site in 1:dim(Mod)[3]) {
    # 1. Organize data ---------------------------
    
    M = Mod[,,site,]
    S = Syn[,,,site,]
    
    Mfull = c()
    Sfull = c()
    for (d in 1:length(Decades)) {
      MD = M[,,d]; MDmelt = melt(MD); MDmelt$Decade = Decades[d]
      
      MDsort = MDmelt[order(MDmelt$value, decreasing = TRUE),]
      MDsort$Realization = 0
      MDsort$P = ((1:nrow(MDsort))-0.5)/(nrow(MDsort)+1)
      
      SD = S[,,d]
      for (rlz in 1:num_rlz) {
        SDR = SD[rlz,]
        SDRmat = matrix(SDR, ncol = 12, byrow = TRUE)
        SDRmelt = melt(SDRmat)
        SDRmelt$Decade = Decades[d]
        SDRmelt$Realization = rlz
        
        SDRsort = SDRmelt[order(SDRmelt$value, decreasing = TRUE),]
        SDRsort$P = ((1:nrow(SDRsort))-0.5)/(nrow(SDRsort)+1)
        
        Sfull = rbind(Sfull, SDRsort)
      }
      
      Mfull = rbind(Mfull, MDsort)
    }
    Mfull$Site = SOURCES[site]; Mfull$type = "B. RHESSys"
    Sfull$Site = SOURCES[site]; Sfull$type = "A. Synthetic"
    
    temp = rbind(Mfull, Sfull)
    BigSet = rbind(BigSet, temp)
  }
  
  # 2. Plot exceedence curves ------------------
  BSstream = BigSet[which(BigSet$Site %in% c("Haw River", "Little River", "Flat River", "Cane Creek")),]
  BSevap = BigSet[which(BigSet$Site %in% c("Jordan Evap.", "Michie Evap.", "CCR Evap.")),]
  
  BigSet1 = BSstream[which(BSstream$type == "B. RHESSys"),]
  BigSet2 = BSstream[which(BSstream$type != "B. RHESSys"),]
  
  FDC = ggplot() + 
    geom_line(data = BigSet2, aes(x = P, y = value, group = Realization), color = "grey70", size = 2, alpha = 1) + 
    geom_line(data = BigSet1, aes(x = P, y = value, group = Realization), color = "red", size = 1.5, alpha = 1) +
    facet_grid(Site ~ Decade, scales = "free_y") + scale_y_log10() +
    ggtitle("Decadal flow duration curves (of monthly flows) for each site") + xlab("Exceedence Probability") +
    scale_x_continuous(limits = c(0,1), breaks = c(0,1), expand = c(0.01,0.01)) + ylab("Flow (MG)") + theme_light() 
  
  ggsave(paste("transientgeneratoroutput/scenario", scen, "/FDCsMonthly_Streamflow.png", sep = ""),
         units = "in", width = 10, height = 8)
  
  BSevap$value[which(BSevap$type == "B. RHESSys")] = log(BSevap$value[which(BSevap$type == "B. RHESSys")])
  BigSet1 = BSevap[which(BSevap$type == "B. RHESSys"),]
  BigSet2 = BSevap[which(BSevap$type != "B. RHESSys"),]
  
  FDC = ggplot() + 
    geom_line(data = BigSet2, aes(x = P, y = value, group = Realization), color = "grey70", size = 2, alpha = 1) + 
    geom_line(data = BigSet1, aes(x = P, y = value, group = Realization), color = "red", size = 1.5, alpha = 1) +
    facet_grid(Site ~ Decade, scales = "free_y") +
    ggtitle("Decadal flow duration curves (of monthly flows) for each site") + xlab("Exceedence Probability") +
    scale_x_continuous(limits = c(0,1), breaks = c(0,1), expand = c(0.01,0.01)) + ylab("Flow (MG)") + theme_light() 
  
  ggsave(paste("transientgeneratoroutput/scenario", scen, "/FDCsMonthly_Evap.png", sep = ""),
         units = "in", width = 10, height = 6)
}

## ------------------ weekly flow duration curves ---------------------------------------

weeklyexceedenceCurves = function(Syn, Mod, num_syn_yrs = 10, num_rlz = 10, scen = 1, setting = "weekly")
{
  Decades = seq(2010,2060,10)
  SOURCES = c("Haw River", "Little River", "Flat River", "Cane Creek", "Jordan Evap.", "Michie Evap.", "CCR Evap.")
  #Syn = SyntheticWeekly
  #Mod = WeeklyArray
  
  BigSet = c()
  for (site in 1:dim(Mod)[3]) {
    # 1. Organize data ---------------------------
    
    M = Mod[,,site,]
    S = Syn[,,,site,]
    
    Mfull = c()
    Sfull = c()
    for (d in 1:length(Decades)) {
      MD = M[,,d]; MDmelt = melt(MD); MDmelt$Decade = Decades[d]
      
      MDsort = MDmelt[order(MDmelt$value, decreasing = TRUE),]
      MDsort$Realization = 0
      MDsort$P = ((1:nrow(MDsort))-0.5)/(nrow(MDsort)+1)
      
      SD = S[,,,d]
      for (rlz in 1:num_rlz) {
        SDRmat = SD[rlz,,]
        SDRmelt = melt(SDRmat)
        SDRmelt$Decade = Decades[d]
        SDRmelt$Realization = rlz
        
        SDRsort = SDRmelt[order(SDRmelt$value, decreasing = TRUE),]
        SDRsort$P = ((1:nrow(SDRsort))-0.5)/(nrow(SDRsort)+1)
        
        Sfull = rbind(Sfull, SDRsort)
      }
      
      Mfull = rbind(Mfull, MDsort)
    }
    Mfull$Site = SOURCES[site]; Mfull$type = "1. RHEESys"
    Sfull$Site = SOURCES[site]; Sfull$type = "2. Synthetic"
    
    temp = rbind(Mfull, Sfull)
    BigSet = rbind(BigSet, temp)
  }
  
  # 2. Plot exceedence curves ------------------
  BSstream = BigSet[which(BigSet$Site %in% c("Haw River", "Little River", "Flat River", "Cane Creek")),]
  BSevap = BigSet[which(BigSet$Site %in% c("Jordan Evap.", "Michie Evap.", "CCR Evap.")),]
  
  BigSet1 = BSstream[which(BSstream$type == "1. RHEESys"),]
  BigSet2 = BSstream[which(BSstream$type != "1. RHEESys"),]
  
  FDC = ggplot() + 
    geom_line(data = BigSet2, aes(x = P, y = value, group = Realization), color = "grey70", size = 2, alpha = 1) + 
    geom_line(data = BigSet1, aes(x = P, y = value, group = Realization), color = "red", size = 1.5, alpha = 1) +
    facet_grid(Site ~ Decade, scales = "free_y") + scale_y_log10() +
    ggtitle("Decadal flow duration curves (of weekly flows) for each site") + xlab("Exceedence Probability") +
    scale_x_continuous(limits = c(0,1), breaks = c(0,1), expand = c(0.01,0.01)) + ylab("Flow (MG)") + theme_light() 
  
  ggsave(paste("transientgeneratoroutput/scenario", scen, "/FDCsWeekly_Streamflow.png", sep = ""),
         units = "in", width = 10, height = 8)
  
  #BSevap$value[which(BSevap$type == "1. RHEESys")] = log(BSevap$value[which(BSevap$type == "1. RHEESys")])
  BigSet1 = BSevap[which(BSevap$type == "1. RHEESys"),]
  BigSet2 = BSevap[which(BSevap$type != "1. RHEESys"),]
  
  FDC = ggplot() + 
    geom_line(data = BigSet2, aes(x = P, y = value, group = Realization), color = "grey70", size = 2, alpha = 1) + 
    geom_line(data = BigSet1, aes(x = P, y = value, group = Realization), color = "red", size = 1.5, alpha = 1) +
    facet_grid(Site ~ Decade, scales = "free_y") +
    ggtitle("Decadal flow duration curves (of weekly flows) for each site") + xlab("Exceedence Probability") +
    scale_x_continuous(limits = c(0,1), breaks = c(0,1), expand = c(0.01,0.01)) + ylab("Flow (MG)") + theme_light() 
  
  ggsave(paste("transientgeneratoroutput/scenario", scen, "/FDCsWeekly_Evap.png", sep = ""),
         units = "in", width = 10, height = 6)
}

## -------------- daily flow duration curves ------------------------------

dailyexceedenceCurves = function(Syn, Mod, num_syn_yrs = 10, num_rlz = 10, scen = 1, setting = "daily")
{
  Decades = seq(2010,2060,10)
  SOURCES = c("Haw River", "Little River", "Flat River", "Cane Creek", "Jordan Evap.", "Michie Evap.", "CCR Evap.")
  #Syn = SyntheticDailyArray
  #Mod = DailyArray
  
  BigSet = c()
  for (site in 1:dim(Mod)[3]) {
    # 1. Organize data ---------------------------
    
    M = Mod[,,site]
    S = Syn[,,site,]
    
    Mfull = c()
    Sfull = c()
    for (d in 1:length(Decades)) {
      MD = matrix(M[,d], ncol = 365, byrow = TRUE); MDmelt = melt(MD); MDmelt$Decade = Decades[d]
      
      MDsort = MDmelt[order(MDmelt$value, decreasing = TRUE),]
      MDsort$Realization = 0
      MDsort$P = ((1:nrow(MDsort))-0.5)/(nrow(MDsort)+1)
      
      SD = S[,,d]
      for (rlz in 1:num_rlz) {
        SDR = SD[rlz,]
        SDRmat = matrix(SDR, ncol = 365, byrow = TRUE)
        SDRmelt = melt(SDRmat)
        SDRmelt$Decade = Decades[d]
        SDRmelt$Realization = rlz
        
        SDRsort = SDRmelt[order(SDRmelt$value, decreasing = TRUE),]
        SDRsort$P = ((1:nrow(SDRsort))-0.5)/(nrow(SDRsort)+1)
        
        Sfull = rbind(Sfull, SDRsort)
      }
      
      Mfull = rbind(Mfull, MDsort)
    }
    Mfull$Site = SOURCES[site]; Mfull$type = "1. RHEESys"
    Sfull$Site = SOURCES[site]; Sfull$type = "2. Synthetic"
    
    temp = rbind(Mfull, Sfull)
    BigSet = rbind(BigSet, temp)
  }
  
  # 2. Plot exceedence curves ------------------
  BSstream = BigSet[which(BigSet$Site %in% c("Haw River", "Little River", "Flat River", "Cane Creek")),]
  BSevap = BigSet[which(BigSet$Site %in% c("Jordan Evap.", "Michie Evap.", "CCR Evap.")),]
  
  BigSet1 = BSstream[which(BSstream$type == "1. RHEESys"),]
  BigSet2 = BSstream[which(BSstream$type != "1. RHEESys"),]
  
  FDC = ggplot() + 
    geom_line(data = BigSet2, aes(x = P, y = value, group = Realization), color = "grey70", size = 2, alpha = 1) + 
    geom_line(data = BigSet1, aes(x = P, y = value, group = Realization), color = "red", size = 1.5, alpha = 1) +
    facet_grid(Site ~ Decade, scales = "free_y") + scale_y_log10() +
    ggtitle("Decadal flow duration curves (of daily flows) for each site") + xlab("Exceedence Probability") +
    scale_x_continuous(limits = c(0,1), breaks = c(0,1), expand = c(0.01,0.01)) + ylab("Flow (MG)") + theme_light() 
  
  ggsave(paste("transientgeneratoroutput/scenario", scen, "/FDCsDaily_Streamflow.png", sep = ""),
         units = "in", width = 10, height = 8)
  
  #BSevap = BSevap[!is.na(BSevap$Site),]
  BigSet1 = BSevap[which(BSstream$type == "1. RHEESys"),]
  BigSet2 = BSevap[which(BSstream$type != "1. RHEESys"),]
  
  FDC = ggplot() + 
    geom_line(data = BigSet2, aes(x = P, y = value, group = Realization), color = "grey70", size = 2, alpha = 1) + 
    geom_line(data = BigSet1, aes(x = P, y = value, group = Realization), color = "red", size = 1.5, alpha = 1) +
    facet_grid(Site ~ Decade, scales = "free_y") +
    ggtitle("Decadal flow duration curves (of daily flows) for each site") + xlab("Exceedence Probability") +
    scale_x_continuous(limits = c(0,1), breaks = c(0,1), expand = c(0.01,0.01)) + ylab("Flow (MG)") + theme_light() 
  
  ggsave(paste("transientgeneratoroutput/scenario", scen, "/FDCsDaily_Evap.png", sep = ""),
         units = "in", width = 10, height = 6)
}

## --------------- calculate temporal autocorrelations -----------------------

temporalACF = function(Syn = SyntheticArrayHolder, Mod = ModeledArray, 
                       nlags = 12, scen = 1, setting = "monthly", flwst = 6)
{
  library(ggplot2); library(reshape2)
  SOURCES = c("Haw River", "Little River", "Flat River", "Cane Creek", "Jordan Evap.", "Michie Evap.", "CCR Evap.")
  DECS = seq(2010,2060,10); tempall = c()
  for (s in 1:dim(Syn)[length(dim(Syn))-1]) {
    temps = c()
    if (setting == "monthly") {SS = Syn[,,,s,]; MS = Mod[,,s,]}
    if (setting == "weekly") {SS = Syn[,,,s,]; MS = Mod[,,s,]}
    if (setting == "daily") {SS = Syn[,,s,]; MS = Mod[,,s]}
    for (d in 1:dim(Syn)[length(dim(Syn))]) {
      SSDRacrSet = c()
      if (setting == "monthly") {SSD = SS[,,d]; MSD = MS[,,d]}
      if (setting == "weekly") {SSD = SS[,,,d]; MSD = MS[,,d]}
      if (setting == "daily") {SSD = SS[,,d]; MSD = MS[,d]}
      for (r in 1:dim(SSD)[1]) {
        if (setting == "weekly") {SSDR = c(t(SSD[r,,]))} else {SSDR = SSD[r,]}
        SSDR_AC = acf(SSDR, lag.max = nlags, plot = FALSE); SSDRacrSet = rbind(SSDRacrSet, as.vector(SSDR_AC$acf))
      }
      SSDRacrMeans = c(1); SSDRacrCI = c(1,1)
      for (l in 2:(nlags+1)) {
        SSDRacrMeans[l] = t.test(as.vector(SSDRacrSet[,l]))$estimate
        SSDRacrCI = rbind(SSDRacrCI, t.test(as.vector(SSDRacrSet[,l]))$conf.int[1:2])
      }
      MSDacr = as.data.frame(as.vector(acf(c(t(MSD)), lag.max = nlags, plot = FALSE)$acf))
      MSDacr$LB = MSDacr[,1]; MSDacr$UB = MSDacr[,1]; MSDacr$Type = "RHESSys"; colnames(MSDacr)[1] = "Mean Autocorrelation"
      MSDacr$Lag = 0:nlags
      
      SSDacr = cbind(SSDRacrMeans, SSDRacrCI)
      colnames(SSDacr) = c("Mean Autocorrelation", "LB", "UB"); rownames(SSDacr) = 1:nrow(SSDacr)
      SSDacr = as.data.frame(SSDacr); SSDacr$Type = "Synthetic"; SSDacr$Lag = 0:nlags
      
      tempd = as.data.frame(rbind(MSDacr, SSDacr)); tempd$Decade = DECS[d]
      temps = rbind(temps, tempd)
    }
    temps$Site = SOURCES[s]
    tempall = rbind(tempall, temps)
  }
  SynSet = tempall[which(tempall$Type == "Synthetic"),]
  ModSet = tempall[which(tempall$Type != "Synthetic"),]
  ACFplot = ggplot() +
    geom_errorbar(data = SynSet, aes(x = Lag, ymin = LB, ymax = UB), 
                  color = "blue", size = 1) + 
    geom_line(data = SynSet, aes(x = Lag, y = `Mean Autocorrelation`), 
              color = "blue", size = 1) + 
    geom_line(data = ModSet, aes(x = Lag, y = `Mean Autocorrelation`), 
              color = "red", size = 0.75) + 
    facet_grid(Site ~ Decade) +
    ggtitle(paste(setting, " Temporal Autocorrelation for Scenario ", scen, sep = ""))
  ggsave(paste("transientgeneratoroutput/scenario", scen, "/", setting, "flowset", flwst, "TemporalAutocorrelation.png", sep = ""),
         units = "in", width = 8, height = 6)
}

## ------------ calculate spatial correlation --------------------

spatialCor = function(Syn = SyntheticArrayHolder, Mod = ModeledArray, 
                      setting = "monthly", scen = 1, flwst = 6)
{
  library(ggplot2); library(reshape2)
  DECS = seq(2010,2060,10); tempall = c()
  SOURCES = c("Haw River", "Little River", "Flat River", "Cane Creek", "Jordan Evap.", "Michie Evap.", "CCR Evap.")
  for (d in 1:dim(Syn)[length(dim(Syn))]) {
    tempd = c()
    if (setting == "monthly") {SD = Syn[,,,,d]; MD = Mod[,,,d]}
    if (setting == "weekly") {SD = Syn[,,,,d]; MD = Mod[,,,d]}
    if (setting == "daily") {SD = Syn[,,,d]; MD = Mod[,d,]}
    MDlong = c()
    if (setting == "daily") {
      MDlong = MD
    } else {
      for (s in 1:dim(Syn)[length(dim(Syn))-1]) {
        MDlong = cbind(MDlong, c(t(MD[,,s])))
      }
    }
    MDcorr = as.vector(t(cor(MDlong)))
    SDcorr = c()
    for (r in 1:dim(SD)[1]) {
      if (setting == "weekly") {
        SDR = c()
        for (s in 1:dim(Syn)[length(dim(Syn))-1]) {
          SDRS = SD[r,,,s]; SDR = cbind(SDR, c(t(SDRS)))
        }
        temp = as.vector(t(cor(SDR)))
        SDcorr = rbind(SDcorr, temp)
      } else {
        SDR = SD[r,,]
        temp = as.vector(t(cor(SDR)))
        SDcorr = rbind(SDcorr, temp)
      }
    }
    MDcorr = as.matrix(t(MDcorr)); SDcorr = as.matrix(SDcorr)
    
    # for better plotting, take average of synthetic correlations
    # and report variance/sd 
    SDcorrAVG = apply(SDcorr,2,mean)
    #SDcorrSD  = apply(SDcorr,2,sd)
    
    # reorganize correlations into matrix form, but maybe don't use these?
    MDmat = matrix(MDcorr, nrow = 7, byrow = TRUE)
    SDmat = matrix(SDcorrAVG, nrow = 7, byrow = TRUE)
    CorrDiffs = MDmat - SDmat; colnames(CorrDiffs) = SOURCES; rownames(CorrDiffs) = SOURCES
    
    CDm = melt(CorrDiffs)
    CDm$Decade = DECS[d]
    
    tempd = CDm
    tempall = rbind(tempall,tempd)
  }

  SPCORRplot = ggplot() +
    geom_tile(data = tempall, aes(x = Var1, y = Var2, fill = value, group = Decade)) + 
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1)) +
    facet_wrap( ~ Decade, ncol = 3) + xlab("Sites") + ylab("Sites") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1),
                         name = "Difference in\nRHESSys and\nSynthetic\nCorrelation\n") +
    ggtitle(paste(setting, " Spatial Correlation for Scenario ", scen, sep = "")) + coord_fixed()
  ggsave(paste("transientgeneratoroutput/scenario", scen, "/", setting, "flowset", flwst, "SpatialAutocorrelation.png", sep = ""),
         units = "in", width = 7, height = 5, dpi = 1200)
  
  # output difference matrix for this flowset, to be used to create final spatial correlation set of all flowsets
  to_return = tempall; to_return$Flowset = flwst
  return(to_return)
}

## ----------------- boxplots of monthly synthetic and RHESSys comparisons ------------

compareRHESSysHistoricMonthlyFlows = function(Mmonthly, Hmonthly, scen = 1,
                                              sitenames = c("Jordan Lake", "Little River", "Flat River", "Jordan Evap", "Michie Evap"),
                                              printtype = "bysite")
{
  library(ggplot2); library(reshape2)
  
  Hmonthly = HistoricArray
  Mmonthly = ModeledArray
  
  decade = 1 # compare 2010 state
  
  BaseMatrix = array(NA, dim = c((dim(Mmonthly)[1]*dim(Mmonthly)[2] + dim(Hmonthly)[1]*dim(Hmonthly)[2])*dim(Mmonthly)[3], 
                                 4)) # row for every flow observation, col for calendar week, site, decade, model type, obs
  colnames(BaseMatrix) = c("Week", "Site", "Model", "Observations")
  
  Mvec = rep(1:dim(Mmonthly)[2],(dim(Mmonthly)[1]*dim(Mmonthly)[2] + dim(Hmonthly)[1]*dim(Hmonthly)[2])/dim(Mmonthly)[2])
  
  i = 0
  for (site in 1:dim(Hmonthly)[3]) {
    indices = 1:length(Mvec) + length(Mvec)*i
    
    H = c(t(Hmonthly[,,site]))
    M = c(t(Mmonthly[,,site,decade]))
    
    BaseMatrix[indices,4] = as.numeric(as.character(c(H,M))) # observations
    BaseMatrix[indices,3] = c(rep("Historic", length(H)), rep("RHESSys", length(M))) # model type
    BaseMatrix[indices,2] = rep(sitenames[site], length(Mvec)) # site
    BaseMatrix[indices,1] = Mvec # calendar week
    
    i = i + 1
  }
  
  BM = as.data.frame(BaseMatrix); BM$Observations = as.numeric(as.character(BM$Observations)); BM$Week = as.character(BM$Week)
  
  P = ggplot(BM, aes(x = Week, y = Observations, fill = Model, colour = Model), alpha = 0.3) + 
    geom_boxplot(outlier.shape = NA) + scale_x_discrete(limits = 1:12) +
    ylab("Flow (MG)") + facet_wrap(~ Site, nrow = 5, scales = "free_y") + theme_light() + xlab("Month")
  ggsave(paste("transientgeneratoroutput/scenario", scen, "/", 
               if(dim(Mmonthly)[2] == 52) {"Weekly"} else {"Monthly"},
               "BoxplotCompsHistoricVsRHESSys.png", sep = ""),
         units = "in", width = 5, height = 12)
}

## ----------------- boxplots of synthetic and RHESSys comparisons ------------

calendarTrends = function(Mdata, Sdata, scen = 1, 
                          sitenames = c("Jordan Lake", "Little River", "Flat River", "Jordan Evap", "Michie Evap"),
                          printtype = "bysite") # can be "bysite" or "facets"
{
  library(ggplot2); library(reshape2)
  print(paste("Building boxplot sets for ", 
              dim(Mdata)[3], " sites over ", 
              dim(Mdata)[4], " decades at a ", 
              if(dim(Mdata)[2] == 52) {"weekly"} else {"monthly"}, 
              " timestep.", sep = ""))
  
  # ------------ re-organize -------------
  
  BaseMatrix = array(NA, dim = c((dim(Mdata)[1]*dim(Mdata)[2] + dim(Sdata)[1]*dim(Sdata)[2]*dim(Sdata)[3])*dim(Mdata)[3]*dim(Mdata)[4], 
                                 5)) # row for every flow observation, col for calendar week, site, decade, model type, obs
  colnames(BaseMatrix) = c("Week", "Site", "Decade", "Model", "Observations")
  Decades = c(2010,2020,2030,2040,2050,2060)
  
  Wvec = rep(1:dim(Mdata)[2],(dim(Mdata)[1]*dim(Mdata)[2] + dim(Sdata)[1]*dim(Sdata)[2]*dim(Sdata)[3])/dim(Mdata)[2])
  i = 0
  for (site in 1:dim(Mdata)[3]) {
    for (dec in 1:dim(Mdata)[4]) {
      indices = 1:length(Wvec) + length(Wvec)*i
      
      D = c(t(Mdata[,,site,dec])) # timeseries vector of rhessys
      SS = Sdata[,,,site,dec]
      if (dim(SS)[1] > 1) {
        S = c()
        for (n in 1:dim(SS)[1]) {
          SSS = c(t(SS[n,,]))
          S = c(S,SSS)
        }
      } else {
        S = SS # timeseries vector of synthetic
      }
      
      BaseMatrix[indices,5] = as.numeric(as.character(c(D,S))) # observations
      BaseMatrix[indices,4] = c(rep("RHESSys", length(D)), rep("Synthetic", length(S))) # model type
      BaseMatrix[indices,3] = as.character(rep(Decades[dec], length(Wvec))) # decade
      BaseMatrix[indices,2] = rep(sitenames[site], length(Wvec)) # site
      BaseMatrix[indices,1] = Wvec # calendar week
      
      i = i + 1
    }
  }
  BM = as.data.frame(BaseMatrix); BM$Observations = as.numeric(as.character(BM$Observations)); BM$Week = as.character(BM$Week)
  
  if (printtype == "bysite") {
    for (site in 1:length(unique(BM$Site))) {
      #      for (decade in 1:length(unique(BM$Decade))) {
      BMsub = BM[which(BM$Site == unique(BM$Site)[site]),]
      P = ggplot(BMsub, aes(x = Week, y = Observations, fill = Model, colour = Model), alpha = 0.3) + 
        geom_boxplot(outlier.shape = NA) + scale_x_discrete(limits = 1:52) +
        ylab("Flow (MG)") + facet_wrap(Site ~ Decade, nrow = 1) + theme_light() + 
        coord_cartesian(ylim = quantile(BMsub$Observations, c(0.01, 0.99)))
      ggsave(paste("transientgeneratoroutput/scenario", scen, "/", 
                   if(dim(Mdata)[2] == 52) {"Weekly"} else {"Monthly"}, unique(BM$Site)[site],
                   "BoxplotComps.png", sep = ""),
             units = "in", width = 20, height = 5)
      #      }
    }
  } else {
    P = ggplot(BM, aes(x = Week, y = Observations, fill = Model, colour = Model), alpha = 0.3) + geom_boxplot(outlier.size = 0.5) + 
      ylab("Flow (MG)") + facet_grid(Site ~ Decade, scales = "free_y") + theme_light() + scale_x_discrete(limits = 1:52)
    ggsave(paste("transientgeneratoroutput/scenario", scen, "/", 
                 if(dim(Mdata)[2] == 52) {"Weekly"} else {"Monthly"}, 
                 "BoxplotCompsFacets.png", sep = ""),
           units = "in", width = 15, height = 12)
  }
  
}


## ------------------ master ACF data collector ---------------------

masterACFs = function(synpath = "C:/Users/David/OneDrive - University of North Carolina at Chapel Hill/UNC/Research/WSC/Coding/Dynamic Streamflow generator/transientgeneratoroutput",
                      GI = "C:/Users/David/OneDrive - University of North Carolina at Chapel Hill/UNC/Research/WSC/Coding/Dynamic Streamflow generator/inflow-data_updated",
                      type = "monthly", s = 1) {
  if (type == "monthly") {nc = 120} else if (type == "weekly") {nc = 520} else {nc = 3650}
  bigset = matrix(NA, nrow = (1000)*1*6*5, ncol = nc + 4)
  i = 1
  
  if (s %in% c(1,3)) {mp = 9; nr = 250} else {mp = 21; nr = 70}
  if (s == 1) {
    EVAPSOURCES = c("jordan", "michielittle") # JORDAN FOR SCENARIOS 1 and 2
    SOURCES     = c("TotalJordan", "little", "Flat", EVAPSOURCES)
  } else if (s == 2) {
    EVAPSOURCES = c("jordan", "michielittle") # JORDAN FOR SCENARIOS 1 and 2
    SOURCES     = c("totaljordan", "little", "Flat", EVAPSOURCES) # scenario 2
  } else {
    EVAPSOURCES = c("jordanevap", "michielittle") # JORDANEVAP FOR SCENARIOS 3 and 4
    SOURCES     = c("Jordan", "little", "Flat", EVAPSOURCES) # scenario 3 and 4
  }
  
  D = collectModeledInputFlows(type = "monthly", SCENARIO = s, SOURCES = SOURCES, num_timestep = nc/10,
                               DECADES = seq(2010,2060,10), GeneratorInputsFolder = GI, longform = TRUE)
  
  for (l in 1:length(SOURCES)) {
    for (d in 1:6) {
      for (set in 6:mp) {
        temp = read.csv(paste(synpath, "/scenario", s, "/", type, SOURCES[l], "_decade", d, "_set", set, ".csv", sep = ""), 
                        header = FALSE)
        bigset[i:(i+nrow(temp)-1),] = as.matrix(cbind(rep(s, nrow(temp)),
                                                      rep(l, nrow(temp)),
                                                      rep(d, nrow(temp)),
                                                      rep(set, nrow(temp)),
                                                      temp))
        i = i + nrow(temp)
      }
    }
  }
  
  ## ---------- get autocorrelations -----------
  
  AllMods = matrix(NA, nrow = nc/10*1*6*5*length(6:mp), ncol = 6); j = 1
  AllSyns = matrix(NA, nrow = nc/10*1*6*5*nr*length(6:mp), ncol = 7); k = 1
  for (l in 1:length(SOURCES)) {
    for (d in 1:6) {
      temp = D[,l,d,]
      temp3 = c()
      temp4 = c()
      for (set in 6:mp) {
        things = matrix(NA, nrow = )
        for (r in 1:nr) {
          temp2 = bigset[which(bigset[,1] == s &
                                 bigset[,2] == l &
                                 bigset[,3] == d &
                                 bigset[,4] == set),5:ncol(bigset)]
          
          
          temp3 = rbind(temp3, as.vector(acf(temp2[r,], lag.max = nc/10 - 1, plot = FALSE)$acf))
        }
        
        temp4 = rbind(temp4, as.vector(acf(temp[,set-5], lag.max = nc/10 - 1, plot = FALSE)$acf))
        
      }
      
      SynacrMeans = c(1); SynacrCI = c(1,1); ModacrMeans = c(1); ModacrCI = c(1,1)
      for (l in 2:(nc/10)) {
        SynacrMeans[l] = t.test(as.vector(temp3[,l]))$estimate
        SynacrCI = rbind(SynacrCI, t.test(as.vector(temp3[,l]))$conf.int[1:2])
        
        ModacrMeans[l] = t.test(as.vector(temp4[,l]))$estimate
        ModacrCI = rbind(ModacrCI, t.test(as.vector(temp4[,l]))$conf.int[1:2])
      }
      
      MSDacr = as.data.frame(as.vector(acf(c(t(MSD)), lag.max = nlags, plot = FALSE)$acf))
      MSDacr$LB = MSDacr[,1]; MSDacr$UB = MSDacr[,1]; MSDacr$Type = "RHESSys"; colnames(MSDacr)[1] = "Mean Autocorrelation"
      MSDacr$Lag = 0:nlags
    }
  }
  
  # -------------- build the plot... --------------
  
  cn = c("Scenario", "Site", "Decade", "Lag", "Autocorrelation", "Type")
  
  Mod = as.data.frame(AllMods[,c(1:3,5:6)])
  Mod$Type = "RHESSys"
  
  Syn = as.data.frame(AllSyns[,c(1:3,6:7)])
  Syn$Type = "Synthetic"
  
  library(reshape2); library(ggplot2)
  allD = rbind(Mod, Syn); colnames(allD) = cn
  
  MSDacr = as.data.frame(as.vector(acf(c(t(MSD)), lag.max = nlags, plot = FALSE)$acf))
  MSDacr$LB = MSDacr[,1]; MSDacr$UB = MSDacr[,1]; MSDacr$Type = "RHESSys"; colnames(MSDacr)[1] = "Mean Autocorrelation"
  MSDacr$Lag = 0:nlags
  
  ACFplot = ggplot() +
    geom_errorbar(data = SynSet, aes(x = Lag, ymin = LB, ymax = UB), 
                  color = "blue", size = 1) + 
    geom_line(data = SynSet, aes(x = Lag, y = `Mean Autocorrelation`), 
              color = "blue", size = 1) + 
    geom_line(data = ModSet, aes(x = Lag, y = `Mean Autocorrelation`), 
              color = "red", size = 0.75) + 
    facet_grid(Site ~ Decade) +
    ggtitle(paste(setting, " Temporal Autocorrelation for Scenario ", scen, sep = ""))
  ggsave(paste("transientgeneratoroutput/scenario", scen, "/", setting, "flowset", flwst, "TemporalAutocorrelation.png", sep = ""),
         units = "in", width = 8, height = 6)
}
