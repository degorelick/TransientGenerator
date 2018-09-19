## Required functions for transient streamflow generator
##  D Gorelick (Oct 2017)
## -------------------------------------------------------------

## ------------------ collect all input flows --------------------------------

collectModeledInputFlows = function(num_years = 59, num_timestep = 52, type = "weekly", longform = TRUE,
                                      SCENARIO = 1, SOURCES = c("Jordan"), DECADES = 2010, GeneratorInputsFolder)
{
  library(lubridate)
  
  if (SCENARIO %in% c(1,3)) {maxflowsets = 9} else {maxflowsets = 21}
  
  if (longform) {
    DataHolder = array(NA, dim = c(num_years*num_timestep,length(SOURCES),length(DECADES), length(6:maxflowsets))) 
    # 60 years * 52 weeks by 3 sources by 6 decades by number of RHESSys records
  } else {
    DataHolder = array(NA, dim = c(num_years,num_timestep,length(SOURCES),length(DECADES), length(6:maxflowsets))) 
      # 60 years by 52 weeks by 3 sources by 6 decades by number of RHESSys records
  }
  
  n_source = 1
  for (SOURCE in SOURCES) {
    
    n_decade = 1
    for (DECADE in DECADES) {
      
      ## Reading data from source, based on scenario and decade
      ## if the source is an evaporation dataset, do not convert units to MG as it already is in MG/acre
      convert.units = TRUE
      if (SOURCE %in% c("jordanevap", "michielittle", "owasaevap")) {
        D = read.csv(paste(GeneratorInputsFolder, "/scenario", SCENARIO, "/evap/", SOURCE, "/netevap_", DECADE, ".csv", sep = ""))
        convert.units = FALSE
      } else {
        D = read.csv(paste(GeneratorInputsFolder, "/scenario", SCENARIO, "/", SOURCE, "_", DECADE, ".csv", sep = ""))
      }
      
      OrderedDate = ymd(paste(D$year, D$month, D$day))

      if (type == "weekly") {
        WeekOfYear = week(OrderedDate)

        # set week of year for each daily observation, set to 52 weeks per year
        WeeklyD = aggregate(D, by = list(D$year, WeekOfYear), FUN = sum)
        WeeklyD = WeeklyD[which(WeeklyD$Group.2 != 53),] # remove week 53 records
        
        MD = WeeklyD[which(WeeklyD$Group.1 != min(WeeklyD$Group.1) & WeeklyD$Group.1 != max(WeeklyD$Group.1)),]
        
        if (convert.units) {
          MDflowsonly = MD[order(MD$Group.1),6:maxflowsets] * 264.172/1000000 # convert to MGperMonth, sorted by year THEN MONTH
        } else {
          MDflowsonly = MD[order(MD$Group.1),6:maxflowsets] # sorted by year THEN MONTH
        }
        
        for (fst in 1:ncol(MDflowsonly)) {
          if (longform) {
            WeeklyMatrix = MDflowsonly[1:(dim(DataHolder)[1]),fst] # data organized into vector of weekly flows
            DataHolder[, n_source, n_decade, fst] = WeeklyMatrix # add to data holder, trim last row for even 60 year record 
          } else {
            WeeklyMatrix = matrix(MDflowsonly[,fst], ncol = num_timestep, byrow = F) # data organized into 60x52 matrix of weekly flows
            DataHolder[, , n_source, n_decade, fst] = WeeklyMatrix[1:num_years,] # add to data holder, trim last row for even 60 year record 
          }
        }
      } else if (type == "monthly") {
        MonthOfYear = month(OrderedDate)
        MonthOfYear[MonthOfYear > num_timestep] = num_timestep # set month of year for each daily observation, set to 12 months per year
        MonthlyD = aggregate(D, by = list(D$year, MonthOfYear), FUN = sum) # aggregate monthly 
        
        MD = MonthlyD[which(MonthlyD$Group.1 != min(MonthlyD$Group.1) & MonthlyD$Group.1 != max(MonthlyD$Group.1)),]
        if (convert.units) {
          MDflowsonly = MD[order(MD$Group.1),6:maxflowsets] * 264.172/1000000 # convert to MGperMonth, sorted by year THEN MONTH
        } else {
          MDflowsonly = MD[order(MD$Group.1),6:maxflowsets] # sorted by year THEN MONTH
        }
        
        for (fst in 1:ncol(MDflowsonly)) {
          if (longform) {
            WeeklyMatrix = MDflowsonly[1:(dim(DataHolder)[1]),fst] # data organized into vector of weekly flows
            DataHolder[, n_source, n_decade, fst] = WeeklyMatrix # add to data holder, trim last row for even 60 year record 
          } else {
            WeeklyMatrix = matrix(MDflowsonly[,fst], ncol = num_timestep, byrow = F) # data organized into 60x52 matrix of weekly flows
            DataHolder[, , n_source, n_decade, fst] = WeeklyMatrix[1:num_years,] # add to data holder, trim last row for even 60 year record 
          }
        }
      } else {
        DayOfYear = day(OrderedDate)

        WeeklyD = aggregate(D, by = list(D$year, DayOfYear), FUN = sum)
        WeeklyD = WeeklyD[which(WeeklyD$Group.2 > 365),] 
        
        MD = WeeklyD[which(WeeklyD$Group.1 != min(WeeklyD$Group.1) & WeeklyD$Group.1 != max(WeeklyD$Group.1)),]
        
        if (convert.units) {
          MDflowsonly = MD[order(MD$Group.1),6:maxflowsets] * 264.172/1000000 # convert to MGperMonth, sorted by year THEN MONTH
        } else {
          MDflowsonly = MD[order(MD$Group.1),6:maxflowsets] # sorted by year THEN MONTH
        }
        
        for (fst in 1:ncol(MDflowsonly)) {
          if (longform) {
            WeeklyMatrix = MDflowsonly[1:(dim(DataHolder)[1]),fst] # data organized into vector of weekly flows
            DataHolder[, n_source, n_decade, fst] = WeeklyMatrix # add to data holder, trim last row for even 60 year record 
          } else {
            WeeklyMatrix = matrix(MDflowsonly[,fst], ncol = num_timestep, byrow = F) # data organized into 60x52 matrix of weekly flows
            DataHolder[, , n_source, n_decade, fst] = WeeklyMatrix[1:num_years,] # add to data holder, trim last row for even 60 year record 
          }
        }
      }
      
      n_decade = n_decade + 1
    }
    
    n_source = n_source + 1
  }
  
  return(DataHolder)
}

## ------------------ aggregate input flows ----------------------------------

aggregateModeledInputFlows = function(num_years = 59, num_timestep = 12, aggregation = "months", 
                                      SCENARIO = 1, SOURCES = c("TotalJordan"), DECADES = 2010, GeneratorInputsFolder,
                                      chosenflows = 6)
{
  DataHolder = array(NA, dim = c(num_years,num_timestep,length(SOURCES),length(DECADES))) 
  # 60 years by 52 weeks by 3 sources by 6 decades
  DailyRecordHolder = array(NA, dim = c(365*num_years,length(DECADES),length(SOURCES))) 
  # 60 years by 52 weeks by 3 sources by 6 decades
  
  SourceIndexHolder = list()
  SourceShiftedDataHolder = list() 
  
  n_source = 1
  for (SOURCE in SOURCES) {
    
    DecadalIndexHolder = list()
    DecadalShiftedDataHolder = list()
    
    n_decade = 1
    for (DECADE in DECADES) {
      
      ## Reading data from different scenarios
      convert.units = TRUE
      if (SOURCE %in% c("jordanevap", "michielittle", "owasaevap")) {
        D = read.csv(paste(GeneratorInputsFolder, "/scenario", SCENARIO, "/evap/", SOURCE, "/netevap_", DECADE, ".csv", sep = ""))
        convert.units = FALSE
      } else {
        D = read.csv(paste(GeneratorInputsFolder, "/scenario", SCENARIO, "/", SOURCE, "_", DECADE, ".csv", sep = ""))
      }
      
      OrderedDate = ymd(paste(D$year, D$month, D$day))
      if (aggregation == "weeks") {
        WeekOfYear = week(OrderedDate)
        WeekOfYear[WeekOfYear > num_timestep] = num_timestep # set week of year for each daily observation, set to 52 weeks per year
        if (convert.units) {
          WeeklyD = aggregate(D, by = list(D$year, WeekOfYear), FUN = sum)*264.172/1000000
            # aggregate weekly and convert to MG
        } else {
          WeeklyD = aggregate(D, by = list(D$year, WeekOfYear), FUN = sum)
            # aggregate weekly
        }
        WeeklyMatrix = matrix(WeeklyD[,chosenflows], ncol = num_weeks, byrow = F) # data organized into 60x52 matrix of weekly flows
        DataHolder[, , n_source, n_decade] = WeeklyMatrix[1:num_years,] # add to data holder, trim last row for even 60 year record each time
      } else if (aggregation == "months") {
        # THIS IS JUST SET UP FOR ONE RHESSYS SIMULATION TO BE READ
        # remember RHESSys goes by water year, need to cut first two months and last 10 of record to have Jan-Dec order
        MonthOfYear = month(OrderedDate)
        MonthOfYear[MonthOfYear > num_timestep] = num_timestep # set month of year for each daily observation, set to 12 months per year
        MonthlyD = aggregate(D, by = list(D$year, MonthOfYear), FUN = sum) # aggregate monthly 
        
        ## BECAUSE RECORD GOES FROM OCT-OCT, REMOVE FIRST YEAR AND LAST YEAR RECORDS TO SHIFT TO JAN-JAN
        MD = MonthlyD[which(MonthlyD$Group.1 != min(MonthlyD$Group.1) & MonthlyD$Group.1 != max(MonthlyD$Group.1)),]
        if (convert.units) {
          MDflowsonly = MD[order(MD$Group.1),chosenflows] * 264.172/1000000 # convert to MGperMonth, sorted by year THEN MONTH
        } else {
          MDflowsonly = MD[order(MD$Group.1),chosenflows] # sorted by year THEN MONTH
        }
        
        MonthlyMatrix = matrix(MDflowsonly, ncol = num_timestep, byrow = TRUE) # data organized into 60x12 matrix of weekly flows
        
        if(length(MDflowsonly)/num_timestep < num_years) {stop(paste("Modeled record only has ", 
                                                                     length(MDflowsonly)/num_timestep, " years of data", sep = ""))}
        
        if (!convert.units) {MonthlyMatrix = exp(MonthlyMatrix)} # exponentiate evap data
        DataHolder[, , n_source, n_decade] = MonthlyMatrix[1:num_years,] # add to data holder, trim last row for even 60 year record each time
        
        ## create shifted record for kNN calculations later ##
        #   THIS IS FOR A SINGLE SITE AND DECADE             #
        ## ------------------------------------------------ ##
        
        FlowRecord = D[order(D$year, D$month, D$day),c(1:3,chosenflows-2)] # single timeseries of daily flows
        FlowRecord = FlowRecord[which(FlowRecord$year != min(FlowRecord$year) & FlowRecord$year != max(FlowRecord$year)),]
          # remove end years to convert timeseries to Jan-Dec form
        if (convert.units) {
          Days365_Record = FlowRecord[-which(D$day == 29 & D$month == 2),4]*264.172/1000000
        } else {
          Days365_Record = FlowRecord[-which(D$day == 29 & D$month == 2),4]
        }
        Days365_Record = Days365_Record[1:(365*num_years)] # ensure proper number of years of data
        # matrix of single timeseries w/o leap days
        
        #if (!convert.units) {Days365_Record = exp(Days365_Record)} # exponentiate evap data
        DailyRecordHolder[,n_decade,n_source] = Days365_Record # store for use later
        
        ExtendedRecord = c(Days365_Record[(length(Days365_Record)-7):length(Days365_Record)], Days365_Record, Days365_Record[1:8])
        # add at least 7 days to each end of record
        
        IndexHolder = list()
        ShiftedDataSet = list()
        for (calendarmonth in 1:12) {

          nTotals = (floor(length(Days365_Record)/365)) * 15
          
          Qmonthly_shifted = matrix(NA, nrow = nTotals, ncol = 1)
          Indices = matrix(NA, nrow = nTotals, 2)
          for (shift in 1:15) {
            shiftedRecord = ExtendedRecord[shift:(length(Days365_Record) + shift - 1)]
            shiftedmonthlyRecord = convert_data_to_monthly(shiftedRecord)
            
            if (nrow(shiftedmonthlyRecord) != num_years) {
              print(paste("Calendar Month: ", calendarmonth, ", Shift: ", shift, sep = ""))
              print(paste("Rows in Shifted Monthly Record: ", nrow(shiftedmonthlyRecord), sep = ""))
              stop("Mismatch in shifted record years and desired years")
            }
            
            shiftedmonthlyRecord = shiftedmonthlyRecord[,calendarmonth] # should be number of years in length
            Qmonthly_shifted[((shift-1)*length(shiftedmonthlyRecord)+1):(length(shiftedmonthlyRecord)*shift),1] = shiftedmonthlyRecord
            
            Indices[((shift-1)*length(shiftedmonthlyRecord)+1):(length(shiftedmonthlyRecord)*shift),1] = 1:length(shiftedmonthlyRecord)
            Indices[((shift-1)*length(shiftedmonthlyRecord)+1):(length(shiftedmonthlyRecord)*shift),2] = shift
            
          }
          IndexHolder[[calendarmonth]] = Indices
          ShiftedDataSet[[calendarmonth]] = Qmonthly_shifted
        }
        DecadalIndexHolder[[n_decade]] = IndexHolder
        DecadalShiftedDataHolder[[n_decade]] = ShiftedDataSet
      }
      n_decade = n_decade + 1
    }
    SourceIndexHolder[[n_source]] = DecadalIndexHolder
    SourceShiftedDataHolder[[n_source]] = DecadalShiftedDataHolder
    
    n_source = n_source + 1
  }
  
  Output = list(DataHolder, SourceShiftedDataHolder, SourceIndexHolder, DailyRecordHolder)
  return(Output)
}


## ------------------ convert daily data from 365 days to monthly ------------

convert_data_to_monthly = function(Qt) 
{
  Nyears = length(Qt)/365; Nmonths = 12; DaysPerMonth = c(31,28,31,30,31,30,31,31,30,31,30,31)
  Qmonthly = matrix(0, nrow = Nyears, ncol = Nmonths)
  
  for (year in 1:Nyears) {
    for (month in 1:Nmonths) {
      if (month == 1) {
        start = (year-1)*365 + 1
      } else {
        start = (year-1)*365 + sum(DaysPerMonth[1:(month-1)])+1
      }
      Qmonthly[year,month] = sum(Qt[start:(start+DaysPerMonth[month]-1)])
    }
  }
  
  return(Qmonthly)
}

## ------------------ aggregate observed flows -------------------------------

aggregateObservedFlows = function(num_years = 80, num_timesteps = 12, SOURCES, aggregation = "months",
                                  GeneratorInputsFolder, ObservedRecordFolder)
{
  DataHolder = array(NA, dim = c(num_years,num_timesteps,length(SOURCES))) 
  # 80 years by (52 weeks or 12 months) by 3 sources
  n_source = 1
  for (SOURCE in SOURCES) {
    filename = paste(GeneratorInputsFolder, ObservedRecordFolder, "/updated", SOURCE, "Inflow", ".csv", sep = "")
    if (file.exists(filename)) {
      D = read.csv(paste(GeneratorInputsFolder, ObservedRecordFolder, "/updated", SOURCE, "Inflow", ".csv", sep = ""), header = F)
      evapdata = FALSE
    } else {
      D = read.csv(paste(GeneratorInputsFolder, ObservedRecordFolder, "/", SOURCE, "Evap", ".csv", sep = ""), header = F)
      evapdata = TRUE
        # assume evaporation data
    }
        # historic data in weeks
    if (aggregation == "weeks") {
      FINALMATRIX = D
      #if (evapdata) {FINALMATRIX = exp(FINALMATRIX)}
    } else if (aggregation == "months") {
      library(reshape2); D$YR = 1926:(1926+nrow(D)-1); colnames(D)[1:52] = 1:52
      weeklytimeseries = melt(D, id = "YR"); colnames(weeklytimeseries)[2:3] = c("WK", "FLOW")
      weeklytimeseries = weeklytimeseries[order(weeklytimeseries$YR),]
      dailytimeseries = c(); DAY = c(); YEAR = c()
      for (yr in sort(unique(weeklytimeseries$YR))) 
      {
        DAY  = c(DAY,1:365)
        YEAR = c(YEAR, rep(yr,365))
        for (wk in 1:52) 
        {
          if (wk == 52) {
            dailytimeseries = c(dailytimeseries, rep(weeklytimeseries[which(weeklytimeseries$YR == yr & weeklytimeseries$WK == wk),"FLOW"]/8,8))
          } else {
            dailytimeseries = c(dailytimeseries, rep(weeklytimeseries[which(weeklytimeseries$YR == yr & weeklytimeseries$WK == wk),"FLOW"]/7,7))
          }
        }
      }
      #if (evapdata) {dailytimeseries = exp(dailytimeseries)}
      
      Date = strptime(paste(YEAR, DAY), format = "%Y %j")
      MonthOfYear = month(Date)
      monthlytimeseries = aggregate(dailytimeseries, by = list(YEAR, MonthOfYear), sum); colnames(monthlytimeseries) = c("YR", "MON", "FLOW")
      
      FINALMATRIX = matrix(monthlytimeseries$FLOW, nrow = nrow(D)) # aggregated flows by calendar month
      if (evapdata) {FINALMATRIX = exp(FINALMATRIX)}
    }
    
    DataHolder[, , n_source] = FINALMATRIX[1:num_years,] # add to data holder, trim last row for even 80 year record each time
    n_source = n_source + 1
  }
  
  return(DataHolder)
}

## ------------------ k nearest neighbor identification ----------------------

## Disaggregate monthly flows to weekly ##
#   this approach uses a k-nearest       #
#   neighbor method to map monthly       #
#   trends to member weeks               #
## ------------------------------------ ##

kNNidentification =  function(SyntheticTimeseries, ShiftedTimeseries, month, decade, syntheticmonth, realization = 1)
  # month and decade input arguments are values 1-12 and 1-6, respectively.
  # identifies kNNs for SINGLE SYNTHETIC MONTH AT EACH SITE, SyntheticTimeseries is total set,
  #   SM is what is indexed to extract only a single synthetic month's flows at each site
{
  SM = SyntheticTimeseries[realization,,syntheticmonth,,decade] 
    # synthetic monthly flows at each site for given synthetic month 
    # if evaporation, should be exponentiated at this point
  HT = c() # observed shifted flows for a calendar month at every site
  for (site in 1:length(SM)) {
    SiteData = ShiftedTimeseries[[site]][[decade]][[month]] 
      # observed flows from given calendar month at one site
      # these flows should be raw values, not exponentiated
    HT = cbind(HT, SiteData)
  }
  
  ## 1. Determine k -------------------------------------------------------------------------------
  k = floor(sqrt(nrow(HT)))
  
  ## 2. Identify k nearest neighbors for each month, across each site -----------------------------
  kNNdist = vector("numeric", length = nrow(HT))
  for (hist_month in 1:nrow(HT)) {
    for (site in 1:length(SM)) {
      kNNdist[hist_month] = kNNdist[hist_month] + (HT[hist_month,site] - SM[site])^2
        # month of record with lowest cumulative "distance" between historic and synthetic data at each site
    }
  }
  kNNdist = sqrt(kNNdist)
  
  kNNdistMat = data.frame(kNNid = 1:nrow(HT), kNNdist)
  OrderedkNN = kNNdistMat[order(kNNdistMat$kNNdist),]
  
  kNNids = OrderedkNN$kNNid[1:k] # k nearest neighbor index IDs
  
  ## 3. Identify ID weights -----------------------------------------------------------------------
  invW = 1/c(1:k)
  W = invW/sum(invW)
  
  ## 4. Return kNN IDs and weights ----------------------------------------------------------------
  Output = cbind(kNNids, W)
  return(Output)
}

## ------------------ k nearest neighbor selection ---------------------------

kNNselection = function(identificationOutput, indicesList, HistoricalDailyArray, month, decade, syntheticmonth)
{
  #identificationOutput = kNNinfo
  #indicesList = ShiftedIndex
  #HistoricalDailyArray = DailyArray
  
  # 1. Locate required inputs ---------------------------------------------------------------------
  
  kNNids = identificationOutput[,1]; CumWeights = cumsum(identificationOutput[,2])
  Indices = c() 
  for (site in 1:length(indicesList)) {
    SiteData = indicesList[[site]][[decade]][[month]] 
    Indices = cbind(Indices, SiteData)
  }
  DailyFlows = HistoricalDailyArray[,decade,]
    # if evaporation, this is NOT exponentiated
  
  # 2. Select one kNN -----------------------------------------------------------------------------
  
  r = runif(1) # random number between 0 and 1
  Weights = c(0,CumWeights)
  for (i in 1:length(CumWeights)) {if (r > Weights[i] & r <= Weights[i+1]) {kNNid = i}}
  kNNyearID = kNNids[kNNid] # year ID of kNN selected
  #print(paste("Chosen k nearest neighbor year: ", kNNyearID, sep = ""))
  
  # 3. Locate daily data from selected NN year ----------------------------------------------------
  
  ExtendedRecord = rbind(DailyFlows[(nrow(DailyFlows)-7):nrow(DailyFlows),], DailyFlows, DailyFlows[1:8,])
  # add at least 7 days to each end of record
  YearSelected  = Indices[kNNyearID,1]
  ShiftSelected = Indices[kNNyearID,2]
  
  shifted_DailyFlows = ExtendedRecord[ShiftSelected:(ShiftSelected+nrow(DailyFlows)-1),]
  
  DaysPerMonth = c(31,28,31,30,31,30,31,31,30,31,30,31)
  SynYear = YearSelected - 1
  CalMonth = (syntheticmonth-1) %% 12 + 1
  DayRange = (SynYear*365 + sum(DaysPerMonth[1:CalMonth]) - DaysPerMonth[CalMonth] + 1):(SynYear*365 + sum(DaysPerMonth[1:CalMonth]))
  
  #DayRange = ((YearSelected-1)*365 + sum(DaysPerMonth[1:month]) - DaysPerMonth[month] + 1):
  #  ((YearSelected-1)*365 + sum(DaysPerMonth[1:month]))
  
  #start = 365*(YearSelected-1) + sum(DaysPerMonth[1:(month-1)])+1
  dailyFlows = shifted_DailyFlows[DayRange,]
  #dailyFlows[,4:5] = exp(dailyFlows[,4:5]) # exponentiate evap
  
  # 4. Normalize daily data to proportion of total monthly flow -----------------------------------
  
   dailySign = dailyFlows
   dailySign[dailySign < 0] = -1; dailySign[dailySign >= 0] = 1
   
   dailyFlows = abs(dailyFlows) 
    # in order to deal with evaporation, take abs value before determining proportions
    # this should not impact streamflows, as they have no negative values
    # for evaporation, should the sign be flipped after disaggregation? yes
    # otherwise negative values (rain) is turned into evap, and vice versa
  
  DailyProportions = t(t(dailyFlows*dailySign)/colSums(dailyFlows)) # chosen month's daily flow proportions for each site
  return(DailyProportions)
}

## ------------------ build transient generator ------------------------------

transienttrajectorygenerator = function(ModeledArray, HistoricArray, num_synthetic_years = 59, IS.CONCURRENT = TRUE,
                                        aggregation = "months", num_synthetic_realizations = 1, num_generated_years = 10,
                                        use.modeled.corr = FALSE) 
{
  ## Initialize necessary data holders, etc. ##
  ## --------------------------------------- ##
  
  if (aggregation == "weeks") {
    dim_timestep = 52
  } else if (aggregation == "months") {
    dim_timestep = 12
  } else {
    stop("Data aggregation not recognized - should be weeks or months")
  }
  
  ModeledResiduals  = as.array(ModeledArray)
  HistoricResiduals = as.array(HistoricArray)
  
  LogMeanHolder = array(NA, dim = c(dim(ModeledArray)[3], dim(ModeledArray)[4], dim_timestep))
  LogSDHolder   = array(NA, dim = c(dim(ModeledArray)[3], dim(ModeledArray)[4], dim_timestep))
  
  ## Create vector of random numbers that determine the bootstrap sample index to be taken in each week  ##
  #   THIS RANDOM NUMBER STAYS THE SAME ACROSS ALL SITES FOR EACH DECADE DURING A GIVEN                   #
  #   REALIZATION OF THE GENERATOR                                                                        #
  # ---------------------------------------------------------------------------------------------------- ##
  
  if (IS.CONCURRENT) {RandomSeeds = matrix(floor(runif(ncol(ModeledArray[,,1,1])*num_synthetic_realizations*(num_generated_years+1), 
                                                       min = 1, max = nrow(ModeledArray[,,1,1]))), 
                                           nrow = (num_generated_years+1)*num_synthetic_realizations)}
  
  ResidualDraws = array(NA, dim = c(nrow(RandomSeeds), ncol(RandomSeeds), dim(ModeledArray)[3], dim(ModeledArray)[4]))
  SyntheticFlowTimeseries = array(NA, dim = c(num_synthetic_realizations, num_generated_years*dim_timestep, 
                                              dim(ModeledArray)[3], dim(ModeledArray)[4]))
  
  if (dim_timestep != ncol(RandomSeeds)) {stop("Mismatch in dimensions of random sampling matrix")}
  
  ## Remove seasonal trends and generate residual records ##
  #   assumes input arrays are aggregated to the same      #
  #   temporal scale (i.e. weeks, months, etc.)            #
  ## ---------------------------------------------------- ##
  
  for (site in 1:dim(HistoricArray)[3]) {
    for (dec in 1:dim(ModeledArray)[4]) {
      
      # Build dataset of residuals and store values for later back-transformation
      LogModeled  = log(ModeledArray[,,site,dec])
      LogHistoric = log(HistoricArray[,,site])
      
      LogModeledMeans  = apply(LogModeled, MARGIN = 2, mean)
      LogModeledSD     = apply(LogModeled, MARGIN = 2, sd)
      LogHistoricMeans = apply(LogHistoric, MARGIN = 2, mean)
      LogHistoricSD    = apply(LogHistoric, MARGIN = 2, sd)
      
      ModeledResiduals[,,site,dec] = t((t(LogModeled) - LogModeledMeans) / LogModeledSD)
      HistoricResiduals[,,site] = t((t(LogHistoric) - LogHistoricMeans) / LogHistoricSD)
      
      temp = as.data.frame(ModeledResiduals[,,site,dec]) # build temporary data frame of modeled data
      ResidualDraws[,,site,dec] = matrix(temp[cbind(c(t(RandomSeeds)),
                                                    rep(1:dim_timestep,(num_generated_years+1)*num_synthetic_realizations))],
                                         nrow = (num_generated_years+1)*num_synthetic_realizations, byrow = TRUE)
      
      temp = as.vector(t(ResidualDraws[,,site,dec])) # timeseries vector of modeled draws
      ShiftedResidualDraws = matrix(temp[ceiling(dim_timestep/2+1):(length(temp)-floor(dim_timestep/2))], 
                                    ncol = 12, byrow = TRUE)
      
      if (use.modeled.corr) {
        temp = as.vector(t(ModeledResiduals[,,site,dec]))
        ShiftedModeledResiduals = matrix(temp[ceiling(dim_timestep/2+1):(length(temp)-floor(dim_timestep/2))], 
                                         ncol = 12, byrow = TRUE)
        
        ModeledResidualCorrelation = cor(ModeledResiduals[,,site,dec])
        ShiftedModeledResidualCorrelation = cor(ShiftedModeledResiduals)
        CholeskyDecomp  = chol(ModeledResidualCorrelation)
        CholeskyDecompShift = chol(ShiftedModeledResidualCorrelation) # cholesky decomposition of correlation
      } else {
        temp = as.vector(t(HistoricResiduals[,,site])) # timeseries vector of historic residuals
        ShiftedHistoricResiduals = matrix(temp[ceiling(dim_timestep/2+1):(length(temp)-floor(dim_timestep/2))], 
                                          ncol = 12, byrow = TRUE)
        
        HistoricResidualCorrelation = cor(HistoricResiduals[,,site])
        ShiftedHistoricResidualCorrelation = cor(ShiftedHistoricResiduals)
        CholeskyDecomp  = chol(HistoricResidualCorrelation)
        CholeskyDecompShift = chol(ShiftedHistoricResidualCorrelation) # cholesky decomposition of correlation
      }
      
      ModeledAutoCorrelatedSyntheticResiduals = ResidualDraws[,,site,dec] %*% CholeskyDecomp
      ShiftedModeledAutoCorrelatedSyntheticResiduals = ShiftedResidualDraws %*% CholeskyDecompShift
      
      FinalSyntheticResiduals = cbind(ShiftedModeledAutoCorrelatedSyntheticResiduals[1:(num_generated_years*num_synthetic_realizations),
                                                                                     ceiling(dim_timestep/2+1):dim_timestep],
                                      ModeledAutoCorrelatedSyntheticResiduals[2:(num_generated_years*num_synthetic_realizations+1),
                                                                              ceiling(dim_timestep/2+1):dim_timestep])
        # concatenate last 6 months of columns from each auto-correlated, synthetic residual record
      
      ## back-convert residuals into real streamflow, separate into 10-year vectors for each synthetic realization
      SyntheticFlowTimeseries[,,site,dec] = matrix(exp(t(FinalSyntheticResiduals) * LogModeledSD + LogModeledMeans), 
                                                   nrow = num_synthetic_realizations)
      if (site > 4) {
        SyntheticFlowTimeseries[,,site,dec] = log(SyntheticFlowTimeseries[,,site,dec])
          # transform again, back into net evap units from exponentiated ones
      }
    }
  }
  return(SyntheticFlowTimeseries)
}

## ----------------- merge sets of synthetic records created ----------------

makemastersets = function(scen = 1, sets = 4, startset = 6, SRCE = "Little", 
                          realstotake = 1, keep.hist.yrs = 45, keep.syn.yrs = 51) {
  allD = c()
  for (s in 1:sets) {
    D = read.csv(paste("transientgeneratoroutput/", "scenario", scen, "/", 
                       SRCE, "set", startset+s-1, ".csv", sep = ""),
                 header = FALSE)
    D = D[,1:(keep.syn.yrs*52)]
    
    if (SRCE == "Little") {
      histD = read.csv(paste("inflow-data_updated/historic/updatedLittleRiverInflow.csv", sep = ""),
                       header = FALSE)
    } else if (SRCE == "Flat") {
      histD = read.csv(paste("inflow-data_updated/historic/updatedMichieInflow.csv", sep = ""),
                       header = FALSE)
    } else if (SRCE == "Jordan") {
      histD = read.csv(paste("inflow-data_updated/historic/updatedJordanLakeInflow.csv", sep = ""),
                       header = FALSE)
    } else {
      histD = read.csv(paste("inflow-data_updated/historic/updatedOWASAInflow.csv", sep = ""),
                       header = FALSE)
    }
    histD = as.matrix(histD); histD = c(t(histD))
    
    keepweeks = keep.hist.yrs * 52 + 7 # account for waterpaths start week 
    
    histD = histD[(length(histD)-keepweeks+1):length(histD)]
    
    histD = t(replicate(realstotake,histD))
    D = cbind(histD,D[1:realstotake,]) # realstotake must be > 0
    
    allD = rbind(allD,D)
  }
  
  write.table(allD, sep = ",", 
              paste("transientgeneratoroutput/", "scenario", scen, "/", SRCE, "_fullset.csv", sep = ""), 
              row.names = FALSE, col.names = FALSE)
  write.table(allD, sep = ",", 
              paste("C:/Users/David/OneDrive - University of North Carolina at Chapel Hill/UNC/Research/WSC/Coding/WaterPaths/TestFiles/inflows/", 
                    "scenario", scen, "/", SRCE, "_fullset.csv", sep = ""), 
              row.names = FALSE, col.names = FALSE)
}

## ----------------- merge sets of synthetic records created ----------------

makemasterevapsets = function(scen = 1, sets = 4, startset = 6, SRCE = "little", 
                              realstotake = 1, keep.hist.yrs = 45, keep.syn.yrs = 51) {
  allD = c()
  for (s in 1:sets) {
    D = read.csv(paste("transientgeneratoroutput/", "scenario", scen, "/", 
                       SRCE, "set", startset+s-1, ".csv", sep = ""),
                 header = FALSE)
    D = D[,1:(keep.syn.yrs*52)]
    
    if (SRCE == "michielittle") {
      histD = read.csv(paste("inflow-data_updated/historic/michielittleEvap.csv", sep = ""),
                       header = FALSE)
    } else if (SRCE == "jordanevap") {
      histD = read.csv(paste("inflow-data_updated/historic/jordanevapEvap.csv", sep = ""),
                       header = FALSE)
    } else {
      histD = read.csv(paste("inflow-data_updated/historic/owasaevapEvap.csv", sep = ""),
                       header = FALSE)
    }
    histD = as.matrix(histD); histD = c(t(histD))
    
    keepweeks = keep.hist.yrs * 52 + 7
    histD = histD[(length(histD)-keepweeks+1):length(histD)]
    
    histD = t(replicate(realstotake,histD))
    D = cbind(histD,D[1:realstotake,])
    
    allD = rbind(allD,D)
  }
  
  write.table(allD, sep = ",", 
              paste("transientgeneratoroutput/", "scenario", scen, "/", SRCE, "_netEvap_fullset.csv", sep = ""), 
              row.names = FALSE, col.names = FALSE)
  write.table(allD, sep = ",", 
              paste("C:/Users/David/OneDrive - University of North Carolina at Chapel Hill/UNC/Research/WSC/Coding/WaterPaths/TestFiles/evaporation/", 
                    "scenario", scen, "/", SRCE, "_netEvap_fullset.csv", sep = ""), 
              row.names = FALSE, col.names = FALSE)
}
