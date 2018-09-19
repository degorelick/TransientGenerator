## Transient Synthetic Streamflow Generator
## David Gorelick (Oct 2017)
##
##  Accept 60-year records representing decadal flow regimes from 2010-2060
##  from RHESSys and concatenate synthetic records decadally to form 
##  synthetic records from 2010 to 2060
## ---------------------------------------------------------------------------

## ------------------ clear memory and set paths -----------------------------

rm(list=ls())
setwd("C:/Users/David/OneDrive - University of North Carolina at Chapel Hill/UNC/Research/WSC/Coding/Dynamic Streamflow generator")
source("transientgeneratorfunctions.R"); source("generatorvalidationandplottingfunctions.R"); library(lubridate)

## ------------------ run generator ------------------------------------------

flowset  = 6; startset = 6
num_flowsets_s1 = 4; num_flowsets_s2 = 16; num_flowsets_s3 = 4; num_flowsets_s4 = 16

GeneratorInputsFolder = "inflow-data_updated"; ObservedRecordFolder  = "/historic"

DECADES = seq(from = 2010, to = 2060, by = 10)
num_weeks = 52; num_months = 12; num_years = 59; num_obs_years = 80; num_syn_years = 10

EVAPSOURCES = c("jordanevap", "michielittle", "owasaevap")
SOURCES     = c("Jordan", "Little", "Flat", "OWASA", EVAPSOURCES)
OBSSOURCES  = c("JordanLake", "LittleRiver", "Michie", "OWASA", EVAPSOURCES)

for (scenario in 3:4) {

  print(paste("Building records for scenario ", scenario, sep = ""))
  
  if (scenario %in% c(1,3)) {
    a = 1000
    num_realizations = a/num_flowsets_s1 # FOR SCENARIO 1 and 3
    maxflowset = 9
  } else {
    a = 1120
    num_realizations = a/num_flowsets_s2 # FOR SCENARIO 2 and 4
    maxflowset = 21
  } 
  
  # num_realizations = 5 #TEMPORARY
  # maxflowset = 7 #TEMPORARY
  
  ## ------------------ load input data ----------------------------------------
  
  print("Loading historic data...")
  HistoricArray = aggregateObservedFlows(num_obs_years, num_months, OBSSOURCES, 
                                         aggregation = "months", GeneratorInputsFolder, 
                                         ObservedRecordFolder)
  print("Done reading historic data.")
    # array of 80x12 matrices for each site
  
  for (flowset in startset:maxflowset) {
    print(paste("Flowset ", flowset, " will have ", num_realizations, " realizations.", sep = ""))
    
    print("Reading RHESSys output flows...")
    ModeledInfo   = aggregateModeledInputFlows(num_years, num_months, 
                                               aggregation = "months", SCENARIO = scenario, 
                                               SOURCES, DECADES, 
                                               GeneratorInputsFolder, chosenflows = flowset)
    print("Done reading modeled flows.")

      # list containing:
      #   1. Modeled Data (array of 60x12 matrices for each decade at each site)
    ModeledArray = ModeledInfo[[1]]
      #   2. Shifted Data (list containing:)
      #       a. List of shifted sets for each site (contains lists)
      #           a. List of shifted sets for each decade (contains lists)
      #               a. List of shifted timeseries for each calendar month
      #                   a. Timeseries of a single calendar month's shifted flows over modeled record years
    ShiftedLists = ModeledInfo[[2]]
      #   3. Shifted Data indices (list with same organization as Shifted Data)
    ShiftedIndex = ModeledInfo[[3]]
      #   4. Daily Modeled Data (array of [365*60]x6 matrices for each site (for 6 decades) with daily data)
    DailyArray   = ModeledInfo[[4]]
    
    rm(ModeledInfo)
    
    DaysPerMonth = c(31,28,31,30,31,30,31,31,30,31,30,31)
    SyntheticDailyArray = array(NA, dim = c(num_realizations, 365*num_syn_years, length(SOURCES), length(DECADES)))
    SyntheticWeekly = array(NA, dim = c(dim(SyntheticDailyArray)[1], dim(SyntheticDailyArray)[2]/365, 52, length(SOURCES), length(DECADES)))
    TransientWeeklyTimeseries = array(NA, dim = c(dim(SyntheticDailyArray)[1], dim(SyntheticDailyArray)[2]/365*52*length(DECADES), length(SOURCES)))
    SyntheticArrayHolder = array(NA, dim = c(num_realizations, 1, num_syn_years * num_months, length(SOURCES), length(DECADES)))
    
    for (realization in 1:num_realizations) {
      
      print(paste("Generating realization ", realization, " at ", 
                  length(SOURCES), " sites, over ", length(DECADES), " decades.", 
                  sep = ""))
      
      ## ------------------ run generator ------------------------------------------
      
      SyntheticArrayHolder[realization,,,,] = transienttrajectorygenerator(ModeledArray, 
                                                                           HistoricArray, 
                                                                           use.modeled.corr = TRUE)
      
      ## ------------------ disaggregate DAILY -------------------------------------
      #   select nearest neighbor month and disaggregate synthetic month flows    
      
      i = 1
      for (dec in 1:length(DECADES)) {
        calmonth = 1
        for (synmonth in 1:dim(SyntheticArrayHolder)[3]) {
          kNNinfo = kNNidentification(SyntheticTimeseries = SyntheticArrayHolder, ShiftedTimeseries = ShiftedLists, 
                                      month = calmonth, decade = dec, syntheticmonth = synmonth, 
                                      realization = realization)
          
          DProp   = kNNselection(identificationOutput = kNNinfo, 
                                 indicesList = ShiftedIndex, 
                                 HistoricalDailyArray = DailyArray, 
                                 month = calmonth, 
                                 decade = dec, 
                                 syntheticmonth = synmonth)
            # this doesnt work well for data distributed around 0 (like net evaporation...)
          
          invsites = c(5:7)
          for (q in invsites) {
            if (SyntheticArrayHolder[realization,,synmonth,q,dec] < 0) {
              DProp[,q] = DProp[,q]*-1
                # only invert the evap proportions if the monthly value is negative
                # this will ensure that rainfall days stay as such and aren't flipped
            }
          }
          
          SynYear     = floor((synmonth-1)/12)
          SynCalMonth = (synmonth-1) %% 12 + 1
          SynDayRange = (SynYear*365 + sum(DaysPerMonth[1:SynCalMonth]) - 
                           DaysPerMonth[SynCalMonth] + 1):(SynYear*365 + sum(DaysPerMonth[1:SynCalMonth]))
          
          SyntheticDailyArray[realization,SynDayRange,,dec] = t(t(DProp) * SyntheticArrayHolder[realization,,synmonth,,dec])
          
          calmonth = calmonth + 1
          if (calmonth > 12) {calmonth = 1}
          i = i + 1
        }
      }
      rm(i, calmonth, SynYear, SynCalMonth, SynDayRange, invsites, q, dec, synmonth)
      
      ## ------------------ aggregate WEEKLY ---------------------------------------
      
      for (site in 1:length(SOURCES)) {
        for (decade in 1:length(DECADES)) {
          OneRealizationSyn = matrix(SyntheticDailyArray[realization,,site,decade],nrow = dim(SyntheticDailyArray)[2]/365, 
                                     byrow = T)
          for (wk in 1:52) {
            if (wk < 52) {
              SyntheticWeekly[realization,,wk,site,decade] = apply(OneRealizationSyn[,((wk-1)*7+1):(wk*7)], MARGIN = 1, sum)
            } else {
              SyntheticWeekly[realization,,wk,site,decade] = apply(OneRealizationSyn[,((wk-1)*7+1):(wk*7)], MARGIN = 1, sum)
                # only take 7 days of data
              # SyntheticWeekly[realization,,wk,site,decade] = apply(OneRealizationSyn[,((wk-1)*7+1):ncol(OneRealizationSyn)], 
              #                                                      MARGIN = 1, sum)
                # take the rest of the year (8-9 days)
            }
          }
        }
      }
      rm(site, decade, wk)
      
      for (site in 1:length(SOURCES)) {
        temp = c()
        for (decade in 1:length(DECADES)) {
          temp = c(temp, c(t(SyntheticWeekly[realization,,,site,decade])))
        }
        TransientWeeklyTimeseries[realization,,site] = temp
      }
      rm(site, temp, decade)
    }
    rm(realization)
    
    ## ------------------- export synthetic records for each site ----------------
    
    for (sc in 1:length(SOURCES)) {
      print(paste("Printing ", SOURCES[sc], " synthetic flows...", sep = ""))
      write.table(TransientWeeklyTimeseries[,,sc], sep = ",", 
                paste("transientgeneratoroutput/", "scenario", scenario, "/", SOURCES[sc], "set", flowset, ".csv", sep = ""), 
                row.names = FALSE, col.names = FALSE)
      for (d in 1:6) {
        write.table(SyntheticDailyArray[,,sc,d], sep = ",", 
                    paste("transientgeneratoroutput/", "scenario", scenario, "/daily", SOURCES[sc], "_decade", d, "_set", flowset, ".csv", sep = ""), 
                    row.names = FALSE, col.names = FALSE)
        temp = matrix(NA, nrow = num_realizations, ncol = 520)
        for (r in 1:num_realizations) {temp[r,] = c(t(SyntheticWeekly[r,,,sc,d]))}; rm(r)
        write.table(temp, sep = ",", 
                    paste("transientgeneratoroutput/", "scenario", scenario, "/weekly", SOURCES[sc], "_decade", d, "_set", flowset, ".csv", sep = ""), 
                    row.names = FALSE, col.names = FALSE)
        write.table(SyntheticArrayHolder[,,,sc,d], sep = ",", 
                    paste("transientgeneratoroutput/", "scenario", scenario, "/monthly", SOURCES[sc], "_decade", d, "_set", flowset, ".csv", sep = ""), 
                    row.names = FALSE, col.names = FALSE)
      }
      rm(d, temp)
    }
    rm(sc)

    ## ---------- autocorrelation ----------------
     
    temporalACF(Syn = SyntheticArrayHolder, Mod = ModeledArray, nlags = 12, setting = "monthly", scen = scenario, flwst = flowset)
    #temporalACF(Syn = SyntheticWeekly, Mod = WeeklyArray, nlags = 28, setting = "weekly", scen = scenario, flwst = flowset)
    temporalACF(Syn = SyntheticDailyArray, Mod = DailyArray, nlags = 31, setting = "daily", scen = scenario, flwst = flowset)
    
    ## ---------- spatial correlation ------------
    
    a = spatialCor(Syn = SyntheticArrayHolder, Mod = ModeledArray, setting = "monthly", scen = scenario, flwst = flowset)
    #b = spatialCor(Syn = SyntheticWeekly, Mod = WeeklyArray, setting = "weekly", scen = scenario, flwst = flowset)
    c = spatialCor(Syn = SyntheticDailyArray, Mod = DailyArray, setting = "daily", scen = scenario, flwst = flowset)
    
    rm(a,b,c)
  }
  
  ## ---------- flow duration curves -----------
  # NEED TO LOG TRANSFORM RHESSYS FLOWS BEFORE PLOTTING
  WeeklyInfo = aggregateModeledInputFlows(num_years, 52, aggregation = "weeks", 
                                          SCENARIO = scenario, SOURCES, DECADES, GeneratorInputsFolder)
  WeeklyArray = WeeklyInfo[[1]]
  weeklyexceedenceCurves(Syn = SyntheticWeekly, Mod = WeeklyArray, num_rlz = num_realizations, scen = scenario)
  #dailyexceedenceCurves(SyntheticDailyArray, DailyArray, num_rlz = num_realizations, scen = scenario)
  monthlyexceedenceCurves(Syn = SyntheticArrayHolder, Mod = ModeledArray, num_rlz = num_realizations, scen = scenario)
  
  rm(WeeklyInfo, WeeklyArray)
}

## -------------------- build full synthetic sets -----------------------------

for (s in 2:4) {
  if (s %in% c(1,3)) {num_keep = 250} else {num_keep = 70}
  
  # num_keep = 1 #TEMPORARY
  
  flowsources = c("Little", "Flat", "Jordan", "OWASA")
  evapsources = c("jordanevap", "michielittle", "owasaevap")
  
  if (s == 1) {nsets = num_flowsets_s1} else if (s == 2) {nsets = num_flowsets_s2} else if (s == 3) {nsets = num_flowsets_s3} else {nsets = num_flowsets_s4}
  
  #nsets = 1 #TEMPORARY
  
  for (l in flowsources) {
    makemastersets(scen = s, sets = nsets, startset = 6, SRCE = l, 
                   realstotake = num_keep, keep.hist.yrs = 45, keep.syn.yrs = 51)
  }
  
  for (l in evapsources) {
    makemasterevapsets(scen = s, sets = nsets, startset = 6, SRCE = l, 
                       realstotake = num_keep, keep.hist.yrs = 45, keep.syn.yrs = 51)
  }
}

