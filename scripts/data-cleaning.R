#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
################################################################################
################################################################################
## Title: Data cleaning and imputation
## Author: Steve Lane
## Date: Wednesday, 08 March 2017
## Synopsis: Cleans data for manuscript and model fitting, and performs
## imputation on the vessel level.
## Time-stamp: <2017-11-02 21:35:38 (overlordR)>
################################################################################
################################################################################
if(!(length(args) %in% 0:1)){
    stop("One argument can be supplied: numMI.\nnumMI, the number of multiply imputed datasets is optional, and defaults to 50.\nRscript data-cleaning.R numMI=numMI",
         call. = FALSE)
} else {
    if(length(args) == 0){
        ## Default if option not specified
        numMI <- 50
    }
    hasOpt <- grepl("=", args)
    argLocal <- strsplit(args[hasOpt], "=")
    for(i in seq_along(argLocal)){
        value <- NA
        tryCatch(value <- as.double(argLocal[[i]][2]), warning = function(e){})
        if(!is.na(value)){
            ## Assume int/double
            assign(argLocal[[i]][1], value, inherits = TRUE)
        } else {
            assign(argLocal[[i]][1], argLocal[[i]][2], inherits = TRUE)
        }
    }
}
## Add github packages using gitname/reponame format
source("../scripts/imputation-functions.R")
packages <- c("dplyr", "mice")
ipak(packages)
## Turn off messages
options(warn = -1, verbose = FALSE)
samplesdata <- read.csv("../data-raw/samples.csv")
## Bring in vessel data as well.
vessels <- read.csv("../data-raw/vessel.csv")
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Filter out locations not used and mutate data as necessary.
################################################################################
################################################################################
## First prepare data.
data <- samplesdata %>% filter(LocID %in% c("HA", "PJ", "HP")) %>%
    select(-boatType, -boatCode, -paintRating, -Location, -LocCode, -samLab,
           -wetWeight1, -wetWeight2)
vessels <- vessels %>%
    distinct(boatID, .keep_all = TRUE) %>%
    select(boatID, samLoc, boatType, ApproxHullSA) %>%
    filter(!is.na(boatID))
## Merge on vessel data
data <- left_join(
    data,
    vessels,
    by = "boatID"
) %>%
    mutate(
        samLoc = recode(samLoc,
                        "SaNAringham Yacht Club" = "Sandringham Yacht Club",
                        "Hobsons Bay Yacht Club " = "Hobsons Bay Yacht Club"),
        boatType = recode(boatType,
                          "Fishing vessel (Lobster/scallop)" = "Fishing",
                          "Fishing vessel (Abalone mothership)" = "Fishing",
                          "Fishing vessel (Long line)" = "Fishing",
                          "Motor cruiser" = "Motor cruiser/Other",
                          "Ferry" = "Motor cruiser/Other",
                          "Tug" = "Motor cruiser/Other"),
        boatType = as.factor(boatType),
        LocID = recode(LocID,
                       "HA" = "Hull", "HP" = "Keel", "PJ" = "Rudder"),
        LocID = factor(LocID),
        cens = ifelse(wetWeight < 1.5, 1, 0),
        paintType = as.factor(paintType)
    )
## Boat ID 24 has ridiculously high wet weights, as does 36 in the rudder
## locations. Remove 24 completely, and the bad 36 observations. Relabel boatID
## for future use.
data <- data %>%
    filter(wetWeight <= 1000) %>%
    mutate(boatIDOld = boatID,
           boatID = ifelse(boatID > 24, boatID - 1, boatID))
## Loop to create multiple imputations - give a loop as I want to start each
## imputation off with a random draw from a U(0, 1.5) for the censored data just
## to inject a little randomness into it.
## Create numMI imputations, join to full data, and also create stan data.
## Furthermore, I need indicators for categorical values, so probably easiest to
## create a series of lookup tables.
locLookup <- data_frame(
    LocID = levels(data$LocID),
    LocIDInt = seq_len(length(LocID))
)
boatLookup <- data_frame(
    boatType = levels(data$boatType),
    boatTypeInt = seq_len(length(boatType))
)
paintLookup <- data_frame(
    paintType = levels(data$paintType),
    paintTypeInt = seq_len(length(paintType))
)
## Data for imputations/modelling
impData <- data %>% select(-samLoc, -cens, -LocID)
## Level 1 data
lvl1Data <- data %>%
    left_join(., locLookup, by = "LocID") %>%
    select(boatID, wetWeight, LocIDInt, cens)
set.seed(13)
impList <- lapply(1:numMI, function(i){
    imp <- lvl2Imp(impData) %>%
        left_join(., boatLookup, by = "boatType") %>%
        left_join(., paintLookup, by = "paintType") %>%
        select(-boatType, -paintType)
    stanData <- createStanData(lvl1Data, imp)
    list(lvl2 = imp, stanData = stanData)
})
data <- data %>%
    left_join(., locLookup, by = "LocID") %>%
    left_join(., boatLookup, by = "boatType") %>%
    left_join(., paintLookup, by = "paintType")
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Create newdata for predictions
################################################################################
################################################################################
## Days2 at 0, 3, and 6 months.
subData <- data %>% select(boatID, days2) %>% distinct(boatID, .keep_all = TRUE)
newData <- expand.grid(
    days2 = c(0, (365/4) / sd(subData$days2), (365/2) / sd(subData$days2)),
    locIDInt = 1:3, paintTypeInt = 1:3, boatTypeInt = 1:3) %>%
    mutate(days1 = 0, midTrips = 0)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Save as rds for further use.
################################################################################
################################################################################
if(!dir.exists("../data/")) dir.create("../data/")
saveRDS(data, file = "../data/biofouling.rds")
saveRDS(impList, file = "../data/imputations.rds")
saveRDS(impList[1:10], file = "../data/imputations-short.rds")
saveRDS(newData, file = "../data/newData.rds")
################################################################################
################################################################################
