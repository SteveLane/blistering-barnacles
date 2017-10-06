#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
################################################################################
################################################################################
## Title: Data cleaning and imputation
## Author: Steve Lane
## Date: Wednesday, 08 March 2017
## Synopsis: Cleans data for manuscript and model fitting, and performs
## imputation on the vessel level.
## Time-stamp: <2017-10-06 14:44:24 (slane)>
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
packages <- c("dplyr", "mice", "parallel")
ipak(packages)
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
        LocID = as.factor(LocID),
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
## Data for imputations/modelling
impData <- data %>% select(-samLoc, -cens, -LocID)
## Level 1 data
lvl1Data <- data %>% select(boatID, wetWeight, LocID, cens)
## Loop to create multiple imputations - give a loop as I want to start each
## imputation off with a random draw from a U(0, 1.5) for the censored data just
## to inject a little randomness into it.
## Create numMI imputations, join to full data, and also create stan data
set.seed(787, "L'Ecuyer")
impList <- mclapply(1:numMI, function(i){
    imp <- lvl2Imp(impData)
    stanData <- createStanData(lvl1Data, imp)
    list(lvl2 = imp, stanData = stanData)
})
## Save as rds for further use.
if(!dir.exists("../data/")) dir.create("../data/")
saveRDS(data, file = "../data/biofouling.rds")
saveRDS(impList, file = "../data/imputations.rds")
################################################################################
################################################################################
