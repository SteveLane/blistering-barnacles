################################################################################
################################################################################
## Title: Imputation functions
## Author: Steve Lane
## Date: Wednesday, 29 March 2017
## Synopsis: Functions to run the censored regression imputation models
## Time-stamp: <2017-10-11 01:36:59 (overlordR)>
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Function to load packages (and install if required).
################################################################################
################################################################################
ipak <- function(pkg){
    ## Check for github packages (throw away github username)
    chk.git <- gsub(".*/", "", pkg)    
    new.pkg <- pkg[!(chk.git %in% installed.packages()[, "Package"])]
    if(!(length(new.pkg) == 0)){
        git.ind <- grep("/", new.pkg)
        if(length(git.ind) == 0){
            install.packages(new.pkg, dependencies = TRUE,
                             repos = "https://cran.csiro.au/")
        } else {
            devtools::install_github(new.pkg[git.ind])
        }
    }
    sapply(chk.git, require, character.only = TRUE)
}
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Perform level 2 imputation.
################################################################################
################################################################################
## Function takes in the whole dataset, summarises within-boat wet weights by
## the median, performs one iteration of random forest imputation using mice,
## then passes back the imputed data.
lvl2Imp <- function(data){
    data <- data %>%
        mutate(wwImp = ifelse(wetWeight < 1.5, runif(n(), 0, 1.5), wetWeight))
    lvl2 <- data %>%
        group_by(boatID, days1, days2, midTrips, paintType, boatType,
                 ApproxHullSA) %>%
        summarise(m = median(wwImp)) %>% ungroup() %>%
        mutate(
            days1 = as.numeric(scale(days1)),
            days2 = as.numeric(scale(days2)),
            midTrips = as.numeric(scale(midTrips)),
            ApproxHullSA = as.numeric(scale(ApproxHullSA))
        )
    pred.mat <- matrix(1 - diag(ncol(lvl2)), ncol(lvl2))
    pred.mat[1, ] <- pred.mat[, 1] <- 0
    miLvl2 <- mice(lvl2, m = 1, method = "rf", predictorMatrix = pred.mat,
                   maxit = 20, printFlag = FALSE)
    lvl2 <- complete(miLvl2) %>% select(-m)
    return(lvl2)
}
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Create stan data function
################################################################################
################################################################################
createStanData <- function(lvl1, lvl2){
    ## Can pass list with too many variables to stan, so just add all as it
    ## makes it easier in the long run.
    obsData <- lvl1 %>% filter(cens == 0)
    censData <- lvl1 %>% filter(cens == 1)
    stanObs <- with(
        obsData,
        list(Y = wetWeight, N = nrow(obsData), numLoc = max(LocIDInt),
             locID = LocIDInt, numBoat = max(boatID), boatID = boatID)
    )
    stanMiss <- with(
        censData,
        list(nCens = nrow(censData), U = 1.5, locIDCens = LocIDInt,
             boatIDCens = boatID)
    )
    stanLvl2 <- with(
        lvl2,
        list(days1 = days1, days2 = days2, midTrips = midTrips,
             hullSA = ApproxHullSA, numPaint = max(paintTypeInt),
             paintType = paintTypeInt, numType = max(boatTypeInt),
             boatType = boatTypeInt)
    )
    stanData <- c(stanObs, stanMiss, stanLvl2)
    return(stanData)
}
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Function to create coefficient summary data
################################################################################
################################################################################
sumMC <- function(draws, qnts){
    if(is.na(ncol(draws))){
        summ <- quantile(draws, qnts, names = FALSE)
        summ <- tibble::as_tibble(t(summ))
        names(summ) <- names(qnts)
    } else {
        summ <- apply(draws, 2, quantile, probs = qnts, names = FALSE)
        summ <- tibble::as_tibble(t(summ))
        names(summ) <- names(qnts)
    }
    summ
}
summRename <- function(summList){
    summ <- lapply(names(summList), function(nm){
        dat <- summList[[nm]]
        nr <- nrow(dat)
        if(nr > 1){
            dat <- dat %>%
                mutate(coef = paste0(nm, seq_len(nr)))
        } else {
            dat <- dat %>%
                mutate(coef = nm)
        }
    })
    bind_rows(summ)
}
################################################################################
################################################################################
