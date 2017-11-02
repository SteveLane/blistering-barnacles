#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
################################################################################
################################################################################
## Title: Fit model (VB)
## Author: Steve Lane
## Date: Friday, 21 April 2017
## Time-stamp: <2017-11-02 00:43:46 (overlordR)>
## Synopsis: Script that drives the censored regression model. Designed to be
## called from the Makefile, it requires the model name, a seed for rng, and
## number of iterations to be set on the command line, or prior to sourcing the
## script.
## This particular version runs a variational bayes sampler for speed.
################################################################################
################################################################################
if(!(length(args) %in% 2:3)){
    stop("Three arguments must be supplied: model name, myseed, and iter.\niter (for HMC) is optional, and defaults to 2000.\nRscript fit-model.R mname=model-name myseed=my-seed iter=num-iter",
         call. = FALSE)
} else {
    if(length(args) == 2){
        ## Default if option not specified
        iter <- 1000
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
packages <- c("rstan", "parallel")
ipak(packages)
## Load stan model (it should already be compiled from the compile-model.R
## script)
rstan_options(auto_write = TRUE)
## Want cores to be one, we're only running one chain, then combining. Each
## imputation will be sent out via mclapply.
options(mc.cores = parallel::detectCores()/2)
model <- stan_model(paste0("../stan/", mname, ".stan"))
## Load data
impList <- readRDS("../data/imputations.rds")
newData <- readRDS("../data/newData.rds")
newStan <- with(
    newData,
    list(newN = nrow(newData), days2New = days2, locIDNew = locIDInt,
         paintTypeNew = paintTypeInt, boatTypeNew = boatTypeInt,
         days1New = days1, midTripsNew = midTrips)
)
set.seed(myseed)
out <- mclapply(impList, function(dat){
    locMod <- try(vb(model, data = c(dat$stanData, newStan), elbo_samples = 250,
                     adapt_iter = 100))
    locMod
})
keep <- sapply(out, function(x) {
    if(inherits(x, "try-error")) {
        return(FALSE)
    } else {
        return(TRUE)
    }
})
outname <- paste0("../data/", mname, "-var-bayes.rds")
output <- sflist2stanfit(out[keep])
saveRDS(output, file = outname)
################################################################################
################################################################################
