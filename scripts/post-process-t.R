################################################################################
################################################################################
## Title: Post processing, robust model
## Author: Steve Lane
## Date: Thursday, 04 May 2017
## Synopsis: Post process the output from the regression models
## Time-stamp: <2017-10-13 00:16:47 (overlordR)>
################################################################################
################################################################################
## Add github packages using gitname/reponame format
source("../scripts/imputation-functions.R")
packages <- c("tidyr", "dplyr", "tibble", "rstan", "loo", "ggplot2",
              "RColorBrewer")
ipak(packages)
## This is set for readable text when included at half page width.
theme_set(theme_bw())
m3 <- readRDS("../data/censored-mle-m3-t.rds")
imps <- readRDS("../data/imputations.rds")
biofoul <- readRDS("../data/biofouling.rds")
vessels <- biofoul %>% distinct(boatID, .keep_all = TRUE)
bLookup <- biofoul %>%
    select(boatType, paintType, boatTypeInt, paintTypeInt) %>%
    distinct() %>%
    na.omit()
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Observed weight figures
################################################################################
################################################################################
histData <- biofoul %>%
    filter(wetWeight >= 1.5) %>%
    mutate(wwLog = log(wetWeight)) %>%
    select(LocID, `Weight (gm)` = wetWeight,
           `Weight (gm), log-scale` = wwLog) %>%
    gather("logged", "ww", 2:3)
plHist <- ggplot(histData, aes(x = ww)) + geom_histogram(bins = 11) +
    facet_grid(LocID ~ logged, scales = "free") +
    xlab("") +
    ylab("Count") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
ggsave("../graphics/obs-hist.pdf", plHist)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Imputation figures
################################################################################
################################################################################
lvl2Imp <- lapply(imps, function(x) x$lvl2) %>% bind_rows() %>%
    mutate(nummi = rep(seq_along(imps), each = nrow(vessels)))
vessImps <- vessels %>%
    mutate(
        days1 = as.numeric(scale(days1)),
        days2 = as.numeric(scale(days2)),
        midTrips = as.numeric(scale(midTrips)),
        ApproxHullSA = as.numeric(scale(ApproxHullSA)),
        nummi = 0
    )
vessImps <- bind_rows(vessImps, lvl2Imp)
plDays1 <- ggplot(vessImps, aes(x = factor(nummi), y = days1)) +
    geom_boxplot(outlier.size = 0.5) +
    xlab("Imputation number") +
    ylab("Standardised value") +
    theme_bw(base_size = 7.7) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
ggsave("../graphics/imp-days1.pdf", plDays1, width = 4.9, height = 4.9)
plMidTrips <- ggplot(vessImps, aes(x = factor(nummi), y = midTrips)) +
    geom_boxplot(outlier.size = 0.5) +
    xlab("Imputation number") +
    ylab("Standardised value") +
    theme_bw(base_size = 7.7) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
ggsave("../graphics/imp-trips.pdf", plMidTrips, width = 4.9, height = 4.9)
## Adjust based on length of imps
set.seed(76)
if(length(imps) < 15){
    impSelect <- seq_along(imps)
} else {
    impSelect <- sample(seq_along(imps), 15)
}
plPaint <- ggplot(vessImps %>% filter(nummi %in% c(0, impSelect),
                                      !is.na(paintType)),
                  aes(x = paintType)) +
    geom_bar() +
    facet_wrap(~ factor(nummi)) +
    xlab("Anti-fouling paint type") +
    theme_bw(base_size = 7.7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("../graphics/imp-paint.pdf", plPaint, width = 4.9, height = 4.9)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: M3 figures
################################################################################
################################################################################
set.seed(13)
lvl2 <- imps[[sample(seq_along(imps), 1)]]$lvl2 %>%
    left_join(., bLookup) %>%
    select(-boatTypeInt, -paintTypeInt)
a3 <- extract(m3, "alphaBoat")$alphaBoat
a3Sum <- t(apply(a3, 2, quantile, probs = c(0.1, 0.5, 0.9)))
a3Dat <- tibble(low = a3Sum[,1], mid = a3Sum[,2], high = a3Sum[,3],
                boatID = 1:nrow(a3Sum))
a3Dat <- left_join(a3Dat, lvl2) %>%
    select(-ApproxHullSA, -days2) %>%
    gather(type, value, -low, -mid, -high, -boatID, -paintType, -boatType) %>%
    mutate(type = recode(type,
                         days1 = "Days since last used",
                         midTrips = "Median number of trips")) %>%
    rename(`Vessel type` = boatType, `Paint type` = paintType)
slopes <- extract(m3, c("betaDays1", "betaDays2", "betaMidTrips",
                        "betaPaint", "betaType", "betaDaysType",
                        "betaTripsType", "betaTripsPaint"))
intType <- apply(slopes$betaType, 2, median)
slopeDaysType <- c(median(slopes$betaDays1 + slopes$betaDaysType[,1]),
                   median(slopes$betaDays1 + slopes$betaDaysType[,2]),
                   median(slopes$betaDays1 + slopes$betaDaysType[,3]))
slopeTripsType <- c(median(slopes$betaMidTrips + slopes$betaTripsType[,1]),
                    median(slopes$betaMidTrips + slopes$betaTripsType[,2]),
                    median(slopes$betaMidTrips + slopes$betaTripsType[,3]))
parsType <- tibble(intercept = rep(intType, 2),
                   slope = c(slopeDaysType, slopeTripsType),
                   boatTypeInt = rep(1:3, 2),
                   type = rep(c("Days since last used",
                                "Median number of trips"), each = 3)) %>%
    left_join(., bLookup %>% select(contains("boat")) %>% distinct(),
              by = "boatTypeInt") %>%
    rename(`Vessel type` = boatType)
plM3boat <- ggplot(a3Dat, aes(x = value, y = mid, ymin = low, ymax = high,
                              colour = `Vessel type`)) +
    geom_pointrange(fatten = 0.5) +
    facet_wrap(~ type, ncol = 2) +
    geom_abline(aes(slope = slope, intercept = intercept,
                    colour = `Vessel type`),
                data = parsType) +
    ylab("Vessel-level intercept") +
    xlab("Scaled value") +
    scale_colour_brewer(palette = "Dark2") +
    theme(legend.position = "bottom")
ggsave("../graphics/plM3boat-robust.pdf", plM3boat, height = 3.5)
intPaint <- apply(slopes$betaPaint, 2, median)
slopeTripsPaint <- c(median(slopes$betaMidTrips + slopes$betaTripsPaint[,1]),
                     median(slopes$betaMidTrips + slopes$betaTripsPaint[,2]),
                     median(slopes$betaMidTrips + slopes$betaTripsPaint[,3]))
parsType <- tibble(intercept = intPaint,
                   slope = slopeTripsPaint,
                   paintTypeInt = 1:3,
                   type = rep("Median number of trips", each = 3)) %>%
    left_join(., bLookup %>% select(contains("paint")) %>% distinct(),
              by = "paintTypeInt") %>%
    rename(`Paint type` = paintType)
plM3paint <- ggplot(a3Dat %>% filter(type == "Median number of trips"),
                    aes(x = value, y = mid, ymin = low, ymax = high,
                        colour = `Paint type`)) +
    geom_pointrange(fatten = 1) +
    facet_wrap(~ type) +
    geom_abline(aes(slope = slope, intercept = intercept,
                    colour = `Paint type`),
                data = parsType) +
    ylab("Vessel-level intercept") +
    xlab("Scaled value") +
    scale_colour_brewer(palette = "Dark2") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom")
ggsave("../graphics/plM3paint-robust.pdf", plM3paint, width = 3.5, height = 3.5)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Regression coefficients
################################################################################
################################################################################
coef3 <- extract(m3, pars = c("nu", "mu", "betaLoc", "betaDays1", "betaDays2",
                              "betaMidTrips", "betaPaint",
                              "betaType", "betaDaysType", "betaTripsType",
                              "betaTripsPaint", "sigma_alphaBoat", "sigma"))
m3Summary <- lapply(coef3, sumMC, qnts = qnts) %>% summRename %>%
    mutate(model = "M3")
ords <- unique(m3Summary$coef)
inds <- which(ords == "betaMidTrips")
ords <- c(ords[1:inds], "betaHullSA", ords[(inds+1):length(ords)])
m3Summary <- m3Summary %>%
    mutate(coef = factor(coef, levels = rev(ords)))
labs <- c("nu" = expression(nu),
          "mu" = expression(mu),
          "betaLoc1" = expression(beta[1]^{L}),
          "betaLoc2" = expression(beta[2]^{L}),
          "betaLoc3" = expression(beta[3]^{L}),
          "betaDays1" = expression(beta^{d1}),
          "betaDays2" = expression(beta^{d2}),
          "betaMidTrips" = expression(beta^{m}),
          "betaHullSA" = expression(beta^{h}),
          "betaPaint1" = expression(beta[1]^{P}),
          "betaPaint2" = expression(beta[2]^{P}),
          "betaPaint3" = expression(beta[3]^{P}),
          "betaType1" = expression(beta[1]^{T}),
          "betaType2" = expression(beta[2]^{T}),
          "betaType3" = expression(beta[3]^{T}),
          "betaDaysType1" = expression(beta[1]^{DT}),
          "betaDaysType2" = expression(beta[2]^{DT}),
          "betaDaysType3" = expression(beta[3]^{DT}),
          "betaTripsType1" = expression(beta[1]^{MT}),
          "betaTripsType2" = expression(beta[2]^{MT}),
          "betaTripsType3" = expression(beta[3]^{MT}),
          "betaTripsPaint1" = expression(beta[1]^{MP}),
          "betaTripsPaint2" = expression(beta[2]^{MP}),
          "betaTripsPaint3" = expression(beta[3]^{MP}),
          "sigma_alphaBoat" = expression(sigma[alpha]),
          "sigma" = expression(sigma))
plSummary <- ggplot(m3Summary, aes(x = coef, y = med, ymin = low5,
                                   ymax = high95)) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_pointrange(aes(ymin = low25, ymax = high75),
                    position = position_dodge(width = 1.0), fatten = 0.75,
                    size = 1, show.legend = FALSE) +
    geom_pointrange(position = position_dodge(width = 1.0), fatten = 0.5) +
    coord_flip() +
    xlab("Variable") +
    ylab("Value") +
    scale_x_discrete(labels = labs)
ggsave("../graphics/plSummary-robust.pdf", plSummary, width = 3.5)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Posterior predictive comparisons
################################################################################
################################################################################
## Yacht, with mean days1, mean days2 -> 3 or 6 months longer?, measured at the
## hull, and ablative paint.
## Data is scaled, so figure out the scaling.
## I also need to calculate the actual predictions and difference them, given
## we're doing posterior predictions.
## This is probably easy to do via matrix calcs.
## First the coefficients for alpha hat
## aHatCoefs <- with(coef3, cbind(betaDays1, betaDays2, betaMidTrips, betaPaint,
##                                betaType, betaDaysType, betaTripsType,
##                                betaTripsPaint))

## newData <- data.frame(
##     days1 = 0,
##     days2 = c(0, (365/4) / sd(vessels$days2), (365/2) / sd(vessels$days2)),
##     midTrips = 0,
##     paintType = factor("Ablative", levels = levels(vessels$paintType)),
##     boatType = factor("Yacht", levels(vessels$boatType))
## )
## newMat <- model.matrix(
##     ~ days1 + days2 + midTrips + paintType + boatType +
##         days1:boatType + midTrips:boatType + midTrips:paintType,
##     data = newData)
## aHat <- aHatCoefs %*% t(newMat[, -1])

## set.seed(37)
## days2 <- c((365/4) / sd(vessels$days2), (365/2) / sd(vessels$days2))
## aHat <- t(days2 %*% t(coef3$betaDays2)) +
##     cbind(coef3$betaType[, 2], coef3$betaType[, 2])
## alphaHat <- t(sapply(seq_len(nrow(aHat)), function(i){
##     a0 <- rcauchy(1, aHat[i, 1], coef3$sigma_alphaBoat[i])
##     a1 <- rcauchy(1, aHat[i, 2], coef3$sigma_alphaBoat[i])
##     a2 <- rcauchy(1, aHat[i, 3], coef3$sigma_alphaBoat[i])
##     c(a0, a1, a2)
## }))
## muHat <- alphaHat + cbind(coef3$mu, coef3$mu, coef3$mu)
## tStd <- t(sapply(seq_len(nrow(aHat)), function(i){
##     coef3$sigma[i] * rt(3, coef3$nu[i])
## }))
## lY <- muHat + tStd
## diffDays2 <- apply(exp(lY), 2, quantile, probs = 0.2, names = FALSE)
## ## Difference between yachts and motor cruisers/fishing. Base case is motor
## ## cruisers. At mean values (so they're zero). Ablative paint (although it
## ## shouldn't matter).
## alphaHat <- t(sapply(seq_len(nrow(aHat)), function(i){
##     aHat1 <- 
##     a1 <- rcauchy(1, coef3$betaType[i, 2], coef3$sigma_alphaBoat[i])
##     a2 <- rcauchy(1, coef3$betaType[i, 2] - coef3$betaType[i, 1],
##                   coef3$sigma_alphaBoat[i])
##     c(a1, a2)
## }))
## muHat <- alphaHat + cbind(coef3$mu, coef3$mu)
## tStd <- t(sapply(seq_len(nrow(aHat)), function(i){
##     coef3$sigma[i] * rt(2, coef3$nu[i])
## }))
## lY <- muHat + tStd
## diffType <- apply(exp(lY), 2, quantile, probs = 0.2, names = FALSE)
## saveRDS(list(diffDays2 = diffDays2, diffType = diffType),
##         "../data/diffs.rds")
################################################################################
################################################################################
