################################################################################
################################################################################
## Title: Post processing, robust model
## Author: Steve Lane
## Date: Thursday, 04 May 2017
## Synopsis: Post process the output from the regression models
## Time-stamp: <2017-11-08 22:07:49 (overlordR)>
################################################################################
################################################################################
## Add github packages using gitname/reponame format
source("../scripts/imputation-functions.R")
packages <- c("tidyr", "dplyr", "tibble", "rstan", "loo", "ggplot2",
              "RColorBrewer")
ipak(packages)
options(loo.cores = 1)
## This is set for readable text when included at half page width.
theme_set(theme_bw())
m0 <- readRDS("../data/censored-mle-m0-t.rds")
m1 <- readRDS("../data/censored-mle-m1-t.rds")
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
    mutate(nummi = rep(seq_along(imps), each = nrow(vessels))) %>%
    left_join(., bLookup)
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
paintDat <- vessImps %>% filter(nummi %in% c(0, impSelect),
                                !is.na(paintType)) %>%
    mutate(numMI = factor(nummi))
plPaint <- ggplot(paintDat, aes(x = paintType)) +
    geom_bar() +
    facet_wrap(~ numMI) +
    xlab("Anti-fouling paint type") +
    ylab("Count") +
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
    left_join(., bLookup, by = c("boatTypeInt", "paintTypeInt")) %>%
    select(-paintTypeInt, -boatTypeInt)
a3 <- extract(m3, "alphaHat")$alphaHat
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
coefs <- extract(m3, c("betaDays1", "betaDays2", "betaMidTrips", "betaHullSA",
                       "betaPaint", "betaType", "betaDaysType",
                       "betaTripsType", "betaTripsPaint"))
intType <- apply(coefs$betaType, 2, median)
coefDaysType <- c(median(coefs$betaDays1 + coefs$betaDaysType[,1]),
                  median(coefs$betaDays1 + coefs$betaDaysType[,2]),
                  median(coefs$betaDays1 + coefs$betaDaysType[,3]))
coefTripsType <- c(median(coefs$betaMidTrips + coefs$betaTripsType[,1]),
                   median(coefs$betaMidTrips + coefs$betaTripsType[,2]),
                   median(coefs$betaMidTrips + coefs$betaTripsType[,3]))
parsType <- tibble(intercept = rep(intType, 2),
                   slope = c(coefDaysType, coefTripsType),
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
ggsave("../graphics/plM3boat-t.pdf", plM3boat, height = 3.5)
intPaint <- apply(coefs$betaPaint, 2, median)
coefTripsPaint <- c(median(coefs$betaMidTrips + coefs$betaTripsPaint[,1]),
                    median(coefs$betaMidTrips + coefs$betaTripsPaint[,2]),
                    median(coefs$betaMidTrips + coefs$betaTripsPaint[,3]))
parsType <- tibble(intercept = intPaint,
                   slope = coefTripsPaint,
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
ggsave("../graphics/plM3paint-t.pdf", plM3paint, width = 3.5, height = 3.5)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Comparison of regression coefficients
################################################################################
################################################################################
qnts <- c("low5" = 0.05, "low25" = 0.25, "med" = 0.5, "high75" = 0.75,
          "high95" = 0.95)
coef0 <- extract(m0, pars = c("nu", "mu", "betaLoc", "sigma_alphaBoat",
                              "sigma"))
m0Summary <- lapply(coef0, sumMC, qnts = qnts) %>% summRename %>%
    mutate(model = "M0")
coef1 <- extract(m1, pars = c("nu", "mu", "betaLoc", "betaDays1", "betaDays2",
                              "betaMidTrips", "betaHullSA", "betaPaint",
                              "betaType", "sigma_alphaBoat", "sigma"))
m1Summary <- lapply(coef1, sumMC, qnts = qnts) %>% summRename %>%
    mutate(model = "M1")
coef3 <- extract(m3, pars = c("nu", "mu", "betaLoc", "betaDays1", "betaDays2",
                              "betaMidTrips", "betaHullSA", "betaPaint",
                              "betaType", "betaDaysType", "betaTripsType",
                              "betaTripsPaint", "sigma_alphaBoat", "sigma"))
m3Summary <- lapply(coef3, sumMC, qnts = qnts) %>% summRename %>%
    mutate(model = "M2")
allSummary <- bind_rows(m0Summary, m1Summary, m3Summary)
ords <- unique(m3Summary$coef)
allSummary <- allSummary %>%
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
plSummary <- ggplot(allSummary, aes(x = coef, y = med, ymin = low5,
                                    ymax = high95, colour = model)) +
    geom_pointrange(aes(ymin = low25, ymax = high75),
                    position = position_dodge(width = 1.0), fatten = 0.75,
                    size = 1, show.legend = FALSE) +
    geom_pointrange(position = position_dodge(width = 1.0), fatten = 0.5) +
    coord_flip() +
    xlab("Variable") +
    ylab("Value") +
    scale_x_discrete(labels = labs) +
    scale_colour_brewer(palette = "Dark2", name = "Model") +
    theme(legend.position = "bottom")
ggsave("../graphics/plSummary-t.pdf", plSummary, width = 3.5)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Posterior predictive comparisons
################################################################################
################################################################################
## New data exists.
newData <- readRDS("../data/newData.rds") %>%
    left_join(., bLookup)
yNew1 <- extract(m3, "yNew1")$yNew1
## Obs 1, 2, and 3 contain the varying days2 data.
diffDays2 <- c(
    mean(exp(yNew1[, 2]) - exp(yNew1[, 1]) > 0),
    mean(exp(yNew1[, 3]) - exp(yNew1[, 1]) > 0),
    quantile(exp(yNew1[, 2]) - exp(yNew1[, 1]), probs = 0.5),
    quantile(exp(yNew1[, 3]) - exp(yNew1[, 1]), probs = 0.5)
)
## Obs 55, 28, and 1 contain yachts, fishing and motor cruisers with all other
## variables the same.
diffType <- c(
    mean(exp(yNew1[, 55]) - exp(yNew1[, 28]) > 0),
    mean(exp(yNew1[, 55]) - exp(yNew1[, 1]) > 0),
    quantile(exp(yNew1[, 55]) - exp(yNew1[, 28]), probs = 0.5),
    quantile(exp(yNew1[, 55]) - exp(yNew1[, 1]), probs = 0.5)
)
## Differences in location, everything else held same. Hull (row 1) and Rudder
## (row 7) vs keel (row 4).
diffLoc <- c(
    mean(exp(yNew1[, 1]) - exp(yNew1[, 4]) > 0),
    mean(exp(yNew1[, 7]) - exp(yNew1[, 4]) > 0),
    quantile(exp(yNew1[, 1]) - exp(yNew1[, 4]), probs = 0.5),
    quantile(exp(yNew1[, 7]) - exp(yNew1[, 4]), probs = 0.5)
)
saveRDS(list(diffDays2 = diffDays2, diffType = diffType, diffLoc = diffLoc),
        "../data/diffs.rds")
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Create data for looic table.
################################################################################
################################################################################
m0ll <- extract_log_lik(m0)
m0loo <- loo(m0ll)
m1ll <- extract_log_lik(m1)
m1loo <- loo(m1ll)
m3ll <- extract_log_lik(m3)
m3loo <- loo(m3ll)
looTab <- as.data.frame(compare(m0loo, m1loo, m3loo)) %>%
    mutate(mID = rownames(.))
names(looTab) <- c("LOOIC", "se(LOOIC)", "ELPD", "se(ELPD)", "Eff. P",
                   "se(Eff. P)", "mID")
looLookup <- tibble(
    mID = paste0("m", 0:3, "loo"),
    Model = paste0("M", 0:3)
)
looTab <- left_join(looTab, looLookup) %>%
    select(-mID)
## Model 7 has the lowest looic/elpd, but not more so than model 5:
diffs <- rbind(
    rep(NA, 2),
    compare(m0loo, m3loo),
    compare(m1loo, m3loo)
) %>%
    as_tibble() %>%
    mutate(Model = c("M2", "M0", "M1"))
## Put differences on LOOIC scale (LOOIC = -2*ELPD)
looTab <- left_join(looTab, diffs) %>%
    mutate(elpd_diff = 2 * elpd_diff,
           se = sqrt(2) * se) %>%
    rename(`$\\Delta$LOOIC` = elpd_diff,
           `se($\\Delta$LOOIC)` = se) %>%
    select(Model, everything())
saveRDS(looTab, "../data/looic-t.rds")
################################################################################
################################################################################
