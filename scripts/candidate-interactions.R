################################################################################
################################################################################
## Title: Candidate Interactions
## Author: Steve Lane
## Date: Monday, 06 November 2017
## Synopsis: Plot marginal random effects to look at interactions.
## Time-stamp: <2017-11-27 04:08:08 (overlordR)>
################################################################################
################################################################################
rm(list = ls())
library(dplyr)
library(tidyr)
library(rstan)
library(ggplot2)
theme_set(theme_bw())
bLookup <- readRDS("../data/biofouling.rds") %>%
    select(boatType, paintType, boatTypeInt, paintTypeInt) %>%
    distinct() %>%
    na.omit()
m1 <- readRDS("../data/censored-mle-m1-t.rds")

################################################################################
################################################################################
## Begin Section: Vessel type
################################################################################
################################################################################
a1 <- extract(m1, "alphaHat")$alphaHat
a1Sum <- t(apply(a1, 2, quantile, probs = c(0.1, 0.5, 0.9)))
a1Dat <- tibble(low = a1Sum[,1], mid = a1Sum[,2], high = a1Sum[,3],
                boatID = 1:nrow(a1Sum))
imps <- readRDS("../data/imputations.rds")
set.seed(13)
lvl2 <- imps[[sample(seq_along(imps), 1)]]$lvl2 %>%
    left_join(., bLookup, by = c("boatTypeInt", "paintTypeInt")) %>%
    select(-paintTypeInt, -boatTypeInt)
a1Dat <- left_join(a1Dat, lvl2, by = "boatID") %>%
    gather(type, value, -low, -mid, -high, -boatID, -paintType, -boatType) %>%
    mutate(type = recode(type,
                         days1 = "Days since last used",
                         days2 = "Days since last cleaned",
                         midTrips = "Median number of trips",
                         ApproxHullSA = "Hull surface area")) %>%
    rename(`Vessel type` = boatType, `Paint type` = paintType)
slopes <- extract(m1, c("betaDays1", "betaDays2", "betaMidTrips",
                        "betaHullSA"))
slopes <- sapply(slopes, quantile, probs = 0.5)
ints <- extract(m1, "betaType")$betaType
ints <- apply(ints, 2, quantile, probs = 0.5)
coefs <- tibble(slope = rep(slopes, 3),
                type = rep(unique(a1Dat$type), 3),
                intercept = rep(ints, each = 4),
                boatTypeInt = rep(1:3, each = 4)) %>%
    left_join(., bLookup %>% select(contains("boat")) %>% distinct(),
              by = "boatTypeInt") %>%
    rename(`Vessel type` = boatType)
pl <- ggplot(a1Dat, aes(x = value, y = mid, ymin = low, ymax = high,
                        colour = `Vessel type`, shape = `Vessel type`)) +
    ## geom_pointrange(fatten = 0.5) +
    geom_linerange(size = 0.5, alpha = 0.5) +
    geom_point(size = 1.5) +
    facet_wrap(~ type) +
    geom_abline(aes(slope = slope, intercept = intercept,
                    colour = `Vessel type`), data = coefs,
                size = 0.5, alpha = 0.5) +
    ylab(expression(paste("Vessel-level intercept, ",
                          gamma[j],"*"))) +
    xlab("Scaled value") +
    scale_colour_brewer(palette = "Dark2") +
    theme_bw(base_size = 7.7) +
    theme(legend.position = "bottom")
ggsave("../graphics/plM1boat-full.pdf", pl, width = 4.9, height = 4.9)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Paint type
################################################################################
################################################################################
ints <- extract(m1, "betaPaint")$betaPaint
ints <- apply(ints, 2, quantile, probs = 0.5)
coefs <- tibble(slope = rep(slopes, 3),
                type = rep(unique(a1Dat$type), 3),
                intercept = rep(ints, each = 4),
                paintTypeInt = rep(1:3, each = 4)) %>%
    left_join(., bLookup %>% select(contains("paint")) %>% distinct(),
              by = "paintTypeInt") %>%
    rename(`Paint type` = paintType)
pl <- ggplot(a1Dat, aes(x = value, y = mid, ymin = low, ymax = high,
                        colour = `Paint type`, shape = `Paint type`)) +
    ## geom_pointrange(fatten = 0.5) +
    geom_linerange(size = 0.5, alpha = 0.5) +
    geom_point(size = 1.5) +
    facet_wrap(~ type) +
    geom_abline(aes(slope = slope, intercept = intercept,
                    colour = `Paint type`), data = coefs,
                size = 0.5, alpha = 0.5) +
    ylab(expression(paste("Vessel-level intercept, ",
                          gamma[j],"*"))) +
    xlab("Scaled value") +
    scale_colour_brewer(palette = "Dark2") +
    theme_bw(base_size = 7.7) +
    theme(legend.position = "bottom")
ggsave("../graphics/plM1paint-full.pdf", pl, width = 4.9, height = 4.9)
################################################################################
################################################################################
