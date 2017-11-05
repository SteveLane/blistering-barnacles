################################################################################
################################################################################
## Title: Post-processing comparison
## Author: Steve Lane
## Date: Thursday, 11 May 2017
## Synopsis: Performs comparison between outcome models.
## Time-stamp: <2017-11-05 23:07:37 (overlordR)>
################################################################################
################################################################################
## Add github packages using gitname/reponame format
source("../scripts/imputation-functions.R")
packages <- c("tidyr", "dplyr", "tibble", "rstan", "loo", "ggplot2",
              "RColorBrewer")
options(mc.cores = 1)
options(loo.cores = 1)
ipak(packages)
## This is set for readable text when included at half page width.
theme_set(theme_bw())
biofoul <- readRDS("../data/biofouling.rds")
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: looic table for t model
################################################################################
################################################################################
m0 <- readRDS("../data/censored-mle-m0-t.rds")
m0ll <- loo::extract_log_lik(m0)
m0loo <- loo::loo(m0ll)
message("m0 finished.")
rm(m0ll)
m1 <- readRDS("../data/censored-mle-m1-t.rds")
m1ll <- loo::extract_log_lik(m1)
m1loo <- loo::loo(m1ll)
message("m1 finished.")
rm(m1, m1ll)
m2 <- readRDS("../data/censored-mle-m2-t.rds")
m2ll <- loo::extract_log_lik(m2)
m2loo <- loo::loo(m2ll)
message("m2 finished.")
rm(m2, m2ll)
m3 <- readRDS("../data/censored-mle-m3-t.rds")
m3ll <- loo::extract_log_lik(m3)
m3loo <- loo::loo(m3ll)
message("m3 finished.")
rm(m3ll)
looTab <- as.data.frame(compare(m0loo, m1loo, m2loo, m3loo)) %>%
    mutate(mID = rownames(.))
names(looTab) <- c("LOOIC", "se(LOOIC)", "ELPD", "se(ELPD)", "Eff. P",
                   "se(Eff. P)", "mID")
looLookup <- tibble(
    mID = paste0("m", 0:3, "loo"),
    Model = paste0("M", 0:3)
)
looTab <- left_join(looTab, looLookup) %>%
    select(-mID)
message("Finished t section.")
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: looic for normal model
################################################################################
################################################################################
m0N <- readRDS("../data/censored-mle-m0.rds")
m0Nll <- loo::extract_log_lik(m0N)
m0Nloo <- loo::loo(m0Nll)
message("m0 finished.")
rm(m0Nll)
m1N <- readRDS("../data/censored-mle-m1.rds")
m1Nll <- loo::extract_log_lik(m1N)
m1Nloo <- loo::loo(m1Nll)
message("m1 finished.")
rm(m1N, m1Nll)
m2N <- readRDS("../data/censored-mle-m2.rds")
m2Nll <- loo::extract_log_lik(m2N)
m2Nloo <- loo::loo(m2Nll)
message("m2 finished.")
rm(m2N, m2Nll)
m3N <- readRDS("../data/censored-mle-m3.rds")
m3Nll <- loo::extract_log_lik(m3N)
m3Nloo <- loo::loo(m3Nll)
message("m3 finished.")
rm(m3Nll)
looTabN <- as.data.frame(compare(m0Nloo, m1Nloo, m2Nloo, m3Nloo)) %>%
    mutate(mID = rownames(.))
names(looTabN) <- c("LOOIC", "se(LOOIC)", "ELPD", "se(ELPD)", "Eff. P",
                    "se(Eff. P)", "mID")
looLookupN <- tibble(
    mID = paste0("m", 0:3, "Nloo"),
    Model = paste0("M", 0:3)
)
looTabN <- left_join(looTabN, looLookupN) %>%
    select(-mID)
message("Finished normal section.")
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: Join up for output.
################################################################################
################################################################################
compTab <- left_join(looTab, looTabN, by = "Model")
diffs <- rbind(
    compare(m3Nloo, m3loo),
    compare(m2Nloo, m2loo),
    compare(m1Nloo, m1loo),
    compare(m0Nloo, m0loo)
) %>%
    as_tibble() %>%
    mutate(Model = paste0("M", 3:0))
compTab <- left_join(compTab, diffs)
compTab$elpd_diff <- 2*compTab$elpd_diff
compTab$se <- sqrt(2)*compTab$se
saveRDS(compTab, "../data/looic-compare-full.rds")
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: PPCs for M0: prop < cens; median; iqr
################################################################################
################################################################################
yPPC <- extract(m0, "y_ppc")$y_ppc
yPPCExp <- exp(yPPC)
ppc1 <- rowMeans(yPPC < log(1.5))
ppc2 <- apply(yPPCExp, 1, quantile, probs = 0.5, names = FALSE)
ppc3 <- apply(yPPCExp, 1, IQR)
ppc <- tibble(
    `Prop(hat(Y)<1.5)` = ppc1,
    `Median(hat(Y))` = ppc2,
    `IQR(hat(Y))` = ppc3,
    model = "O2"
)
yPPCN <- extract(m0N, "y_ppc")$y_ppc
yPPCNExp <- exp(yPPCN)
ppc1 <- rowMeans(yPPCN < log(1.5))
ppc2 <- apply(yPPCNExp, 1, quantile, probs = 0.5, names = FALSE)
ppc3 <- apply(yPPCNExp, 1, IQR)
ppcN <- tibble(
    `Prop(hat(Y)<1.5)` = ppc1,
    `Median(hat(Y))` = ppc2,
    `IQR(hat(Y))` = ppc3,
    model = "O1"
)
allPPC <- bind_rows(ppc, ppcN) %>%
    gather(ppc, value, -model)
## Observed data
obs <- tibble(ppc = c("Prop(hat(Y)<1.5)", "Median(hat(Y))", "IQR(hat(Y))"),
              value = with(biofoul,
                           c(mean(wetWeight < 1.5), median(wetWeight),
                             IQR(wetWeight))))
plPPC <- ggplot(allPPC, aes(x = value)) +
    geom_histogram() +
    facet_grid(model ~ ppc, scales = "free", labeller = label_parsed) +
    geom_vline(aes(xintercept = value), data = obs) +
    xlab("") +
    ylab("Count") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
ggsave("../graphics/ppc-compare-m0-full.pdf", plPPC)
################################################################################
################################################################################

################################################################################
################################################################################
## Begin Section: PPCs for M3: prop < cens; median; iqr
################################################################################
################################################################################
yPPC <- extract(m3, "y_ppc")$y_ppc
yPPCExp <- exp(yPPC)
ppc1 <- rowMeans(yPPC < log(1.5))
ppc2 <- apply(yPPCExp, 1, quantile, probs = 0.5, names = FALSE)
ppc3 <- apply(yPPCExp, 1, IQR)
ppc <- tibble(
    `Prop(hat(Y)<1.5)` = ppc1,
    `Median(hat(Y))` = ppc2,
    `IQR(hat(Y))` = ppc3,
    model = "O2"
)
yPPCN <- extract(m3N, "y_ppc")$y_ppc
yPPCNExp <- exp(yPPCN)
ppc1 <- rowMeans(yPPCN < log(1.5))
ppc2 <- apply(yPPCNExp, 1, quantile, probs = 0.5, names = FALSE)
ppc3 <- apply(yPPCNExp, 1, IQR)
ppcN <- tibble(
    `Prop(hat(Y)<1.5)` = ppc1,
    `Median(hat(Y))` = ppc2,
    `IQR(hat(Y))` = ppc3,
    model = "O1"
)
allPPC <- bind_rows(ppc, ppcN) %>%
    gather(ppc, value, -model)
plPPC <- ggplot(allPPC, aes(x = value)) +
    geom_histogram() +
    facet_grid(model ~ ppc, scales = "free", labeller = label_parsed) +
    geom_vline(aes(xintercept = value), data = obs) +
    xlab("") +
    ylab("Count") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
ggsave("../graphics/ppc-compare-m3-full.pdf", plPPC)
################################################################################
################################################################################
