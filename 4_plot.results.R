library(gridExtra)
library(dplyr)
library(ggplot2)
library(grid)
library(svglite)
source("functions.R")

load("APC_df.RData") 

## create simulated df 
sim <- expand.grid(AGEY_M=sort(unique(APC_df$AGEY_M)), 
                   period=sort(unique(APC_df$period)), 
                   birthcohort2=sort(unique(APC_df$birthcohort2)))
#repeat 10 times
sim <- sim[rep(seq_len(nrow(sim)), each = 10), ]


# specify bootstrap iterations
bs_it=499


# Plot main results -------------------------------------------------------

## load aggregated output for each health behavior
load("Mediationresults_aggregated_2023-03-07.RData")

##calculate relative difference as (cf/nc)-1
## smoking
BS.diff.SMOK <- (SMOK.BS.cf2[,4:ncol(SMOK.BS.cf2)] / SMOK.BS.nc2[,4:ncol(SMOK.BS.nc2)])-1
yintercept_temp <- 0

BS.diff.SMOK2 <- BS.diff.SMOK
BS.diff.SMOK2$lb <- apply(BS.diff.SMOK, 1, quantile, probs = c(0.025))
BS.diff.SMOK2$ub <- apply(BS.diff.SMOK, 1, quantile, probs = c(0.975))
BS.diff.SMOK2$diff <- rowMeans(BS.diff.SMOK)

reldiff.SMOK <- cbind(SMOK.BS.nc2[,1:3], BS.diff.SMOK2)

#set age and period to the reference
reldiff.SMOK1 <- subset(reldiff.SMOK, reldiff.SMOK$AGEY_M == 50 &
                          reldiff.SMOK$period == 1996)

#keep only aggregated results
reldiff.SMOK2 <- reldiff.SMOK1[,c(1:3,(bs_it+4):(bs_it+6))]

## alcohol consumption
BS.diff.DRINK <- (DRINK.BS.cf2[,4:ncol(DRINK.BS.cf2)] / DRINK.BS.nc2[,4:ncol(DRINK.BS.nc2)])-1

BS.diff.DRINK2 <- BS.diff.DRINK
BS.diff.DRINK2$lb <- apply(BS.diff.DRINK, 1, quantile, probs = c(0.025))
BS.diff.DRINK2$ub <- apply(BS.diff.DRINK, 1, quantile, probs = c(0.975))
BS.diff.DRINK2$diff <- rowMeans(BS.diff.DRINK)

#keep only aggregated results
reldiff.DRINK <- cbind(DRINK.BS.nc2[,1:3], BS.diff.DRINK2)

#set age and period to the reference
reldiff.DRINK1 <- subset(reldiff.DRINK, reldiff.DRINK$AGEY_M == 50 &
                          reldiff.DRINK$period == 1996)
reldiff.DRINK2 <- reldiff.DRINK1[,c(1:3,(bs_it+4):(bs_it+6))]

## physical activity
BS.diff.PHYSACT <- (PHYSACT.BS.cf2[,4:ncol(PHYSACT.BS.cf2)] / PHYSACT.BS.nc2[,4:ncol(PHYSACT.BS.nc2)])-1

BS.diff.PHYSACT2 <- BS.diff.PHYSACT
BS.diff.PHYSACT2$lb <- apply(BS.diff.PHYSACT, 1, quantile, probs = c(0.025))
BS.diff.PHYSACT2$ub <- apply(BS.diff.PHYSACT, 1, quantile, probs = c(0.975))
BS.diff.PHYSACT2$diff <- rowMeans(BS.diff.PHYSACT)

reldiff.PHYSACT <- cbind(PHYSACT.BS.nc2[,1:3], BS.diff.PHYSACT2)
#set age and period to the reference
reldiff.PHYSACT1 <- subset(reldiff.PHYSACT, reldiff.PHYSACT$AGEY_M == 50 &
                           reldiff.PHYSACT$period == 1996)
#keep only aggregated results
reldiff.PHYSACT2 <- reldiff.PHYSACT1[,c(1:3,(bs_it+4):(bs_it+6))]

## BMI
BS.diff.BMI <- (BMI.BS.cf2[,4:ncol(BMI.BS.cf2)] / BMI.BS.nc2[,4:ncol(BMI.BS.nc2)])-1

BS.diff.BMI2 <- BS.diff.BMI
BS.diff.BMI2$lb <- apply(BS.diff.BMI, 1, quantile, probs = c(0.025))
BS.diff.BMI2$ub <- apply(BS.diff.BMI, 1, quantile, probs = c(0.975))
BS.diff.BMI2$diff <- rowMeans(BS.diff.BMI)

reldiff.BMI <- cbind(BMI.BS.nc2[,1:3], BS.diff.BMI2)
#set age and period to the reference
reldiff.BMI1 <- subset(reldiff.BMI, reldiff.BMI$AGEY_M == 50 &
                           reldiff.BMI$period == 1996)
#keep only aggregated results
reldiff.BMI2 <- reldiff.BMI1[,c(1:3,(bs_it+4):(bs_it+6))]

## merge datasets together and plot
reldiff.DRINK2$mediator <- "Alcohol Consumption" 
reldiff.SMOK2$mediator <- "Smoking" 
reldiff.PHYSACT2$mediator <- "Physical Activity" 
reldiff.BMI2$mediator <- "BMI" 

reldiff.plot <- rbind(reldiff.DRINK2, reldiff.SMOK2, reldiff.PHYSACT2, reldiff.BMI2)
reldiff.plot$mediator <- factor(reldiff.plot$mediator,
                                   levels=c("Alcohol Consumption",
                                            "Smoking",
                                            "Physical Activity",
                                            "BMI"))
reldiff.plot.p <- reldiff.plot %>%
  ggplot(aes(x = birthcohort2, y = diff)) +
  geom_line() + 
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=13),
        title = element_text(size=13),
        strip.text = element_text(size=13)) +
  geom_hline(yintercept=yintercept_temp) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = .15) +
  labs(title = "Relative difference in depression probability", y = "Relative difference", x="Birth year") +
  facet_wrap(~mediator)

ggsave(file="RelDiff.ALLHB.svg", reldiff.plot.p)




# plot results by gender --------------------------------------------------
load("Mediationresults_bysex_aggregated.RData")

BS.diff.SMOK <- (SMOK.BS.cf2[,5:ncol(SMOK.BS.cf2)] / SMOK.BS.nc2[,5:ncol(SMOK.BS.nc2)])-1

BS.diff.SMOK2 <- BS.diff.SMOK
BS.diff.SMOK2$lb <- apply(BS.diff.SMOK, 1, quantile, probs = c(0.025))
BS.diff.SMOK2$ub <- apply(BS.diff.SMOK, 1, quantile, probs = c(0.975))
BS.diff.SMOK2$diff <- rowMeans(BS.diff.SMOK)

reldiff.SMOK <- cbind(SMOK.BS.nc2[,1:4], BS.diff.SMOK2)
reldiff.SMOK1 <- subset(reldiff.SMOK, reldiff.SMOK$AGEY_M == 50 &
                          reldiff.SMOK$period == 1996)
reldiff.SMOK2 <- reldiff.SMOK1[,c(1:4,(bs_it+5):(bs_it+7))]

BS.diff.DRINK <- (DRINK.BS.cf2[,5:ncol(DRINK.BS.cf2)] / DRINK.BS.nc2[,5:ncol(DRINK.BS.nc2)])-1

BS.diff.DRINK2 <- BS.diff.DRINK
BS.diff.DRINK2$lb <- apply(BS.diff.DRINK, 1, quantile, probs = c(0.025))
BS.diff.DRINK2$ub <- apply(BS.diff.DRINK, 1, quantile, probs = c(0.975))
BS.diff.DRINK2$diff <- rowMeans(BS.diff.DRINK)

reldiff.DRINK <- cbind(DRINK.BS.nc2[,1:4], BS.diff.DRINK2)
reldiff.DRINK1 <- subset(reldiff.DRINK, reldiff.DRINK$AGEY_M == 50 &
                          reldiff.DRINK$period == 1996)
reldiff.DRINK2 <- reldiff.DRINK1[,c(1:4,(bs_it+5):(bs_it+7))]

BS.diff.PHYSACT <- (PHYSACT.BS.cf2[,5:ncol(PHYSACT.BS.cf2)] / PHYSACT.BS.nc2[,5:ncol(PHYSACT.BS.nc2)])-1

BS.diff.PHYSACT2 <- BS.diff.PHYSACT
BS.diff.PHYSACT2$lb <- apply(BS.diff.PHYSACT, 1, quantile, probs = c(0.025))
BS.diff.PHYSACT2$ub <- apply(BS.diff.PHYSACT, 1, quantile, probs = c(0.975))
BS.diff.PHYSACT2$diff <- rowMeans(BS.diff.PHYSACT)

reldiff.PHYSACT <- cbind(PHYSACT.BS.nc2[,1:4], BS.diff.PHYSACT2)
reldiff.PHYSACT1 <- subset(reldiff.PHYSACT, reldiff.PHYSACT$AGEY_M == 50 &
                          reldiff.PHYSACT$period == 1996)
reldiff.PHYSACT2 <- reldiff.PHYSACT1[,c(1:4,(bs_it+5):(bs_it+7))]

BS.diff.BMI <- (BMI.BS.cf2[,5:ncol(BMI.BS.cf2)] / BMI.BS.nc2[,5:ncol(BMI.BS.nc2)])-1

BS.diff.BMI2 <- BS.diff.BMI
BS.diff.BMI2$lb <- apply(BS.diff.BMI, 1, quantile, probs = c(0.025))
BS.diff.BMI2$ub <- apply(BS.diff.BMI, 1, quantile, probs = c(0.975))
BS.diff.BMI2$diff <- rowMeans(BS.diff.BMI)

reldiff.BMI <- cbind(BMI.BS.nc2[,1:4], BS.diff.BMI2)
reldiff.BMI1 <- subset(reldiff.BMI, reldiff.BMI$AGEY_M == 50 &
                          reldiff.BMI$period == 1996)
reldiff.BMI2 <- reldiff.BMI1[,c(1:4,(bs_it+5):(bs_it+7))]

## merge datasets together and plot
reldiff.DRINK2$mediator <- "Alcohol Consumption" 
reldiff.SMOK2$mediator <- "Smoking" 
reldiff.PHYSACT2$mediator <- "Physical Activity" 
reldiff.BMI2$mediator <- "BMI" 

reldiff.plot <- rbind(reldiff.DRINK2, reldiff.SMOK2, reldiff.PHYSACT2, reldiff.BMI2)
reldiff.plot$mediator <- factor(reldiff.plot$mediator,
                                levels=c("Alcohol Consumption",
                                         "Smoking",
                                         "Physical Activity",
                                         "BMI"))

yintercept_temp <- 0
reldiff.plot.sex.p <- reldiff.plot %>%
  ggplot(aes(x = birthcohort2, y = diff, color=sex)) +
  geom_line() + 
  theme(legend.position = "bottom", legend.text = element_text(size = 13), 
        legend.title = element_blank(),
        axis.text = element_text(size=13),
        axis.title = element_text(size=13),
        title = element_text(size=13),
        strip.text = element_text(size=13)) +
  geom_hline(yintercept=yintercept_temp) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = .15) +
  labs(title = "Relative difference in depression probability (cf/nc-1) by sex", y = "Relative difference", x="Birth year") +
  facet_wrap(~mediator)

ggsave(file="RelDiff.ALLHB.bysex.svg", reldiff.plot.sex.p)



# plot results by race/ethnicity ------------------------------------------
load("Mediationresults_byethn_aggregated.RData")

BS.diff.SMOK <- (SMOK.BS.cf2[,5:ncol(SMOK.BS.cf2)] / SMOK.BS.nc2[,5:ncol(SMOK.BS.nc2)])-1

BS.diff.SMOK2 <- BS.diff.SMOK
BS.diff.SMOK2$lb <- apply(BS.diff.SMOK, 1, quantile, probs = c(0.025))
BS.diff.SMOK2$ub <- apply(BS.diff.SMOK, 1, quantile, probs = c(0.975))
BS.diff.SMOK2$diff <- rowMeans(BS.diff.SMOK)

reldiff.SMOK <- cbind(SMOK.BS.nc2[,1:4], BS.diff.SMOK2)
reldiff.SMOK1 <- subset(reldiff.SMOK, reldiff.SMOK$AGEY_M == 50 &
                          reldiff.SMOK$period == 1996)
reldiff.SMOK2 <- reldiff.SMOK1[,c(1:4,(bs_it+5):(bs_it+7))]

BS.diff.DRINK <- (DRINK.BS.cf2[,5:ncol(DRINK.BS.cf2)] / DRINK.BS.nc2[,5:ncol(DRINK.BS.nc2)])-1

BS.diff.DRINK2 <- BS.diff.DRINK
BS.diff.DRINK2$lb <- apply(BS.diff.DRINK, 1, quantile, probs = c(0.025))
BS.diff.DRINK2$ub <- apply(BS.diff.DRINK, 1, quantile, probs = c(0.975))
BS.diff.DRINK2$diff <- rowMeans(BS.diff.DRINK)

reldiff.DRINK <- cbind(DRINK.BS.nc2[,1:4], BS.diff.DRINK2)
reldiff.DRINK1 <- subset(reldiff.DRINK, reldiff.DRINK$AGEY_M == 50 &
                           reldiff.DRINK$period == 1996)
reldiff.DRINK2 <- reldiff.DRINK1[,c(1:4,(bs_it+5):(bs_it+7))]

BS.diff.PHYSACT <- (PHYSACT.BS.cf2[,5:ncol(PHYSACT.BS.cf2)] / PHYSACT.BS.nc2[,5:ncol(PHYSACT.BS.nc2)])-1

BS.diff.PHYSACT2 <- BS.diff.PHYSACT
BS.diff.PHYSACT2$lb <- apply(BS.diff.PHYSACT, 1, quantile, probs = c(0.025))
BS.diff.PHYSACT2$ub <- apply(BS.diff.PHYSACT, 1, quantile, probs = c(0.975))
BS.diff.PHYSACT2$diff <- rowMeans(BS.diff.PHYSACT)

reldiff.PHYSACT <- cbind(PHYSACT.BS.nc2[,1:4], BS.diff.PHYSACT2)
reldiff.PHYSACT1 <- subset(reldiff.PHYSACT, reldiff.PHYSACT$AGEY_M == 50 &
                             reldiff.PHYSACT$period == 1996)
reldiff.PHYSACT2 <- reldiff.PHYSACT1[,c(1:4,(bs_it+5):(bs_it+7))]

BS.diff.BMI <- (BMI.BS.cf2[,5:ncol(BMI.BS.cf2)] / BMI.BS.nc2[,5:ncol(BMI.BS.nc2)])-1

BS.diff.BMI2 <- BS.diff.BMI
BS.diff.BMI2$lb <- apply(BS.diff.BMI, 1, quantile, probs = c(0.025))
BS.diff.BMI2$ub <- apply(BS.diff.BMI, 1, quantile, probs = c(0.975))
BS.diff.BMI2$diff <- rowMeans(BS.diff.BMI)

reldiff.BMI <- cbind(BMI.BS.nc2[,1:4], BS.diff.BMI2)
reldiff.BMI1 <- subset(reldiff.BMI, reldiff.BMI$AGEY_M == 50 &
                         reldiff.BMI$period == 1996)
reldiff.BMI2 <- reldiff.BMI1[,c(1:4,(bs_it+5):(bs_it+7))]

## merge datasets together and plot
reldiff.DRINK2$mediator <- "Alcohol Consumption" 
reldiff.SMOK2$mediator <- "Smoking" 
reldiff.PHYSACT2$mediator <- "Physical Activity" 
reldiff.BMI2$mediator <- "BMI" 

reldiff.plot <- rbind(reldiff.DRINK2, reldiff.SMOK2, reldiff.PHYSACT2, reldiff.BMI2)
reldiff.plot$mediator <- factor(reldiff.plot$mediator,
                                levels=c("Alcohol Consumption",
                                         "Smoking",
                                         "Physical Activity",
                                         "BMI"))

cbPalette <- c("#999999", "#E69F00", "#009E73")
reldiff.plot$ethn2 <- factor(reldiff.plot$ethn, labels=c("Black", "Hispanic","White"))

yintercept_temp <- 0

reldiff.plot.ethn.p <- reldiff.plot %>%
  ggplot(aes(x = birthcohort2, y = diff, color=ethn2)) +
  geom_line() + scale_color_manual(values = cbPalette) +
  theme(legend.position = "bottom", legend.text = element_text(size = 13), 
        legend.title = element_blank(),
        axis.text = element_text(size=13),
        axis.title = element_text(size=13),
        title = element_text(size=13),
        strip.text = element_text(size=13)) +
  geom_hline(yintercept=yintercept_temp) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = .15) +
  labs(title = "Relative difference in depression probability (cf/nc-1) by race/ethnicity", y = "Relative difference", x="Birth year") +
  facet_wrap(~mediator)

ggsave(file="RelDiff.ALLHB.byethn.svg", reldiff.plot.ethn.p)
