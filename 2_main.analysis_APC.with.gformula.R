library(foreach)
library("doParallel")
library("foreign")
library(dplyr)
library(tidyselect)
library("lattice")
library(Epi)
library(ggplot2)
library(nnet)
library(MASS)
library(Hmisc)
library(cfdecomp)
Sys.setenv(LANG="EN") 

load("APC_df.RData") 

## load MEDIATION ANALYSIS function
##function works for every mediator
source("functions.R")

# create simulated dataframe ----------------------------------------------

## 10 rows for each combination of age period and cohort for prediction later
sim <- expand.grid(AGEY_M=sort(unique(APC_df$AGEY_M)), 
                   period=sort(unique(APC_df$period)), 
                   birthcohort2=sort(unique(APC_df$birthcohort2)))
#repeat 10 times
sim <- sim[rep(seq_len(nrow(sim)), each = 10), ]

## allow the confounders to vary by birth cohort
## find % males and females by birthcohort; % education; % ethnicity 

prop.sex  <- t(prop.table(table(APC_df$RAGENDER,APC_df$birthcohort2), margin=2))
prop.educ <- t(prop.table(table(APC_df$RAEDUC,APC_df$birthcohort2), margin=2))
prop.ethn <- t(prop.table(table(APC_df$ethn,APC_df$birthcohort2), margin=2))

## rep 51 times 51 is the number of unique birth cohorts we have
RAGENDER <- RAEDUC <- ethn <- rep(1,51)

## randomly sample confounder values and assign them to birth cohort rows
set.seed(123)
for (i in 1:51){
  for (j in seq(1,510,10)){
    RAGENDER[j:(j+9)] <- rbinom(n = 10, size = 1, prob = prop.sex[i,])
    
    sim.educ    <- t(rmultinom(n = 10, size = 1, prob = prop.educ[i,]))
    RAEDUC[j:(j+9)] <- colnames(sim.educ)[max.col(sim.educ)]
    
    sim.ethn     <- t(rmultinom(n = 10, size = 1, prob = prop.ethn[i,]))
    ethn[j:(j+9)] <- colnames(sim.ethn)[max.col(sim.ethn)]
  }}

birthcohort <- rep(sort(unique(APC_df$birthcohort2)),each=10)

sim1 <- as.data.frame(cbind(birthcohort, RAGENDER, RAEDUC, ethn))
sim2 <- sim1[rep(seq_len(nrow(sim1)), times = 341), ]
sim2 <- sim2[order(sim2$birthcohort),]

sim3 <- cbind(sim, sim2[,-1])

sim3$RAGENDER <- as.factor(sim3$RAGENDER)
levels(sim3$RAGENDER) <- c("1.Male", "2.Female")
sim3$RAEDUC <- as.factor(sim3$RAEDUC)
sim3$ethn <- as.factor(sim3$ethn)

## copy simulated dataframe and actual dataframe
DRINK.sim <- SMOK.sim <- BMI.sim <- PHYSACT.sim <- sim3
df.DRINK <- df.SMOK <- df.BMI <- df.PHYSACT <- APC_df


# run main analysis --------------------------------------------------
## run in parallel 
cl <- makeCluster(10) #use 10 cores
registerDoParallel(cl) 

## specify bootstrap and monte carlo iterations
bs_it = 499
mc_it = 50

set.seed(123)

## run mediation analysis for smoking
start.time <- Sys.time()

SMOK <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("cfdecomp","Epi")) %dopar% {
  mediation_analysis(mediator="SMOKN", dataframe=df.SMOK, sim.df=SMOK.sim, mc_iteration=mc_it, reference_cohort = 1945) #default ref.cohort = 1945
}
end.time <- Sys.time()
end.time-start.time

SMOK.nc <- SMOK[1:(nrow(SMOK)/2),] #first half is natural course, second half is counterfactual
SMOK.cf <- SMOK[(nrow(SMOK)/2+1):nrow(SMOK),]

## DRINK
start.time <- Sys.time()
DRINK <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("cfdecomp","nnet","Epi","Hmisc")) %dopar% {
  mediation_analysis(mediator="DRINK_cat", dataframe=df.DRINK, sim.df=DRINK.sim, mc_iteration=mc_it, reference_cohort = 1945) #default ref.cohort = 1945
}
end.time <- Sys.time()
end.time-start.time

DRINK.nc <- DRINK[1:(nrow(DRINK)/2),] #first half is nc, second half is cf
DRINK.cf <- DRINK[(nrow(DRINK)/2+1):nrow(DRINK),]

## BMI
start.time <- Sys.time()
BMI <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("cfdecomp","nnet","Epi","Hmisc")) %dopar% {
  mediation_analysis(mediator="BMI_cat", dataframe=df.BMI, sim.df=BMI.sim, mc_iteration=mc_it, reference_cohort = 1945) #default ref.cohort = 1945
}
end.time <- Sys.time()
end.time-start.time

BMI.nc <- BMI[1:(nrow(BMI)/2),] #first half is nc, second half is cf
BMI.cf <- BMI[(nrow(BMI)/2+1):nrow(BMI),]

## PHYSICAL ACTIVITY
start.time <- Sys.time()
PHYSACT <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("cfdecomp","Epi","Hmisc")) %dopar% {
  mediation_analysis(mediator="phys_act", dataframe=df.PHYSACT, sim.df=PHYSACT.sim, mc_iteration=mc_it, reference_cohort = 1945) #default ref.cohort = 1945
}
end.time <- Sys.time()
end.time-start.time

PHYSACT.nc <- PHYSACT[1:(nrow(PHYSACT)/2),] #first half is nc, second half is cf
PHYSACT.cf <- PHYSACT[(nrow(PHYSACT)/2+1):nrow(PHYSACT),]

stopCluster(cl) # stop cluster

save(list=c("SMOK.nc", "SMOK.cf",
            "DRINK.nc", "DRINK.cf",
            "BMI.nc", "BMI.cf",
            "PHYSACT.nc", "PHYSACT.cf"),
     file = paste0("Mediationresults_BS499.MC50",Sys.Date(),".Rdata"))

rm(SMOK, DRINK, PHYSACT, BMI, SMOK.nc, SMOK.cf, DRINK.nc, DRINK.cf, BMI.nc, BMI.cf, PHYSACT.nc, PHYSACT.cf)


# Stratified by gender ----------------------------------------------------

bs_it = 499
mc_it = 50

cl <- makeCluster(15) 
registerDoParallel(cl) 
set.seed(123)

#SMOK
start.time <- Sys.time()

SMOK <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("dplyr","cfdecomp","Epi")) %dopar% {
  mediation_analysis(mediator="SMOKN", dataframe=df.SMOK, sim.df=SMOK.sim, mc_iteration=mc_it, reference_cohort = 1945, strata = "sex") #default ref.cohort = 1945
}
end.time <- Sys.time()
end.time-start.time

## save results for females and males for natural course and counterfactual
SMOK.F.nc <- SMOK[1:(nrow(SMOK)/4),]
SMOK.F.cf <- SMOK[(nrow(SMOK)/4+1):(nrow(SMOK)/4*2),]
SMOK.M.nc <- SMOK[(nrow(SMOK)/4*2+1):(nrow(SMOK)/4*3),]
SMOK.M.cf <- SMOK[(nrow(SMOK)/4*3+1):nrow(SMOK),]

#DRINK
start.time <- Sys.time()

DRINK <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("dplyr","cfdecomp","nnet","Epi","Hmisc")) %dopar% {
  mediation_analysis(mediator="DRINK_cat", dataframe=df.DRINK, sim.df=DRINK.sim, mc_iteration=mc_it, reference_cohort = 1945, strata = "sex") #default ref.cohort = 1945
}
end.time <- Sys.time()
end.time-start.time

DRINK.F.nc <- DRINK[1:(nrow(DRINK)/4),]
DRINK.F.cf <- DRINK[(nrow(DRINK)/4+1):(nrow(DRINK)/4*2),]
DRINK.M.nc <- DRINK[(nrow(DRINK)/4*2+1):(nrow(DRINK)/4*3),]
DRINK.M.cf <- DRINK[(nrow(DRINK)/4*3+1):nrow(DRINK),]

#BMI
start.time <- Sys.time()

BMI <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("dplyr","cfdecomp","nnet","Epi","Hmisc")) %dopar% {
  mediation_analysis(mediator="BMI_cat", dataframe=df.BMI, sim.df=BMI.sim, mc_iteration=mc_it, reference_cohort = 1945, strata = "sex") #default ref.cohort = 1945
}
end.time <- Sys.time()
end.time-start.time

BMI.F.nc <- BMI[1:(nrow(BMI)/4),]
BMI.F.cf <- BMI[(nrow(BMI)/4+1):(nrow(BMI)/4*2),]
BMI.M.nc <- BMI[(nrow(BMI)/4*2+1):(nrow(BMI)/4*3),]
BMI.M.cf <- BMI[(nrow(BMI)/4*3+1):nrow(BMI),]

#PHYSICAL ACTIVITY
start.time <- Sys.time()

PHYSACT <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("dplyr","cfdecomp","Epi","Hmisc")) %dopar% {
  mediation_analysis(mediator="phys_act", dataframe=df.PHYSACT, sim.df=PHYSACT.sim, mc_iteration=mc_it, reference_cohort = 1945, strata = "sex") #default ref.cohort = 1945
}
end.time <- Sys.time()
end.time-start.time

PHYSACT.F.nc <- PHYSACT[1:(nrow(PHYSACT)/4),]
PHYSACT.F.cf <- PHYSACT[(nrow(PHYSACT)/4+1):(nrow(PHYSACT)/4*2),]
PHYSACT.M.nc <- PHYSACT[(nrow(PHYSACT)/4*2+1):(nrow(PHYSACT)/4*3),]
PHYSACT.M.cf <- PHYSACT[(nrow(PHYSACT)/4*3+1):nrow(PHYSACT),]

stopCluster(cl)

save(list=c("SMOK.F.nc", "SMOK.F.cf", "SMOK.M.nc", "SMOK.M.cf",
            "DRINK.F.nc", "DRINK.F.cf", "DRINK.M.nc", "DRINK.M.cf",
            "BMI.F.nc", "BMI.F.cf", "BMI.M.nc", "BMI.M.cf",
            "PHYSACT.F.nc", "PHYSACT.F.cf", "PHYSACT.M.nc", "PHYSACT.M.cf"),
            file = paste0("Mediationresults_bysex_BS499.MC50",Sys.Date(),".Rdata"))

rm(SMOK, DRINK, PHYSACT, BMI, SMOK.F.nc, SMOK.F.cf, SMOK.M.nc, SMOK.M.cf,
   DRINK.F.nc, DRINK.F.cf, DRINK.M.nc, DRINK.M.cf,
   BMI.F.nc, BMI.F.cf, BMI.M.nc, BMI.M.cf,
   PHYSACT.F.nc, PHYSACT.F.cf, PHYSACT.M.nc, PHYSACT.M.cf)


# Stratified by ethnicity -------------------------------------------------

bs_it = 499
mc_it = 50

cl <- makeCluster(15) 
registerDoParallel(cl) 

#SMOK
set.seed(123)

start.time <- Sys.time()
SMOK <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("dplyr","cfdecomp","Epi")) %dopar% {
  mediation_analysis(mediator="SMOKN", dataframe=df.SMOK, sim.df=SMOK.sim, mc_iteration=mc_it, reference_cohort = 1945, strata = "ethnicity") #default ref.cohort = 1945
}
end.time <- Sys.time()
end.time-start.time

## save results for each race/ethn group for natural course and counterfactual
SMOK.White.nc <- SMOK[1:(nrow(SMOK)/8),]
SMOK.White.cf <- SMOK[(nrow(SMOK)/8+1):(nrow(SMOK)/8*2),]
SMOK.Black.nc <- SMOK[(nrow(SMOK)/8*2+1):(nrow(SMOK)/8*3),]
SMOK.Black.cf <- SMOK[(nrow(SMOK)/8*3+1):(nrow(SMOK)/8*4),]
SMOK.Hisp.nc <- SMOK[(nrow(SMOK)/8*4+1):(nrow(SMOK)/8*5),]
SMOK.Hisp.cf <- SMOK[(nrow(SMOK)/8*5+1):(nrow(SMOK)/8*6),]
SMOK.other.nc <- SMOK[(nrow(SMOK)/8*6+1):(nrow(SMOK)/8*7),]
SMOK.other.cf <- SMOK[(nrow(SMOK)/8*7+1):(nrow(SMOK)),]

#DRINK
start.time <- Sys.time()

DRINK <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("dplyr","cfdecomp","nnet","Epi","Hmisc")) %dopar% {
  mediation_analysis(mediator="DRINK_cat", dataframe=df.DRINK, sim.df=DRINK.sim, mc_iteration=mc_it, reference_cohort = 1945, strata = "ethnicity") #default ref.cohort = 1945
}
end.time <- Sys.time()
end.time-start.time

DRINK.White.nc <- DRINK[1:(nrow(DRINK)/8),]
DRINK.White.cf <- DRINK[(nrow(DRINK)/8+1):(nrow(DRINK)/8*2),]
DRINK.Black.nc <- DRINK[(nrow(DRINK)/8*2+1):(nrow(DRINK)/8*3),]
DRINK.Black.cf <- DRINK[(nrow(DRINK)/8*3+1):(nrow(DRINK)/8*4),]
DRINK.Hisp.nc <- DRINK[(nrow(DRINK)/8*4+1):(nrow(DRINK)/8*5),]
DRINK.Hisp.cf <- DRINK[(nrow(DRINK)/8*5+1):(nrow(DRINK)/8*6),]
DRINK.other.nc <- DRINK[(nrow(DRINK)/8*6+1):(nrow(DRINK)/8*7),]
DRINK.other.cf <- DRINK[(nrow(DRINK)/8*7+1):(nrow(DRINK)),]


#BMI
start.time <- Sys.time()

BMI <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("dplyr","cfdecomp","nnet","Epi","Hmisc")) %dopar% {
  mediation_analysis(mediator="BMI_cat", dataframe=df.BMI, sim.df=BMI.sim, mc_iteration=mc_it, reference_cohort = 1945, strata = "ethnicity") #default ref.cohort = 1945
}
end.time <- Sys.time()
end.time-start.time

BMI.White.nc <- BMI[1:(nrow(BMI)/8),]
BMI.White.cf <- BMI[(nrow(BMI)/8+1):(nrow(BMI)/8*2),]
BMI.Black.nc <- BMI[(nrow(BMI)/8*2+1):(nrow(BMI)/8*3),]
BMI.Black.cf <- BMI[(nrow(BMI)/8*3+1):(nrow(BMI)/8*4),]
BMI.Hisp.nc <- BMI[(nrow(BMI)/8*4+1):(nrow(BMI)/8*5),]
BMI.Hisp.cf <- BMI[(nrow(BMI)/8*5+1):(nrow(BMI)/8*6),]
BMI.other.nc <- BMI[(nrow(BMI)/8*6+1):(nrow(BMI)/8*7),]
BMI.other.cf <- BMI[(nrow(BMI)/8*7+1):(nrow(BMI)),]


#PHYSICAL ACTIVITY
start.time <- Sys.time()

PHYSACT <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("dplyr","cfdecomp","Epi","Hmisc")) %dopar% {
  mediation_analysis(mediator="phys_act", dataframe=df.PHYSACT, sim.df=PHYSACT.sim, mc_iteration=mc_it, reference_cohort = 1945, strata = "ethnicity") #default ref.cohort = 1945
}
end.time <- Sys.time()
end.time-start.time

PHYSACT.White.nc <- PHYSACT[1:(nrow(PHYSACT)/8),]
PHYSACT.White.cf <- PHYSACT[(nrow(PHYSACT)/8+1):(nrow(PHYSACT)/8*2),]
PHYSACT.Black.nc <- PHYSACT[(nrow(PHYSACT)/8*2+1):(nrow(PHYSACT)/8*3),]
PHYSACT.Black.cf <- PHYSACT[(nrow(PHYSACT)/8*3+1):(nrow(PHYSACT)/8*4),]
PHYSACT.Hisp.nc <- PHYSACT[(nrow(PHYSACT)/8*4+1):(nrow(PHYSACT)/8*5),]
PHYSACT.Hisp.cf <- PHYSACT[(nrow(PHYSACT)/8*5+1):(nrow(PHYSACT)/8*6),]
PHYSACT.other.nc <- PHYSACT[(nrow(PHYSACT)/8*6+1):(nrow(PHYSACT)/8*7),]
PHYSACT.other.cf <- PHYSACT[(nrow(PHYSACT)/8*7+1):(nrow(PHYSACT)),]


stopCluster(cl)

save(list=c("SMOK.White.nc", "SMOK.White.cf", "SMOK.Black.nc", "SMOK.Black.cf",
            "SMOK.Hisp.nc", "SMOK.Hisp.cf", "SMOK.other.nc", "SMOK.other.cf",
            "DRINK.White.nc", "DRINK.White.cf", "DRINK.Black.nc", "DRINK.Black.cf",
            "DRINK.Hisp.nc", "DRINK.Hisp.cf", "DRINK.other.nc", "DRINK.other.cf",
            "BMI.White.nc", "BMI.White.cf", "BMI.Black.nc", "BMI.Black.cf",
            "BMI.Hisp.nc", "BMI.Hisp.cf", "BMI.other.nc", "BMI.other.cf",
            "PHYSACT.White.nc", "PHYSACT.White.cf", "PHYSACT.Black.nc", "PHYSACT.Black.cf",
            "PHYSACT.Hisp.nc", "PHYSACT.Hisp.cf", "PHYSACT.other.nc", "PHYSACT.other.cf"),
     file = paste0("Mediationresults_by.ethn_BS.499.MC50",Sys.Date(),".Rdata"))

rm(SMOK, DRINK, PHYSACT, BMI, SMOK.White.nc, SMOK.White.cf, SMOK.Black.nc, SMOK.Black.cf,
   SMOK.Hisp.nc, SMOK.Hisp.cf, SMOK.other.nc, SMOK.other.cf,
   DRINK.White.nc, DRINK.White.cf, DRINK.Black.nc, DRINK.Black.cf,
   DRINK.Hisp.nc, DRINK.Hisp.cf, DRINK.other.nc, DRINK.other.cf,
   BMI.White.nc, BMI.White.cf, BMI.Black.nc, BMI.Black.cf,
   BMI.Hisp.nc, BMI.Hisp.cf, BMI.other.nc, BMI.other.cf,
   PHYSACT.White.nc, PHYSACT.White.cf, PHYSACT.Black.nc, PHYSACT.Black.cf,
   PHYSACT.Hisp.nc, PHYSACT.Hisp.cf, PHYSACT.other.nc, PHYSACT.other.cf)
