library(dplyr)

load("APC_df.RData") 

source("functions.R")

sim <- expand.grid(AGEY_M=sort(unique(APC_df$AGEY_M)), 
                   period=sort(unique(APC_df$period)), 
                   birthcohort2=sort(unique(APC_df$birthcohort2)))
#repeat 10 times
sim <- sim[rep(seq_len(nrow(sim)), each = 10), ]

###main analysis----
load("Mediationresults_aggregated.RData")

## make tables
SMOK <- decomp_table(SMOK.BS.nc2, SMOK.BS.cf2)
SMOK2 <- rel.diff_contr_table(SMOK.BS.nc2, SMOK.BS.cf2)

BMI <- decomp_table(BMI.BS.nc2, BMI.BS.cf2)
BMI2 <- rel.diff_contr_table(BMI.BS.nc2, BMI.BS.cf2)

DRINK <- decomp_table(DRINK.BS.nc2, DRINK.BS.cf2)
DRINK2 <- rel.diff_contr_table(DRINK.BS.nc2, DRINK.BS.cf2)

PHYSACT <- decomp_table(PHYSACT.BS.nc2, PHYSACT.BS.cf2)
PHYSACT2 <- rel.diff_contr_table(PHYSACT.BS.nc2, PHYSACT.BS.cf2)

###BY SEX----
load("Mediationresults_bysex_aggregated.RData")

## make tables
SMOK <- decomp_table(SMOK.BS.nc2, SMOK.BS.cf2, strata="sex")
SMOK2 <- rel.diff_contr_table(SMOK.BS.nc2, SMOK.BS.cf2, strata="sex")

BMI <- decomp_table(BMI.BS.nc2, BMI.BS.cf2, strata="sex")
BMI2 <- rel.diff_contr_table(BMI.BS.nc2, BMI.BS.cf2, strata="sex")

DRINK <- decomp_table(DRINK.BS.nc2, DRINK.BS.cf2, strata="sex")
DRINK2 <- rel.diff_contr_table(DRINK.BS.nc2, DRINK.BS.cf2, strata="sex")

PHYSACT <- decomp_table(PHYSACT.BS.nc2, PHYSACT.BS.cf2, strata="sex")
PHYSACT2 <- rel.diff_contr_table(PHYSACT.BS.nc2, PHYSACT.BS.cf2, strata="sex")

####by ethnicity----
load("Mediationresults_byethn_aggregated.RData")

## make tables
SMOK <- decomp_table(SMOK.BS.nc2, SMOK.BS.cf2, strata="ethn")
SMOK2 <- rel.diff_contr_table(SMOK.BS.nc2, SMOK.BS.cf2, strata="ethn")

BMI <- decomp_table(BMI.BS.nc2, BMI.BS.cf2, strata="ethn")
BMI2 <- rel.diff_contr_table(BMI.BS.nc2, BMI.BS.cf2, strata="ethn")

DRINK <- decomp_table(DRINK.BS.nc2, DRINK.BS.cf2, strata="ethn")
DRINK2 <- rel.diff_contr_table(DRINK.BS.nc2, DRINK.BS.cf2, strata="ethn")

PHYSACT <- decomp_table(PHYSACT.BS.nc2, PHYSACT.BS.cf2, strata="ethn")
PHYSACT2 <- rel.diff_contr_table(PHYSACT.BS.nc2, PHYSACT.BS.cf2, strata="ethn")

