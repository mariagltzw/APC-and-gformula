#########################################
########### DATA HANDLING FOR ########### 
########### CONTRIBUTION OF HEALTH ######
####BEHAVIORS TO COHORT DIFFERENCES #####
########## IN DEPRESSION RISK ########### 
########### GUELTZOW ET AL    ########### 
#########################################

# > sessionInfo()
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] compiler_4.1.1 tools_4.1.1     

####load packages####
library("foreign")
library(reshape2)
library(dplyr)
library(tidyselect)
library("lattice")


# EXTRACT RELEVANT VARIABLES FROM HRS DATASET ----------------------------------------------------------------------

## load HRS RAND 1922 to 2016 v2 file
df <- read.spss("randhrs1992_2016v2.sav", to.data.frame = TRUE)


## only keep respondents data
## we are not interested in spouse data 
df_resp   <- df %>%
           select(vars_select(names(df), !starts_with('S')))

## find all variables that are constant across time 
idnames   <- (df_resp[,grep("[1,2,3,4,5,6,7,8,9,10,11,12,13]", names(df_resp),
                          value=TRUE, invert=TRUE)])

## create a new dataset with variables of interest
base      <- df_resp[,c("HHIDPN","RAGENDER", "RACOHBYR", "RAHISPAN", "RARACEM",
                       "RABMONTH","RABYEAR", "RABDATE", "RABPLACE", "RADMONTH","RADYEAR", "RADDATE",
                       "RAFEDUC","RAMEDUC", "RAEDUC", "RAEDYRS", "RAEDEGRM", "RARELIG")]

## whether proxy interview
proxy <- df_resp[,grep("([P][R][O][X][Y])", names(df_resp), value=TRUE)]
proxy <- proxy[,1:13]

## response status
response_stat  <- df_resp[,grep("[I][W][S][T][A][T]", names(df_resp), value=TRUE)] #3 variables, drinks y/n, N days/ week, drinks per day

## wave
time        <- df_resp[,grep("[I][N][W]", names(df_resp), value=TRUE)]

## age
age       <- df_resp[,grep("([A][G][E][M]|[A][G][E][Y])", names(df_resp), value=TRUE)]
age[,79:ncol(age)] <- NULL

#### outcome ####
## depressive symptoms measured with CESD
depr      <- df_resp[,grep("[C][E][S][D]", names(df_resp), value=TRUE)] #R1CESD & R13CESDM is missing
## was not measured in wave 1, add in wave 1 as NA for consistency
depr$R1CESD <- as.numeric(NA)
depr        <- depr[,c(26,1:25)]

#### mediator ####
## BMI
## PMBMI-physical measurement, BMI - self reported
## we will exclude physical measurements of BMI, only available from wave 8
bmi         <- df_resp[,grep("[B][M][I]", names(df_resp), value=TRUE)]

##alcohol consumption
alc         <- df_resp[,grep("[D][R][I][N][K]", names(df_resp), value=TRUE)] #3 variables, drinks y/n, N days/ week, drinks per day
#DRINKR wave 1-2 differs from DRINKN and DRINKD
#code DRINKR 1-2 as DRINKN 1-2, code DRINKD 1-2 as NA
#colnames(alc)[colnames(alc) == 'R1DRINKR'] <- 'R1DRINKN'
#colnames(alc)[colnames(alc) == 'R2DRINKR'] <- 'R2DRINKN'
alc$R1DRINKN <- as.numeric(NA)
alc$R2DRINKN <- as.numeric(NA)
alc$R1DRINKD <- as.numeric(NA)
alc$R2DRINKD <- as.numeric(NA)
alc         <- alc[,c(1:13,40,41,16:26,38,39,27:37)]

## smoking
smok        <- df_resp[,grep("[S][M][O][K][E]", names(df_resp), value=TRUE)] #smoke ever, smoke now -> recode to 1 variable

## physical activity
act         <- df_resp[,grep("[V][G][A][C][T]|[V][I][G][A][C][T]", names(df_resp), value=TRUE)] 
act$R1VGACTF <- NULL
act$R2VGACTN <- NULL
act$R2VGACTP <- NULL
act$R7VIGACT <- as.numeric(NA)
act$R8VIGACT <- as.numeric(NA)
act$R9VIGACT <- as.numeric(NA)
act$R10VIGACT <- as.numeric(NA)
act$R11VIGACT <- as.numeric(NA)
act$R12VIGACT <- as.numeric(NA)
act$R13VIGACT <- as.numeric(NA)

act$R1VIGACTX <- as.numeric(NA)
act$R2VIGACTX <- as.numeric(NA)
act$R3VIGACTX <- as.numeric(NA)
act$R4VIGACTX <- as.numeric(NA)
act$R5VIGACTX <- as.numeric(NA)
act$R6VIGACTX <- as.numeric(NA)
act         <- act[,c(1:6,14:20,21:26,7:13)]


APC_df      <- cbind(base, proxy, response_stat, time, age, depr, bmi[,1:13], alc, smok, act) 

#transform from wide to long
APC_long    <-  reshape(APC_df,
                direction = "long",
                varying = do.call(list,lapply(seq(19,ncol(APC_df),by=13),
                                              function(i) {
                                              (names(APC_df)[i:(i+12)])
                                               })),
                v.names = c("proxy","response_stat","INW","AGEM_B","AGEM_E","AGEM_M","AGEY_B","AGEY_E","AGEY_M",
                            "CESD", "CESDM", "BMI", "DRINK", "DRINKD", "DRINKN", "SMOKV", "SMOKN", "VIGACT", "VIGACTX"),
                idvar = "HHIDPN",
                timevar = "wave",
                times = 1:13)

rm(df);rm(df_resp);rm(idnames);rm(APC_df);(base);rm(proxy);rm(response_stat);rm(time);rm(age);rm(depr);rm(bmi);rm(alc);rm(smok);rm(act)

# Exclude non eligible respondents ----------------------------------------
## keep only alive respondents
APC_alive <- APC_long[APC_long$response_stat=="1.Resp,alive",]

## exclude wave 1 because there are no depression measurements
APC_excl.w1 <- subset(APC_alive, APC_alive$wave != "1") 

## exclude anyone older than 80 or younger than 50
APC_age_eligible <- APC_excl.w1[APC_excl.w1$AGEY_M>=50 & APC_excl.w1$AGEY_M<=80,]


# Recode some variables ---------------------------------------------------

## make year/time period variable (same as wave)
APC_age_eligible$period <- factor(APC_age_eligible$wave, labels=c("1994","1996","1998","2000","2002","2004","2006","2008","2010","2012","2014","2016"))
APC_age_eligible$period <- as.numeric(as.character(APC_age_eligible$period))

## create new birth cohort covariate in case there are discrepancies with the measured birth year
APC_age_eligible$birthcohort2 <- APC_age_eligible$period - APC_age_eligible$AGEY_M

## create elevated depressive symptoms variable (CESD > 3 or higher = elevated depressive symptoms)
APC_age_eligible$depress <- ifelse(APC_age_eligible$CESD >= "3",1, 
                                   ifelse(APC_age_eligible$CESD < "3",0, NA))

## recode alcohol consumption
APC_age_eligible$DRINKN[APC_age_eligible$DRINKN=="0.0 or doesnt drink"] <- 0
APC_age_eligible$DRINKN[APC_age_eligible$DRINKN=="99.All day"] <- 99
APC_age_eligible$DRINKN <- as.numeric(APC_age_eligible$DRINKN)

## recode drinks per week based on guidelines for women and men
APC_age_eligible$DRINK_cat <- ifelse(APC_age_eligible$DRINKN == 0, "0.or doesnt drink",
                                     ifelse(APC_age_eligible$DRINKN == 1, "1.moderate drinker",
                                            ifelse(APC_age_eligible$DRINKN == 2 & APC_age_eligible$RAGENDER == "1.Male",   "1.moderate drinker",   
                                                   ifelse(APC_age_eligible$DRINKN == 2 & APC_age_eligible$RAGENDER == "2.Female", "2.heavy drinker",
                                                          ifelse(APC_age_eligible$DRINKN ==3, "2.heavy drinker", 
                                                                 ifelse(APC_age_eligible$DRINKN == 4 & APC_age_eligible$RAGENDER == "1.Male", "2.heavy drinker",
                                                                        ifelse(APC_age_eligible$DRINKN >=4 & APC_age_eligible$RAGENDER == "2.Female","3.excessive drinker",
                                                                               ifelse(APC_age_eligible$DRINKN >=5 & APC_age_eligible$RAGENDER == "1.Male", "3.excessive drinker",NA))))))))
APC_age_eligible$DRINK_cat <- as.factor(APC_age_eligible$DRINK_cat)

#group BMI 
APC_age_eligible$BMI_cat <- ifelse(APC_age_eligible$BMI < 18, "1.underweight",
                         ifelse(APC_age_eligible$BMI >= 18 & APC_age_eligible$BMI < 25, "2.normal",
                                ifelse(APC_age_eligible$BMI >= 25 & APC_age_eligible$BMI < 30, "3.overweight",
                                       ifelse(APC_age_eligible$BMI >= 30, "4.obese",NA))))
APC_age_eligible$BMI_cat <- as.factor(APC_age_eligible$BMI_cat)


## physical activity
## there are some discrepancies here with VGACT 1-6 as nominal and VGACTX 7-13 as categories
## we will first recode physical activity in wave 7 to 13 to be binary
APC_age_eligible$VIGACTX <- as.factor(APC_age_eligible$VIGACTX)
APC_age_eligible$vigactx_bin <- ifelse(APC_age_eligible$VIGACTX == "5.Never", 0,
                      ifelse(is.na(APC_age_eligible$VIGACTX),NA,1))
APC_age_eligible$vigactx_bin <- as.factor(APC_age_eligible$vigactx_bin)

## then make new binary physical activity variables that combines all waves
APC_age_eligible$phys_act <- ifelse(APC_age_eligible$wave <= 6, APC_age_eligible$VIGACT,
                               ifelse(APC_age_eligible$wave > 6, APC_age_eligible$vigactx_bin, NA))
APC_age_eligible$phys_act<- as.factor(APC_age_eligible$phys_act)
levels(APC_age_eligible$phys_act) <- c("no vig. phys.activity", "vig. phy. activity")


## combine race and ethnicity into one
APC_age_eligible$ethn <- ifelse(APC_age_eligible$RARACEM == "1.White/Caucasian", "1.White/Caucasian",
                            ifelse(APC_age_eligible$RARACEM == "2.Black/African American", "2.Black/African American",
                                   ifelse(APC_age_eligible$RARACEM == "3.Other" & APC_age_eligible$RAHISPAN == "1.Hispanic", "3.Hispanic", 
                                          ifelse(APC_age_eligible$RARACEM == "3.Other" & APC_age_eligible$RAHISPAN == "0.Not Hispanic", "4.Other", NA))))
APC_age_eligible$ethn <- as.factor(APC_age_eligible$ethn)

## keep only necessary variables
APC_df <- APC_age_eligible[,(c("HHIDPN","RAEDUC",
                                  "RAGENDER","ethn",
                                  "SMOKN","DRINK_cat",
                                  "BMI_cat", "phys_act", "depress",
                                  "AGEY_M","period", "birthcohort2"))]

## remove missing values
APC_df <- na.omit(APC_df)

## save
save(APC_df, file="APC_df.RData")
