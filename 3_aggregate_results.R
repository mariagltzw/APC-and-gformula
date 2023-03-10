library(dplyr)


load("APC_df.RData")

## create simulated df
sim <- expand.grid(AGEY_M=sort(unique(APC_df$AGEY_M)), 
                   period=sort(unique(APC_df$period)), 
                   birthcohort2=sort(unique(APC_df$birthcohort2)))
#repeat 10 times
sim <- sim[rep(seq_len(nrow(sim)), each = 10), ]


# aggregate bootstrap results for main analysis -------------------------------------
## load mediation analysis output
load("Mediationresults_BS499.MC50.Rdata")

#specify bootstrap iterations
bs_it =499

# smoking
SMOK.BS.nc2 <- cbind(sim, (SMOK.nc))
SMOK.BS.nc2 <- SMOK.BS.nc2 %>% 
  group_by(AGEY_M, period, birthcohort2) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

SMOK.BS.cf2 <- cbind(sim, (SMOK.cf))
SMOK.BS.cf2 <- SMOK.BS.cf2 %>% 
  group_by(AGEY_M, period, birthcohort2) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

# alcohol consumption
DRINK.BS.nc2 <- cbind(sim, (DRINK.nc))
DRINK.BS.nc2 <- DRINK.BS.nc2 %>% 
  group_by(AGEY_M, period, birthcohort2) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

DRINK.BS.cf2 <- cbind(sim, (DRINK.cf))
DRINK.BS.cf2 <- DRINK.BS.cf2 %>% 
  group_by(AGEY_M, period, birthcohort2) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

# BMI
BMI.BS.nc2 <- cbind(sim, (BMI.nc))
BMI.BS.nc2 <- BMI.BS.nc2 %>% 
  group_by(AGEY_M, period, birthcohort2) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

BMI.BS.cf2 <- cbind(sim, (BMI.cf))
BMI.BS.cf2 <- BMI.BS.cf2 %>% 
  group_by(AGEY_M, period, birthcohort2) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

# physical activity
PHYSACT.BS.nc2 <- cbind(sim, (PHYSACT.nc))
PHYSACT.BS.nc2 <- PHYSACT.BS.nc2 %>% 
  group_by(AGEY_M, period, birthcohort2) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

PHYSACT.BS.cf2 <- cbind(sim, (PHYSACT.cf))
PHYSACT.BS.cf2 <- PHYSACT.BS.cf2 %>% 
  group_by(AGEY_M, period, birthcohort2) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

save(list=c("SMOK.BS.nc2", "SMOK.BS.cf2",
                       "BMI.BS.nc2", "BMI.BS.cf2",
                       "DRINK.BS.nc2", "DRINK.BS.cf2",
                       "PHYSACT.BS.nc2", "PHYSACT.BS.cf2"),
                file = paste0("Mediationresults_aggregated_",Sys.Date(),".Rdata"))



#BY SEX----
load("Mediationresults_bysex_bs_it499.MC.Rdata")

SMOK.BS.nc2.F <- cbind(sim, (SMOK.F.nc))
SMOK.BS.nc2.F$sex <- "female"
SMOK.BS.nc2.M <- cbind(sim, (SMOK.M.nc))
SMOK.BS.nc2.M$sex <- "male"

SMOK.BS.nc2 <- rbind(SMOK.BS.nc2.F, SMOK.BS.nc2.M)
SMOK.BS.nc2 <- SMOK.BS.nc2 %>%
  group_by(AGEY_M, period, birthcohort2, sex) %>%
  summarise_at(.vars = names(.) [4:(bs_it+3)],
               .funs = c(mean="mean"))

SMOK.BS.cf2.F <- cbind(sim, (SMOK.F.cf))
SMOK.BS.cf2.M <- cbind(sim, (SMOK.M.cf))
SMOK.BS.cf2.F$sex <- "female"
SMOK.BS.cf2.M$sex <- "male"

SMOK.BS.cf2 <- rbind(SMOK.BS.cf2.F, SMOK.BS.cf2.M)
SMOK.BS.cf2 <- SMOK.BS.cf2 %>%
  group_by(AGEY_M, period, birthcohort2, sex) %>%
  summarise_at(.vars = names(.) [4:(bs_it+3)],
               .funs = c(mean="mean"))

DRINK.BS.nc2.F <- cbind(sim, (DRINK.F.nc))
DRINK.BS.nc2.F$sex <- "female"
DRINK.BS.nc2.M <- cbind(sim, (DRINK.M.nc))
DRINK.BS.nc2.M$sex <- "male"

DRINK.BS.nc2 <- rbind(DRINK.BS.nc2.F, DRINK.BS.nc2.M)
DRINK.BS.nc2 <- DRINK.BS.nc2 %>%
  group_by(AGEY_M, period, birthcohort2, sex) %>%
  summarise_at(.vars = names(.) [4:(bs_it+3)],
               .funs = c(mean="mean"))

DRINK.BS.cf2.F <- cbind(sim, (DRINK.F.cf))
DRINK.BS.cf2.M <- cbind(sim, (DRINK.M.cf))
DRINK.BS.cf2.F$sex <- "female"
DRINK.BS.cf2.M$sex <- "male"

DRINK.BS.cf2 <- rbind(DRINK.BS.cf2.F, DRINK.BS.cf2.M)
DRINK.BS.cf2 <- DRINK.BS.cf2 %>%
  group_by(AGEY_M, period, birthcohort2, sex) %>%
  summarise_at(.vars = names(.) [4:(bs_it+3)],
               .funs = c(mean="mean"))

BMI.BS.nc2.F <- cbind(sim, (BMI.F.nc))
BMI.BS.nc2.F$sex <- "female"
BMI.BS.nc2.M <- cbind(sim, (BMI.M.nc))
BMI.BS.nc2.M$sex <- "male"

BMI.BS.nc2 <- rbind(BMI.BS.nc2.F, BMI.BS.nc2.M)
BMI.BS.nc2 <- BMI.BS.nc2 %>%
  group_by(AGEY_M, period, birthcohort2, sex) %>%
  summarise_at(.vars = names(.) [4:(bs_it+3)],
               .funs = c(mean="mean"))

BMI.BS.cf2.F <- cbind(sim, (BMI.F.cf))
BMI.BS.cf2.M <- cbind(sim, (BMI.M.cf))
BMI.BS.cf2.F$sex <- "female"
BMI.BS.cf2.M$sex <- "male"

BMI.BS.cf2 <- rbind(BMI.BS.cf2.F, BMI.BS.cf2.M)
BMI.BS.cf2 <- BMI.BS.cf2 %>%
  group_by(AGEY_M, period, birthcohort2, sex) %>%
  summarise_at(.vars = names(.) [4:(bs_it+3)],
               .funs = c(mean="mean"))

PHYSACT.BS.nc2.F <- cbind(sim, (PHYSACT.F.nc))
PHYSACT.BS.nc2.F$sex <- "female"
PHYSACT.BS.nc2.M <- cbind(sim, (PHYSACT.M.nc))
PHYSACT.BS.nc2.M$sex <- "male"

PHYSACT.BS.nc2 <- rbind(PHYSACT.BS.nc2.F, PHYSACT.BS.nc2.M)
PHYSACT.BS.nc2 <- PHYSACT.BS.nc2 %>%
  group_by(AGEY_M, period, birthcohort2, sex) %>%
  summarise_at(.vars = names(.) [4:(bs_it+3)],
               .funs = c(mean="mean"))

PHYSACT.BS.cf2.F <- cbind(sim, (PHYSACT.F.cf))
PHYSACT.BS.cf2.M <- cbind(sim, (PHYSACT.M.cf))
PHYSACT.BS.cf2.F$sex <- "female"
PHYSACT.BS.cf2.M$sex <- "male"

PHYSACT.BS.cf2 <- rbind(PHYSACT.BS.cf2.F, PHYSACT.BS.cf2.M)
PHYSACT.BS.cf2 <- PHYSACT.BS.cf2 %>%
  group_by(AGEY_M, period, birthcohort2, sex) %>%
  summarise_at(.vars = names(.) [4:(bs_it+3)],
               .funs = c(mean="mean"))

save(list=c("SMOK.BS.nc2", "SMOK.BS.cf2",
            "BMI.BS.nc2", "BMI.BS.cf2",
            "DRINK.BS.nc2", "DRINK.BS.cf2",
            "PHYSACT.BS.nc2", "PHYSACT.BS.cf2"),
     file = paste0("Mediationresults_bysex_aggregated.",Sys.Date(),".Rdata"))

#BY ETHNICITY----
# because the category "other" is quite small and heterogeneous we exclude it from the analysis

load("Mediationresults_by.ethn_BS.499.MC502023-03-09.Rdata")

SMOK.BS.nc2.White <- cbind(sim, (SMOK.White.nc))
SMOK.BS.nc2.White$ethn <- "White/Caucasian"
SMOK.BS.nc2.Black <- cbind(sim, (SMOK.Black.nc))
SMOK.BS.nc2.Black$ethn <- "Black/African American"
SMOK.BS.nc2.Hisp <- cbind(sim, (SMOK.Hisp.nc))
SMOK.BS.nc2.Hisp$ethn <- "Hispanic"
SMOK.BS.nc2.other <- cbind(sim, (SMOK.other.nc))
SMOK.BS.nc2.other$ethn <- "Other"

SMOK.BS.nc2 <- rbind(SMOK.BS.nc2.Black,
                     SMOK.BS.nc2.Hisp,SMOK.BS.nc2.White)#, SMOK.BS.nc2.other)
SMOK.BS.nc2 <- SMOK.BS.nc2 %>%
  group_by(AGEY_M, period, birthcohort2, ethn) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

SMOK.BS.cf2.White <- cbind(sim, (SMOK.White.cf))
SMOK.BS.cf2.White$ethn <- "White/Caucasian"
SMOK.BS.cf2.Black <- cbind(sim, (SMOK.Black.cf))
SMOK.BS.cf2.Black$ethn <- "Black/African American"
SMOK.BS.cf2.Hisp <- cbind(sim, (SMOK.Hisp.cf))
SMOK.BS.cf2.Hisp$ethn <- "Hispanic"
SMOK.BS.cf2.other <- cbind(sim, (SMOK.other.cf))
SMOK.BS.cf2.other$ethn <- "Other"

SMOK.BS.cf2 <- rbind(SMOK.BS.cf2.Black,
                     SMOK.BS.cf2.Hisp, SMOK.BS.cf2.White)#,SMOK.BS.cf2.other)
SMOK.BS.cf2 <- SMOK.BS.cf2 %>%
  group_by(AGEY_M, period, birthcohort2, ethn) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

DRINK.BS.nc2.White <- cbind(sim, (DRINK.White.nc))
DRINK.BS.nc2.White$ethn <- "White/Caucasian"
DRINK.BS.nc2.Black <- cbind(sim, (DRINK.Black.nc))
DRINK.BS.nc2.Black$ethn <- "Black/African American"
DRINK.BS.nc2.Hisp <- cbind(sim, (DRINK.Hisp.nc))
DRINK.BS.nc2.Hisp$ethn <- "Hispanic"
DRINK.BS.nc2.other <- cbind(sim, (DRINK.other.nc))
DRINK.BS.nc2.other$ethn <- "Other"

DRINK.BS.nc2 <- rbind(DRINK.BS.nc2.Black,
                      DRINK.BS.nc2.Hisp, DRINK.BS.nc2.White)#,DRINK.BS.nc2.other)
DRINK.BS.nc2 <- DRINK.BS.nc2 %>%
  group_by(AGEY_M, period, birthcohort2, ethn) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

DRINK.BS.cf2.White <- cbind(sim, (DRINK.White.cf))
DRINK.BS.cf2.White$ethn <- "White/Caucasian"
DRINK.BS.cf2.Black <- cbind(sim, (DRINK.Black.cf))
DRINK.BS.cf2.Black$ethn <- "Black/African American"
DRINK.BS.cf2.Hisp <- cbind(sim, (DRINK.Hisp.cf))
DRINK.BS.cf2.Hisp$ethn <- "Hispanic"
DRINK.BS.cf2.other <- cbind(sim, (DRINK.other.cf))
DRINK.BS.cf2.other$ethn <- "Other"

DRINK.BS.cf2 <- rbind(DRINK.BS.cf2.Black,
                      DRINK.BS.cf2.Hisp, DRINK.BS.cf2.White)#,DRINK.BS.cf2.other)
DRINK.BS.cf2 <- DRINK.BS.cf2 %>%
  group_by(AGEY_M, period, birthcohort2, ethn) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

BMI.BS.nc2.White <- cbind(sim, (BMI.White.nc))
BMI.BS.nc2.White$ethn <- "White/Caucasian"
BMI.BS.nc2.Black <- cbind(sim, (BMI.Black.nc))
BMI.BS.nc2.Black$ethn <- "Black/African American"
BMI.BS.nc2.Hisp <- cbind(sim, (BMI.Hisp.nc))
BMI.BS.nc2.Hisp$ethn <- "Hispanic"
BMI.BS.nc2.other <- cbind(sim, (BMI.other.nc))
BMI.BS.nc2.other$ethn <- "Other"

BMI.BS.nc2 <- rbind(BMI.BS.nc2.Black,
                    BMI.BS.nc2.Hisp, BMI.BS.nc2.White)#,BMI.BS.nc2.other)
BMI.BS.nc2 <- BMI.BS.nc2 %>%
  group_by(AGEY_M, period, birthcohort2, ethn) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

BMI.BS.cf2.White <- cbind(sim, (BMI.White.cf))
BMI.BS.cf2.White$ethn <- "White/Caucasian"
BMI.BS.cf2.Black <- cbind(sim, (BMI.Black.cf))
BMI.BS.cf2.Black$ethn <- "Black/African American"
BMI.BS.cf2.Hisp <- cbind(sim, (BMI.Hisp.cf))
BMI.BS.cf2.Hisp$ethn <- "Hispanic"
BMI.BS.cf2.other <- cbind(sim, (BMI.other.cf))
BMI.BS.cf2.other$ethn <- "Other"

BMI.BS.cf2 <- rbind(BMI.BS.cf2.Black,
                    BMI.BS.cf2.Hisp, BMI.BS.cf2.White)#,BMI.BS.cf2.other)
BMI.BS.cf2 <- BMI.BS.cf2 %>%
  group_by(AGEY_M, period, birthcohort2, ethn) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

PHYSACT.BS.nc2.White <- cbind(sim, (PHYSACT.White.nc))
PHYSACT.BS.nc2.White$ethn <- "White/Caucasian"
PHYSACT.BS.nc2.Black <- cbind(sim, (PHYSACT.Black.nc))
PHYSACT.BS.nc2.Black$ethn <- "Black/African American"
PHYSACT.BS.nc2.Hisp <- cbind(sim, (PHYSACT.Hisp.nc))
PHYSACT.BS.nc2.Hisp$ethn <- "Hispanic"
PHYSACT.BS.nc2.other <- cbind(sim, (PHYSACT.other.nc))
PHYSACT.BS.nc2.other$ethn <- "Other"

PHYSACT.BS.nc2 <- rbind(PHYSACT.BS.nc2.Black,
                        PHYSACT.BS.nc2.Hisp, PHYSACT.BS.nc2.White)#,PHYSACT.BS.nc2.other)
PHYSACT.BS.nc2 <- PHYSACT.BS.nc2 %>%
  group_by(AGEY_M, period, birthcohort2, ethn) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

PHYSACT.BS.cf2.White <- cbind(sim, (PHYSACT.White.cf))
PHYSACT.BS.cf2.White$ethn <- "White/Caucasian"
PHYSACT.BS.cf2.Black <- cbind(sim, (PHYSACT.Black.cf))
PHYSACT.BS.cf2.Black$ethn <- "Black/African American"
PHYSACT.BS.cf2.Hisp <- cbind(sim, (PHYSACT.Hisp.cf))
PHYSACT.BS.cf2.Hisp$ethn <- "Hispanic"
PHYSACT.BS.cf2.other <- cbind(sim, (PHYSACT.other.cf))
PHYSACT.BS.cf2.other$ethn <- "Other"

PHYSACT.BS.cf2 <- rbind(PHYSACT.BS.cf2.Black,
                        PHYSACT.BS.cf2.Hisp, PHYSACT.BS.cf2.White)#,PHYSACT.BS.cf2.other)
PHYSACT.BS.cf2 <- PHYSACT.BS.cf2 %>%
  group_by(AGEY_M, period, birthcohort2, ethn) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)],
               .funs = c(mean="mean"))

save(list=c("SMOK.BS.nc2", "SMOK.BS.cf2",
            "BMI.BS.nc2", "BMI.BS.cf2",
            "DRINK.BS.nc2", "DRINK.BS.cf2",
            "PHYSACT.BS.nc2", "PHYSACT.BS.cf2"),
     file = paste0("Mediationresults_byethn_aggregated_",Sys.Date(),".Rdata"))
