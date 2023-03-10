## function below runs the mediation analysis for each mediator
## mediator = name of mediator, e.g. "BMI_cat"
## dataframe = df to use for mediation
## sim.df = simulated dataframe to use for predictions
## mc_iteration = iteration for monte carlo error reduction loop
## reference cohort = default is 1945
## strata = default is FALSE, options: "sex", "ethnicity"


#output will be c(nc results, cf results)
#output stratified by sex will be c(FEMALES (c(nc results, cf results)), MALES (c(nc results, cf results)))
#output stratified by ethn will be c(white (c(nc results, cf results)), black (c(nc results, cf results))), hisp (c(nc results, cf results))), other (c(nc results, cf results)))

mediation_analysis <- function(mediator, dataframe, sim.df, mc_iteration, reference_cohort = 1945, strata = FALSE){
  
  ## check to make sure the correct dataset is being used
  if (nrow(dataframe) != 163760) {return(paste0("CHECK IF DATASET IS CORRECT! Current dataset has ", nrow(dataframe)," rows"))}
  
  if (strata =="FALSE"){
    
    ## sample from dataframe
    bootstrap_sample  <- cluster.resample(dataframe, cluster.name = "HHIDPN", 
                                          size = length(unique(dataframe[,"HHIDPN"])))
    
    ## fit outcome model with APC specifications (Carstensen approach)
    ## check the paper for more information on the APC approach 
    outcome.model     <- glm(paste("depress ~ Ns(AGEY_M, df=6) + Ns(period, df=5,detrend=TRUE) + Ns(birthcohort2, df=10) +", mediator, "+ RAGENDER + RAEDUC + ethn",sep=" "),
                             family ="binomial",
                             data = bootstrap_sample)
    
    ## copy simulated df
    sim.df.nc <- sim.df.cf <- data.frame(1:nrow(sim.df))
    
    ## fit mediator model
    if (class(dataframe[,mediator]) == "factor" & length(unique(dataframe[,mediator])) == 2){
      mediator.model <- glm(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=5,detrend=TRUE) + Ns(birthcohort2, df=10) + RAGENDER + RAEDUC + ethn",sep=""), 
                            family ="binomial",
                            data = bootstrap_sample)
      probs <- predict(mediator.model, type = "response", newdata = sim.df)
      
      ## start of Monte Carlo loop
      for (m in 1:mc_iteration){
        ## natural course
        ## predict mediator values
        sim.df[,mediator] <- as.factor(rbinom(n = nrow(sim.df), size = 1,
                                              prob = probs))
        levels(sim.df[,mediator]) <- levels(dataframe[,mediator])
        sim.df.nc[,m]   <- predict(outcome.model, type="response", newdata = sim.df) 
        
        ## counterfactual
        ## our counterfactual is the health behavior factor distribution of the 1945 cohort
        ## so we predict just with that cohort
        sim.df[,mediator] <- sim.df[,mediator][sim.df$birthcohort2==reference_cohort]
        sim.df.cf[,m]  <- predict(outcome.model, type="response", newdata = sim.df) 
      }
    }   else if (class(dataframe[,mediator]) == "factor" & length(unique(dataframe[,mediator])) > 2){
      mediator.model <- multinom(bootstrap_sample[,mediator] ~ Ns(AGEY_M, df=6) + 
                                   Ns(period, df=5,detrend=TRUE) + 
                                   Ns(birthcohort2, df=10) + RAGENDER + RAEDUC + ethn, 
                                 data = bootstrap_sample)
      simulated_mediator <- predict(mediator.model, type="probs",newdata = sim.df)
      
      for (m in 1:mc_iteration){
        #natural course
        #predict mediator values
        sim.df[,mediator] <- rMultinom(m = 1, probs = simulated_mediator)
        sim.df.nc[,m]   <- predict(outcome.model, type="response", newdata = sim.df) 
        
        #counterfactual
        sim.df[,mediator] <- sim.df[,mediator][sim.df$birthcohort2==reference_cohort]
        sim.df.cf[,m]  <- predict(outcome.model, type="response", newdata = sim.df) 
      }
      
    } else if (class(dataframe[,mediator]) == "numeric"){
      mediator.model <- lm(bootstrap_sample[,mediator] ~ Ns(AGEY_M, df=6) + 
                             Ns(period, df=5,detrend=TRUE) + 
                             Ns(birthcohort2, df=10) + RAGENDER + RAEDUC + ethn, 
                           data = bootstrap_sample)
      simulated_mediator  <- predict(mediator.model, type = "response", newdata = sim.df)
      
      for (m in 1:mc_iteration){
        #natural course
        #predict mediator values
        sim.df[,mediator]     <- rnorm(n=nrow(sim.df), mean = mean(simulated_mediator), sd=sd(mediator.model$residuals))
        sim.df.nc[,m]   <- predict(outcome.model, type="response", newdata = sim.df) 
        
        #counterfactual
        sim.df[,mediator] <- sim.df[,mediator][sim.df$birthcohort2==reference_cohort]
        sim.df.cf[,m]  <- predict(outcome.model, type="response", newdata = sim.df) 
      }
    } else {
      return("Mediator has the wrong class")
    }
    
    ## aggregate MC results for natural course and counterfactual
    df.BS.nc <- rowMeans(sim.df.nc)
    df.BS.cf <- rowMeans(sim.df.cf)
    c(df.BS.nc,df.BS.cf)
    }
  
  ## strata for sex
  else if(strata == "sex"){
    df.F <- subset(dataframe, dataframe[,"RAGENDER"]=="2.Female")
    df.M <- subset(dataframe, dataframe[,"RAGENDER"]=="1.Male")
    bootstrap_sample.F  <- sample_n(df.F, size = nrow(dataframe), replace = TRUE)
    bootstrap_sample.M  <- sample_n(df.M, size = nrow(dataframe), replace = TRUE)
    
    outcome.model.F     <- glm(paste("depress ~ Ns(AGEY_M, df=6) + Ns(period, df=5,detrend=TRUE) + Ns(birthcohort2, df=10) +", mediator, " + RAEDUC + ethn",sep=" "),
                               family ="binomial",
                               data = bootstrap_sample.F)
    outcome.model.M     <- glm(paste("depress ~ Ns(AGEY_M, df=6) + Ns(period, df=5,detrend=TRUE) + Ns(birthcohort2, df=10) +", mediator, " + RAEDUC + ethn",sep=" "),
                               family ="binomial",
                               data = bootstrap_sample.M)
    
    sim.df.F.nc <- sim.df.F.cf <- sim.df.M.nc <- sim.df.M.cf <- data.frame(1:nrow(sim.df))
    sim.df.F <- sim.df.M <- sim.df
    
    if (class(dataframe[,mediator]) == "factor" & length(unique(dataframe[,mediator])) == 2){
      mediator.model.F <- glm(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=5,detrend=TRUE) + Ns(birthcohort2, df=10) + RAEDUC + ethn",sep=""), 
                              family ="binomial",
                              data = bootstrap_sample.F)
      mediator.model.M <- glm(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=5,detrend=TRUE) + Ns(birthcohort2, df=10) + RAEDUC + ethn",sep=""), 
                              family ="binomial",
                              data = bootstrap_sample.M)
      probs.F <- predict(mediator.model.F, type = "response", newdata = sim.df.F)
      probs.M <- predict(mediator.model.M, type = "response", newdata = sim.df.M)
      
      for (m in 1:mc_iteration){
        #FEMALES
        #natural course
        #predict mediator values
        sim.df.F[,mediator] <- as.factor(rbinom(n = nrow(sim.df.F), size = 1,
                                                prob = probs.F))
        levels(sim.df.F[,mediator]) <- levels(dataframe[,mediator])
        sim.df.F.nc[,m]   <- predict(outcome.model.F, type="response", newdata = sim.df.F) 
        
        #counterfactual
        sim.df.F[,mediator] <- sim.df.F[,mediator][sim.df.F$birthcohort2==reference_cohort]
        sim.df.F.cf[,m]  <- predict(outcome.model.F, type="response", newdata = sim.df.F) 
        
        #MALES
        #natural course
        #predict mediator values
        sim.df.M[,mediator] <- as.factor(rbinom(n = nrow(sim.df.M), size = 1,
                                                prob = probs.M))
        levels(sim.df.M[,mediator]) <- levels(dataframe[,mediator])
        sim.df.M.nc[,m]   <- predict(outcome.model.M, type="response", newdata = sim.df.M) 
        
        #counterfactual
        sim.df.M[,mediator] <- sim.df.M[,mediator][sim.df.M$birthcohort2==reference_cohort]
        sim.df.M.cf[,m]  <- predict(outcome.model.M, type="response", newdata = sim.df.M) 
        
      }
      
    }   else if (class(dataframe[,mediator]) == "factor" & length(unique(dataframe[,mediator])) > 2){
      mediator.model.F <- multinom(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=5,detrend=TRUE) + Ns(birthcohort2, df=10) + RAEDUC + ethn",sep=""), 
                                   data = bootstrap_sample.F)
      mediator.model.M <- multinom(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=5,detrend=TRUE) + Ns(birthcohort2, df=10) + RAEDUC + ethn",sep=""), 
                                   data = bootstrap_sample.M)
      
      simulated_mediator.F <- predict(mediator.model.F, type="probs",newdata = sim.df.F)
      simulated_mediator.M <- predict(mediator.model.M, type="probs",newdata = sim.df.M)
      
       for (m in 1:mc_iteration){
        #FEMALES
        #natural course
        #predict mediator values
        sim.df.F[,mediator] <- rMultinom(m = 1, probs = simulated_mediator.F)
        sim.df.F.nc[,m]   <- predict(outcome.model.F, type="response", newdata = sim.df.F) 
        #counterfactual
        sim.df.F[,mediator] <- sim.df.F[,mediator][sim.df.F$birthcohort2==reference_cohort]
        sim.df.F.cf[,m]  <- predict(outcome.model.F, type="response", newdata = sim.df.F) 
        
        #MALES
        #natural course
        #predict mediator values
        sim.df.M[,mediator] <- rMultinom(m = 1, probs = simulated_mediator.M)
        sim.df.M.nc[,m]   <- predict(outcome.model.M, type="response", newdata = sim.df.M) 
        #counterfactual
        sim.df.M[,mediator] <- sim.df.M[,mediator][sim.df.M$birthcohort2==reference_cohort]
        sim.df.M.cf[,m]  <- predict(outcome.model.M, type="response", newdata = sim.df.M) 
      }
      
    } else if (class(dataframe[,mediator]) == "numeric"){
      mediator.model.F <- lm(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=5,detrend=TRUE) + Ns(birthcohort2, df=10) + RAEDUC + ethn",sep=""), 
                             family ="binomial",
                             data = bootstrap_sample.F)
      mediator.model.M <- lm(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=5,detrend=TRUE) + Ns(birthcohort2, df=10) + RAEDUC + ethn",sep=""), 
                             family ="binomial",
                             data = bootstrap_sample.M)
      
      simulated_mediator.F    <- predict(mediator.model.F, type = "response", newdata = sim.df.F)
      simulated_mediator.M    <- predict(mediator.model.M, type = "response", newdata = sim.df.M)
      
      for (m in 1:mc_iteration){
        #FEMALES
        #natural course
        #predict mediator values
        sim.df.F[,mediator]     <- rnorm(n=nrow(sim.df.F), mean = mean(simulated_mediator.F), sd=sd(mediator.model.F$residuals))
        sim.df.F.nc[,m]   <- predict(outcome.model.F, type="response", newdata = sim.df.F) 
        
        #counterfactual
        sim.df.F[,mediator] <- sim.df.F[,mediator][sim.df.F$birthcohort2==reference_cohort]
        sim.df.F.cf[,m]  <- predict(outcome.model.F, type="response", newdata = sim.df.F) 
        
        #MALES
        #natural course
        #predict mediator values
        sim.df.M[,mediator]     <- rnorm(n=nrow(sim.df.M), mean = mean(simulated_mediator.M), sd=sd(mediator.model.M$residuals))
        sim.df.M.nc[,m]   <- predict(outcome.model.M, type="response", newdata = sim.df.M) 
        
        #counterfactual
        sim.df.M[,mediator] <- sim.df.M[,mediator][sim.df.M$birthcohort2==reference_cohort]
        sim.df.M.cf[,m]  <- predict(outcome.model.M, type="response", newdata = sim.df.M) 
      }
    } else {
      return("Mediator has the wrong class")
    }
    df.BS.F.nc <- rowMeans(sim.df.F.nc)
    df.BS.F.cf <- rowMeans(sim.df.F.cf)
    results_FEMALES <- c(df.BS.F.nc, df.BS.F.cf)
    
    df.BS.M.nc <- rowMeans(sim.df.M.nc)
    df.BS.M.cf <- rowMeans(sim.df.M.cf)
    results_MALES <- c(df.BS.M.nc, df.BS.M.cf)

    c(results_FEMALES,results_MALES)
  }
  
  ## strate for race/ethnicity
  else if(strata == "ethnicity"){
    df.1 <- subset(dataframe, dataframe[,"ethn"]=="1.White/Caucasian")
    df.2 <- subset(dataframe, dataframe[,"ethn"]=="2.Black/African American")
    df.3 <- subset(dataframe, dataframe[,"ethn"]=="3.Hispanic")
    df.4 <- subset(dataframe, dataframe[,"ethn"]=="4.Other")
    
    bootstrap_sample.1  <- sample_n(df.1, size = nrow(dataframe), replace = TRUE)
    bootstrap_sample.2  <- sample_n(df.2, size = nrow(dataframe), replace = TRUE)
    bootstrap_sample.3  <- sample_n(df.3, size = nrow(dataframe), replace = TRUE)
    bootstrap_sample.4  <- sample_n(df.4, size = nrow(dataframe), replace = TRUE)
    
    outcome.model.1     <- glm(paste0("depress ~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) +", mediator, " + RAGENDER + RAEDUC"),
                               family ="binomial",
                               data = bootstrap_sample.1)
    outcome.model.2     <- glm(paste0("depress ~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) +", mediator, " + RAGENDER + RAEDUC"),
                               family ="binomial",
                               data = bootstrap_sample.2)
    outcome.model.3     <- glm(paste0("depress ~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) +", mediator, " + RAGENDER + RAEDUC"),
                               family ="binomial",
                               data = bootstrap_sample.3)
    outcome.model.4     <- glm(paste0("depress ~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) +", mediator, " + RAGENDER + RAEDUC"),
                               family ="binomial",
                               data = bootstrap_sample.4)
    
    sim.df.1.nc <- sim.df.1.cf <- sim.df.2.nc <- sim.df.2.cf <- sim.df.3.nc <- sim.df.3.cf <- sim.df.4.nc <- sim.df.4.cf <- data.frame(1:nrow(sim.df))
    sim.df.1 <- sim.df.2 <- sim.df.3 <- sim.df.4 <- sim.df
    
    if (class(dataframe[,mediator]) == "factor" & length(unique(dataframe[,mediator])) == 2){
      mediator.model.1 <- glm(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) + RAGENDER + RAEDUC",sep=""), 
                              family ="binomial",
                              data = bootstrap_sample.1)
      mediator.model.2 <- glm(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) + RAGENDER + RAEDUC",sep=""), 
                              family ="binomial",
                              data = bootstrap_sample.2)
      mediator.model.3 <-  glm(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) + RAGENDER + RAEDUC",sep=""), 
                               family ="binomial",
                               data = bootstrap_sample.3)
      mediator.model.4 <- glm(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) + RAGENDER + RAEDUC",sep=""), 
                              family ="binomial",
                              data = bootstrap_sample.4)
      
      probs.1 <- predict(mediator.model.1, type = "response", newdata = sim.df.1)
      probs.2 <- predict(mediator.model.2, type = "response", newdata = sim.df.2)
      probs.3 <- predict(mediator.model.3, type = "response", newdata = sim.df.3)
      probs.4 <- predict(mediator.model.4, type = "response", newdata = sim.df.4)
      
      for (m in 1:mc_iteration){
        #1
        #natural course
        #predict mediator values
        sim.df.1[,mediator] <- as.factor(rbinom(n = nrow(sim.df.1), size = 1,
                                                prob = probs.1))
        levels(sim.df.1[,mediator]) <- levels(dataframe[,mediator])
        sim.df.1.nc[,m]   <- predict(outcome.model.1, type="response", newdata = sim.df.1) 
        
        #counterfactual
        sim.df.1[,mediator] <- sim.df.1[,mediator][sim.df.1$birthcohort2==reference_cohort]
        sim.df.1.cf[,m]  <- predict(outcome.model.1, type="response", newdata = sim.df.1) 
        
        #2
        #natural course
        #predict mediator values
        sim.df.2[,mediator] <- as.factor(rbinom(n = nrow(sim.df.2), size = 1,
                                                prob = probs.2))
        levels(sim.df.2[,mediator]) <- levels(dataframe[,mediator])
        sim.df.2.nc[,m]   <- predict(outcome.model.2, type="response", newdata = sim.df.2) 
        
        #counterfactual
        sim.df.2[,mediator] <- sim.df.2[,mediator][sim.df.2$birthcohort2==reference_cohort]
        sim.df.2.cf[,m]  <- predict(outcome.model.2, type="response", newdata = sim.df.2) 
        
        #3
        #natural course
        #predict mediator values
        sim.df.3[,mediator] <- as.factor(rbinom(n = nrow(sim.df.3), size = 1,
                                                prob = probs.3))
        levels(sim.df.3[,mediator]) <- levels(dataframe[,mediator])
        sim.df.3.nc[,m]   <- predict(outcome.model.3, type="response", newdata = sim.df.3) 
        
        #counterfactual
        sim.df.3[,mediator] <- sim.df.3[,mediator][sim.df.3$birthcohort2==reference_cohort]
        sim.df.3.cf[,m]  <- predict(outcome.model.3, type="response", newdata = sim.df.3) 
        
        #4
        #natural course
        #predict mediator values
        sim.df.4[,mediator] <- as.factor(rbinom(n = nrow(sim.df.4), size = 1,
                                                prob = probs.4))
        levels(sim.df.4[,mediator]) <- levels(dataframe[,mediator])
        sim.df.4.nc[,m]   <- predict(outcome.model.4, type="response", newdata = sim.df.4) 
        
        #counterfactual
        sim.df.4[,mediator] <- sim.df.4[,mediator][sim.df.4$birthcohort2==reference_cohort]
        sim.df.4.cf[,m]  <- predict(outcome.model.4, type="response", newdata = sim.df.4) 
      }
      
    }   else if (class(dataframe[,mediator]) == "factor" & length(unique(dataframe[,mediator])) > 2){
      mediator.model.1 <- multinom(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) + RAGENDER + RAEDUC",sep=""), 
                                   data = bootstrap_sample.1)
      mediator.model.2 <- multinom(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) + RAGENDER + RAEDUC",sep=""), 
                                   data = bootstrap_sample.2)
      mediator.model.3 <- multinom(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) + RAGENDER + RAEDUC",sep=""), 
                                   data = bootstrap_sample.3)
      mediator.model.4 <- multinom(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) + RAGENDER + RAEDUC",sep=""), 
                                   data = bootstrap_sample.4)
      
      simulated_mediator.1 <- predict(mediator.model.1, type="probs",newdata = sim.df.1)
      simulated_mediator.2 <- predict(mediator.model.2, type="probs",newdata = sim.df.2)
      simulated_mediator.3 <- predict(mediator.model.3, type="probs",newdata = sim.df.3)
      simulated_mediator.4 <- predict(mediator.model.4, type="probs",newdata = sim.df.4)
      
      for (m in 1:mc_iteration){
        #1
        #natural course
        #predict mediator values
        sim.df.1[,mediator] <- rMultinom(m = 1, probs = simulated_mediator.1)
        sim.df.1.nc[,m]   <- predict(outcome.model.1, type="response", newdata = sim.df.1) 
        #counterfactual
        sim.df.1[,mediator] <- sim.df.1[,mediator][sim.df.1$birthcohort2==reference_cohort]
        sim.df.1.cf[,m]  <- predict(outcome.model.1, type="response", newdata = sim.df.1) 
        
        #2
        #natural course
        #predict mediator values
        sim.df.2[,mediator] <- rMultinom(m = 1, probs = simulated_mediator.2)
        sim.df.2.nc[,m]   <- predict(outcome.model.2, type="response", newdata = sim.df.2) 
        #counterfactual
        sim.df.2[,mediator] <- sim.df.2[,mediator][sim.df.2$birthcohort2==reference_cohort]
        sim.df.2.cf[,m]  <- predict(outcome.model.2, type="response", newdata = sim.df.2) 
        
        #3
        #natural course
        #predict mediator values
        sim.df.3[,mediator] <- rMultinom(m = 1, probs = simulated_mediator.3)
        sim.df.3.nc[,m]   <- predict(outcome.model.3, type="response", newdata = sim.df.3) 
        #counterfactual
        sim.df.3[,mediator] <- sim.df.3[,mediator][sim.df.3$birthcohort2==reference_cohort]
        sim.df.3.cf[,m]  <- predict(outcome.model.3, type="response", newdata = sim.df.3) 
        
        #4
        #natural course
        #predict mediator values
        sim.df.4[,mediator] <- rMultinom(m = 1, probs = simulated_mediator.4)
        sim.df.4.nc[,m]   <- predict(outcome.model.4, type="response", newdata = sim.df.4) 
        #counterfactual
        sim.df.4[,mediator] <- sim.df.4[,mediator][sim.df.4$birthcohort2==reference_cohort]
        sim.df.4.cf[,m]  <- predict(outcome.model.4, type="response", newdata = sim.df.4) 
      }
      
    } else if (class(dataframe[,mediator]) == "numeric"){
      mediator.model.1 <- lm(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) + RAGENDER + RAEDUC",sep=""), 
                             family ="binomial",
                             data = bootstrap_sample.1)
      mediator.model.2 <- lm(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) + RAGENDER + RAEDUC",sep=""), 
                             family ="binomial",
                             data = bootstrap_sample.2)
      mediator.model.3 <- lm(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) + RAGENDER + RAEDUC",sep=""), 
                             family ="binomial",
                             data = bootstrap_sample.3)
      mediator.model.4 <- lm(paste(mediator, "~ Ns(AGEY_M, df=6) + Ns(period, df=4,detrend=TRUE) + Ns(birthcohort2, df=10) + RAGENDER + RAEDUC",sep=""), 
                             family ="binomial",
                             data = bootstrap_sample.4)
      
      simulated_mediator.1    <- predict(mediator.model.1, type = "response", newdata = sim.df.1)
      simulated_mediator.2    <- predict(mediator.model.2, type = "response", newdata = sim.df.2)
      simulated_mediator.3    <- predict(mediator.model.3, type = "response", newdata = sim.df.3)
      simulated_mediator.4    <- predict(mediator.model.4, type = "response", newdata = sim.df.4)
      
      for (m in 1:mc_iteration){
        #1
        #natural course
        #predict mediator values
        sim.df.1[,mediator]     <- rnorm(n=nrow(sim.df.1), mean = mean(simulated_mediator.1), sd=sd(mediator.model.1$residuals))
        sim.df.1.nc[,m]   <- predict(outcome.model.1, type="response", newdata = sim.df.1) 
        
        #counterfactual
        sim.df.1[,mediator] <- sim.df.1[,mediator][sim.df.1$birthcohort2==reference_cohort]
        sim.df.1.cf[,m]  <- predict(outcome.model.1, type="response", newdata = sim.df.1) 
        
        #2
        #natural course
        #predict mediator values
        sim.df.2[,mediator]     <- rnorm(n=nrow(sim.df.2), mean = mean(simulated_mediator.2), sd=sd(mediator.model.2$residuals))
        sim.df.2.nc[,m]   <- predict(outcome.model.2, type="response", newdata = sim.df.2) 
        
        #counterfactual
        sim.df.2[,mediator] <- sim.df.2[,mediator][sim.df.2$birthcohort2==reference_cohort]
        sim.df.2.cf[,m]  <- predict(outcome.model.2, type="response", newdata = sim.df.2) 
        #4
        #natural course
        #predict mediator values
        sim.df.3[,mediator]     <- rnorm(n=nrow(sim.df.3), mean = mean(simulated_mediator.3), sd=sd(mediator.model.3$residuals))
        sim.df.3.nc[,m]   <- predict(outcome.model.3, type="response", newdata = sim.df.3) 
        
        #counterfactual
        sim.df.3[,mediator] <- sim.df.3[,mediator][sim.df.3$birthcohort2==reference_cohort]
        sim.df.3.cf[,m]  <- predict(outcome.model.3, type="response", newdata = sim.df.3) 
        
        #4
        #natural course
        #predict mediator values
        sim.df.4[,mediator]     <- rnorm(n=nrow(sim.df.4), mean = mean(simulated_mediator.4), sd=sd(mediator.model.4$residuals))
        sim.df.4.nc[,m]   <- predict(outcome.model.4, type="response", newdata = sim.df.4) 
        
        #counterfactual
        sim.df.4[,mediator] <- sim.df.4[,mediator][sim.df.4$birthcohort2==reference_cohort]
        sim.df.4.cf[,m]  <- predict(outcome.model.4, type="response", newdata = sim.df.4) 
      }
    } else {
      return("Mediator has the wrong class")
    }
    df.BS.1.nc <- rowMeans(sim.df.1.nc)
    df.BS.1.cf <- rowMeans(sim.df.1.cf)
    results_white <- c(df.BS.1.nc, df.BS.1.cf)
    
    df.BS.2.nc <- rowMeans(sim.df.2.nc)
    df.BS.2.cf <- rowMeans(sim.df.2.cf)
    results_black <- c(df.BS.2.nc, df.BS.2.cf)
    
    df.BS.3.nc <- rowMeans(sim.df.3.nc)
    df.BS.3.cf <- rowMeans(sim.df.3.cf)
    results_hisp <- c(df.BS.3.nc, df.BS.3.cf)
    
    df.BS.4.nc <- rowMeans(sim.df.4.nc)
    df.BS.4.cf <- rowMeans(sim.df.4.cf)
    results_other <- c(df.BS.4.nc, df.BS.4.cf)
    
    c(results_white,results_black,results_hisp, results_other)
  } else {
    return("Specify strata as FALSE, sex or ethnicity")
  }
}

## decomposition table function
## needs results_nc and results_cf income
## strata=FALSE -> main analysis, can also accommodate strata="sex" and strata="ethnicity"
## bs_it are the bootstrap iteration. If not specified 499 iterations are assumed.
decomp_table <- function(results_nc, results_cf, bs_it=499, strata=FALSE){
  x <- ncol(results_nc)
  if(strata=="FALSE"){
    results_nc$lb <- round(apply(results_nc[,4:x], 1, quantile, probs = c(0.025))*100,2)
    results_nc$ub <- round(apply(results_nc[,4:x], 1, quantile, probs = c(0.975))*100,2)
    results_nc$prob.depress <- round(rowMeans(results_nc[,4:x])*100,2)
    results_nc$scenario <- "natural course"
    #remove unaggregated columns
    results_nc <- results_nc[,c(1:3, (bs_it+4):(bs_it+7))]
    
    results_cf$lb <- round(apply(results_cf[,4:x], 1, quantile, probs = c(0.025))*100,2)
    results_cf$ub <- round(apply(results_cf[,4:x], 1, quantile, probs = c(0.975))*100,2)
    results_cf$prob.depress <- round(rowMeans(results_cf[,4:x])*100,2)
    results_cf$scenario <- "counterfactual"
    #remove unaggregated columns
    results_cf <- results_cf[,c(1:3, (bs_it+4):(bs_it+7))]
    
    plot.results  <- merge(results_nc, results_cf, by = c("birthcohort2", "period", "AGEY_M"))
    plot.results.c <- subset(plot.results, plot.results$AGEY_M == 50 &
                               plot.results$period == 1996)
    #remove age and period columns
    plot.results.c <- plot.results.c[,c(1,4:11)]
    plot.results.c <- as.data.frame(plot.results.c)
    plot.results.c$naturalcourse <- paste0(plot.results.c$prob.depress.x,"% (",plot.results.c$lb.x,"-",plot.results.c$ub.x,")")
    plot.results.c$counterfactual <- paste0(plot.results.c$prob.depress.y,"% (",plot.results.c$lb.y,"-",plot.results.c$ub.y,")")
    table <- plot.results.c[,c(1,10,11)] 
  } 
  else if (strata=="sex"){
    results_nc$lb <- round(apply(results_nc[,5:x], 1, quantile, probs = c(0.025))*100,2)
    results_nc$ub <- round(apply(results_nc[,5:x], 1, quantile, probs = c(0.975))*100,2)
    results_nc$prob.depress <- round(rowMeans(results_nc[,5:x])*100,2)
    results_nc$scenario <- "natural course"
    #remove unaggregated columns
    results_nc <- results_nc[,c(1:4, (bs_it+5):(bs_it+8))]
    
    results_cf$lb <- round(apply(results_cf[,5:x], 1, quantile, probs = c(0.025))*100,2)
    results_cf$ub <- round(apply(results_cf[,5:x], 1, quantile, probs = c(0.975))*100,2)
    results_cf$prob.depress <- round(rowMeans(results_cf[,5:x])*100,2)
    results_cf$scenario <- "counterfactual"
    #remove unaggregated columns
    results_cf <- results_cf[,c(1:4, (bs_it+5):(bs_it+8))]
    
    plot.results  <- merge(results_nc, results_cf, by = c("birthcohort2", "period", "AGEY_M","sex"))
    plot.results.c <- subset(plot.results, plot.results$AGEY_M == 50 &
                               plot.results$period == 1996)
    plot.results.c <- plot.results.c[,c(1,4,5:12)]
    
    plot.results.c <- as.data.frame(plot.results.c[order(plot.results.c$sex),])
    plot.results.c$naturalcourse <- paste0(plot.results.c$prob.depress.x,"% (",plot.results.c$lb.x,"-",plot.results.c$ub.x,")")
    plot.results.c$counterfactual <- paste0(plot.results.c$prob.depress.y,"% (",plot.results.c$lb.y,"-",plot.results.c$ub.y,")")
    table <- plot.results.c[,c(1,2,11,12)] 
  } 
  else if (strata=="ethn"){
    results_nc$lb <- round(apply(results_nc[,5:x], 1, quantile, probs = c(0.025))*100,2)
    results_nc$ub <- round(apply(results_nc[,5:x], 1, quantile, probs = c(0.975))*100,2)
    results_nc$prob.depress <- round(rowMeans(results_nc[,5:x])*100,2)
    results_nc$scenario <- "natural course"
    results_nc <- results_nc[,c(1:4, (bs_it+5):(bs_it+8))]
    
    results_cf$lb <- round(apply(results_cf[,5:x], 1, quantile, probs = c(0.025))*100,2)
    results_cf$ub <- round(apply(results_cf[,5:x], 1, quantile, probs = c(0.975))*100,2)
    results_cf$prob.depress <- round(rowMeans(results_cf[,5:x])*100,2)
    results_cf$scenario <- "counterfactual"
    results_cf <- results_cf[,c(1:4, (bs_it+5):(bs_it+8))]
    
    plot.results  <- merge(results_nc, results_cf, by = c("birthcohort2", "period", "AGEY_M","ethn"))
    plot.results.c <- subset(plot.results, plot.results$AGEY_M == 50 &
                               plot.results$period == 1996)
    plot.results.c <- plot.results.c[,c(1,4,5:12)]
    plot.results.c <- as.data.frame(plot.results.c[order(plot.results.c$ethn),])
    plot.results.c$naturalcourse <- paste0(plot.results.c$prob.depress.x,"% (",plot.results.c$lb.x,"-",plot.results.c$ub.x,")")
    plot.results.c$counterfactual <- paste0(plot.results.c$prob.depress.y,"% (",plot.results.c$lb.y,"-",plot.results.c$ub.y,")")
    table <- plot.results.c[,c(1,2,11,12)] 
  }
  return(table)
}

## make relative difference and contribution function
## needs results_nc and results_cf income
## strata=FALSE -> main analysis, can also accommodate strata="sex" and strata="ethnicity"
## bs_it are the bootstrap iteration. If not specified 499 iterations are assumed.
rel.diff_contr_table <- create_mediationplot_abs.rel.diff <- function(results_nc, results_cf, bs_it=499, strata =FALSE){
  results_nc <- subset(results_nc, results_nc$AGEY_M == 50 &
                         results_nc$period == 1996)
  results_cf <- subset(results_cf, results_cf$AGEY_M == 50 &
                         results_cf$period == 1996)
  x <- ncol(results_nc)
  if (strata=="FALSE"){
    #relative difference
    BS.diff <- (results_cf[,4:x] / results_nc[,4:x])-1
    
    BS.diff2 <- BS.diff
    BS.diff2$lb <- round(apply(BS.diff, 1, quantile, probs = c(0.025))*100,2)
    BS.diff2$ub <- round(apply(BS.diff, 1, quantile, probs = c(0.975))*100,2)
    BS.diff2$diff <- round(rowMeans(BS.diff)*100,2)
    
    #contribution
    #ref.cf prob.
    #ref cohort prob same for cf and nc
    p.depr.refcohort <- rowMeans(results_nc[results_nc$birthcohort2==1945,4:x])
    perc.med <- 1-((results_cf[,4:x] - p.depr.refcohort)/((results_nc[,4:x]) - p.depr.refcohort))
    
    perc.med2 <- perc.med
    perc.med2$contr <- round(apply(perc.med, 1, median)*100,2)
    perc.med2$lb_c <- round(apply(perc.med, 1, quantile, probs = c(0.025))*100,2)
    perc.med2$ub_c <- round(apply(perc.med, 1, quantile, probs = c(0.975))*100,2)
    
    make.table <- cbind(results_nc[,1:3], BS.diff2[,(ncol(BS.diff2)-2):ncol(BS.diff2)], perc.med2[,(ncol(perc.med2)-2):ncol(perc.med2)])
    make.table <- as.data.frame(make.table)
    make.table$reldiff <- paste0(make.table$diff,"% (",make.table$lb,"-",make.table$ub,")")
    make.table$contr <- paste0(make.table$contr,"% (",make.table$lb_c,"-",make.table$ub_c,")")
    return(make.table[,c(1:3,10,7)])
  }
  else if (strata=="sex"){
    #relative difference
    BS.diff <- (results_cf[,5:x] / results_nc[,5:x])-1
    
    BS.diff2 <- BS.diff
    BS.diff2$lb <- round(apply(BS.diff, 1, quantile, probs = c(0.025))*100,2)
    BS.diff2$ub <- round(apply(BS.diff, 1, quantile, probs = c(0.975))*100,2)
    BS.diff2$diff <- round(rowMeans(BS.diff)*100,2)
    
    #contribution
    p.depr.refcohort <- ifelse(results_nc$sex == "female",
                               rowMeans(results_nc[results_nc$birthcohort2==1945 & results_nc$sex=="female",5:x]),
                               ifelse(results_nc$sex == "male",
                                      rowMeans(results_nc[results_nc$birthcohort2==1945 & results_nc$sex=="male",5:x]), NA))
    perc.med <- 1-((results_cf[,5:x] - p.depr.refcohort)/((results_nc[,5:x]) - p.depr.refcohort))
    
    perc.med2 <- perc.med
    perc.med2$contr <- round(apply(perc.med, 1, median)*100,2)
    perc.med2$lb_c <- round(apply(perc.med, 1, quantile, probs = c(0.025))*100,2)
    perc.med2$ub_c <- round(apply(perc.med, 1, quantile, probs = c(0.975))*100,2)
    
    make.table <- cbind(results_nc[,1:4], BS.diff2[,(ncol(BS.diff2)-2):ncol(BS.diff2)], perc.med2[,(ncol(perc.med2)-2):ncol(perc.med2)])
    make.table <- as.data.frame(make.table[order(make.table$sex),])
    make.table$reldiff <- paste0(make.table$diff,"% (",make.table$lb,"-",make.table$ub,")")
    make.table$contr <- paste0(make.table$contr,"% (",make.table$lb_c,"-",make.table$ub_c,")")
    return(make.table[,c(1:4,11,8)])
  }
  else if (strata=="ethn"){
    #relative difference
    BS.diff <- (results_cf[,5:x] / results_nc[,5:x])-1
    
    BS.diff2 <- BS.diff
    BS.diff2$lb <- round(apply(BS.diff, 1, quantile, probs = c(0.025))*100,2)
    BS.diff2$ub <- round(apply(BS.diff, 1, quantile, probs = c(0.975))*100,2)
    BS.diff2$diff <- round(rowMeans(BS.diff)*100,2)
    
    #contribution
    p.depr.refcohort <- ifelse(results_nc$ethn == "White/Caucasian",
                               rowMeans(results_nc[results_nc$birthcohort2==1945 & results_nc$ethn=="White/Caucasian",5:x]),
                               ifelse(results_nc$ethn == "Black/African American",
                                      rowMeans(results_nc[results_nc$birthcohort2==1945 & results_nc$ethn=="Black/African American",5:x]),
                                      ifelse(results_nc$ethn == "Hispanic",
                                             rowMeans(results_nc[results_nc$birthcohort2==1945 & results_nc$ethn=="Hispanic",5:x]), NA)))
    
    perc.med <- 1-((results_cf[,5:x] - p.depr.refcohort)/((results_nc[,5:x]) - p.depr.refcohort))
    
    perc.med2 <- perc.med
    perc.med2$contr <- round(apply(perc.med, 1, median)*100,2)
    perc.med2$lb_c <- round(apply(perc.med, 1, quantile, probs = c(0.025))*100,2)
    perc.med2$ub_c <- round(apply(perc.med, 1, quantile, probs = c(0.975))*100,2)
    
    make.table <- cbind(results_nc[,1:4], BS.diff2[,(ncol(BS.diff2)-2):ncol(BS.diff2)], perc.med2[,(ncol(perc.med2)-2):ncol(perc.med2)])
    make.table <- as.data.frame(make.table[order(make.table$ethn),])
    make.table$reldiff <- paste0(make.table$diff,"% (",make.table$lb,"-",make.table$ub,")")
    make.table$contr <- paste0(make.table$contr,"% (",make.table$lb_c,"-",make.table$ub_c,")")
    return(make.table[,c(1:4,11,8)])
  } 
}

