#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# 0. clean environment and load packages -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
#setwd("C:/Users/ytid13aw/Documents/Forschung/Projekt 3 Multilevel Imputation/multilevel_imputation_simluation/multilevel_imputation_simluation/LiDO")


#rm(list=ls())



required.packages <- c("lme4", "missRanger", "lmerTest", 
                       "foreach", "doParallel", "dplyr", "fungible", "merTools",
                       "broom.mixed", "sandwich",  "gamm4", "missMethods", "mixgb")


missing.packages <- required.packages[!required.packages %in% installed.packages()]

# Install missing packages
#if (length(missing.packages) > 0) {
 #   install.packages(missing.packages)
#}

# Load all the packages
lapply(required.packages, library, character.only = TRUE)
library("mice")
# load the auxiliary functions
#source("01_helper.R")

#date <- "20231016"


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# 1. Simulation Design (settings) and Function -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  


# Simulation function

single_case <- function(j, Seed){
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# i. set the design -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

    imputationen <- 5
    
    Design <- expand.grid(
        J = c(25, 50), #Cluster. Both are "few clusters". SE's will always be difficult to be estimated in the number of clusters is small: Cameron, R. C., & Miller, D. L. (2015). A Practitioner’s Guide to Cluster-Robust Inference. Journal of Human Resources, 50, 317–372. https://doi.org/10.3368/jhr.50.2.317
        # With 20 and 50, there is one case with moderate-few (25) and moderate-large (50). (KG: was 5, 20 changed to 25, 50) 
        N = c(1000), # Number of observations. I would increase it to 1000. With 1000 there are 40 and 20 individuals per cluster. (KG: was 500)
        misrate = c(0.1, 0.2 ,0.3, 0.4, 0.5), # I think 0.05 is "boring". I would set the minimum to 0.1 (was 0.05, 0.5 changed to 0.1, 0.5)
        model = c("random intercept", "random intercept and random slope"),
        mechanism = c("MAR","MNAR", "MCAR")
    )  
    
    # number of clusters
    J <- Design$J[j]
    # sample size
    N <- Design$N[j]
    # missigness rate
    misrate <- Design$misrate[j]
    # imputation and analysis model
    model <- as.character(Design$model[j])
    # missingness mechanism
    mechanism <- as.character(Design$mechanism[j])

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# ii. data generation and complete data model -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    #generate data with dataGen from helper.R
    
    Daten <- dataGen(N, J, model, Seed)
    DataGen <- Daten$DataGen
    mySeed <- Daten$mySeed

    #########################################old################################################################################ 
    # # there will be two variables at Level 2
    # # one ordinal (values 1-J) and one normally distributed with mu = 100, sd = 15
    # Level2 <- data.frame(l1 = rep(sample(J), each = N/J),
    #                      l2 = rep(rnorm(J, 100, 15), each= N/J))
    # 
    # # generating clustered data with 4 predictor variables at level 1 with indicator validities
    # # from R package "fungible" (simulates clustered data with user-defined properties)
    # multilevelData <- monte(seed = 20500, nvar = 4, nclus = J, clus.size = rep(N/J, J), eta2 = runif(4)) #JS: Why 500/?  N
    # 
    # # JS: Question: Is there any correlation between the Level2 variables and the multilevelData? Or should one generate to additional variables in multilevelData and average them on the second level? For imputation to work, variables should be correlated with each other.    
    # #KG: as of now Level 2 Variables are not correlated with multilevelData. Since it is important to have correlations I can generate two additional variables and aggregate tgem
    # 
    # # make the data set with only positive values (##K: imputation problems with negative values) #JS: What problems because of negative values?
    # #KG: mice doesn't work, problems with Cholesky decomposition because the matrix seems to be not positive semidefinit
    # 
    # monte(seed = 20500, nvar = 6, nclus = J, clus.size = rep(N/J, J), eta2 = runif(4))
    # 
    # 
    # 
    # 
    # Data <- cbind(abs(multilevelData$data), Level2)
    # 
    # colnames(Data) <- c("ID", "x1", "x2", "x3", "x4", "L21", "L22")
    # 
    # # fitting model with pre-defined properties for an outcome variable y
    # DataGen <- as.data.frame(cbind(Data, y = model.fit(J, model, N, Data))) #JS: How do you set the the beta-values here? Or how wo we know the true coefficients? --> MyModel shows you the true coefficients? Just in case: the true coefficients should not vary with each simulation, I would say.           
    # #KG: the values are  set in helper.R and fixed coefficients are set constant for all the simulations and do not vary (beta <- c(0.3, 0.5, 1, 1.5, 1.3, 1.4, 0.07)) 
    # #if that's what you mean
    # 
    # colnames(DataGen) <- c("ID", "x1", "x2", "x3", "x4", "L21", "L22", "y")#, "L23")
    ##############################################################################################################################
    
    
    # analysis model for the complete data
    if(model == "random intercept"){
        
        MyModel <- lmer(y ~ x1 + x2 + x3 + x4 + x5 + x6 + L21 + L22 +L23 + L24 + L25 + L26 + (1|ID), data = DataGen)
        
    }else if(model == "random intercept and random slope"){
        # random intercept and random slope
        MyModel  <- lmer(y ~  x1 + x2 + x3 + x4 + x5 + x6 + L21 + L22 +L23 + L24 + L25 + L26 + (1 + x1 + x2|ID), data = DataGen) 
        #JS: Why the random slope for x1 and x2? How about we expand it to plm(, model = "random") and plm(, model = "within")?
        #KG: what does this mean? that all covariates are random slopes? I only chose 2 for my thesis to not make it complicated
        #should I add another model? is plm better than lmer?
    }else{
        print("Error: undefined model fit")
    }
    
    #OLS
    #OLSData <- DataGen
    #OLSData$ID <- as.factor(OLSData$ID)
    #MyModelOLS <- lm(y ~ . , data = OLSData)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# iii. amputation -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    #amputation functions from 'helper.R'
    
    indices <- seq(from = 1, to = N, by = N/J)
    L2.val <- as.data.frame(cbind(
        L21 = DataGen$L21[indices], 
        L22 = DataGen$L22[indices], 
        L23 = DataGen$L23[indices], 
        L24 = DataGen$L24[indices], 
        L25 = DataGen$L25[indices], 
        L26 = DataGen$L26[indices]
    ))

    if(mechanism == "MAR"){ 
        #KG: according ro markus pauly ampute does not work well with MAR
        #mar_missingness has problems with mice
      col. <- sample(1:6, size = 1)
      MisDat <- mar_normal(DataGen[, -c(1,8:13)], misrate, column = col.)
      L2.mis <- mcar_missingness(L2.val, misrate, nrow(L2.val))
        
    }else if(mechanism == "MCAR"){
        MisDat <- mcar_missingness(DataGen[, -c(1,8:13)], misrate, N)
        L2.mis <- mcar_missingness(L2.val, misrate, nrow(L2.val))
        
    }else {
        #KG: I'm not sure how to go about MNAR, so I'll just leave ampute for now
        MisDat <- ampute(DataGen[ , -c(1,8:13)], misrate, mech = mechanism)$amp
        L2.mis <- ampute(L2.val, misrate, mech = mechanism)$amp
    }

    #Level 2 variables have problemc with mice with 50%, ampute has too low missing values so I use generateNA and reduce misrate
    #Error in .imputation.level2(y = y, ry = ry, x = x, type = type, wy = wy,  : 
    #No class variable
    
    # L2.mis <- L2.val
    # if(misrate == 0.5){
    #         L2.mis$L21 <- generateNA(L2.val$L21, misrate-0.1)
    #         L2.mis$L22 <- generateNA(L2.val$L22, misrate-0.1)
    # }else{
    #     L2.mis$L21 <- generateNA(L2.val$L21, misrate)
    #     L2.mis$L22 <- generateNA(L2.val$L22, misrate)
    # }

     #L2.mis <- ampute(L2.val, misrate )$amp
    
     # L2.mis <- L2.val
     # 
     # L2.mis[as.logical(rbinom(nrow(L2.val), 1, misrate-0.1)), 1] <- NA
     # L2.mis[as.logical(rbinom(nrow(L2.val), 1, misrate-0.1)), 2] <- NA
     # 
    L2.var <- as.data.frame(cbind(rep(L2.mis$L21, each = N/J), rep(L2.mis$L22, each = N/J),  rep(L2.mis$L23, each = N/J),
                                   rep(L2.mis$L24, each = N/J) ,  rep(L2.mis$L25, each = N/J),  rep(L2.mis$L26, each = N/J)))
    DataAmp <- cbind("ID" = DataGen$ID, MisDat, L2.var)
    colnames(DataAmp) <- c("ID", "x1", "x2", "x3", "x4", "x5", "x6", "y", "L21", "L22", "L23", "L24", "L25", "L26")

    is_l21_miss <- anyNA(DataAmp$L21)
    is_l22_miss <- anyNA(DataAmp$L22)
    
    ###########################################old#######################################
    # # Amputation with ampute from mice
    # 
    # # ID variable without missing data
    # Amp <- ampute(DataGen[ , -c(1,6:7)], misrate, mech = mechanism ) # Should we use ampute? Markus Pauly said that this method is not the best with MAR, I think? Why do variables 6 and 7 do not get missings? Because they are level2? Do they get missings later in the code? Yes.
    # #KG:then should the amputation be implemented manually for all methods or just MAR?
    # 
    # # ampute some clusters at level 2 for all the observations in this cluster
    # # indices of cluster values
    # indices <- seq(from = 1, to = N, by = N/J)
    # 
    # # for 5 clusters and 50% missingness reduce missing rate otherwise it is possible that all the clusters are missing
    # ##K: even less missing?
    # if(misrate == 0.5 & J == 5){ # If this still a problem if we set J to at least 25? #KG: it should not be a problem anymore
    #     
    #     lev2mis <- 0.25
    #     
    # }else{
    #     
    #     lev2mis <- misrate
    #     
    # }
    # 
    # Lev2 <- ampute(as.data.frame(cbind(DataGen$L21[indices], DataGen$L22[indices])), lev2mis, mech = mechanism)$amp
    # Level2 <- as.data.frame(cbind(rep(Lev2$V1, each = N/J), rep(Lev2$V2, each = N/J)))
    # 
    # # amputed dataset
    # DataAmp <- cbind("ID" = DataGen$ID, Amp$amp, Level2)
    # colnames(DataAmp) <- c("ID", "x1", "x2", "x3", "x4", "y", "L21", "L22")
    #######################################################################################
    
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# iv. imputations -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# iv.a. MICE -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    #JS: with MICE method? Just name it here once so that one does not have to open the helper
    # prediction matrix and imputation methods predefined in helper.R for both models
    #KG: 2l.norm for Level1 variables, 2lonly.pmm for Level2 variables
    PredMar <-  predMet(DataAmp, model)
    # prediction matrix
    predmat <- PredMar$predmat
    # prediction methods
    impMethod <- PredMar$impMethod
    # imputation settings 
    Start.mice <- Sys.time()
    imp <- mice(DataAmp, method = impMethod, predictorMatrix = predmat, maxit = 10, ##K:10
                m = imputationen) 
    time.mice <- Sys.time() - Start.mice
    
    dataMice <- as.list(1:imputationen)
    
    for(i in 1:imputationen){ 
        
        dataMice[[i]] <- complete(imp, action=i)
        
    }
 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# iv.a.2 MICE without considering the multilevel structure
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    PredMar2 <-  predMet(DataAmp, model)
    # prediction matrix
    predmat <- PredMar2$predmat
    # imputation settings 
    Start.mice <- Sys.time()
    imp <- mice(DataAmp, method = "pmm", predictorMatrix = predmat, maxit = 10, ##K:10
                m = imputationen) 
    time.mice2 <- Sys.time() - Start.mice
    
    dataMice2 <- as.list(1:imputationen)
    
    for(i in 1:imputationen){ 
        
        dataMice2[[i]] <- complete(imp, action=i)
        
    }
    
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# iv.b. missRanger standard -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    Start.mR <- Sys.time()
    rangerImp.list <- list()
    # run a multiple imputation process for missRanger
    for (ii in 1:imputationen) { 
        rangerImp.list[[ii]] <- missRanger::missRanger(DataAmp,
                                                       formula = . ~ . - ID , 
                                                       pmm.k = 0, # Why N/J? And not just 3? I would fix it, because now the variation is not just N and J, but also that the methods is performed differently. If you have a good reason for it, I am happy to hear it. But it could also be an expansion. Once with a fixed pmm.k and once with your N/J. 
                                                       splitrule = "variance",
                                                       num.trees = 300,
        )
    }
    
    time.mr <- Sys.time() - Start.mR
    
    # aggregating the level two variables: all observations within the same cluster should have
    # same values for level 2 variables
    
    #for (ii in 1:imputationen){ 
     #   
      #  rangerImp.list[[ii]]$L21 <- round(rep((rangerImp.list[[ii]] %>%
       #                                            dplyr::group_by(ID) %>%
        #                                           summarise_at(vars(L21), list(name = mean)))$name, each = N/J), 2)
        
    #}
    #
    #for (ii in 1:imputationen){ 
     #   
      #  rangerImp.list[[ii]]$L22 <- round(rep((rangerImp.list[[ii]] %>%
       #                                            dplyr:: group_by(ID) %>%
        #                                           summarise_at(vars(L22), list(name = mean)))$name, each = N/J), 2)
        
    #}
    #
    #for (ii in 1:imputationen){ 
     #   
      #  rangerImp.list[[ii]]$L23 <- round(rep((rangerImp.list[[ii]] %>%
       #                                            dplyr:: group_by(ID) %>%
        #                                           summarise_at(vars(L23), list(name = mean)))$name, each = N/J), 2)
        
    #}

    
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# iv.b. missRanger adjusted with dummies -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
#+

    #Psych needed for dummy.code
    # dummy variables indicating the corresponding cluster 
    Dummies <- as.data.frame(psych::dummy.code(as.factor(DataGen$ID)))
    
    # name and convert to factors for imputation model
    for(i in 1:J){
        
        colnames(Dummies)[i] <- paste0("School", i) 
        
    }
    
    Dummies[, 1:J] <- lapply(Dummies[, 1:J], as.factor)
    
    DataAmpDummies <- cbind(DataAmp, Dummies)
    

    Start.dmR <- Sys.time()
    
    # run a multiple imputation process for missRanger taking into consideration
    # dummy factor variables for the imputation model
    levelRangerImp <- list()
    for (ii in 1:imputationen) { 
        
        levelRangerImp[[ii]]<-  missRanger(DataAmpDummies, 
                                           formula = .  ~ . -ID  ,
                                           pmm.k = 0, 
                                           splitrule = "variance",
                                           num.trees = 300,
        )
    }
    
    time.dmR <- Sys.time() - Start.dmR
    
    # aggregating the level two variables: all observations within the same cluster should have
    # same values for level 2 variables
    
   # for (ii in 1:imputationen){ 
        
    #    levelRangerImp[[ii]]$L21 <- round(rep((levelRangerImp[[ii]] %>%
      #                                             dplyr:: group_by(ID) %>%
     #                                              summarise_at(vars(L21), list(name = mean)))$name, each = N/J), 2)
    #    
    #}
    #
    #for (ii in 1:imputationen){ 
     #   
     #   levelRangerImp[[ii]]$L22 <- round(rep((levelRangerImp[[ii]] %>%
     #                                              dplyr:: group_by(ID) %>%
     #                                              summarise_at(vars(L22), list(name = mean)))$name, each = N/J), 2)
     #   
    #}
    #
    #for (ii in 1:imputationen){ 
     #   
      #  levelRangerImp[[ii]]$L23 <- round(rep((levelRangerImp[[ii]] %>%
      #                                             dplyr:: group_by(ID) %>%
       #                                            summarise_at(vars(L23), list(name = mean)))$name, each = N/J), 2)
        
#    }

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# iv.c. missRanger with pmm = 5 -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    
    Start.mR5 <- Sys.time()
    rangerImp5.list <- list()
    # run a multiple imputation process for missRanger
    for (ii in 1:imputationen) { 
        rangerImp5.list[[ii]] <- missRanger::missRanger(DataAmp,
                                                       formula = . ~ . - ID , 
                                                       pmm.k = 5, # Why N/J? And not just 3? I would fix it, because now the variation is not just N and J, but also that the methods is performed differently. If you have a good reason for it, I am happy to hear it. But it could also be an expansion. Once with a fixed pmm.k and once with your N/J. 
                                                       splitrule = "variance",
                                                       num.trees = 300,
        )
    }
    
    time.mr5 <- Sys.time() - Start.mR5
    
    # aggregating the level two variables: all observations within the same cluster should have
    # same values for level 2 variables
    
   # for (ii in 1:imputationen){ 
    #    
     #   rangerImp5.list[[ii]]$L21 <- round(rep((rangerImp5.list[[ii]] %>%
     #                                              dplyr::group_by(ID) %>%
     #                                              summarise_at(vars(L21), list(name = mean)))$name, each = N/J), 2)
        
    #}
    
    #for (ii in 1:imputationen){ 
     #   
      #  rangerImp5.list[[ii]]$L22 <- round(rep((rangerImp5.list[[ii]] %>%
      #                                             dplyr:: group_by(ID) %>%
       #                                            summarise_at(vars(L22), list(name = mean)))$name, each = N/J), 2)
        #
    #}
    
    #for (ii in 1:imputationen){ 
     #   
      #  rangerImp5.list[[ii]]$L23 <- round(rep((rangerImp5.list[[ii]] %>%
       #                                             dplyr:: group_by(ID) %>%
        #                                            summarise_at(vars(L22), list(name = mean)))$name, each = N/J), 2)
        
    #}
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# iv.d. missRanger adjusted with dummies and 5 PMM -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
#
    Start.dmR5 <- Sys.time()
    
    # run a multiple imputation process for missRanger taking into consideration
    # dummy factor variables for the imputation model
    levelRangerImp5 <- list()
    for (ii in 1:imputationen) { 
        
        levelRangerImp5[[ii]]<-  missRanger(DataAmpDummies, 
                                           formula = .  ~ . -ID  ,
                                           pmm.k = 5, 
                                           splitrule = "variance",
                                           num.trees = 300,
        )
    }
    
    time.dmR5 <- Sys.time() - Start.dmR
    
    # aggregating the level two variables: all observations within the same cluster should have
    # same values for level 2 variables
    
   # for (ii in 1:imputationen){ 
    #    
     #   levelRangerImp5[[ii]]$L21 <- round(rep((levelRangerImp5[[ii]] %>%
      #                                             dplyr:: group_by(ID) %>%
       #                                            summarise_at(vars(L21), list(name = mean)))$name, each = N/J), 2)
        
    #}
    #
    #for (ii in 1:imputationen){ 
     #   
      #  levelRangerImp5[[ii]]$L22 <- round(rep((levelRangerImp5[[ii]] %>%
       #                                            dplyr:: group_by(ID) %>%
        #                                          summarise_at(vars(L22), list(name = mean)))$name, each = N/J), 2)
        
    #}
    
    #for (ii in 1:imputationen){ 
     #   
      #  levelRangerImp5[[ii]]$L23 <- round(rep((levelRangerImp5[[ii]] %>%
        #                                            dplyr:: group_by(ID) %>%
       #                                             summarise_at(vars(L23), list(name = mean)))$name, each = N/J), 2)
        
    #}
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
 #iv.e. mixgb -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
  
   startgb <- Sys.time() 
   MXGB <-  mixgb(DataAmp, m = 5, maxit = 5)
    timegb <-  Sys.time()  - startgb  
 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# iv.f. mixgb with dummies -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    
    startgb.d <- Sys.time() 
    MXGB.dummies <-  mixgb(DataAmpDummies, m = 5, maxit = 5)
    timegb.d <-  Sys.time()  - startgb.d  
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# v. model analysis -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    # fit the models with the pre defined function from helper.R
    
    imp.fit.Mice <- imp.fit(dataMice, model)
    imp.fit.Mice.no.cluster <- imp.fit(dataMice2, model)
    imp.fit.missRanger <- imp.fit(rangerImp.list, model)
    imp.fit.LevelmissRanger <- imp.fit(levelRangerImp, model)
    imp.fit.missRanger5 <- imp.fit(rangerImp5.list, model)
    imp.fit.LevelmissRanger5 <- imp.fit(levelRangerImp5, model)
    imp.fit.boost <- imp.fit(MXGB, model)
    imp.fit.boost.dummies <- imp.fit(MXGB.dummies, model)
    # OLS.mice <- ols.fit(dataMice)
    # OLS.ranger <- ols.fit(rangerImp.list)
    # OLS.Dummyranger <- ols.fit(levelRangerImp)
    # simple pooling of the fixed effects of the linear mixed effects models
    pooled.Mice <- pool(imp.fit.Mice)
    pooled.Mice.no.cluster <- pool(imp.fit.Mice.no.cluster)
    pooled.ranger <- pool(imp.fit.missRanger)
    pooled.levelRanger <- pool(imp.fit.LevelmissRanger)
    pooled.ranger5 <- pool(imp.fit.missRanger5)
    pooled.levelRanger5 <- pool(imp.fit.LevelmissRanger5)
    pooled.boost <- pool(imp.fit.boost)
    pooled.boost.dummies <- pool(imp.fit.boost.dummies)
    
    # OLS.pooled.mice <- pool(OLS.mice)
    # OLS.pooled.ranger <- pool(OLS.ranger)
    # OLS.pooled.Dummyranger <- pool(OLS.Dummyranger)
    # model summary tables
    Summary.model <- as.data.frame(summary(MyModel)$coefficients)
    Summary.Mice <- summary(pooled.Mice)
    Summary.Mice.no.cluster <- summary(pooled.Mice.no.cluster)
    Summary.ranger <- summary(pooled.ranger)
    Summary.levelRanger <-  summary(pooled.levelRanger)
    Summary.ranger5 <- summary(pooled.ranger5)
    Summary.levelRanger5 <-  summary(pooled.levelRanger5)
    Summary.boost <- summary(pooled.boost)
    Summary.boost.dummies <- summary(pooled.boost.dummies)
    
    # Summary.OLS <- as.data.frame(summary(MyModelOLS)$coefficients)
    # Summary.OLS.mice <- summary(OLS.pooled.mice)
    # Summary.OLS.ranger <- summary(OLS.pooled.ranger)
    # Summary.OLS.Dummyranger <- summary(OLS.pooled.Dummyranger)
    
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# v.a. additional: goodness of fit measures -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
 
    # goodness of fit measures with predefined functions from helper.R
    
   # LogLiks <- rbind(c(logLik(MyModel, REML = FALSE)), LogLs(imp.fit.Mice), 
    #                 LogLs(imp.fit.missRanger), LogLs(imp.fit.LevelmissRanger), LogLs(imp.fit.missRanger5),
     #                LogLs(imp.fit.LevelmissRanger5), LogLs(imp.fit.boost), LogLs(imp.fit.boost.dummies))
    
    #Devs <- -2*LogLiks
    
    #AICS <- rbind(AIC(MyModel), (-2*LogLiks[2, 1] + 2*7), (-2*LogLiks[3, 1] + 2*7), (-2*LogLiks[4, 1] + 2*7),
    #              (-2*LogLiks[5, 1] + 2*7), (-2*LogLiks[6, 1] + 2*7), (-2*LogLiks[7, 1] + 2*7),  (-2*LogLiks[8, 1] + 2*7))
    
    #BICS <- rbind(BIC(MyModel), (-2*LogLiks[2, 1] + 7*log(1000)), (-2*LogLiks[3, 1] + 7*log(1000)),
     #             (-2*LogLiks[4, 1] + 7*log(1000)), (-2*LogLiks[5, 1] + 7*log(1000)),
      #            (-2*LogLiks[6, 1] + 7*log(1000)), (-2*LogLiks[7, 1] + 7*log(1000)), (-2*LogLiks[8, 1] + 7*log(1000))) 
    
    #GoF <- cbind(LogLiks, Devs, AICS, BICS)
    #colnames(GoF) <- c("LogPooled", "DevPooled", "AICPooled", "BICPooled")
    #rownames(GoF) <- c("Data", "Mice", "standard missRanger", "adapted missRanger",
       #                "missRanger 5PMM", "adapted missRanger 5PMM", "mixgboost", "mixgboost dummies")
    

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# vi. additional: imputed level 2 variable values -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
     
    # true cluster values
    Level2Data <- as.data.frame(cbind(DataGen$L21[indices], DataGen$L22[indices]))
    
    # cluster values for imputation with mice
    Level2Mice <- as.data.frame(cbind( rowMeans(matrix(unlist(lapply(dataMice, function(df) df$L21[indices])), nrow = J, byrow = FALSE)),
                                       rowMeans(matrix(unlist(lapply(dataMice, function(df) df$L22[indices])), nrow = J, byrow = FALSE))))
    
    # cluster values for imputation with mice without considering clustered structure
    Level2Mice.no.cluster <- as.data.frame(cbind( rowMeans(matrix(unlist(lapply(dataMice2, function(df) df$L21[indices])), nrow = J, byrow = FALSE)),
                                       rowMeans(matrix(unlist(lapply(dataMice2, function(df) df$L22[indices])), nrow = J, byrow = FALSE))))

    # cluster values for imputation with standard missRanger
    Level2Ranger <- as.data.frame(cbind( rowMeans(matrix(unlist(lapply(rangerImp.list, function(df) df$L21[indices])), nrow = J, byrow = FALSE)),
                                         rowMeans(matrix(unlist(lapply(rangerImp.list, function(df) df$L22[indices])), nrow = J, byrow = FALSE))))

    # cluster values for imputation with adapted missRanger
    Level2LevelRanger <- as.data.frame(cbind( rowMeans(matrix(unlist(lapply(levelRangerImp, function(df) df$L21[indices])), nrow = J, byrow = FALSE)),
                                              rowMeans(matrix(unlist(lapply(levelRangerImp, function(df) df$L22[indices])), nrow = J, byrow = FALSE))))
    
    Level2Ranger5 <- as.data.frame(cbind( rowMeans(matrix(unlist(lapply(rangerImp5.list, function(df) df$L21[indices])), nrow = J, byrow = FALSE)),
                                         rowMeans(matrix(unlist(lapply(rangerImp5.list, function(df) df$L22[indices])), nrow = J, byrow = FALSE))))
    
    Level2LevelRanger5 <- as.data.frame(cbind( rowMeans(matrix(unlist(lapply(levelRangerImp5, function(df) df$L21[indices])), nrow = J, byrow = FALSE)),
                                              rowMeans(matrix(unlist(lapply(levelRangerImp5, function(df) df$L22[indices])), nrow = J, byrow = FALSE))))
    
    Level2boost <- as.data.frame(cbind( rowMeans(matrix(unlist(lapply(MXGB, function(df) df$L21[indices])), nrow = J, byrow = FALSE)),
                                        rowMeans(matrix(unlist(lapply(MXGB, function(df) df$L22[indices])), nrow = J, byrow = FALSE))))
    
    Level2boost.dummies <- as.data.frame(cbind( rowMeans(matrix(unlist(lapply(MXGB.dummies, function(df) df$L21[indices])), nrow = J, byrow = FALSE)),
                                                rowMeans(matrix(unlist(lapply(MXGB.dummies, function(df) df$L22[indices])), nrow = J, byrow = FALSE))))
    
    
    DesignM <- as.data.frame(cbind(Design[j, ], L21NA = is_l21_miss, L22NA = is_l22_miss))
    Times <- c(time.mice, time.mice2, time.mr, time.mr5, timegb )    
    names(Times) <- c("mice",  "missranger",  "missRanger5", "mixgboost", "mixgboost.dummies")
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# vii. list of results -----
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
   
    myRes <- list(Design = DesignM, Summary.model = Summary.model, Summary.Mice = Summary.Mice,
                  Summary.Mice.no.cluster = Summary.Mice.no.cluster,
                  Summary.ranger = Summary.ranger, Summary.levelRanger = Summary.levelRanger,
                  Summary.ranger5 = Summary.ranger5, Summary.levelRanger5 = Summary.levelRanger5,
                  Summary.boost = Summary.boost, Summary.boost.dummies = Summary.boost.dummies,
                  #Summary.OLS = Summary.OLS, Summary.OLS.mice = Summary.OLS.mice, 
                  #Summary.OLS.ranger = Summary.OLS.ranger, Summary.OLS.Dummyranger = Summary.OLS.Dummyranger, 
                  Level2Data = Level2Data, Level2Mice = Level2Mice, Level2Mice.no.cluster = Level2Mice.no.cluster,
                  Level2Ranger= Level2Ranger, 
                  Level2LevelRanger = Level2LevelRanger, Level2Ranger5 = Level2Ranger5,
                  Level2LevelRanger5 = Level2LevelRanger5, Level2boost = Level2boost, 
                  Level2boost.dummies = Level2boost.dummies,
                  #GoF = GoF,
                  runtimes = Times, mySeed=mySeed )
    
    return(myRes)
    
    
}


