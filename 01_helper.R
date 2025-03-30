#this document just creates some auxiliary functions for the simulation



# to fit a model and create an outcome variable
 model.fit <- function(J, fit.model, N, data){
   # coefficients
   beta <- c(0.3, 0.5, 1, 1.5, 0, 0, 0, 0.5, 1, 1.5, 0, 0, 0)
   # errors
   eps <- rnorm(N, 0, 1)  
    
   # data matrix
   X <- as.matrix(cbind(1, data[, 2:13]))

   if(fit.model == "random intercept"){
     # random intercept
     b0 <- rnorm(J, 0, 1)
     # generate y
     y <-  b0[data$ID] + X%*%beta + eps #JS: "b0[data$ID]" is the individual specific effect which is constant over time but different for each indidivual?             
     #b0[data$ID] is the cluster specific effect, constant for all individuals in the same cluster and different for each clusters
     
     data$y <- as.numeric(y)

   } else if(fit.model == "random intercept and random slope"){ #JS: I would suggest to also include an if here and add an else with "print(Something is not defined.)"
     # random intercept  
     b0 <- rnorm(J, 0, 1)
     # random slope x1 and x2
     b1 <- rnorm(J, 0, 1)
     b2 <- rnorm(J, 0, 1)
     # generate y
     y <- b0[data$ID] + b1[data$ID]*data$x1 + b2[data$ID]*data$x2 + X%*%beta + eps
     data$y <- as.numeric(y)

   }else{
       print("Error: undefined model fit")
   }
   
   return(y = data$y)
 }

 
 #Data Generation
 dataGen <- function(N, J, model, Seed){
      #KG: for each different N x J x model design different data but same for replications or without seed
      #op <- options(digit.secs = 6)
      #Seed <- sample(1000, 1)
      #gsub("[^[:digit:]]", "",Seed)
     Seed
      multilevelData <- fungible::monte(seed = Seed, nvar = 12, nclus = J, clus.size = rep(N/J, J), eta2 = runif(12, max = 0.85), sortMeans = FALSE)
      Data <- data.frame(multilevelData$data)
      mySeed <- multilevelData$seed
      colnames(Data) <- c("ID", "x1", "x2", "x3", "x4", "x5", "x6", "L21", "L22", "L23", "L24", "L25", "L26")
      Data$ID <- rep(sample(J), each = N/J)
      #aggregate 2 variables at level 2 (level 2 variables are correlated to level 1 variables)
      Data$L21 <- rep(round(aggregate(Data$L21, by = list(Data$ID), FUN = mean)$x*10, 2), each = N/J)
      Data$L22 <- rep(round(aggregate(Data$L22, by = list(Data$ID), FUN = mean)$x*10, 2), each = N/J)
      Data$L23 <- rep(round(aggregate(Data$L23, by = list(Data$ID), FUN = mean)$x*10, 2), each = N/J)
      Data$L24 <- rep(round(aggregate(Data$L24, by = list(Data$ID), FUN = mean)$x*10, 2), each = N/J)
      Data$L25 <- rep(round(aggregate(Data$L25, by = list(Data$ID), FUN = mean)$x*10, 2), each = N/J)
      Data$L26 <- rep(round(aggregate(Data$L26, by = list(Data$ID), FUN = mean)$x*10, 2), each = N/J)
      DataGen <- as.data.frame(cbind(Data, y = model.fit(J, model, N, Data)))
      colnames(DataGen) <- c("ID", "x1", "x2", "x3", "x4", "x5", "x6", "L21", "L22", "L23", "L24", "L25", "L26", "y")
      return(list(DataGen = DataGen, mySeed = mySeed))
 }


 #missingness
 #nach Xijuan Zhang, Tutorial: How to Generate Missing Data For Simulation Studies Seite 13-14
 #Error in .imputation.level2(y = y, ry = ry, x = x, type = type, wy = wy,  : 
 #No class variable
 
 # mcar_missingness <- function(data, missingness_prob ) {
 #     n <- nrow(data)
 #     p <- ncol(data)
 #     
 #     missing_data <- data
 #     
 #     for (row_idx in 1:n) {
 #         for (col_idx in 1:p) {
 #             if (runif(1) < missingness_prob) {
 #                 missing_data[row_idx, col_idx] <- NA
 #             }
 #         }
 #     }
 #     
 #     return(missing_data)
 # }
 # 
 
mcar_missingness <- function(data, missing.percentage, N){
    
    for(i in 1:ncol(data)){
     ind <- as.logical(rbinom(N, 1, missing.percentage))
     data[ind,i] <- NA
    }
    return(data)
} 
 
              
 
 
 
 # 7. MAR Mechanism ----
 # Algorithm described in Thurow et al. (2021) (Goodness (of fit) of imputation accuracy: The GoodImpact analysis)
 #' MAR mechanism
 #'
 #' Returns the given dataframe with MAR entries according to the missingess rate.
 #' 
 #' @param data The dataframe to be modified by adding NAs
 #' @param miss_rate Value between 0 and 1 corresponding to the expected missingness rate.
 #' 
 #' @return The initial dataframe after adding missing values.
 #' 
 #' c:andrés with some adjustments for numerical data
 #' 
mar_normal <- function(mat, miss_rate, column = 2){
    # Checking the missingness rate
    if(miss_rate > 1 | miss_rate < 0){
        stop("The missingness rate is not between 0 and 1")
    }else{ 

        # Making a vector of the desired missing rates (equal to all columns)
        miss_rates <- rep(miss_rate,ncol(mat))
        # Adding missing values to the first column in the dataframe
        
        if(miss_rate == 0.5){
          mat[,column] <- mcar_missingness(data.frame(mat[,column]), 0.3, nrow(mat)) 
        }else{
           mat[,column] <- mcar_missingness(data.frame(mat[,column]), miss_rate, nrow(mat)) 
        }
        
        #first mat[, column] to grouped factor
        allvar <- as.factor(cut(mat[, column], 25))
        
        # Indices in the chosen column for which an NA has been placed
        mis_1 <- !complete.cases(mat[,column])
        # initializing m: the matrix with the missing values (1 is NA, 0 is not NA)
        m <- matrix(0,nrow(mat),ncol(mat)) 
        # Adding the missingness to the chosen column
        m[,column] <- mis_1+0
        
        # vec1' is a dataframe with the information for the MAR mechanism.
        # It starts with 2 columns, the first one 'a1' has the unique non-NA values
        # and the second one 'H1' the absolute frequencies of each non-Na value.
        
        vec1 <- data.frame(table(allvar[!mis_1]))
        names(vec1) <- c("a1","H1")
        #vec1$a1 <- allvar[!mis_1]
        #vec1$a1 <- as.factor(vec1$a1)
        vec1 <- vec1[order(vec1$H1, decreasing = FALSE),]
        row.names(vec1) <- NULL
        
        
        # 'pa' is a vector of random uniform values between 0 and 1 used to then add 
        # randomness to the other columns of the dataframe
        vec1["pa"] <- sort(runif(nrow(vec1),0,1)) 
        
        # Denominator for normalizing the probabilities in 'pa'
        den <- sum(vec1$H1*vec1$pa)
        # 'pi' has the normalized probabilities for each non-NA value
        vec1["pi"] <- vec1$H1*vec1$pa/den
        # Number of missing entries per column
        n_mis <- nrow(m)*miss_rate # Number of missing entries per column
        
        # Chosen entries
        vec1["cho"] <- round(vec1$pi*n_mis)
        # Adjusting the chosen entries to meet the expected number after rounding
        vec1$cho[1] <- n_mis - sum(vec1$cho[2:nrow(vec1)])
        
        # Adjustment for Missingness Rate = 0.5
        if(miss_rate >=0.5){
            vec1$cho <- vec1$H1
        }
        
        #
        
        # Loop that adds NA's to the other columns according to 'vec1'
        for (c in 1:length(miss_rates)){
            if(c != column){
                mis_count <- 0
                for (i in 1:nrow(vec1)){
                    if(vec1$cho[i] > 0 & mis_count < n_mis){
                        #browser()
                        if(vec1$a1[i] %in% allvar){
                             I <- which(allvar == vec1$a1[i])
                        } else{
                            I <- sample(nrow(mat), vec1$cho[i])
                        }
                       
                        mis_count <- mis_count + vec1$cho[i]
                        if(vec1$cho[1] <= length(I)){
                         I_final <- sample(x = I, size = vec1$cho[i]) 
                        }else{
                            I_final <- I 
                        }
                         
                        m[I_final,c] <- 1
                    }
                }
            }
        }
        
        mat[m==1] <- NA
        mat
    }
}
 
 
 
 
 
 
 
 
 
 # mice pred matrix and methods
 # pred 1 for level 2 variables
 # 3 for random intercept all others
 # 4 for random slopes
 
 predMet <- function(data, model){
     ini <- mice(data, maxit = 0)
     impMethod <- character(ncol(data))
     names(impMethod) <- colnames(data)
     impMethod[2:8] <- "2l.norm"
     impMethod[9:14] <- "2lonly.pmm"
     
     if (model == "random intercept") {
         predmat <- ini$predictorMatrix
         # predmat["ID", ] <- 0
         predmat["x1", ] <- c(-2, 0, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x2", ] <- c(-2, 3, 0, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x3", ] <- c(-2, 3, 3, 0, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x4", ] <- c(-2, 3, 3, 3, 0, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x5", ] <- c(-2, 3, 3, 3, 3, 0, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x6", ] <- c(-2, 3, 3, 3, 3, 3, 0, 3, 1, 1, 1, 1, 1, 1)
         predmat["y", ]  <- c(-2, 3, 3, 3, 3, 3, 3, 0, 1, 1, 1, 1, 1, 1)
         predmat["L21", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1)
         predmat["L22", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1)
         predmat["L23", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1)
         predmat["L24", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1)
         predmat["L25", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1)
         predmat["L26", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)
     } else {
         predmat <- ini$predictorMatrix
         # predmat["ID", ] <- 0
         predmat["x1", ] <- c(-2, 0, 4, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x2", ] <- c(-2, 4, 0, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x3", ] <- c(-2, 4, 4, 0, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x4", ] <- c(-2, 4, 4, 3, 0, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x5", ] <- c(-2, 4, 4, 3, 3, 0, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x6", ] <- c(-2, 4, 4, 3, 3, 3, 0, 3, 1, 1, 1, 1, 1, 1)
         predmat["y", ]  <- c(-2, 4, 4, 3, 3, 3, 3, 0, 1, 1, 1, 1, 1, 1)
         predmat["L21", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1)
         predmat["L22", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1)
         predmat["L23", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1)
         predmat["L24", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1)
         predmat["L25", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1)
         predmat["L26", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)
     }
     
     PredMat <- list(predmat = predmat, impMethod = impMethod)
     
 }
     

 
 predMet2 <- function(data, model){
     ini <- mice(data, maxit = 0)
     impMethod <- character(ncol(data))
     names(impMethod) <- colnames(data)
     impMethod[2:8] <- "2l.norm"
     impMethod[9:14] <- "2lonly.pmm"
     
     if (model == "random intercept") {
         predmat <- ini$predictorMatrix
         # predmat["ID", ] <- 0
         predmat["x1", ] <- c(-2, 0, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x2", ] <- c(-2, 3, 0, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x3", ] <- c(-2, 3, 3, 0, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x4", ] <- c(-2, 3, 3, 3, 0, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x5", ] <- c(-2, 3, 3, 3, 3, 0, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x6", ] <- c(-2, 3, 3, 3, 3, 3, 0, 3, 1, 1, 1, 1, 1, 1)
         predmat["y", ]  <- c(-2, 3, 3, 3, 3, 3, 3, 0, 1, 1, 1, 1, 1, 1)
         predmat["L21", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1)
         predmat["L22", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1)
         predmat["L23", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1)
         predmat["L24", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1)
         predmat["L25", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1)
         predmat["L26", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)
     } else {
         predmat <- ini$predictorMatrix
         # predmat["ID", ] <- 0
         predmat["x1", ] <- c(-2, 0, 4, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x2", ] <- c(-2, 4, 0, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x3", ] <- c(-2, 4, 4, 0, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x4", ] <- c(-2, 4, 4, 3, 0, 3, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x5", ] <- c(-2, 4, 4, 3, 3, 0, 3, 3, 1, 1, 1, 1, 1, 1)
         predmat["x6", ] <- c(-2, 4, 4, 3, 3, 3, 0, 3, 1, 1, 1, 1, 1, 1)
         predmat["y", ]  <- c(-2, 4, 4, 3, 3, 3, 3, 0, 1, 1, 1, 1, 1, 1)
         predmat["L21", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1)
         predmat["L22", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1)
         predmat["L23", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1)
         predmat["L24", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1)
         predmat["L25", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1)
         predmat["L26", ]<- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)
     }
     
     PredMat <- list(predmat = predmat, impMethod = impMethod)
     
 }
 

 
 # fit a model on imputed list
 imp.fit <- function(data, model){
     if(model == "random intercept" ) {

         imp.fit <- lapply(data, FUN=function(xx){
             lmer(y ~  x1 + x2 + x3 + x4 + x5 + x6 + L21 + L22 +L23 + L24 + L25 + L26  +(1|ID), data=xx)
         })
     }else{
         imp.fit <- lapply(data, FUN=function(xx){
             lmer(y ~  x1 + x2 + x3 + x4 + x5 + x6 + L21 + L22 +L23 + L24 + L25 + L26 + (1 + x1 + x2 |ID), data = xx)
         })

     }
     return(imp.fit)
 }



 ols.fit <- function(data){
         
         ols.fit <- lapply(data, FUN=function(xx){
             OLSData <- xx
             OLSData$ID <- as.factor(OLSData$ID)
             MyModelOLS <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + L21 + L22 +L23 + L24 + L25 + L26 + ID, data = OLSData)
         })
     return(ols.fit)
         
 }











 # pool goodness of measure
 # loglikelihood
 LogLs <- function(impfit){

     m <- length(impfit)
     LL <- sapply(impfit, FUN = function(xx) logLik(xx, REML = FALSE))
     # between variance
     #VB <- var(LL)
     # within variance
     #VW <- mean((LL-mean(LL))^2)
     # pooled LL
     pooled <- mean(LL)
     # standard error
     #SE <- sqrt(VB/m)
     return(pooled)

 }

 #deviance and SE
 # Ds <- function(impfit){
 #
 #     m <- length(impfit)
 #     DD <- sapply(impfit, FUN = function(xx) deviance(xx, REML = FALSE))
 #     # between variance
 #     VB <- var(DD)
 #     # within variance
 #     #VW <- mean((DD-mean(DD))^2)
 #     # pooled deviance
 #     pooled <- mean(DD)
 #     # standard error
 #     SE <- sqrt(VB/m)
 #     cbind(pooled, SE)
 #
 # }
