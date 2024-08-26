#This file generates the dataset in the correct format for the simulation study.

##Load the necessary libraries
library(Rcpp)
library(mvtnorm)
library(mipfp)
###

##
#install.packages("Rcpp", dependencies = TRUE)
#install.packages("mvtnorm", dependencies = TRUE)
#install.packages("mipfp", dependencies = TRUE)
##

##Load the necessary C++ script file
sourceCpp("likelihoodwithAGH.cpp")

#The baseline hazards functions are:
#$h_0(t) = 5+0.2t$
#$r_0(t) = 8+0.2t$

#Global variables
rho=0.25 #correlation coefficient
b=0 
n=500 #sample size
B=200 #number of Monte Carlo replications
Z_pn=10 #Number of covariates
Z_pncont <- ceiling(Z_pn/2) # Number of continuous covariates
Z_pnbin <- Z_pn - Z_pncont # Number of binary covariates
beta_1true <- c(1,rep(0,Z_pn-2),-1) #true values of beta_1
beta_2true <- c(1,-0.5,rep(0,Z_pn-2)) #true values of beta_2
gamma=1 
phi=1
niStore <- matrix(0, nrow = B, ncol = n)
cens <- numeric(B)

#Set seed
set.seed(2023)

while(b < B){
  
  recevent_datacontandbin <- function(beta1, beta2, n, gamma, psi, rho, cont, bin) {
    
    if(cont > 0){
      mu_Z1 <- rep(0, cont)
      Sigma_Z1 <- matrix(0, nrow = cont, ncol = cont)
      
      for(ii in 1:cont){
        for(jj in 1:cont){
          Sigma_Z1[ii,jj] = (rho)^{abs(ii-jj)}
        }
      }
      
      Zcont <- rmvnorm(n, mu_Z1, Sigma_Z1) #simulate covariates
    }
    
    if(bin > 0){
      
      mu_Z2 <- rep(0.5,bin)
      Sigma_Z2 <- matrix(0, nrow = bin, ncol = bin)
      
      for(ii in 1:bin){
        for(jj in 1:bin){
          Sigma_Z2[ii,jj] = (rho)^{abs(ii-jj)}
        }
      }
      
      multdist <- ObtainMultBinaryDist(corr=Sigma_Z2, marg.probs=mu_Z2)
      Zbin <- RMultBinary(n, mult.bin.dist = multdist)$binary.sequences
    }
    
    Z <- cbind(Zcont,Zbin)
    C <- runif(n, 0, 2) #censoring time
    
    nu <- rnorm(n, 0, sd=psi)
    u1 <- runif(n)
    temp1 = -log(u1)/exp(Z%*%beta2 + gamma*nu) #inverse sampling
    D <- -25 + 5*sqrt(0.4*temp1 + 25)
    X <- numeric(n)
    
    for(kk in 1:n){
      X[kk] <- min(C[kk], D[kk])
    }
    
    delta <- as.numeric(D <= C) 
    T <- matrix(0, nrow = n, ncol = 2000)
    id <- vector()
    count_row = 0
    rec_start = 1
    
    for(i in 1:n){
      #set a seed
      u <- runif(1,0,1)
      temp2 <- as.vector(-log(u)/exp(Z[i,]%*%beta1 + nu[i]))
      tempT <- -40 + 5*sqrt(0.4*temp2 + 64)
      k = 1
      
      while(tempT <= X[i]){
        T[i,k] = tempT
        u = runif(1,0,1)
        temp2 = -log(u)/(exp(nu[i] + Z[i,]%*%beta1)) + 8*tempT + 0.1*tempT^2
        tempT = -40 + 5*sqrt(0.4*temp2 + 64)
        k = k+1
      }
      
      count_row = count_row + 1
      rec_end = rec_start + k -1
      id[rec_start:rec_end] = rep(i, k) 
      rec_start = rec_end + 1  
    }
    
    Zs <- matrix(0, nrow = rec_end, ncol = length(beta1))
    Xs = vector()
    deltas = vector()
    Ddeltas = vector()
    T_s = vector()
    temp_start = 1
    
    for(l in 1:n){
      temp_count <- sum(id == l)
      temp_end <- temp_start+temp_count-1
      Zs[temp_start:temp_end,] <- matrix(rep(Z[l,],temp_count), nrow = temp_count, byrow = TRUE)
      Xs[temp_start:temp_end] <- rep(X[l], temp_count)
      deltas[temp_start:temp_end] <- rep(delta[l], temp_count)
      Ddeltas[temp_start] <- 1
      if(temp_count == 1){
        T_s[temp_start] <- X[l]
      } else {
        T_s[temp_start] <- X[l]
        T_s[(temp_start+1):temp_end] <- T[l, 1:(temp_count-1)]
        Ddeltas[(temp_start+1):temp_end] <- 0
      }
      temp_start = temp_end + 1
    }
    
    N <- rec_end - n
    data_rec <- data.frame(cbind(Zs,Xs,T_s,id,deltas,Ddeltas,rep(N,rec_end)))
    #colnames(data_rec) <- c("Z1","Z2","Z3","Z4","Z5","Z6","X","T","id","delta","Ddelta","N")
    #function returns a single dataset
    return(data_rec)
  }
  
  my_data <- recevent_datacontandbin(beta_1true, beta_2true, n, gamma, phi, rho, Z_pncont, Z_pnbin)
  cat("Total number of recurrent events is:", my_data[1,ncol(my_data)], "\n")
  recur_data <- RecurDataCpp(data.matrix(my_data)) #works
  recur_indicators <- RecurIndCpp(data.matrix(my_data), n)
  Zs <- ZsCpp(data.matrix(my_data), Z_pn, n)
  ni <- niCpp(recur_indicators, recur_data)
  delta <- TermIndCpp(data.matrix(my_data))
  cens[b+1] <- (1-sum(delta)/n)
  
  if(max(ni) > 100){ 
    next;
  }
  
  #Comment: For any given dataset, if the maximum number of recurrent events of 
  #any observation exceeds 100, skip to the next dataset
  niStore[b+1, ] <- ni
  b=b+1

  write.csv(my_data, file=paste("MixedData",b,".csv",sep=""), row.names=FALSE)

}

#Check summary statistics regarding number of recurrent events and censoring
apply(niStore, 1, max)
apply(niStore, 1, sum)
apply(niStore,1,mean)
mean(delta)
mean(cens)
summary(cens)
