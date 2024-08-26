rm(list=ls())
cat("The program starts at", date(), "\n")

##
install.packages("fastGHQuad", dependencies = TRUE)
install.packages("Rcpp", dependencies = TRUE)
install.packages("RcppArmadillo", dependencies = TRUE)
install.packages("mvtnorm", dependencies = TRUE)
##
##Load the necessary libraries##
library(fastGHQuad)
library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)

##

##
###Load files
sourceCpp("MICagq.cpp")
sourceCpp("MICAGQmanualcacl.cpp")

##Global variables##
Z_pn = 10
pn = 2*Z_pn #Total number of covariates 
beta_1true <- c(1,rep(0, (Z_pn-2)),-1)
beta_2true <- c(1,-0.5, rep(0, (Z_pn-2)))
np = 5 #number of pieces
B = 200 #number of datasets  
n = 500
gamma = 1
phi = 1
udiff = 0.025 #The interval of $u_i$ for the Gauss-Hermite quadrature. 
truevalues <- c(beta_1true,beta_2true, rep(5,np), rep(8,np), gamma,phi) 
estimates <- matrix(0, nrow = B, ncol = length(truevalues))
useq = seq(-2.5, 2.5, udiff)
rule10 <- gaussHermiteData(10)
sigma=0.1

data_recurrent <- function(data){
  #function that returns the recurrent events data
  a <- which(data[,ncol(data)-1] == 0)
  data_reduc <- data[a,]
  id <- data_reduc[,ncol(data)-3]
  T_r <- data_reduc[,ncol(data)-4]
  return(cbind(T_r,id))
}

set.seed(2023)

for(b in 1:B){
  
  skip = FALSE
  
  theta_init = truevalues + rnorm(length(truevalues), 0, sigma)
  
  cat("The initial values are:", theta_init, "\n")
  name <- paste("MixedData",b,".csv",sep="")
  my_data <- read.csv(name, header = TRUE)
  cat("Total number of recurrent events is:", my_data[1,ncol(my_data)], "\n")
  
  recur_data <- RecurDataCpp(data.matrix(my_data)) 
  recur_indicators <- RecurIndCpp(data.matrix(my_data), n) 
  Zs <- ZsCpp(data.matrix(my_data), Z_pn, n) 
  delta <- TermIndCpp(data.matrix(my_data)) 
  X <- YCpp(data.matrix(my_data)) 
  ni <- niCpp(recur_indicators, recur_data) 
  
  full_uhatCpp <- function(theta){
    #function to find the Uhat
    beta_r <- theta[1:Z_pn]
    beta_h <- theta[(Z_pn+1):(2*Z_pn)]
    hs <- theta[(2*Z_pn+1):(2*Z_pn+np)]
    rs <- theta[(2*Z_pn+np+1):(2*(Z_pn+np))]
    gamma <- theta[2*(Z_pn+np) + 1]
    psi <- theta[length(theta)]
    
    uhat <- numeric(n)
    for(i in 1:nrow(Zs)){
      
      term_lik <- function(u){
        #returns the likelihood of terminal event per sample id
        ((exp(elemmultsum(beta_h, Zs[i,]) + gamma*u)*elemmultsum(hs, smhrCpp(X,X[i])))^delta[i])*
          (exp(-elemmultsum(hs, caphrCpp(X,X[i]))*exp(elemmultsum(beta_h, Zs[i,]) + gamma*u)))
      }
      
      rec_lik <- function(u){
        #returns the likelihood for recurrent events per sample id
        if(recur_indicators[i] == 1){
          recur_prod <- prodRecurCpp(rs, recur_data, ni[i], i)
          recur <- exp(-exp(elemmultsum(beta_r, Zs[i,]) + u)*(elemmultsum(rs, caphrCpp(X, X[i]))))*recur_prod*(exp(elemmultsum(beta_r, Zs[i,])+u))^(ni[i])
        } else {
          recur <- exp(-exp(elemmultsum(beta_r, Zs[i,]) + u)*(elemmultsum(rs, caphrCpp(X, X[i]))))
        }
        return(recur)
      }
      
      mll1 <- function(u){
        #returns the value inside the integral per sample id
        return(term_lik(u)*rec_lik(u)*dnorm(u, mean = 0, sd = psi))
      }
      
      #compute the value of the 2nd derivative with respect to the mode
      uhat[i] <- useq[which.max(mll1(useq))]
      
    }
    return(uhat)
  }
  
  UHat <- full_uhatCpp(theta_init)
  
  MIC_likelihoodCpp <- function(theta){
    #initialize the parameters
    alpha_r <- theta[1:Z_pn]
    alpha_h <- theta[(Z_pn+1):pn]
    hs <- theta[(pn+1):(pn+np)]
    rs <- theta[(pn+np+1):(pn+2*np)]
    gamma <- theta[length(theta_init)-1]
    psi <- theta[length(theta_init)]
    
    ####the hyperbolic function
    a_n <- n/4
    omega_alphar <- numeric(length(alpha_r))
    omega_alphah <- numeric(length(alpha_h))
    for(i in 1:length(alpha_r)){
      omega_alphar[i] <- (exp(2*a_n*alpha_r[i]^2) - 1)/(exp(2*a_n*alpha_r[i]^2) + 1)
      omega_alphah[i] <- (exp(2*a_n*alpha_h[i]^2) - 1)/(exp(2*a_n*alpha_h[i]^2) + 1)
    }
    
    #Reparameterization
    beta_r <- alpha_r*omega_alphar
    beta_h <- alpha_h*omega_alphah
    
    sdH <- numeric(n)
    adapQuad <- numeric(n)
    log_likid_uc_Cpp <- numeric(n)
    
    for(i in 1:n){
      
      term_lik <- function(u){
        #returns the likelihood of terminal event per sample id
        ((exp(elemmultsum(beta_h,Zs[i,]) + gamma*u)*elemmultsum(hs, smhrCpp(X, X[i])))^delta[i])*
          (exp(-elemmultsum(hs, caphrCpp(X, X[i]))*exp(elemmultsum(beta_h, Zs[i,]) + gamma*u)))
      }
      
      rec_lik <- function(u){
        #returns the likelihood for recurrent events per sample id
        if(recur_indicators[i] == 1){
          recur_prod <- prodRecurCpp(rs, recur_data, ni[i], i)
          recur <- exp(-exp(elemmultsum(beta_r, Zs[i,]) + u)*(elemmultsum(rs, caphrCpp(X,X[i]))))*recur_prod*(exp(elemmultsum(beta_r, Zs[i,])+u))^(ni[i])
        } else {
          recur <- exp(-exp(elemmultsum(beta_r, Zs[i,]) + u)*(elemmultsum(rs, caphrCpp(X,X[i]))))
        }
        return(recur)
      }
      
      mll <- function(u){
        return(term_lik(u)*dnorm(u, mean=0, sd=psi)*rec_lik(u))
      }
      
      #compute the value of the 2nd derivative with respect to the mode
      sdH[i] <- -gamma^2*(elemmultsum(hs, caphrCpp(X,X[i]))*exp(elemmultsum(beta_h, Zs[i,]) + gamma*UHat[i])) - 
        (elemmultsum(rs,caphrCpp(X,X[i]))*exp(elemmultsum(beta_r, Zs[i,]) + UHat[i])) - psi^(-2)
      adapQuad[i] <- aghQuad(g = mll, muHat = UHat[i], sigmaHat = sqrt(-1/sdH[i]), rule = rule10)
      log_likid_uc_Cpp[i] <- log(adapQuad[i])
    }
    #the penalty
    pen <- (log(nrow(recur_data))/2)*(sum(omega_alphah) + sum(omega_alphar))
    return(-sum(log_likid_uc_Cpp) + pen)
  }
  
  estimates[b,] <- tryCatch(nlm(f=MIC_likelihoodCpp, p=theta_init, iterlim=240, steptol = 1e-10)$estimate, error = function(e) { skip <<- TRUE})
  
  if(skip == TRUE){
    cat("Non-finite values in optimization", "\n")
  }
  
  estimates[b,] <- estimates[b,]*(abs(estimates[b,]) > tol)
  
  cat("Iteration", b, "was finished at", date(), "\n \n")
  cat("The estimates are", estimates[b,1:pn], "\n")
  
}