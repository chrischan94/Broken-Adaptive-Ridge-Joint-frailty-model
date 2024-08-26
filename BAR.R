rm(list=ls())
cat("The program starts at", date(), "\n")

##Load the necessary libraries##
library(fastGHQuad)
library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)

###Load files

sourceCpp("likelihoodwithAGH.cpp")
sourceCpp("BARalgorithmwAGQ.cpp")
source("MuHatfunctions.R")

##Global variables##
Z_pn = 10
pn = 2*Z_pn
beta_1true <- c(1,rep(0, (Z_pn-2)),-1)
beta_2true <- c(1,-0.5, rep(0, (Z_pn-2)))
z_pn = pn-4
qn = pn - z_pn
np = 5 #number of pieces
B = 200 #number of datasets to simulate 
C = 60 #maximum number of iterations allowed for the BAR method
n = 500
gamma = 1
phi = 1
tol = 1e-4 #tolerance for the final estimates of the BAR algorithm
udiff = 0.025 #The interval of $u_i$ for the Gauss-Hermite quadrature. 
rtol = 0.001 
d = 10^(-4) #Small positive constant added to the diagonal entries of $D(\widehat{beta})$ 
truevalues <- c(beta_1true, beta_2true, rep(5,np), rep(8,np), gamma,phi)
estimates <- matrix(0, nrow = B, ncol = length(truevalues))
useq = seq(-2.5, 2.5, udiff)
rule10 <- gaussHermiteData(10)
rule50 <- gaussHermiteData(50) #for BAR algorithm
lambda = seq(2.5,4.5,by=0.5)
sigma=0.1

data_recurrent <- function(data){
  #function that returns the recurrent events data
  a <- which(data[,ncol(data)-1] == 0)
  data_reduc <- data[a,]
  id <- data_reduc[,ncol(data)-3]
  T_r <- data_reduc[,ncol(data)-4]
  return(cbind(T_r,id))
}

estimates <- read.csv("sigma10n500p10BARNEGNEW.csv", header=TRUE)
set.seed(2023)

for(b in 1:B){
  
  skip = FALSE
  
  theta_init = truevalues + rnorm(length(truevalues), 0, sigma)
  
  cat("The initial values are:", theta_init, "\n")
  name <- paste("MixedData",b,".csv",sep="")
  my_data <- read.csv(name, header = TRUE)
  cat("Total number of recurrent events is:", my_data[1,ncol(my_data)], "\n")
  
  recur_data <- RecurDataCpp(data.matrix(my_data)) 
  recur_indicators <- RecurIndCpp(data.matrix(my_data), n) #Indicator of recurrent events
  Zs <- ZsCpp(data.matrix(my_data), Z_pn, n) #covariate (design) matrix
  delta <- TermIndCpp(data.matrix(my_data)) #terminal event indicator
  X <- YCpp(data.matrix(my_data)) #follow-up time
  ni <- niCpp(recur_indicators, recur_data) #number of recurrent events
  
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
  
  Uhat <- tryCatch(full_uhatCpp(theta_init), error = function(e) {skip <<- TRUE})
  
  if(skip == TRUE){
    cat("Problem with Uhat", "\n")
    next;
  }
  
  likelihoodCpp <- function(theta){
    #initialize the parameters
    beta_r <- theta[1:Z_pn]
    beta_h <- theta[(Z_pn+1):pn]
    hs <- theta[(pn+1):(pn+np)]
    rs <- theta[(pn+np+1):(pn+2*np)]
    gamma <- theta[length(theta_init)-1]
    psi <- theta[length(theta_init)]
    
    sdH <- numeric(n)
    adapQuad <- numeric(n)
    log_likid_uc_Cpp <- numeric(n)
    
    for(i in 1:n){
      
      term_lik <- function(u){
        #returns the likelihood of terminal event per sample id
        ((exp(elemmultsum(beta_h,Zs[i,]) + gamma*u)*elemmultsum(hs, smhrCpp(X, X[i])))^delta[i])*
          (exp(-elemmultsum(hs, caphrCpp(X, X[i]))*exp(elemmultsum(beta_h,Zs[i,]) + gamma*u)))
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
      
      mll2 <- function(u){
        return(term_lik(u)*dnorm(u, mean=0, sd=psi)*rec_lik(u))
      }
      
      #compute the value of the 2nd derivative with respect to the mode
      sdH[i] <- -gamma^2*(elemmultsum(hs, caphrCpp(X,X[i]))*exp(elemmultsum(beta_h, Zs[i,]) + gamma*Uhat[i])) - 
        (elemmultsum(rs, caphrCpp(X,X[i]))*exp(elemmultsum(beta_r, Zs[i,]) + Uhat[i])) - psi^(-2)
      adapQuad[i] <- aghQuad(g = mll2, muHat = Uhat[i], sigmaHat = sqrt(-1/sdH[i]), rule = rule10)
      log_likid_uc_Cpp[i] <- log(adapQuad[i])
    }
    return(-sum(log_likid_uc_Cpp))
  }
  
  thetahat0 <- nlm(f=likelihoodCpp,p=theta_init, iterlim=150, steptol = 1e-7)$estimate
  
  MuHat1 = full_uhatCpp(thetahat0)
  recur_data <- data.frame(data_recurrent(my_data))
  MuHat2 = full_uhatCpp2(thetahat0)
  MuHat3 = full_uhatCpp3(thetahat0)
  MuHat4 = full_uhatCpp4(thetahat0)
  MuHat5 = full_uhatCpp5(thetahat0)
  MuHat6 = full_uhatCpp6(thetahat0)
  MuHats <- cbind(MuHat1, MuHat2, MuHat3, MuHat4, MuHat5, MuHat6)
  
  recur_data <- RecurDataCpp(data.matrix(my_data))
  sdHs = sigmaHCpp(thetahat0, X,Zs,ni,recur_data, delta, MuHats,Z_pn,Z_pn,np,np)
  
  GCV <- numeric(length(lambda))
  recur_data <- RecurDataCpp(data.matrix(my_data))
  ests <- matrix(0, nrow=length(lambda), ncol=length(theta_init))
  
  for(l in 1:length(lambda)){
    
    maxdiff=2
    skip = FALSE
    
    estloop <- matrix(0, nrow = C+1, ncol = length(theta_init))
    estloop[1, ] <- thetahat0
    
    betar <- thetahat0[1:Z_pn]
    betah <- thetahat0[(Z_pn+1):pn]
    hs <- thetahat0[(pn+1):(pn+5)]
    rs <- thetahat0[(pn+6):(pn+10)]
    gamma <- thetahat0[length(thetahat0)-1]
    phi <- thetahat0[length(thetahat0)]
    
    LikeAGQlp <- matrix(0, nrow = C, ncol = n)
    g3AGQ <- matrix(0, nrow = C, ncol = n)
    g4AGQ <- matrix(0, nrow = C, ncol = n)
    g5AGQ <- matrix(0, nrow = C, ncol = n)
    g6AGQ <- matrix(0, nrow = C, ncol = n)
    g7AGQ <- matrix(0, nrow = C, ncol = n)
    a=1
    
    while(maxdiff > rtol){
      
      par1_beta1 <- matrix(0, nrow = length(beta_1true), ncol = n)
      par1_beta2 <- matrix(0, nrow = length(beta_2true), ncol = n)
      par2_beta11 <- array(0, dim = c(length(beta_1true), length(beta_1true), n))
      par2_beta12 <- array(0, dim = c(length(beta_1true), length(beta_1true), n))
      par2_beta22 <- array(0, dim = c(length(beta_1true), length(beta_1true), n))
      
      for(i in 1:n){
        
        g1 <- function(u){ #corresponds to $g_1(Y_i|u_i)$ 
          exp(-elemmultsum(hs, caphrCpp(X, X[i]))*exp(elemmultsum(betah, Zs[i,]) + gamma*u))*
            ((exp(elemmultsum(betah, Zs[i,]) + gamma*u)*elemmultsum(hs, smhrCpp(X,X[i])))^delta[i])
        }
        
        g2 <- function(u){ #corresponds to $g_2(Y_i|u_i)$
          exp(-elemmultsum(rs, caphrCpp(X,X[i]))*exp(elemmultsum(betar, Zs[i,]) + u))
        }
        
        r0Tik <- function(u){ #corresponds to $\prod^{n_i}_{k=1} \widetilde{r}_0(T_{ik}) \exp(\boldsymbol{\beta}^\top_1 \textbf{Z}_i + u_i)$
          if(ni[i] > 0){
            recur_dataid <- recur_data[which(recur_data[,2] == i),]
            recur_prod <- prodRecurCpp(rs, recur_data, ni[i], i)
            recur <- recur_prod*(exp(elemmultsum(betar, Zs[i,]) + u))^(ni[i])
          } else {
            recur <- 1
          }
          return(recur)
        }
        
        Likelihood <- function(u){
          #corresponds to $L_i$
          g1(u)*g2(u)*r0Tik(u)*dnorm(u, mean = 0, sd = phi)
        }
        
        g3 <- function(u){ 
          #corresponds to $g_3(Y_i)$
          g1(u)*g2(u)*r0Tik(u)*(ni[i] - (elemmultsum(rs,caphrCpp(X,X[i]))*exp(elemmultsum(betar,Zs[i,]) + u)))*dnorm(u, mean = 0, sd = phi)
        }
        
        g4 <- function(u){ 
          #corresponds to $g_4(Y_i)$
          g1(u)*g2(u)*r0Tik(u)*(delta[i] - (elemmultsum(hs,caphrCpp(X,X[i]))*exp(elemmultsum(betah,Zs[i,]) + gamma*u)))*dnorm(u, mean = 0, sd = phi)
        }
        
        g5 <- function(u){
          #corresponds to $m_4(u_i; Y_i)$
          g1(u)*g2(u)*r0Tik(u)*((ni[i] - (elemmultsum(rs,caphrCpp(X,X[i]))*exp(elemmultsum(betar, Zs[i,]) + u)))^2 - (elemmultsum(rs, caphrCpp(X,X[i]))*exp(elemmultsum(betar,Zs[i,]) + u)))*dnorm(u, mean = 0,sd = phi)
        }
        
        g6 <- function(u){
          #corresponds to $m_5(u_i; Y_i)$
          g1(u)*g2(u)*r0Tik(u)*((delta[i] - elemmultsum(hs, caphrCpp(X,X[i]))*exp(elemmultsum(betah, Zs[i,]) + gamma*u))^2 - (elemmultsum(hs,caphrCpp(X,X[i]))*exp(elemmultsum(betah,Zs[i,]) + gamma*u)))*dnorm(u, mean = 0, sd=phi)
        }
        
        g7 <- function(u){
          #corresponds to $m_6(u_i; Y_i)$
          g1(u)*g2(u)*r0Tik(u)*(ni[i] - (elemmultsum(rs,caphrCpp(X,X[i]))*exp(elemmultsum(betar,Zs[i,]) + u)))*(delta[i] - (elemmultsum(hs,caphrCpp(X,X[i]))*exp(elemmultsum(betah,Zs[i,]) + gamma*u)))*dnorm(u, mean = 0, sd = phi)
        }
        
        LikeAGQlp[a,i] <- aghQuad(g = Likelihood, muHat = MuHats[i,1], sigmaHat = sqrt(-1/sdHs[i,1]), rule = rule50)
        g3AGQ[a,i] <- aghQuad(g = g3, muHat = MuHats[i,2], sigmaHat = sqrt(-1/sdHs[i,2]), rule = rule50)
        g4AGQ[a,i] <- aghQuad(g = g4, muHat = MuHats[i,3], sigmaHat = sqrt(-1/sdHs[i,3]), rule = rule50)
        g5AGQ[a,i] <- aghQuad(g = g5, muHat = MuHats[i,4], sigmaHat = sqrt(-1/sdHs[i,4]), rule = rule50)
        g6AGQ[a,i] <- aghQuad(g = g6, muHat = MuHats[i,5], sigmaHat = sqrt(-1/sdHs[i,5]), rule = rule50)
        g7AGQ[a,i] <- aghQuad(g = g7, muHat = MuHats[i,6], sigmaHat = sqrt(-1/sdHs[i,6]), rule = rule50) 
        
        par1_beta1[,i] <- t(t(Zs[i,]))*(g3AGQ[a,i]/LikeAGQlp[a,i])
        par1_beta2[,i] <- t(t(Zs[i,]))*(g4AGQ[a,i]/LikeAGQlp[a,i])
        
        if(is.nan(g3AGQ[a,i]) | is.nan(g4AGQ[a,i]) | is.nan(g5AGQ[a,i]) | is.nan(g6AGQ[a,i]) | is.nan(g7AGQ[a,i])){
          cat("Problems with convergence", "\n")
          break;
        }
        ####
        #Comment: The following ifelse statements are required to prevent the numerator from becoming Infinity
        if(g5AGQ[a,i] > 0){
          par2_beta11[,,i] <- t(t(Zs[i,]))%*%t(Zs[i,])*((min(g5AGQ[a,i]*LikeAGQlp[a,i], 1.7e+308) - min(g3AGQ[a,i]^2, 1.7e+308))/(LikeAGQlp[a,i]^2))
        } else {
          par2_beta11[,,i] <- t(t(Zs[i,]))%*%t(Zs[i,])*((max(max(g5AGQ[a,i]*LikeAGQlp[a,i], -1.7e+308) - min(g3AGQ[a,i]^2, 1.7e+308),-1.7e+308))/(LikeAGQlp[a,i]^2))
        }
        
        if(g6AGQ[a,i] > 0){
          par2_beta22[,,i] <- t(t(Zs[i,]))%*%t(Zs[i,])*((min(g6AGQ[a,i]*LikeAGQlp[a,i], 1.7e+308) - min(g4AGQ[a,i]^2, 1.7e+308))/(LikeAGQlp[a,i]^2))
        } else {
          par2_beta22[,,i] <- t(t(Zs[i,]))%*%t(Zs[i,])*((max(max(g6AGQ[a,i]*LikeAGQlp[a,i], -1.7e+308) - min(g4AGQ[a,i]^2, 1.7e+308),-1.7e+308))/(LikeAGQlp[a,i]^2))
        }
        
        if(g7AGQ[a,i] > 0){
          if(sign(g3AGQ[a,i]) != sign(g4AGQ[a,i])){
            par2_beta12[,,i] <- t(t(Zs[i,]))%*%t(Zs[i,])*((min(min(g7AGQ[a,i]*LikeAGQlp[a,i], 1.7e+308) - max(g3AGQ[a,i]*g4AGQ[a,i], -1.7e+308),1.7e+308))/(LikeAGQlp[a,i]^2))
          } else {
            par2_beta12[,,i] <- t(t(Zs[i,]))%*%t(Zs[i,])*((min(g7AGQ[a,i]*LikeAGQlp[a,i], 1.7e+308) - min(g3AGQ[a,i]*g4AGQ[a,i], 1.7e+308))/(LikeAGQlp[a,i]^2))
          }
        } 
        
        if (g7AGQ[a,i] < 0){
          if(sign(g3AGQ[a,i]) != sign(g4AGQ[a,i])){
            par2_beta12[,,i] <- t(t(Zs[i,]))%*%t(Zs[i,])*((max(g7AGQ[a,i]*LikeAGQlp[a,i], -1.7e+308) - max(g3AGQ[a,i]*g4AGQ[a,i], -1.7e+308))/(LikeAGQlp[a,i]^2))
          } else{
            par2_beta12[,,i] <- t(t(Zs[i,]))%*%t(Zs[i,])*((max(max(g7AGQ[a,i]*LikeAGQlp[a,i], -1.7e+308) - min(g3AGQ[a,i]*g4AGQ[a,i], 1.7e+308),-1.7e+308))/(LikeAGQlp[a,i]^2))
          }
        }
        #####
        
        if(!is.finite(LikeAGQlp[a,i]^2)){
          par2_beta11[,,i] <- matrix(0, nrow=Z_pn, ncol=Z_pn)
          par2_beta22[,,i] <- matrix(0, nrow=Z_pn, ncol=Z_pn)
          par2_beta12[,,i] <- matrix(0, nrow=Z_pn, ncol=Z_pn)
        }
      }
      
      par1_betas <- c(rowSums(par1_beta1), rowSums(par1_beta2)) #vector that stores 
      par2_betas <- matrix(0, nrow = 2*length(beta_1true), ncol = 2*length(beta_2true))
      par2_betas[1:length(beta_1true), 1:length(beta_1true)] <- rowSums(par2_beta11, dim = 2)
      par2_betas[(length(beta_1true)+1):(2*length(beta_1true)), (length(beta_2true)+1):(2*length(beta_2true))] <- rowSums(par2_beta22, dim = 2)
      par2_betas[1:length(beta_1true), (length(beta_2true)+1):(2*length(beta_2true))] <- rowSums(par2_beta12, dim = 2)
      par2_betas[(length(beta_1true)+1):(2*length(beta_1true)), 1:length(beta_2true)] <- rowSums(par2_beta12, dim = 2)
      
      Omega_n <- -par2_betas
      v_n <- par1_betas - par2_betas%*%(c(betar, betah))
      ####
      d = 10^(-4)
      D <- diag(x = 1/(c(betar^2, betah^2) + d^2), pn, pn)
      beta_bar <- tryCatch(solve(Omega_n + lambda[l]*D)%*%v_n, error = function(e) { skip <<- TRUE})
      
      if(skip == TRUE){
        cat("Omega_n was singular couldn't be inverted", "\n")
        break;
      }
      
      betar <- beta_bar[(1:Z_pn)]
      betah <- beta_bar[((Z_pn+1):pn)]
      estloop[a+1,(1:pn)] = t(beta_bar)
      estloop[a+1,-c(1:pn)] = c(hs,rs,gamma,phi)
      
      NuisUpd_Cpp <- function(theta){
        Hs <- theta[1:np]
        Rs <- theta[(np+1):(2*np)]
        Gamma <- theta[2*np + 1]
        Psi <- theta[length(theta)]
        
        sdH <- numeric(n)
        adapQuad <- numeric(n)
        log_likid_uc_Cpp <- numeric(n)
        
        for(i in 1:n){
          
          term_lik <- function(u){
            #returns the likelihood of terminal event per sample id
            ((exp(elemmultsum(betah, Zs[i,]) + Gamma*u)*elemmultsum(Hs, smhrCpp(X, X[i])))^delta[i])*
              (exp(-elemmultsum(Hs, caphrCpp(X, X[i]))*exp(elemmultsum(betah, Zs[i,]) + Gamma*u)))
          }
          
          rec_lik <- function(u){
            #returns the likelihood for recurrent events per sample id
            if(recur_indicators[i] == 1){
              recur_dataid <- recur_data[which(recur_data[,2] == i),]
              nid <- nrow(recur_dataid)
              recur_prod <- prodRecurCpp(Rs, recur_data, ni[i], i)
              recur <- exp(-exp(elemmultsum(betar, Zs[i,]) + u)*(elemmultsum(Rs, caphrCpp(X,X[i]))))*recur_prod*(exp(elemmultsum(betar, Zs[i,]) + u))^(ni[i])
            } else {
              recur <- exp(-exp(elemmultsum(betar, Zs[i,]) + u)*(elemmultsum(Rs, caphrCpp(X,X[i]))))
            }
            return(recur)
          }
          
          mll3 <- function(u){
            return(term_lik(u)*dnorm(u, mean=0, sd=Psi)*rec_lik(u))
          }
          
          #compute the value of the 2nd derivative with respect to the mode
          sdH[i] <- -gamma^2*(elemmultsum(Hs, caphrCpp(X,X[i]))*exp(elemmultsum(betah, Zs[i,]) + Gamma*MuHat1[i])) - 
            (elemmultsum(Rs, caphrCpp(X,X[i]))*exp(elemmultsum(betar, Zs[i,]) + MuHat1[i])) - Psi^(-2)
          adapQuad[i] <- aghQuad(g = mll3, muHat = MuHat1[i], sigmaHat = sqrt(-1/sdH[i]), rule = rule10)
          log_likid_uc_Cpp[i] <- log(adapQuad[i])
        }
        return(-sum(log_likid_uc_Cpp))
      }
      
      #Update the nuisance parameters given $\widehat{\beta}^{(m+1)}$
      nuisnext1 = nlm(f=NuisUpd_Cpp, p=estloop[a,(pn+1):length(theta_init)], iterlim=1)$estimate
      
      hs = nuisnext1[1:np]
      rs = nuisnext1[(np+1):(2*np)]
      gamma = nuisnext1[2*np+1]
      phi = nuisnext1[2*np+2]
      
      estloop[a+1,-c(1:pn)] = nuisnext1
      
      a=a+1
      
      maxdiff = max(abs(estloop[a,1:pn] - estloop[a-1,1:pn]))
      
      if(is.nan(maxdiff)){
        cat("Problems with convergence", "\n")
        break;
      }
      
      cat("The maximum difference for iteration", a-1, "is", maxdiff, "\n")
      
      if(a > C){
        cat("The algorithm has reached its limit", "\n")
        break; 
      }
    }
    
    if(is.nan(maxdiff)){
      cat("Problems with convergence", "\n")
      next;
    }
    
    j = length(estloop[,1][estloop[,1] != 0])
    final_beta <- estloop[j,1:pn]
    hs_GCV <- estloop[j,(pn+1):(pn+np)]
    rs_GCV <- estloop[j,(pn+np):(pn+2*np)]
    gamma_GCV <- estloop[j,pn+2*np+1]
    phi_GCV <- estloop[j,pn+2*np+2]
    b1 = length(final_beta[1:Z_pn][abs(final_beta[1:Z_pn]) > tol ])
    b2 = length(final_beta[(Z_pn+1):pn][abs(final_beta[(Z_pn+1):pn]) > tol ])
    br1 = which(abs(final_beta[1:Z_pn]) > tol )
    br2 = which(abs(final_beta[(Z_pn+1):pn]) > tol)
    
    likelihoodCpp_GCV <- function(theta, b1, b2){
      #initialize the parameters
      beta_r <- theta[1:b1]
      beta_h <- theta[(b1+1):(b1+b2)]
      
      sdH <- numeric(n)
      adapQuad <- numeric(n)
      log_likid_uc_Cpp <- numeric(n)
      
      for(i in 1:n){
        
        term_lik <- function(u){
          #returns the likelihood of terminal event per sample id
          ((exp(sum(beta_h*Zs[i,br2]) + gamma*u)*elemmultsum(hs_GCV, smhrCpp(X, X[i])))^delta[i])*
            (exp(-elemmultsum(hs_GCV, caphrCpp(X, X[i]))*exp(sum(beta_h*Zs[i,br2]) + gamma*u)))
        }
        
        rec_lik <- function(u){
          #returns the likelihood for recurrent events per sample id
          if(recur_indicators[i] == 1){
            recur_prod <- prodRecurCpp(rs_GCV, recur_data, ni[i], i)
            recur <- exp(-exp(sum(beta_r*Zs[i,br1]) + u)*(elemmultsum(rs_GCV, caphrCpp(X,X[i]))))*recur_prod*(exp(sum(beta_r*Zs[i,br1])+u))^(ni[i])
          } else {
            recur <- exp(-exp(sum(beta_r*Zs[i,br1]) + u)*(elemmultsum(rs_GCV, caphrCpp(X,X[i]))))
          }
          return(recur)
        }
        
        mll4 <- function(u){
          return(term_lik(u)*dnorm(u, mean=0, sd=phi_GCV)*rec_lik(u))
        }
        
        #compute the value of the 2nd derivative with respect to the mode
        sdH[i] <- -gamma_GCV^2*(elemmultsum(hs_GCV, caphrCpp(X,X[i]))*exp(sum(beta_h*Zs[i,br2]) + gamma*Uhat[i])) - 
          (elemmultsum(rs_GCV, caphrCpp(X,X[i]))*exp(sum(beta_r*Zs[i,br1]) + Uhat[i])) - phi_GCV^(-2)
        adapQuad[i] <- aghQuad(g = mll4, muHat = Uhat[i], sigmaHat = sqrt(-1/sdH[i]), rule = rule10)
        log_likid_uc_Cpp[i] <- log(adapQuad[i])
      }
      return(-sum(log_likid_uc_Cpp))
    }
    unbias_nz <- nlm(f=likelihoodCpp_GCV, p=final_beta[abs(final_beta) > tol], b1=b1, b2=b2, iterlim=20)$estimate
    
    D_beta <- diag(unbias_nz) #Hessian matrix
    Sigma_lambda <- numeric(length(unbias_nz))
    
    for(i in 1:length(unbias_nz)){
      Sigma_lambda[i] <- (1/(abs(unbias_nz[i])))
    }
    
    r_betahatD <- lambda[l]*diag(Sigma_lambda) 
    s_l <- sum(diag(solve(D_beta - r_betahatD)%*%D_beta))
    GCV[l] <- likelihoodCpp_GCV(unbias_nz,b1,b2)/(n*(1-s_l/n)^2)
    cat("The Generalized CV value is", GCV[l], "\n \n")
    ests[l,] <- estloop[j,]
  }
  
  if(is.nan(maxdiff)){
    cat("Skip to the next dataset", "\n \n")
    next;
  }
  
  opt_lam <- lambda[which.min(GCV)]
  cat("The optimal lambda is", opt_lam, "\n")
  opt_ests <- ests[which.min(GCV), ]
  final_beta <- opt_ests[1:pn]
  final_beta <- final_beta*(abs(final_beta) > tol)
  cat("Final estimates of beta for dataset", b, "is", final_beta, "\n ")
  estimates[b, ] <- c(final_beta, opt_ests[-c(1:pn)])
}