#computes the MuHats for the likelihood
full_uhatCpp <- function(theta){
  #initialize the parameters
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
      ((exp(sum(beta_h*Zs[i,]) + gamma*u)*sum(hs*smhrCpp(X,X[i])))^delta[i])*
        (exp(-sum(hs*caphrCpp(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u)))
    }
    
    rec_lik <- function(u){
      #returns the likelihood for recurrent events per sample id
      if(recur_indicators[i] == 1){
        recur_dataid <- recur_data[which(recur_data$id == i),]
        nid <- nrow(recur_dataid)
        prod_recur <- 1
        for(j in 1:nid){
          prod_recur <- prod_recur*(exp(sum(beta_r*Zs[i,]) + u))*sum(rs*smhrCpp(recur_data[,1], recur_dataid[j,1]))
        }
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp(X, X[i]))))*prod_recur
      } else {
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp(X, X[i]))))
      }
      return(recur)
    }
    
    mll <- function(u){
      #returns the value inside the integral per sample id
      return(term_lik(u)*rec_lik(u)*dnorm(u, mean = 0, sd = psi))
    }
    
    #compute the value of the 2nd derivative with respect to the mode
    uhat[i] <- useq[which.max(mll(useq))]
    
    if(ni[i] > 90){ 
      uhat[i] <- 0 
    }
    
  }
  return(uhat)
}

##Computes muHat for the second function
full_uhatCpp2 <- function(theta){
  #initialize the parameters
  beta_r <- theta[1:Z_pn]
  beta_h <- theta[(Z_pn+1):(2*Z_pn)]
  hs <- theta[(2*Z_pn+1):(2*Z_pn+np)]
  rs <- theta[(2*Z_pn+np+1):(2*(Z_pn+np))]
  gamma <- theta[2*(Z_pn+np) - 1]
  psi <- theta[length(theta)]
  
  uhat <- numeric(n)
  
  for(i in 1:nrow(Zs)){
    
    term_lik <- function(u){
      #returns the likelihood of terminal event per sample id
      ((exp(sum(beta_h*Zs[i,]) + gamma*u)*sum(hs*smhrCpp(X,X[i])))^delta[i])*
        (exp(-sum(hs*caphrCpp(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u)))
    }
    
    rec_lik <- function(u){
      #returns the likelihood for recurrent events per sample id
      if(recur_indicators[i] == 1){
        recur_dataid <- recur_data[which(recur_data$id == i),]
        nid <- nrow(recur_dataid)
        prod_recur <- 1
        for(j in 1:nid){
          prod_recur <- prod_recur*(exp(sum(beta_r*Zs[i,]) + u))*sum(rs*smhrCpp(recur_data[,1], recur_dataid[j,1]))
        }
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp(X, X[i]))))*prod_recur
      } else {
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp(X, X[i]))))
      }
      return(recur)
    }
    g1 <- function(u){
      (ni[i] - sum(rs*caphrCpp(X,X[i]))*exp(sum(beta_r*Zs[i,]) + u))
    }
    
    mll <- function(u){
      #returns the value inside the integral per sample id
      return(term_lik(u)*rec_lik(u)*g1(u)*dnorm(u, mean = 0, sd = psi))
    }
    
    #compute the value of the 2nd derivative with respect to the mode
    uhat[i] <- useq[which.max(mll(useq))]
    
    if(ni[i] > 90){ 
      uhat[i] <- 0 
    }
    
  }
  return(uhat)
}

full_uhatCpp2v2 <- function(theta){
  #initialize the parameters
  beta_r <- theta[1:Z_pn]
  beta_h <- theta[(Z_pn+1):(2*Z_pn)]
  hs <- theta[(2*Z_pn+1):(2*Z_pn+np)]
  rs <- theta[(2*Z_pn+np+1):(2*(Z_pn+np))]
  gamma <- theta[2*(Z_pn+np) - 1]
  psi <- theta[length(theta)]
  
  uhat <- numeric(n)
  
  for(i in 1:nrow(Zs)){
    
    term_lik <- function(u){
      #returns the likelihood of terminal event per sample id
      ((exp(sum(beta_h*Zs[i,]) + gamma*u)*sum(hs*smhrCpp2(X,X[i])))^delta[i])*
        (exp(-sum(hs*caphrCpp2(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u)))
    }
    
    rec_lik <- function(u){
      #returns the likelihood for recurrent events per sample id
      if(recur_indicators[i] == 1){
        recur_dataid <- recur_data[which(recur_data$id == i),]
        nid <- nrow(recur_dataid)
        prod_recur <- 1
        for(j in 1:nid){
          prod_recur <- prod_recur*(exp(sum(beta_r*Zs[i,]) + u))*sum(rs*smhrCpp2(recur_data[,1], recur_dataid[j,1]))
        }
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp2(X, X[i]))))*prod_recur
      } else {
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp2(X, X[i]))))
      }
      return(recur)
    }
    g1 <- function(u){
      (ni[i] - sum(rs*caphrCpp2(X,X[i]))*exp(sum(beta_r*Zs[i,]) + u))
    }
    
    mll <- function(u){
      #returns the value inside the integral per sample id
      return(term_lik(u)*rec_lik(u)*g1(u)*dnorm(u, mean = 0, sd = psi))
    }
    
    #compute the value of the 2nd derivative with respect to the mode
    uhat[i] <- useq[which.max(mll(useq))]
    
    if(ni[i] > 90){ 
      uhat[i] <- 0 
    }
    
  }
  return(uhat)
}

full_uhatCpp3 <- function(theta){
  #initialize the parameters
  beta_r <- theta[1:Z_pn]
  beta_h <- theta[(Z_pn+1):(2*Z_pn)]
  hs <- theta[(2*Z_pn+1):(2*Z_pn+np)]
  rs <- theta[(2*Z_pn+np+1):(2*(Z_pn+np))]
  gamma <- theta[2*(Z_pn+np) - 1]
  psi <- theta[length(theta)]
  
  uhat <- numeric(n)
  
  for(i in 1:nrow(Zs)){
    
    term_lik <- function(u){
      #returns the likelihood of terminal event per sample id
      ((exp(sum(beta_h*Zs[i,]) + gamma*u)*sum(hs*smhrCpp(X,X[i])))^delta[i])*
        (exp(-sum(hs*caphrCpp(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u)))
    }
    
    rec_lik <- function(u){
      #returns the likelihood for recurrent events per sample id
      if(recur_indicators[i] == 1){
        recur_dataid <- recur_data[which(recur_data$id == i),]
        nid <- nrow(recur_dataid)
        prod_recur <- 1
        for(j in 1:nid){
          prod_recur <- prod_recur*(exp(sum(beta_r*Zs[i,]) + u))*sum(rs*smhrCpp(recur_data[,1], recur_dataid[j,1]))
        }
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp(X, X[i]))))*prod_recur
      } else {
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp(X, X[i]))))
      }
      return(recur)
    }
    
    g3 <- function(u){
      (delta[i] - sum(hs*caphrCpp(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u))
    }
    
    mll <- function(u){
      #returns the value inside the integral per sample id
      return(term_lik(u)*rec_lik(u)*g3(u)*dnorm(u, mean = 0, sd = psi))
    }
    
    #compute the value of the 2nd derivative with respect to the mode
    uhat[i] <- useq[which.max(mll(useq))]
    
    if(ni[i] > 90){ 
      uhat[i] <- 0 
    }
    
  }
  return(uhat)
}

full_uhatCpp3v2 <- function(theta){
  #initialize the parameters
  beta_r <- theta[1:Z_pn]
  beta_h <- theta[(Z_pn+1):(2*Z_pn)]
  hs <- theta[(2*Z_pn+1):(2*Z_pn+np)]
  rs <- theta[(2*Z_pn+np+1):(2*(Z_pn+np))]
  gamma <- theta[2*(Z_pn+np) - 1]
  psi <- theta[length(theta)]
  
  uhat <- numeric(n)
  
  for(i in 1:nrow(Zs)){
    
    term_lik <- function(u){
      #returns the likelihood of terminal event per sample id
      ((exp(sum(beta_h*Zs[i,]) + gamma*u)*sum(hs*smhrCpp2(X,X[i])))^delta[i])*
        (exp(-sum(hs*caphrCpp2(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u)))
    }
    
    rec_lik <- function(u){
      #returns the likelihood for recurrent events per sample id
      if(recur_indicators[i] == 1){
        recur_dataid <- recur_data[which(recur_data$id == i),]
        nid <- nrow(recur_dataid)
        prod_recur <- 1
        for(j in 1:nid){
          prod_recur <- prod_recur*(exp(sum(beta_r*Zs[i,]) + u))*sum(rs*smhrCpp2(recur_data[,1], recur_dataid[j,1]))
        }
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp2(X, X[i]))))*prod_recur
      } else {
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp2(X, X[i]))))
      }
      return(recur)
    }
    
    g3 <- function(u){
      (delta[i] - sum(hs*caphrCpp2(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u))
    }
    
    mll <- function(u){
      #returns the value inside the integral per sample id
      return(term_lik(u)*rec_lik(u)*g3(u)*dnorm(u, mean = 0, sd = psi))
    }
    
    #compute the value of the 2nd derivative with respect to the mode
    uhat[i] <- useq[which.max(mll(useq))]
    
    if(ni[i] > 90){ 
      uhat[i] <- 0 
    }
    
  }
  return(uhat)
}

full_uhatCpp4 <- function(theta){
  #initialize the parameters
  beta_r <- theta[1:Z_pn]
  beta_h <- theta[(Z_pn+1):(2*Z_pn)]
  hs <- theta[(2*Z_pn+1):(2*Z_pn+np)]
  rs <- theta[(2*Z_pn+np+1):(2*(Z_pn+np))]
  gamma <- theta[2*(Z_pn+np) - 1]
  psi <- theta[length(theta)]
  
  uhat <- numeric(n)
  
  for(i in 1:nrow(Zs)){
    
    term_lik <- function(u){
      #returns the likelihood of terminal event per sample id
      ((exp(sum(beta_h*Zs[i,]) + gamma*u)*sum(hs*smhrCpp(X,X[i])))^delta[i])*
        (exp(-sum(hs*caphrCpp(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u)))
    }
    
    rec_lik <- function(u){
      #returns the likelihood for recurrent events per sample id
      if(recur_indicators[i] == 1){
        recur_dataid <- recur_data[which(recur_data$id == i),]
        nid <- nrow(recur_dataid)
        prod_recur <- 1
        for(j in 1:nid){
          prod_recur <- prod_recur*(exp(sum(beta_r*Zs[i,]) + u))*sum(rs*smhrCpp(recur_data[,1], recur_dataid[j,1]))
        }
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp(X, X[i]))))*prod_recur
      } else {
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp(X, X[i]))))
      }
      return(recur)
    }
    
    g1 <- function(u){
      ((ni[i] - sum(rs*caphrCpp(X,X[i]))*exp(sum(beta_r*Zs[i,]) + u))^2
       - (sum(rs*caphrCpp(X,X[i]))*exp(sum(beta_r*Zs[i,]) + u)))
    }
    
    mll <- function(u){
      #returns the value inside the integral per sample id
      return(term_lik(u)*rec_lik(u)*g1(u)*dnorm(u, mean = 0, sd = psi))
    }
    
    #compute the value of the 2nd derivative with respect to the mode
    uhat[i] <- useq[which.max(mll(useq))]
    
    if(ni[i] > 90){ 
      uhat[i] <- 0 
    }
    
  }
  return(uhat)
}

full_uhatCpp4v2 <- function(theta){
  #initialize the parameters
  beta_r <- theta[1:Z_pn]
  beta_h <- theta[(Z_pn+1):(2*Z_pn)]
  hs <- theta[(2*Z_pn+1):(2*Z_pn+np)]
  rs <- theta[(2*Z_pn+np+1):(2*(Z_pn+np))]
  gamma <- theta[2*(Z_pn+np) - 1]
  psi <- theta[length(theta)]
  
  uhat <- numeric(n)
  
  for(i in 1:nrow(Zs)){
    
    term_lik <- function(u){
      #returns the likelihood of terminal event per sample id
      ((exp(sum(beta_h*Zs[i,]) + gamma*u)*sum(hs*smhrCpp2(X,X[i])))^delta[i])*
        (exp(-sum(hs*caphrCpp2(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u)))
    }
    
    rec_lik <- function(u){
      #returns the likelihood for recurrent events per sample id
      if(recur_indicators[i] == 1){
        recur_dataid <- recur_data[which(recur_data$id == i),]
        nid <- nrow(recur_dataid)
        prod_recur <- 1
        for(j in 1:nid){
          prod_recur <- prod_recur*(exp(sum(beta_r*Zs[i,]) + u))*sum(rs*smhrCpp2(recur_data[,1], recur_dataid[j,1]))
        }
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp2(X, X[i]))))*prod_recur
      } else {
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp2(X, X[i]))))
      }
      return(recur)
    }
    
    g1 <- function(u){
      ((ni[i] - sum(rs*caphrCpp2(X,X[i]))*exp(sum(beta_r*Zs[i,]) + u))^2
       - (sum(rs*caphrCpp2(X,X[i]))*exp(sum(beta_r*Zs[i,]) + u)))
    }
    
    mll <- function(u){
      #returns the value inside the integral per sample id
      return(term_lik(u)*rec_lik(u)*g1(u)*dnorm(u, mean = 0, sd = psi))
    }
    
    #compute the value of the 2nd derivative with respect to the mode
    uhat[i] <- useq[which.max(mll(useq))]
    
    if(ni[i] > 90){ 
      uhat[i] <- 0 
    }
    
  }
  return(uhat)
}

full_uhatCpp5 <- function(theta){
  #initialize the parameters
  beta_r <- theta[1:Z_pn]
  beta_h <- theta[(Z_pn+1):(2*Z_pn)]
  hs <- theta[(2*Z_pn+1):(2*Z_pn+np)]
  rs <- theta[(2*Z_pn+np+1):(2*(Z_pn+np))]
  gamma <- theta[2*(Z_pn+np) - 1]
  psi <- theta[length(theta)]
  
  uhat <- numeric(n)
  
  for(i in 1:nrow(Zs)){
    
    term_lik <- function(u){
      #returns the likelihood of terminal event per sample id
      ((exp(sum(beta_h*Zs[i,]) + gamma*u)*sum(hs*smhrCpp(X,X[i])))^delta[i])*
        (exp(-sum(hs*caphrCpp(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u)))
    }
    
    rec_lik <- function(u){
      #returns the likelihood for recurrent events per sample id
      if(recur_indicators[i] == 1){
        recur_dataid <- recur_data[which(recur_data$id == i),]
        nid <- nrow(recur_dataid)
        prod_recur <- 1
        for(j in 1:nid){
          prod_recur <- prod_recur*(exp(sum(beta_r*Zs[i,]) + u))*sum(rs*smhrCpp(recur_data[,1], recur_dataid[j,1]))
        }
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp(X, X[i]))))*prod_recur
      } else {
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp(X, X[i]))))
      }
      return(recur)
    }
    
    g1 <- function(u){
      ((delta[i] - sum(hs*caphrCpp(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u))^2
       - (sum(hs*caphrCpp(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u)))
    }
    
    mll <- function(u){
      #returns the value inside the integral per sample id
      return(term_lik(u)*rec_lik(u)*g1(u)*dnorm(u, mean = 0, sd = psi))
    }
    
    #compute the value of the 2nd derivative with respect to the mode
    uhat[i] <- useq[which.max(mll(useq))]
    
    if(ni[i] > 90){ 
      uhat[i] <- 0 
    }
    
  }
  return(uhat)
}

full_uhatCpp5v2 <- function(theta){
  #initialize the parameters
  beta_r <- theta[1:Z_pn]
  beta_h <- theta[(Z_pn+1):(2*Z_pn)]
  hs <- theta[(2*Z_pn+1):(2*Z_pn+np)]
  rs <- theta[(2*Z_pn+np+1):(2*(Z_pn+np))]
  gamma <- theta[2*(Z_pn+np) - 1]
  psi <- theta[length(theta)]
  
  uhat <- numeric(n)
  
  for(i in 1:nrow(Zs)){
    
    term_lik <- function(u){
      #returns the likelihood of terminal event per sample id
      ((exp(sum(beta_h*Zs[i,]) + gamma*u)*sum(hs*smhrCpp2(X,X[i])))^delta[i])*
        (exp(-sum(hs*caphrCpp2(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u)))
    }
    
    rec_lik <- function(u){
      #returns the likelihood for recurrent events per sample id
      if(recur_indicators[i] == 1){
        recur_dataid <- recur_data[which(recur_data$id == i),]
        nid <- nrow(recur_dataid)
        prod_recur <- 1
        for(j in 1:nid){
          prod_recur <- prod_recur*(exp(sum(beta_r*Zs[i,]) + u))*sum(rs*smhrCpp2(recur_data[,1], recur_dataid[j,1]))
        }
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp2(X, X[i]))))*prod_recur
      } else {
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp2(X, X[i]))))
      }
      return(recur)
    }
    
    g1 <- function(u){
      ((delta[i] - sum(hs*caphrCpp2(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u))^2
       - (sum(hs*caphrCpp2(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u)))
    }
    
    mll <- function(u){
      #returns the value inside the integral per sample id
      return(term_lik(u)*rec_lik(u)*g1(u)*dnorm(u, mean = 0, sd = psi))
    }
    
    #compute the value of the 2nd derivative with respect to the mode
    uhat[i] <- useq[which.max(mll(useq))]
    
    if(ni[i] > 90){ 
      uhat[i] <- 0 
    }
    
  }
  return(uhat)
}


full_uhatCpp6 <- function(theta){
  #initialize the parameters
  beta_r <- theta[1:Z_pn]
  beta_h <- theta[(Z_pn+1):(2*Z_pn)]
  hs <- theta[(2*Z_pn+1):(2*Z_pn+np)]
  rs <- theta[(2*Z_pn+np+1):(2*(Z_pn+np))]
  gamma <- theta[2*(Z_pn+np) - 1]
  psi <- theta[length(theta)]
  
  uhat <- numeric(n)
  
  for(i in 1:nrow(Zs)){
    
    term_lik <- function(u){
      #returns the likelihood of terminal event per sample id
      ((exp(sum(beta_h*Zs[i,]) + gamma*u)*sum(hs*smhrCpp(X,X[i])))^delta[i])*
        (exp(-sum(hs*caphrCpp(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u)))
    }
    
    rec_lik <- function(u){
      #returns the likelihood for recurrent events per sample id
      if(recur_indicators[i] == 1){
        recur_dataid <- recur_data[which(recur_data$id == i),]
        nid <- nrow(recur_dataid)
        prod_recur <- 1
        for(j in 1:nid){
          prod_recur <- prod_recur*(exp(sum(beta_r*Zs[i,]) + u))*sum(rs*smhrCpp(recur_data[,1], recur_dataid[j,1]))
        }
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp(X, X[i]))))*prod_recur
      } else {
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp(X, X[i]))))
      }
      return(recur)
    }
    
    g1 <- function(u){
      (ni[i] - sum(rs*caphrCpp(X,X[i]))*exp(sum(beta_r*Zs[i,]) + u))
    }
    
    g2 <- function(u){
      (delta[i] - sum(hs*caphrCpp(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u))
    }
    
    mll <- function(u){
      #returns the value inside the integral per sample id
      return(term_lik(u)*rec_lik(u)*g1(u)*g2(u)*dnorm(u, mean = 0, sd = psi))
    }
    
    #compute the value of the 2nd derivative with respect to the mode
    uhat[i] <- useq[which.max(mll(useq))]
    
    if(ni[i] > 90){ 
      uhat[i] <- 0 
    }
    
  }
  return(uhat)
}

full_uhatCpp6v2 <- function(theta){
  #initialize the parameters
  beta_r <- theta[1:Z_pn]
  beta_h <- theta[(Z_pn+1):(2*Z_pn)]
  hs <- theta[(2*Z_pn+1):(2*Z_pn+np)]
  rs <- theta[(2*Z_pn+np+1):(2*(Z_pn+np))]
  gamma <- theta[2*(Z_pn+np) - 1]
  psi <- theta[length(theta)]
  
  uhat <- numeric(n)
  
  for(i in 1:nrow(Zs)){
    
    term_lik <- function(u){
      #returns the likelihood of terminal event per sample id
      ((exp(sum(beta_h*Zs[i,]) + gamma*u)*sum(hs*smhrCpp2(X,X[i])))^delta[i])*
        (exp(-sum(hs*caphrCpp2(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u)))
    }
    
    rec_lik <- function(u){
      #returns the likelihood for recurrent events per sample id
      if(recur_indicators[i] == 1){
        recur_dataid <- recur_data[which(recur_data$id == i),]
        nid <- nrow(recur_dataid)
        prod_recur <- 1
        for(j in 1:nid){
          prod_recur <- prod_recur*(exp(sum(beta_r*Zs[i,]) + u))*sum(rs*smhrCpp2(recur_data[,1], recur_dataid[j,1]))
        }
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp2(X, X[i]))))*prod_recur
      } else {
        recur <- exp(-exp(sum(beta_r*Zs[i,]) + u)*(sum(rs*caphrCpp2(X, X[i]))))
      }
      return(recur)
    }
    
    g1 <- function(u){
      (ni[i] - sum(rs*caphrCpp2(X,X[i]))*exp(sum(beta_r*Zs[i,]) + u))
    }
    
    g2 <- function(u){
      (delta[i] - sum(hs*caphrCpp2(X,X[i]))*exp(sum(beta_h*Zs[i,]) + gamma*u))
    }
    
    mll <- function(u){
      #returns the value inside the integral per sample id
      return(term_lik(u)*rec_lik(u)*g1(u)*g2(u)*dnorm(u, mean = 0, sd = psi))
    }
    
    #compute the value of the 2nd derivative with respect to the mode
    uhat[i] <- useq[which.max(mll(useq))]
    
    if(ni[i] > 90){ 
      uhat[i] <- 0 
    }
    
  }
  return(uhat)
}
