rm(list=ls())
packages = c("mvtnorm","evd","ggplot2","reshape2","dplyr","ICSNP","expm","doParallel","foreach","ggpubr")
for (package in packages){
  library(package,character.only = TRUE)
}
source("COPES.R")
source("models.R")
set.seed(2022)


one_round <- function(v,model){
  n = 5000
  k_seq <- seq(50,2000,50)
  p = 10
  data = data_generate(n,p,v,model)
  X = data$X
  Y = data$Y
  P = data$P
  d = data$d
  
  mse_COPES_SIR <- rep(NA,length(k_seq))
  mse_COPES_SAVE <- rep(NA,length(k_seq))
  mse_COPES_DR <- rep(NA,length(k_seq))

  for ( i in 1:length(k_seq)){
    k = k_seq[i]
    
    eta = COPES_SIR(k,X_CP,Y,sigma_CP)
    mse_COPES_SIR[i] = norm(eta %*% solve(t(eta)%*%eta) %*%t(eta)-P,"F")^2/(2*d)
    
    eta = COPES_SAVE(k,d,X_CP,Y,sigma_CP)
    mse_COPES_SAVE[i] = norm(eta %*% solve(t(eta)%*%eta) %*%t(eta)-P,"F")^2/(2*d)
    
    eta = COPES_DR(k,d,X_CP,Y,sigma_CP)
    mse_COPES_DR[i] = norm(eta %*% solve(t(eta)%*%eta) %*%t(eta)-P,"F")^2/(2*d)
    
  }
  df = data.frame(k=k_seq,
                  COPES_SIR = mse_COPES_SIR, COPES_SAVE = mse_COPES_SAVE, COPES_DR=mse_COPES_DR)
  return (df)
}
