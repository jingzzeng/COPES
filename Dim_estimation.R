COPES_SAVE_d <- function(k,vX,Y,sigma){
  index = sort(Y,index.return=TRUE,decreasing = TRUE)$ix
  M = matrix(data=0,ncol=ncol(vX),nrow=ncol(vX))
  sigma_inv_sqrt = sqrtm(solve(sigma))
  for ( i in 1:k){
    if (i==1){
      Tm = vX[index[i],] %*% t(vX[index[i],])/k
    }else{
      Tm =Tm+vX[index[i],] %*% t(vX[index[i],]) /k
    }
    tau = median(eigen(sigma_inv_sqrt %*% Tm %*% sigma_inv_sqrt)$values)
    D = Tm -tau*sigma
    M = M + D %*% t(D)/k
  }
  lambda = svd(solve(sigma) %*% M)$d
  lambda = lambda+1e-5
  d_max= 5
  r = rep(NA,d_max)
  for ( i in 1:d_max){
    r[i] = lambda[i]/lambda[i+1]
  }
  return (which.max(r))
# 
#   delta = mean(lambda[1:5])
#   return(which(lambda<delta)[1]-1)
}


COPES_DR_d <- function(k,vX,Y,sigma){
  index = sort(Y,index.return=TRUE,decreasing = TRUE)$ix
  sigma_inv_sqrt = sqrtm(solve(sigma))
  
  K1 = matrix(data=0,ncol=ncol(vX),nrow=ncol(vX))
  K2_inter = matrix(data=0,ncol=ncol(vX),nrow=ncol(vX))
  K3_inter = matrix(data=0,ncol=ncol(vX),nrow=ncol(vX))
  K4_1 = 0
  K5_1 = matrix(data=0,ncol=1,nrow=ncol(vX))
  K5_2 = matrix(data=0,ncol=ncol(vX),nrow=1)
  for ( i in 1:k){
    if (i==1){
      C = vX[index[i],]/k
      Tm = vX[index[i],] %*% t(vX[index[i],])/k
    }else{
      C = C + vX[index[i],]/k
      Tm =Tm+vX[index[i],] %*% t(vX[index[i],]) /k
    }
    tau = median(eigen(sigma_inv_sqrt %*% Tm %*% sigma_inv_sqrt)$values)
    D = Tm -tau*sigma
    K1 = K1+ D%*%t(D)/(3*k)
    K2_inter = K2_inter+i/k^2*D
    K3_inter = C %*% t(C)/k
    K4_1 = K4_1+ as.numeric(t(C) %*% C)/k
    K5_1 = K5_1 + D%*% C/k
    K5_2 = K5_2 + i/k^2*t(C)
  }
  K2 = K2_inter %*% K2_inter
  K3 = K3_inter %*% K3_inter
  K4_2 = K3_inter
  K4 = K4_1*K4_2
  K5 = -K5_1 %*% K5_2
  
  K6_inter = t(K5_2)
  K6 = matrix(data=0,ncol=ncol(vX),nrow=ncol(vX))
  ## the calculation of K6 needs an additional loop
  for ( i in 1:k){
    if (i==1){
      C = vX[index[i],]/k
      Tm = vX[index[i],] %*% t(vX[index[i],])/k
    }else{
      C = C + vX[index[i],]/k
      Tm =Tm+vX[index[i],] %*% t(vX[index[i],]) /k
    }
    tau = median(eigen(sigma_inv_sqrt %*% Tm %*% sigma_inv_sqrt)$values)
    D = Tm -tau*sigma
    K6 = K6 - D %*% K6_inter %*% t(C)/k
  }
  K7 = t(K5)
  K8 = t(K6)
  M = 2*(K1+K2+K3+K4+K5+K6+K7+K8)
  lambda = svd(solve(sigma) %*% M)$d
  lambda = lambda+1e-5
  d_max= 5
  r = rep(NA,d_max)
  for ( i in 1:d_max){
    r[i] = lambda[i]/lambda[i+1]
  }
  return (which.max(r))
  # delta = mean(lambda[1:5])
  # return(which(lambda<delta)[1]-1)
}
# SAVE_d <- function(k,X,Y,sigma){
#   index = sort(Y,index.return=TRUE,decreasing = TRUE)$ix
#   M = matrix(data=0,ncol=ncol(X),nrow=ncol(X))
#   
#   for ( i in 1:k){
#     if (i==1){
#       Tm = (X[index[i],] %*% t(X[index[i],]) - sigma)/k
#     }else{
#       Tm = Tm+(X[index[i],] %*% t(X[index[i],]) -sigma)/k
#     }
#     M = M + Tm %*% t(Tm)/k
#   }
#   lambda = svd(solve(sigma) %*% M)$d
#   d_max= 5
#   r = rep(NA,d_max)
#   for ( i in 1:d_max){
#     r[i] = lambda[i]/lambda[i+1]
#   }
#   return (which.max(r))
# }