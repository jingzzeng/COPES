COPES_SIR <- function(k,vX,Y,sigma){
  index = sort(Y,index.return=TRUE,decreasing = TRUE)$ix
  M = matrix(data=0,ncol=ncol(vX),nrow=ncol(vX))
  for ( i in 1:k){
    if (i==1){
      C = vX[index[i],]/k
    }else{
      C = C + vX[index[i],]/k
    }
    M = M +  C %*% t(C)/k
  }
  eta = svd(solve(sigma) %*% M)$u[,1]
  return (eta)
}

COPES_SAVE <- function(k,d,vX,Y,sigma){
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
  eta = svd(solve(sigma) %*% M)$u[,1:d]
  return (eta)
}




COPES_DR <- function(k,d,vX,Y,sigma){
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
  eta = svd(solve(sigma) %*% M)$u[,1:d]
  return (eta)
}
