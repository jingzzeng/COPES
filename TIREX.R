TIREX_SIR <- function(k,X,Y,sigma){
  index = sort(Y,index.return=TRUE,decreasing = TRUE)$ix
  M = matrix(data=0,ncol=ncol(X),nrow=ncol(X))
  for ( i in 1:k){
    if (i==1){
      C = X[index[i],]/k
    }else{
      C = C + X[index[i],]/k
    }
    M = M +  C %*% t(C)/k
  }
  eta = svd(solve(sigma) %*% M)$u[,1]
  return (eta)
}
TIREX_SAVE <- function(k,d,X,Y,sigma){
  index = sort(Y,index.return=TRUE,decreasing = TRUE)$ix
  M = matrix(data=0,ncol=ncol(X),nrow=ncol(X))
  
  for ( i in 1:k){
    if (i==1){
      Tm = (X[index[i],] %*% t(X[index[i],]) - sigma)/k
    }else{
      Tm = Tm+(X[index[i],] %*% t(X[index[i],]) -sigma)/k
    }
    M = M + Tm %*% t(Tm)/k
  }
 eta = svd(solve(sigma) %*% M)$u[,1:d]
  return (eta)
}

