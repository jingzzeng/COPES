data_generate <- function(n,p,v,model){
    if (v=="u3"){
    X = matrix(0,ncol=p,nrow=n)
    for ( i in 1:p){
      X[,i] = rt(n,df=3)
    }
  }else if (v=="u5"){
    X = matrix(0,ncol=p,nrow=n)
    for ( i in 1:p){
      X[,i] = rt(n,df=5)
    }
  }else if (v=="u2"){
      X = matrix(0,ncol=p,nrow=n)
      for ( i in 1:p){
        X[,i] = rt(n,df=2)
      }
    }else{
    X = rmvt(n,sigma = diag(p),df=v)
    }
  
  
  B = rbinom(n, 1, 0.5)
  xi = rgpd(n,loc=0,scale=1,shape = 1)
  xi2 = rgpd(n,loc=0,scale=1/2,shape = 1/2)
  error = rexp(n)
  if (model=="A"){
    Y = B*sin(X[,1]/2)*xi+(1-B)*sin(X[,2]/2)*xi2
    e = rep(0,p);e[1]=1
    P =  e %*% solve(t(e)%*%e) %*%t(e)
    d=1
  }else if (model=="B"){
    Y = B*(4*X[,1]+(X[,2]+0.25)^2)*xi+(1-B)/((X[,3]+0.25)^2+0.25)*xi2
    e1 = rep(0,p);e1[1]=1
    e2 = rep(0,p); e2[2]=1
    e=matrix(c(e1,e2),ncol=2)
    P =  e %*% solve(t(e)%*%e) %*%t(e)
    d=2
  }else if (model=="C"){
    U = runif(n)
    Z = (1-U)^(-1)-1
    Z2 = (1-U)^(-1/2)-1
    Y = Z*(cos(X[,1]+pi/4)+(U<=0.95)*cos(X[,2]+pi/4))
    e = rep(0,p);e[1]=1
    P =  e %*% solve(t(e)%*%e) %*%t(e)
    d=1
  }else if (model=="D"){
    U = runif(n)
    Z = (1-U)^(-1)-1
    Z2 = (1-U)^(-1/2)-1
    Y = sin(X[,1])/((X[,2]+0.25)^2+0.25)*Z+(U<=0.95)*Z2*sin(X[,3])
    e1 = rep(0,p);e1[1]=1
    e2 = rep(0,p); e2[2]=1
    e=matrix(c(e1,e2),ncol=2)
    P =  e %*% solve(t(e)%*%e) %*%t(e)
    d=2
  }else{
    print("Error: the parameter model is not properly chosen")
  }
  list(X=X,Y=Y,P=P,d=d)
}