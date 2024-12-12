rm(list=ls())
packages = c("mvtnorm","evd","ggplot2","reshape2","dplyr","ICSNP","expm","doParallel","foreach","ggpubr")
for (package in packages){
  library(package,character.only = TRUE)
}
source("TIREX.R")
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
  
  mu_TIREX= apply(X,2,mean)
  sigma_TIREX =cov(X)
  mu_CP = apply(X,2,median)
  sigma_CP = tyler.shape(X,location = mu_CP)
  X_TIREX <- matrix(data=NA,ncol=p,nrow=n)
  X_CP <- matrix(data=NA,ncol=p,nrow=n)
  for ( i in 1:n){
    X_TIREX[i,] =  as.numeric(X[i,]-mu_TIREX)
    X_CP[i,] = as.numeric( X[i,]-mu_CP) / as.numeric(t(X[i,]-mu_CP) %*% solve(sigma_CP) %*% (X[i,]-mu_CP))^(1/2)
  }
  
  mse_TIREX_SIR <- rep(NA,length(k_seq))
  mse_COPES_SIR <- rep(NA,length(k_seq))
  mse_TIREX_SAVE<- rep(NA,length(k_seq))
  mse_COPES_SAVE <- rep(NA,length(k_seq))
  mse_COPES_DR <- rep(NA,length(k_seq))

  for ( i in 1:length(k_seq)){
    k = k_seq[i]
    eta = TIREX_SIR(k,X_TIREX,Y,sigma_TIREX)
    mse_TIREX_SIR[i] = norm(eta %*% solve(t(eta)%*%eta) %*%t(eta)-P,"F")^2/(2*d)
    eta = TIREX_SAVE(k,d,X_TIREX,Y,sigma_TIREX)
    mse_TIREX_SAVE[i] = norm(eta %*% solve(t(eta)%*%eta) %*%t(eta)-P,"F")^2/(2*d)
    
    eta = COPES_SIR(k,X_CP,Y,sigma_CP)
    mse_COPES_SIR[i] = norm(eta %*% solve(t(eta)%*%eta) %*%t(eta)-P,"F")^2/(2*d)
    
    eta = COPES_SAVE(k,d,X_CP,Y,sigma_CP)
    mse_COPES_SAVE[i] = norm(eta %*% solve(t(eta)%*%eta) %*%t(eta)-P,"F")^2/(2*d)
    
    eta = COPES_DR(k,d,X_CP,Y,sigma_CP)
    mse_COPES_DR[i] = norm(eta %*% solve(t(eta)%*%eta) %*%t(eta)-P,"F")^2/(2*d)
    
  }
  df = data.frame(k=k_seq,TIREX_SIR=mse_TIREX_SIR,TIREX_SAVE = mse_TIREX_SAVE,
                  COPES_SIR = mse_COPES_SIR, COPES_SAVE = mse_COPES_SAVE, COPES_DR=mse_COPES_DR)
  return (df)
}



output_model <- function(v,model){
  print(paste0("v = ",v, " model = ", model))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  num_rep = 200
  df=foreach(s = 1:num_rep, 
             .combine=rbind,
             .packages=packages,.export = c("one_round")) %dopar% ({
               source("TIREX.R")
               source("COPES.R")
               source("models.R")
               one_round(v=v,model=model)})
  
  stopCluster(cl)
  df = df %>% group_by(k) %>% summarise(TIREX1 = mean(TIREX_SIR),
                                        TIREX2 = mean(TIREX_SAVE),
                                        COPES_SIR=mean(COPES_SIR),
                                        COPES_SAVE=mean(COPES_SAVE),
                                        COPES_DR=mean(COPES_DR))
  names(df)  =c("k","TIREX1","TIREX2","COPES-SIR","COPES-SAVE","COPES-DR")
  df = melt(df,id.vars = 'k')
  p <- ggplot(df,aes(x=k,y=value,col=variable,linetype=variable))+geom_line()+ylab("MSE")+
   ylim(0,1)+scale_linetype_manual(values=c(2,2,1,1,1))+
    scale_colour_manual(values=c("blue","red","blue","red","black"))+ theme_bw()+theme_grey(base_size = 10)
  p
}

for (model in c("A","B","C","D")){
  p1 <- output_model(v=2,model=model)+theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())+ggtitle(expression(paste(nu, "= ", 2)))
  p2 <- output_model(v=3,model=model)+theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())+ggtitle(expression(paste(nu, "= ", 3)))
  p3 <- output_model(v=5,model=model)+theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())+ggtitle(expression(paste(nu, "= ", 5)))
  p4 <- output_model(v=0,model=model)+theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())+ggtitle("Multivariate Normal")
  ggarrange(p1,p2,p3,p4,common.legend = T,legend="bottom")
  ggsave(filename=paste0("model=" ,model,".pdf"),width = 8, height = 5, dpi = 300, units = "in")
}


for (v in c("u3","u2","u5")){
  p1 <- output_model(v=v,model="A")+theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())+ggtitle("Model (A)")
  p2 <- output_model(v=v,model="B")+theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())+ggtitle("Model (B)")
  p3 <- output_model(v=v,model="C")+theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())+ggtitle("Model (C)")
  p4 <- output_model(v=v,model="D")+theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())+ggtitle("Model (D)")
  ggarrange(p1,p2,p3,p4,common.legend = T,legend="bottom")
  ggsave(filename=paste0("Non_EC_df=" ,v,".pdf"),width = 8, height = 5, dpi = 300, units = "in")
}




