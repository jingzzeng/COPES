rm(list=ls())
packages = c("mvtnorm","evd","ggplot2","reshape2","dplyr","ICSNP","expm","doParallel","foreach")
for (package in packages){
  library(package,character.only = TRUE)
}
source("Dim_estimation.R")
source("models.R")
set.seed(2022)




one_round <- function(v,model){
  df = data.frame()
  n_seq = c(5000,10000,20000)
  p = 10
  for ( i in 1:length(n_seq)){
    n = n_seq[i]
    k_ratio_seq <- c(
                     0.6+log(2)/log(n),
                     2/3,
                     0.7,
                     0.75+log(1/1.5)/log(n))
   
    data = data_generate(n,p,v,model)
    X = data$X
    Y = data$Y
    P = data$P
    d = data$d
    mu_CP = apply(X,2,median)
    sigma_CP = tyler.shape(X,location = mu_CP)
    X_CP <- matrix(data=NA,ncol=p,nrow=n)
    for ( s in 1:n){
      X_CP[s,] = as.numeric( X[s,]-mu_CP) / as.numeric(t(X[s,]-mu_CP) %*% solve(sigma_CP) %*% (X[s,]-mu_CP))^(1/2)
    }
    d_save <- rep(NA,length(k_ratio_seq))
    d_dr <- rep(NA,length(k_ratio_seq))
    for ( j in 1:length(k_ratio_seq)){
      k = floor(n^k_ratio_seq[j])
      d_save[j] = COPES_SAVE_d(k,X_CP,Y,sigma_CP)
      d_dr[j] = COPES_DR_d(k,X_CP,Y,sigma_CP)
    }
    sub_df = data.frame(v = v, n = n ,k_ratio= 1:length(k_ratio_seq),d_save = d_save, d_dr  = d_dr,d=d)
    df  = rbind(df,sub_df)
  }
  return (df)
}



output_model <- function(v,model){
  print(paste0("v = ",v, " model = ", model))
  cl <- makeCluster(10)
  registerDoParallel(cl)
  num_rep = 200
  df=foreach(s = 1:num_rep, 
             .combine=rbind,
             .packages=packages,.export = c("one_round")) %dopar% ({
               source("Dim_estimation.R")
               source("models.R")
               one_round(v=v, model = model)})
  stopCluster(cl)
 
  df$save = (df$d_save==df$d)
  df$dr = (df$d_dr==df$d)
  df = df %>% group_by(v,k_ratio,n) %>% summarise(COPES_SAVE = mean(save), COPES_DR = mean(dr))
  return (df)
}

for (model in c(5,6,7,8)){
  print("model=")
  print(model)
  df_3 = output_model(v=3,model=model)
  print(df_3,n=50)
  # df_2 = output_model(v=2,model=model)
  # print(df_2,n=50)
  # df_0 = output_model(v=0,model=model)
  # print(df_0,n=50)
  # 
  # df_5 = output_model(v=5,model=model)
  # df = rbind(df_0,df_2,df_3,df_5)
  # print(df,n=200)
  # write.csv(df,paste0("M ", model, "_d.csv"),row.names = FALSE)
}



