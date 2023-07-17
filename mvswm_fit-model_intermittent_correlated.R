#Fit correlated auto-logistic regression model to data
setwd('h:/mv_swm')

#Define no. of postprocessing 0/1 matrices to produce
n<-1000

print(paste('start',Sys.time()))

#fixed parameters
ix<-seq(as.Date('1988-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)
locs_full<-c('SHA', 'ORO','YRS','FOL','MKM','NHG','NML','TLG','MRC','MIL','PNF','TRM','SCC','ISB')

#read in data
arr_full<-readRDS('data/arr_full.rds')

#1) Fit autologistic model to data and simulate 0/1 matrices for post-processing

#1a. Define Sequential Autologistic Prediction Function

pred_auto<-function(x,coeffs){
  p<-c()
  p[1]<-exp(coeffs[1] +  x[1,] %*% coeffs[3:length(coeffs)]) / (1 + exp(coeffs[1] + x[1,] %*% coeffs[3:length(coeffs)]))
  prd<-c()
  prd[1]<-rbinom(1,1,p[1])
  for(i in 2:length(x[,1])){
    p[i] <- exp(coeffs[1] + coeffs[2]*prd[i-1] + x[i,] %*% coeffs[3:length(coeffs)]) / (1 + exp(coeffs[1] + coeffs[2]*prd[i-1] + x[i,] %*% coeffs[3:length(coeffs)]))
    prd[i]<-rbinom(1,1,p[i])
  }
  return(prd)
}

#1b. Fit Autologistic model sequentially

alogit_fit_seq_auto<-vector('list',14)

for(i in 1:length(locs_full)){
  idx<-1:length(locs_full)
  idx2<-idx[-c(i)]
  if(i == 1){
    alogit_fit_seq_auto[[i]]<-glm(arr_full[5,,i]~cbind(c(sample(0:1,1),arr_full[5,1:(length(arr_full[5,,i])-1),i]),arr_full[2,,]),family='binomial')} else
    {alogit_fit_seq_auto[[i]]<-glm(arr_full[5,,i]~cbind(c(sample(0:1,1),arr_full[5,1:(length(arr_full[5,,i])-1),i]),arr_full[2,,],arr_full[5,,1:(i-1)]),family='binomial')}
}

saveRDS(alogit_fit_seq_auto,'fit/alogit_fit_seq_auto.rds')

#1c. Create 'n' matrices of 0s and 1s for postprocessing Log Ratio Model simulations

obs_pred_array_seq_auto<-array(0,c(n,dim(arr_full[5,,])))

for(m in 1:n){
  for(i in 1:length(locs_full)){
    if(i == 1){
      obs_pred_array_seq_auto[m,,i]<-pred_auto(arr_full[2,,],matrix(alogit_fit_seq_auto[[i]]$coefficients))} else
      {obs_pred_array_seq_auto[m,,i]<-pred_auto(cbind(arr_full[2,,],obs_pred_array_seq_auto[m,,1:(i-1)]),matrix(alogit_fit_seq_auto[[i]]$coefficients))}
  }
}    

saveRDS(obs_pred_array_seq_auto,'fit/obs_pred_array_seq_auto.rds')

print(paste('end',Sys.time()))

#clean environment
rm(list=ls());gc()

###############################################END################################################################