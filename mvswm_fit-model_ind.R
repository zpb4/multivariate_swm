#Fit Independent AR and SGED models to Differenced errors
setwd('z:/mv_swm')

#Req packages
library(fGarch)
library(bigtime)
library(MASS)
library(robts)

#fixed parameters
ix<-seq(as.Date('1988-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)
locs_full<-c('SHA', 'ORO','YRS','FOL','MKM','NHG','NML','TLG','MRC','MIL','PNF','TRM','SCC','ISB')
ar<-3

#1) Read in data
arr_full<-readRDS('data/arr_full.rds')
nresids<-readRDS('fit/nresids_prinorm-ls-cbias_full.rds')

#2) AR model
#2a. AR model to create uncorrelated residuals

uc_resid <- vector('list',12)
uc_resid_ols <- vector('list',12)
uc_resid_var_ar <- vector('list',12)
ar_coefs<-vector('list',12)
ar_coefs_sub<-vector('list',length(locs_full))
ar_coefs_ols<-vector('list',12)
ar_coefs_sub_ols<-vector('list',length(locs_full))

#Calculate uncorrelated matrices
for(i in 1:12){
  rresid_mat<-nresids[[i]]
  ar_mat<-array(NA,dim(rresid_mat))
  ar_mat_ols<-array(NA,dim(rresid_mat))
  ar_coefs[[i]]<-ar_coefs_sub
  ar_coefs_ols[[i]]<-ar_coefs_sub_ols
  
  for(k in 1:length(locs_full)){
    #if(ar>1){
      #rvar_fit<-sparseVAR(rresid_mat[,k],p=ar,VARpen='HLag',h=1,selection='cv',check_std = F)
      #coefs<-rvar_fit$Phihat
    #}
    #if(ar==1){
      #rvar_fit<-rlm(rresid_mat[,k]~c(rep(0,ar),rresid_mat[1:(length(rresid_mat[,k])-ar),k]),maxit=100)
      #coefs<-rvar_fit$coefficients[2]
    #}
    #mat<-cbind(rresid_mat[,k],rep(0,length(rresid_mat[,k])))
    #mat<-cbind(rresid_mat[,k],rnorm(length(rresid_mat[,k])))
    #b<-varxfit(mat,p=ar,robust=T,constant = F)
    #ar_coef<-b$Bcoef[1,seq(1,(ar*2),2)]
    
    ar_rob<-arrob.regression(rresid_mat[,k],aic=F,order.max = ar,intercept = F)
    #intcpt<-spvar_fit$phi0hat
    #ar_fit<-arima(rresid_mat[,k],order=c(ar,0,0),method='CSS',include.mean = T)
    ar_fit<-arima(rresid_mat[,k],order=c(ar,0,0),fixed=ar_rob$ar,include.mean = F)
    ar_fit_ols<-arima(rresid_mat[,k],order=c(ar,0,0),method='CSS',include.mean = F)
    ar_mat[,k]<-ar_fit$residuals
    ar_mat_ols[,k]<-ar_fit_ols$residuals
    ar_coefs[[i]][[k]]<-ar_fit
    ar_coefs_ols[[i]][[k]]<-ar_fit_ols
  }
  uc_resid[[i]]<-ar_mat
  uc_resid_ols[[i]]<-ar_mat_ols
}

saveRDS(uc_resid,paste('fit/uc_resid_prinorm-ls_ind_pen-cbias_ar',ar,'.rds',sep=''))
saveRDS(ar_coefs,paste('fit/ar_coefs_prinorm-ls_ind_pen-cbias_ar',ar,'.rds',sep=''))
saveRDS(uc_resid_ols,paste('fit/uc_resid_prinorm-ls_ind_ols-cbias_ar',ar,'.rds',sep=''))
saveRDS(ar_coefs_ols,paste('fit/ar_coefs_prinorm-ls_ind_ols-cbias_ar',ar,'.rds',sep=''))

#3) SGED Model (monthly) on independent AR residuals
#3a. Penalized AR SGED

at_lst<-vector('list',12)
gl_par_arr<-array(NA,c(12,length(locs_full),4))
uc_resid<-readRDS(paste('fit/uc_resid_prinorm-ls_ind_pen-cbias_ar',ar,'.rds',sep=''))

for(i in 1:12){
  uc_sim<-uc_resid[[i]]
  seas<-which(ix2$mon==(i-1))
  at_arr<-array(0,c(dim(uc_sim)[1],length(locs_full)))
  for(j in 1:length(locs_full)){
    gl_mle<-sgedFit(uc_sim[,j])
    at<-uc_sim[,j]
    gl_par_arr[i,j,]<-gl_mle$par
    at_arr[,j]<-at
  }
  at_lst[[i]]<-at_arr
}

saveRDS(gl_par_arr,paste('fit/gl_par_arr_prinorm-ls_ind_pen-cbias_ar',ar,'.rds',sep=''))
saveRDS(at_lst,paste('fit/at_lst_prinorm-ls_ind_pen-cbias_ar',ar,'.rds',sep=''))

#3b. OLS AR SGED
at_lst<-vector('list',12)
gl_par_arr<-array(NA,c(12,length(locs_full),4))
uc_resid<-readRDS(paste('fit/uc_resid_prinorm-ls_ind_ols-cbias_ar',ar,'.rds',sep=''))

for(i in 1:12){
  uc_sim<-uc_resid[[i]]
  seas<-which(ix2$mon==(i-1))
  at_arr<-array(0,c(dim(uc_sim)[1],length(locs_full)))
  for(j in 1:length(locs_full)){
    gl_mle<-sgedFit(uc_sim[,j])
    at<-uc_sim[,j]
    gl_par_arr[i,j,]<-gl_mle$par
    at_arr[,j]<-at
  }
  at_lst[[i]]<-at_arr
}

saveRDS(gl_par_arr,paste('fit/gl_par_arr_prinorm-ls_ind_ols-cbias_ar',ar,'.rds',sep=''))
saveRDS(at_lst,paste('fit/at_lst_prinorm-ls_ind_ols-cbias_ar',ar,'.rds',sep=''))


rm(list=ls());gc()

###########################################END################################