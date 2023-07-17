#Fit VAR and SGED models to Differenced errors
setwd('z:/mv_swm')

print(paste('start',Sys.time()))

#Req packages
library(fGarch)
library(BigVAR)
library(rmgarch)

#fixed parameters
ix<-seq(as.Date('1988-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)
locs_full<-c('SHA', 'ORO','YRS','FOL','MKM','NHG','NML','TLG','MRC','MIL','PNF','TRM','SCC','ISB')
ar<-3
sites<-14

#
hyb_loess_fit<-function(obs,pred){
  ls_fit<-loess(pred~obs,span=0.75,degree=2,family='gaussian',control=loess.control(surface='interpolate',statistics = 'none',trace.hat = 'approximate'))
  srt_obs<-sort(obs,index.return=T)
  ls_pred<-predict(ls_fit,srt_obs$x)
  diff_ls_pred<-diff(ls_pred)
  infl_idx<-length(diff_ls_pred[diff_ls_pred>=0])
  slope<-diff_ls_pred[infl_idx]/(srt_obs$x[infl_idx+1]-srt_obs$x[infl_idx])
  y<-ls_pred[infl_idx+1]
  intcpt<-y-slope*srt_obs$x[infl_idx+1]
  cutoff<-srt_obs$x[infl_idx+1]
  return(list(ls_fit,cutoff,slope,intcpt))
}

#this function uses the fitted parameters from hyb_loess_fit to predict conditional mean for new data
hyb_loess_out<-function(loess_mod,cutoff,slope,intcpt,x){
  xout<-sort(x,index.return=T)
  split_idx<-length(which(xout$x<=cutoff))
  if(split_idx==length(x)){
    yout<-predict(loess_mod,x)}
  if(split_idx<length(x)){
    y_loess<-predict(loess_mod,xout$x[1:split_idx])
    y_lm<-lin_mod(intcpt,slope,xout$x[(split_idx+1):length(x)])
    yout<-rep(0,length(x))
    yout[xout$ix[1:split_idx]]<-y_loess
    yout[xout$ix[(split_idx+1):length(x)]]<-y_lm}
  return(yout)
}

#linear model helper functions
lin_mod<-function(intcpt,slope,x){
  y = intcpt + slope * x
  return(y)
}

lin_fit<-function(pars,x,y){
  pred = pars[1] + pars[2] * x
  err = sum((y - pred)^2)
  return(err)
}

#read in data
arr_full<-readRDS('data/arr_full.rds')
obs<-arr_full[1,,]
sim<-arr_full[2,,]



#1a. Fit conditional bias model
sim_cbias<-array(NA,dim(sim)) #array to retain fitted conditional expectation (cexp) estimates
rresids_full<-vector('list',12)

for(i in 1:12){
  seas<-which(ix2$mo==(i-1))
  temp_mat<-array(NA,c(length(seas),sites))
  for(j in 1:sites){
    hy_fit<-hyb_loess_fit(sim[seas,j],obs[seas,j])
    cbias<-hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],sim[seas,j])
    sim_cbias[seas,j]<-cbias
    temp_mat[,j]<-obs[seas,j]-cbias
  }
  rresids_full[[i]]<-temp_mat
}

saveRDS(sim_cbias,'fit/sim_cbias.rds')
saveRDS(rresids_full,'fit/rresids-cbias_full.rds')

#1b. Fit heteroscedastic LOESS normalization model between simulations and errors
lbs<-c(0,0) #lower bounds for intcpt, slope
ubs<-c(1000,10) #upper bounds for intcpt, slope
sts<-c(1,1) #starting parameters for intcpt, slope

norm_fit<-vector('list',12)
norm_fit_sub<-vector('list',length(locs_full))
gl_par_arr<-array(NA,c(12,length(locs_full),4))
sd_arr<-array(NA,c(12,length(locs_full)))
nresids<-vector('list',12)
d_idx<-c(1:14,14)

for(i in 1:12){
  seas<-which(ix2$mon==(i-1))
  obmat<-sim_cbias[seas,]
  rresid_mat<-rresids_full[[i]]
  nresid_mat<-array(NA,dim(rresid_mat))
  norm_fit[[i]]<-norm_fit_sub
  for(k in 1:length(locs_full)){
    ob<-obmat[,d_idx[k]]
    res<-rresid_mat[,d_idx[k]]
    res_sdev<-sqrt(res^2)
    sds<-mean(res_sdev)
    sd_arr[i,k]<-sds
    
    ls_fit<-loess(res_sdev~ob,span=0.75,degree = 2,family='gaussian',control=loess.control(surface='direct'))
    
    norm_vec<-predict(ls_fit,obmat[,d_idx[k]])
    norm_vec[norm_vec<0]<-sds
    
    norm_resids<-rresid_mat[,d_idx[k]] / norm_vec
    
    norm_fit[[i]][[k]]<-ls_fit
    nresid_mat[,k]<-norm_resids
    
  }
  nresids[[i]]<-nresid_mat
}

saveRDS(norm_fit,'fit/norm_fit_prinorm-ls-cbias_full.rds')
saveRDS(nresids,'fit/nresids_prinorm-ls-cbias_full.rds')
saveRDS(sd_arr,'fit/sd_arr_prinorm-ls-cbias_full.rds')

#2) VAR model (monthly) to decorrelate errors

#2a. LASSO Penalized VAR model (BigVAR)
#Calculate uncorrelated matrices
uc_resid <- vector('list',12)
var_coefs<-array(NA,c(12,length(locs_full),(length(locs_full)*ar+1)))

for(i in 1:12){
  rresid_mat<-nresids[[i]]
  m1 = constructModel(rresid_mat, p = ar, struct = "Basic", gran = c(50, 10),IC = F,h=1,rolling_oos = T,
  verbose = F, VARX = list(), separate_lambdas = T,model.controls=list(intercept = F,MN = F))
  m1_res = cv.BigVAR(m1)
  coef<-m1_res@betaPred
  
  resids<-rbind(matrix(0,nrow=ar,ncol=sites),m1_res@resids)
  var_coefs[i,,]<-m1_res@betaPred

  uc_resid[[i]]<-resids
}

saveRDS(uc_resid,paste('fit/uc_resid_prinorm-ls_pen-bvar-cbias_ar',ar,'.rds',sep=''))
saveRDS(var_coefs,paste('fit/var_coefs_prinorm-ls_pen-bvar-cbias_ar',ar,'.rds',sep=''))

#Calculate uncorrelated matrices
uc_resid <- vector('list',12)
var_coefs<-array(NA,c(12,length(locs_full),(length(locs_full)*ar+1)))

for(i in 1:12){
  rresid_mat<-nresids[[i]]
  m1 = constructModel(rresid_mat, p = ar, struct = "Basic", gran = c(50, 10),IC = F,h=1,rolling_oos = T,
                      verbose = F, VARX = list(), separate_lambdas = T,model.controls=list(intercept = T,MN = F))
  m1_res = cv.BigVAR(m1)
  coef<-m1_res@betaPred
  
  resids<-rbind(matrix(0,nrow=ar,ncol=sites),m1_res@resids)
  var_coefs[i,,]<-m1_res@betaPred
  
  uc_resid[[i]]<-resids
}

saveRDS(uc_resid,paste('fit/uc_resid_prinorm-ls_pen-bvar-mn-cbias_ar',ar,'.rds',sep=''))
saveRDS(var_coefs,paste('fit/var_coefs_prinorm-ls_pen-bvar-mn-cbias_ar',ar,'.rds',sep=''))

#VARX model
uc_resid <- vector('list',12)
var_coefs<-array(NA,c(12,length(locs_full),(length(locs_full)*ar+1)))

for(i in 1:12){
  rresid_mat<-nresids[[i]]
  var_fit<-varxfit(rresid_mat,ar,constant=T,robust = T)
  var_resids<-var_fit$xresiduals
  var_resids<-rbind(matrix(0,nrow=ar,ncol=sites),var_resids)

  uc_resid[[i]]<-var_resids
  var_coefs[i,,]<-var_fit$Bcoef
}

saveRDS(uc_resid,paste('fit/uc_resid_prinorm-ls_pen-mn-cbias_ar',ar,'.rds',sep=''))
saveRDS(var_coefs,paste('fit/var_coefs_prinorm-ls_pen-mn-cbias_ar',ar,'.rds',sep=''))

#no mean
uc_resid <- vector('list',12)
var_coefs<-array(NA,c(12,length(locs_full),(length(locs_full)*ar)))

for(i in 1:12){
  rresid_mat<-nresids[[i]]
  var_fit<-varxfit(rresid_mat,ar,constant=F,robust = T)
  var_resids<-var_fit$xresiduals
  var_resids<-rbind(matrix(0,nrow=ar,ncol=sites),var_resids)
  
  uc_resid[[i]]<-var_resids
  var_coefs[i,,]<-var_fit$Bcoef
}

saveRDS(uc_resid,paste('fit/uc_resid_prinorm-ls_pen-cbias_ar',ar,'.rds',sep=''))
saveRDS(var_coefs,paste('fit/var_coefs_prinorm-ls_pen-cbias_ar',ar,'.rds',sep=''))

#2b. OLS VAR model 
uc_resid <- vector('list',12)
var_coefs<-array(NA,c(12,length(locs_full),(length(locs_full)*ar+1)))

for(i in 1:12){
  rresid_mat<-nresids[[i]]
  var_fit<-varxfit(rresid_mat,ar,constant=T,robust = F)
  var_resids<-var_fit$xresiduals
  var_resids<-rbind(matrix(0,nrow=ar,ncol=sites),var_resids)
  
  uc_resid[[i]]<-var_resids
  var_coefs[i,,]<-var_fit$Bcoef

}

saveRDS(uc_resid,paste('fit/uc_resid_prinorm-ls_ols-mn-cbias_ar',ar,'.rds',sep=''))
saveRDS(var_coefs,paste('fit/var_coefs_prinorm-ls_ols-mn-cbias_ar',ar,'.rds',sep=''))

#no-mean
uc_resid <- vector('list',12)
var_coefs<-array(NA,c(12,length(locs_full),(length(locs_full)*ar)))

for(i in 1:12){
  rresid_mat<-nresids[[i]]
  var_fit<-varxfit(rresid_mat,ar,constant=F,robust = F)
  var_resids<-var_fit$xresiduals
  var_resids<-rbind(matrix(0,nrow=ar,ncol=sites),var_resids)
  
  uc_resid[[i]]<-var_resids
  var_coefs[i,,]<-var_fit$Bcoef
  
}

saveRDS(uc_resid,paste('fit/uc_resid_prinorm-ls_ols-cbias_ar',ar,'.rds',sep=''))
saveRDS(var_coefs,paste('fit/var_coefs_prinorm-ls_ols-cbias_ar',ar,'.rds',sep=''))


#3) SGED Model (monthly) for VAR residuals
#3a. Penalized VAR
#BigVAR
at_lst<-vector('list',12)
uc_resid<-readRDS(paste('fit/uc_resid_prinorm-ls_pen-bvar-cbias_ar',ar,'.rds',sep=''))

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

saveRDS(gl_par_arr,paste('fit/gl_par_arr_prinorm-ls_pen-bvar-cbias_ar',ar,'.rds',sep=''))
saveRDS(at_lst,paste('fit/at_lst_prinorm-ls_pen-bvar-cbias_ar',ar,'.rds',sep=''))

#BigVAR+men
at_lst<-vector('list',12)
uc_resid<-readRDS(paste('fit/uc_resid_prinorm-ls_pen-bvar-mn-cbias_ar',ar,'.rds',sep=''))

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

saveRDS(gl_par_arr,paste('fit/gl_par_arr_prinorm-ls_pen-bvar-mn-cbias_ar',ar,'.rds',sep=''))
saveRDS(at_lst,paste('fit/at_lst_prinorm-ls_pen-bvar-mn-cbias_ar',ar,'.rds',sep=''))

#VARX
at_lst<-vector('list',12)
uc_resid<-readRDS(paste('fit/uc_resid_prinorm-ls_pen-mn-cbias_ar',ar,'.rds',sep=''))

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

saveRDS(gl_par_arr,paste('fit/gl_par_arr_prinorm-ls_pen-mn-cbias_ar',ar,'.rds',sep=''))
saveRDS(at_lst,paste('fit/at_lst_prinorm-ls_pen-mn-cbias_ar',ar,'.rds',sep=''))

#no-mean
at_lst<-vector('list',12)
uc_resid<-readRDS(paste('fit/uc_resid_prinorm-ls_pen-cbias_ar',ar,'.rds',sep=''))

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

saveRDS(gl_par_arr,paste('fit/gl_par_arr_prinorm-ls_pen-cbias_ar',ar,'.rds',sep=''))
saveRDS(at_lst,paste('fit/at_lst_prinorm-ls_pen-cbias_ar',ar,'.rds',sep=''))

#3b. OLS VAR model
at_lst<-vector('list',12)
uc_resid<-readRDS(paste('fit/uc_resid_prinorm-ls_ols-mn-cbias_ar',ar,'.rds',sep=''))

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

saveRDS(gl_par_arr,paste('fit/gl_par_arr_prinorm-ls_ols-mn-cbias_ar',ar,'.rds',sep=''))
saveRDS(at_lst,paste('fit/at_lst_prinorm-ls_ols-mn-cbias_ar',ar,'.rds',sep=''))

#no-mean
at_lst<-vector('list',12)
uc_resid<-readRDS(paste('fit/uc_resid_prinorm-ls_ols-cbias_ar',ar,'.rds',sep=''))

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

saveRDS(gl_par_arr,paste('fit/gl_par_arr_prinorm-ls_ols-cbias_ar',ar,'.rds',sep=''))
saveRDS(at_lst,paste('fit/at_lst_prinorm-ls_ols-cbias_ar',ar,'.rds',sep=''))

print(paste('end',Sys.time()))

rm(list=ls());gc()

###########################################END################################