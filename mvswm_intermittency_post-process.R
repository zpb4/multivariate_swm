#setwd('z:/mv_swm/')

ar<-3

arfit<-'ols-cbias' #penalized or ols version of VAR model and residuals
sites <- 14
n<-1000

arr_full<-readRDS('data/arr_full.rds')
obs<-arr_full[1,,]
sim<-arr_full[2,,]
emp_err<-arr_full[3,,]
min_vec1<-apply(sim,2,min)
min_vec2<-apply(obs,2,function(x){x[x==0]<-NA;out<-min(x,na.rm=T);return(out)})

min_vec<-apply(rbind(min_vec1,min_vec2),2,min)
min_mat<-matrix(rep(min_vec,length(sim)),ncol=sites,byrow=T)

#date/time indices
ix<-seq(as.Date('1988-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

syn_flow_corr<-readRDS(paste('out/syn_flow_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))
#syn_flow_corr_nzero<-readRDS(paste('out/syn_flow_prinorm-ls_',arfit,'-nzero_ar',ar,'.rds',sep=''))
obs_pred_array_seq_auto<-readRDS('fit/obs_pred_array_seq_auto.rds')

#syn_flow_out_nozero<-array(NA,dim(syn_flow_corr))
syn_flow_out_samp_int<-array(NA,dim(syn_flow_corr))
syn_flow_out_int<-array(NA,dim(syn_flow_corr))

for (m in 1:n){
  syn_flow<-syn_flow_corr[m,,]
  #syn_flow_sint<-syn_flow_corr_nzero[m,,]
  rand_zero_idx<-which(syn_flow==0)
  syn_flow[rand_zero_idx]<-min_mat[rand_zero_idx]
  #syn_flow_out_nozero[m,,]<-syn_flow
  int_zero_idx<-which(obs_pred_array_seq_auto[m,,]==0)
  syn_flow_int<-syn_flow
  syn_flow_int[int_zero_idx]<-0
  #syn_flow_sint[int_zero_idx]<-0
  syn_flow_out_int[m,,]<-syn_flow_int
  #syn_flow_out_samp_int[m,,]<-syn_flow_sint
}

#saveRDS(syn_flow_out_nozero,paste('out/syn_flow_prinorm-ls_',arfit,'_ar',ar,'_nzero_obsim-min.rds',sep=''))
saveRDS(syn_flow_out_int,paste('out/syn_flow_prinorm-ls_',arfit,'_ar',ar,'_autolog_obsim-min.rds',sep=''))
#saveRDS(syn_flow_out_samp_int,paste('out/syn_flow_prinorm-ls_',arfit,'_ar',ar,'_autolog_samp-min.rds',sep=''))


rm(list=ls());gc()

####################################END###########################################