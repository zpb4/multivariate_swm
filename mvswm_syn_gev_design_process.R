#date
#setwd('z:/mv_swm/')
library(extRemes)
library(doParallel)
library(foreach)
library(abind)

print(paste('start',Sys.time()))

ix<-seq(as.Date('1988-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

sged_mv<-'prinorm-ls_ols-cbias'
sged_ind<-'prinorm-ls_ind_pen-cbias'
ar<-3

#multivariate data files
syn_flow_sged<-readRDS(paste('out/syn_flow_',sged_mv,'_ar',ar,'.rds',sep=''))

#independent data files
syn_flow_sged_ind<-readRDS(paste('out/syn_flow_',sged_ind,'_ar',ar,'.rds',sep=''))

arr_full<-readRDS('data/arr_full.rds')

obs<-arr_full[1,,]
sim<-arr_full[2,,]

syn_flow_sged[syn_flow_sged>10*max(obs)]<-max(obs)
syn_flow_sged_ind[syn_flow_sged_ind>10*max(obs)]<-max(obs)

#parameters
cmb<-2:4
n<-1000
yrs<-88:113


comb_sged_var<-apply(syn_flow_sged[1:n,,cmb],c(1,2),sum)
comb_sged_ind<-apply(syn_flow_sged_ind[1:n,,cmb],c(1,2),sum)

comb_sged_var_annmax<-array(NA,c(n,length(yrs)))
comb_sged_ind_annmax<-array(NA,c(n,length(yrs)))

calfews_sged_var_annmax<-array(NA,c(n,length(yrs),14))
calfews_sged_ind_annmax<-array(NA,c(n,length(yrs),14))

for(i in 1:length(yrs)){
  idx<-which(ix2$year==yrs[i])
  comb_sged_var_annmax[,i]<-apply(comb_sged_var[,idx],1,max)
  comb_sged_ind_annmax[,i]<-apply(comb_sged_ind[,idx],1,max)
  
  calfews_sged_var_annmax[,i,]<-apply(syn_flow_sged[1:n,idx,],c(1,3),max)
  calfews_sged_ind_annmax[,i,]<-apply(syn_flow_sged_ind[1:n,idx,],c(1,3),max)
}

saveRDS(calfews_sged_var_annmax,paste('out/annmax_sged_mv_',sged_mv,'_ar',ar,'_all-sites.rds',sep=''))
saveRDS(calfews_sged_ind_annmax,paste('out/annmax_sged_ind_',sged_ind,'_ar',ar,'_all-sites.rds',sep=''))

saveRDS(comb_sged_var_annmax,paste('out/annmax_sged_mv_',sged_mv,'_ar',ar,'_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))
saveRDS(comb_sged_ind_annmax,paste('out/annmax_sged_ind_',sged_ind,'_ar',ar,'_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))

dsgn_events<-c(10,50,100)

parallel::detectCores()
  n.cores <- parallel::detectCores()-1
  my.cluster<-parallel::makeCluster(n.cores,type='PSOCK')
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()

for(d in 1:length(dsgn_events)){
  dsgn<-dsgn_events[d]
  #Design events
  #obs_est<-readRDS(paste('data/obs_design-est_',dsgn,'-yr-event_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))
  dsgn_prob<-1 - (1/dsgn)

  gev_fun<-function(x){
    fit<-fevd(x,type='GEV')
    dsn<-qevd(dsgn_prob,loc=fit$results$par[1],scale=fit$results$par[2],shape=fit$results$par[3],type='GEV')
    return(dsn)}

  #SGED MV
  out_sged_mv<-foreach(i = 1:n,.combine='rbind',.packages = 'extRemes',.inorder = F) %dopar% {
    dsgn_sged_var<-apply(calfews_sged_var_annmax[i,,],2,gev_fun)
    return(dsgn_sged_var)
  }
  saveRDS(out_sged_mv,paste('out/syn_design-est_sged_mv_',sged_mv,'_ar',ar,'_',dsgn,'-yr-event_all-sites.rds',sep=''))

  comb_out_sged_mv<-foreach(i = 1:n,.combine='c',.packages = 'extRemes',.inorder = F) %dopar% {
    gev_fit<-fevd(comb_sged_var_annmax[i,],type='GEV')
    comb_dsgn_sged_var<-qevd(dsgn_prob,loc=gev_fit$results$par[1],scale =gev_fit$results$par[2],shape = gev_fit$results$par[3],type='GEV')
    return(comb_dsgn_sged_var)
  }
  saveRDS(comb_out_sged_mv,paste('out/syn_design-est_sged_mv_',sged_mv,'_ar',ar,'_',dsgn,'-yr-event_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))

  #SGED IND
  out_sged_ind<-foreach(i = 1:n,.combine='rbind',.packages = 'extRemes',.inorder = F,.errorhandling = 'remove') %dopar% {
    dsgn_sged_ind<-apply(calfews_sged_ind_annmax[i,,],2,gev_fun)
    return(dsgn_sged_ind)
  }
  saveRDS(out_sged_ind,paste('out/syn_design-est_sged_ind_',sged_ind,'_ar',ar,'_',dsgn,'-yr-event_all-sites.rds',sep=''))

  comb_out_sged_ind<-foreach(i = 1:n,.combine='c',.packages = 'extRemes',.inorder = F,.errorhandling = 'remove') %dopar% {
    gev_fit<-fevd(comb_sged_ind_annmax[i,],type='GEV')
    comb_dsgn_sged_ind<-qevd(dsgn_prob,loc=gev_fit$results$par[1],scale =gev_fit$results$par[2],shape = gev_fit$results$par[3],type='GEV')
    return(comb_dsgn_sged_ind)
  }
  saveRDS(comb_out_sged_ind,paste('out/syn_design-est_sged_ind_',sged_ind,'_ar',ar,'_',dsgn,'-yr-event_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))
}

print(paste('end',Sys.time()))

rm(list=ls());gc()

################################################END#########################################################