library(evd)
library(mvtnorm)
library(extRemes)
library(mnormt)
#library(copula)

#date
#setwd('z:/mv_swm/')

print(paste('start',Sys.time()))

ix<-seq(as.Date('1988-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)
locs_full<-c('SHA', 'ORO','YRS','FOL','MKM','NHG','NML','TLG','MRC','MIL','PNF','TRM','SCC','ISB')

sged_mv1<-'prinorm-ls_ols-cbias'
sged_mv2<-'prinorm-ls_ols-cbias'
sged_ind<-'prinorm-ls_ind_pen-cbias'
ar<-3
sites<-14
n<-1000
dsgn<-10
yrs<-88:113
plt_gev<-F
sel_sites<-c('ORO','YRS','FOL')
sel_idx<-c()
for(i in 1:length(sel_sites)){sel_idx[i]<-which(locs_full==sel_sites[i])}
sel_thresh<-c(700,700,500)

arr_full<-readRDS('data/arr_full.rds')

obs<-arr_full[1,,sel_idx]
sim<-arr_full[2,,sel_idx]

#SWM data files
syn_flow_corr<-readRDS(paste('out/syn_flow_',sged_mv1,'_ar',ar,'.rds',sep=''))
#syn_flow_corr_int<-readRDS(paste('out/syn_flow_',sged_mv2,'_ar',ar,'_autolog_samp-min.rds',sep=''))
syn_flow_corr_int<-readRDS(paste('out/syn_flow_',sged_mv2,'_ar',ar,'_autolog_obsim-min.rds',sep=''))
syn_flow_ind<-readRDS(paste('out/syn_flow_',sged_ind,'_ar',ar,'.rds',sep=''))

#functions
gev_fun<-function(x){
  fit<-fevd(x,type='GEV')
  dsn<-qevd(dsgn_prob,loc=fit$results$par[1],scale=fit$results$par[2],shape=fit$results$par[3],type='GEV')
  return(dsn)}

gev_param<-function(x){
  fit<-fevd(x,type='GEV')
  #dsn<-qevd(dsgn_prob,loc=fit$results$par[1],scale=fit$results$par[2],shape=fit$results$par[3],type='GEV')
  return(fit$results$par)}

mov_sum<-function(x,dy){
  idx<-(dy-1)/2
  out<-c()
  for(i in 1:length(x)){
    low_idx<-max(1,(i-idx))
    hi_idx<-min(length(x),(i+idx))
    out[i]<-sum(x[low_idx:hi_idx])
  }
  return(out)
}

mov_avg<-function(x,dy){
  idx<-(dy-1)/2
  out<-c()
  for(i in 1:length(x)){
    low_idx<-max(1,(i-idx))
    hi_idx<-min(length(x),(i+idx))
    out[i]<-mean(x[low_idx:hi_idx])
  }
  return(out)
}


#calculate design event on 7-day flows
wdow<-7
#dsgn_prob<-1 - (1/dsgn)
dsgn_prob<-(1/dsgn)

obs_avg_annmin<-array(NA,c(length(yrs),length(sel_sites)))
obs_avg<-apply(obs,2,function(x){out<-mov_avg(x,wdow);return(out)})

sim_avg_annmin<-array(NA,c(length(yrs),length(sel_sites)))
sim_avg<-apply(sim,2,function(x){out<-mov_avg(x,wdow);return(out)})

for(i in 1:length(yrs)){
  idx<-which(ix2$year==yrs[i])
  obs_avg_annmin[i,]<-apply(obs_avg[idx,],2,min)
  sim_avg_annmin[i,]<-apply(sim_avg[idx,],2,min)
}

obs_avg_param<-apply(obs_avg_annmin,2,gev_param)
sim_avg_param<-apply(sim_avg_annmin,2,gev_param)

obs_avg_est<-apply(obs_avg_annmin,2,gev_fun)
sim_avg_est<-apply(sim_avg_annmin,2,gev_fun)

if(plt_gev==T){
for(i in 1:length(sel_sites)){
  par(mfrow=c(1,2))
  hist(obs_avg_annmin[,i],freq = F)
  lines(0:max(obs_avg_annmin[,i]),devd(0:max(obs_avg_annmin[,i]),loc=obs_avg_param[1,i],scale=obs_avg_param[2,i],shape = obs_avg_param[3,i],type='GEV'),col='red')
  
  hist(sim_avg_annmin[,i],freq = F)
  lines(0:max(sim_avg_annmin[,i]),devd(0:max(sim_avg_annmin[,i]),loc=sim_avg_param[1,i],scale=sim_avg_param[2,i],shape = sim_avg_param[3,i],type='GEV'),col='red')
}
}

#Convert the sim thresholds to u's for each basin 

u_simref_avg_event<-c()
for(i in 1:length(sel_sites)){
  u_simref_avg_event[i]<-pgev(obs_avg_est[i], loc=sim_avg_param[1,i], scale=sim_avg_param[2,i], sim_avg_param[3,i], lower.tail = TRUE)
}

u_obsref_avg_event<-c()
for(i in 1:length(sel_sites)){
  u_obsref_avg_event[i]<-pgev(sim_avg_est[i], loc=obs_avg_param[1,i], scale=obs_avg_param[2,i], obs_avg_param[3,i], lower.tail = TRUE)
}

#Now takes these u's and push them through a standard normal to get Z scores 
z_obs_avg_event<-qnorm(rep(dsgn_prob,length(sel_sites)))
z_simref_avg_event<-qnorm(u_simref_avg_event)
z_obsref_avg_event<-qnorm(u_obsref_avg_event)

#Now we want to fit the copula on our historical flows. First create u's
#obs
u_obs_avg<-array(NA,dim(obs_avg_annmin))
for(i in 1:length(sel_sites)){
  u_obs_avg[,i]<-pgev(obs_avg_annmin[,i], loc=obs_avg_param[1,i], scale=obs_avg_param[2,i], obs_avg_param[3,i], lower.tail = TRUE)
}

#sim
u_sim_avg<-array(NA,dim(sim_avg_annmin))
for(i in 1:length(sel_sites)){
  u_sim_avg[,i]<-pgev(sim_avg_annmin[,i], loc=sim_avg_param[1,i], scale=sim_avg_param[2,i], sim_avg_param[3,i], lower.tail = TRUE)
}

#Then create Z-scores
#obs
#avg
z_obs_avg<-array(NA,dim(obs_avg_annmin))
for(i in 1:length(sel_sites)){
  z_obs_avg[,i]<-qnorm(u_obs_avg[,i])
}

#avg
z_sim_avg<-array(NA,dim(sim_avg_annmin))
for(i in 1:length(sel_sites)){
  z_sim_avg[,i]<-qnorm(u_sim_avg[,i])
}

#Now calculate the historical probability of exceeding the threshold in all basins
#exceedance_obs_avg_annmin=pnorm(z_obs_avg_event[1])-pnorm(z_obs_avg_event[2])-pnorm(z_obs_avg_event[3])+pmnorm(z_obs_avg_event[1:2],mean=rep(0,2),varcov = cor(z_obs_avg[,1:2]))+pmnorm(z_obs_avg_event[c(1,3)],mean=rep(0,2),varcov = cor(z_obs_avg[,c(1,3)]))+pmnorm(z_obs_avg_event[2:3],mean=rep(0,2),varcov = cor(z_obs_avg[,2:3]))-pmnorm(z_obs_avg_event,mean=rep(0,3),varcov = cor(z_obs_avg))
#exceedance_simref_avg_annmin=pnorm(z_obs_avg_event[1])-pnorm(z_obs_avg_event[2])-pnorm(z_obs_avg_event[3])+pmnorm(z_obs_avg_event[1:2],mean=rep(0,2),varcov = cor(z_sim_avg[,1:2]))+pmnorm(z_obs_avg_event[c(1,3)],mean=rep(0,2),varcov = cor(z_sim_avg[,c(1,3)]))+pmnorm(z_obs_avg_event[2:3],mean=rep(0,2),varcov = cor(z_sim_avg[,2:3]))-pmnorm(z_obs_avg_event,mean=rep(0,3),varcov = cor(z_sim_avg))
#exceedance_obsref_avg_annmin=pnorm(z_obsref_avg_event[1])-pnorm(z_obsref_avg_event[2])-pnorm(z_obsref_avg_event[3])+pmnorm(z_obsref_avg_event[1:2],mean=rep(0,2),varcov = cor(z_obs_avg[,1:2]))+pmnorm(z_obsref_avg_event[c(1,3)],mean=rep(0,2),varcov = cor(z_obs_avg[,c(1,3)]))+pmnorm(z_obsref_avg_event[2:3],mean=rep(0,2),varcov = cor(z_obs_avg[,2:3]))-pmnorm(z_obsref_avg_event,mean=rep(0,3),varcov = cor(z_obs_avg))


#saveRDS(c(exceedance_obs_annmax,exceedance_obs_sum_annmax,exceedance_simref_annmax,exceedance_simref_sum_annmax,exceedance_obsref_annmax,exceedance_obsref_sum_annmax),paste('out/trivar-exceedance_obs-sim_',sel_sites[1],'-',sel_sites[2],'-',sel_sites[3],'_',dsgn,'-yr-event.rds',sep=''))

#------------------------------------------------------------------------------------
#calculate empirical occurrences in SWM samples
#daily
noexceed_fun_trivar<-function(x,thresh){
  out<-0
  if(x[1]<thresh[1]&x[2]<thresh[2]&x[3]<thresh[3]){out<-1}
  return(out)
}

#----------------------------------------------------------------------------
#avg empirical occurrences
#obs
obs_avg_occ<-apply(obs_avg,1,function(x){out<-noexceed_fun_trivar(x,sel_thresh);return(out)})
emp_exc_obs_avg<-sum(obs_avg_occ)/length(yrs)

#sim
sim_avg_occ<-apply(sim_avg,1,function(x){out<-noexceed_fun_trivar(x,sel_thresh);return(out)})
emp_exc_sim_avg<-sum(sim_avg_occ)/length(yrs)

#corr
swm_corr<-syn_flow_corr[,,sel_idx]
swm_corr_avg<-array(NA,dim(swm_corr))
for(m in 1:n){
  swm_corr_avg[m,,]<-apply(swm_corr[m,,],2,function(x){out<-mov_avg(x,wdow);return(out)})
}
occurrence_emp_corr_avg<-apply(swm_corr_avg,c(1,2),function(x){out<-noexceed_fun_trivar(x,sel_thresh);return(out)})
emp_exc_corr_avg<-sum(occurrence_emp_corr_avg)/(length(yrs)*n)

#intermittent
swm_corr_int<-syn_flow_corr_int[,,sel_idx]
swm_corr_int_avg<-array(NA,dim(swm_corr))
for(m in 1:n){
  swm_corr_int_avg[m,,]<-apply(swm_corr_int[m,,],2,function(x){out<-mov_avg(x,wdow);return(out)})
}
occurrence_emp_corr_int_avg<-apply(swm_corr_int_avg,c(1,2),function(x){out<-noexceed_fun_trivar(x,sel_thresh);return(out)})
emp_exc_corr_int_avg<-sum(occurrence_emp_corr_int_avg)/(length(yrs)*n)

#independent
swm_ind<-syn_flow_ind[,,sel_idx]
swm_ind_avg<-array(NA,dim(swm_corr))
for(m in 1:n){
  swm_ind_avg[m,,]<-apply(swm_ind[m,,],2,function(x){out<-mov_sum(x,wdow);return(out)})
}
occurrence_emp_ind_avg<-apply(swm_ind_avg,c(1,2),function(x){out<-noexceed_fun_trivar(x,sel_thresh);return(out)})
emp_exc_ind_avg<-sum(occurrence_emp_ind_avg)/(length(yrs)*n)

saveRDS(c(emp_exc_obs_avg,emp_exc_sim_avg,emp_exc_corr_avg,emp_exc_corr_int_avg,emp_exc_ind_avg),paste('out/trivar-loflow_',wdow,'d-avg_',sel_sites[1],'-',sel_sites[2],'-',sel_sites[3],'_',sged_mv2,'_ar',ar,'_chinook-thresh.rds',sep=''))

#-------------------------------------------------------------------------------------------------------
#obs 7Q10
#avg empirical occurrences
#obs
obs_avg_occ<-apply(obs_avg,1,function(x){out<-noexceed_fun_trivar(x,obs_avg_est);return(out)})
emp_exc_obs_avg<-sum(obs_avg_occ)/length(yrs)

#sim
sim_avg_occ<-apply(sim_avg,1,function(x){out<-noexceed_fun_trivar(x,obs_avg_est);return(out)})
emp_exc_sim_avg<-sum(sim_avg_occ)/length(yrs)

#corr
swm_corr<-syn_flow_corr[,,sel_idx]
swm_corr_avg<-array(NA,dim(swm_corr))
for(m in 1:n){
  swm_corr_avg[m,,]<-apply(swm_corr[m,,],2,function(x){out<-mov_avg(x,wdow);return(out)})
}
occurrence_emp_corr_avg<-apply(swm_corr_avg,c(1,2),function(x){out<-noexceed_fun_trivar(x,obs_avg_est);return(out)})
emp_exc_corr_avg<-sum(occurrence_emp_corr_avg)/(length(yrs)*n)

#intermittent
swm_corr_int<-syn_flow_corr_int[,,sel_idx]
swm_corr_int_avg<-array(NA,dim(swm_corr))
for(m in 1:n){
  swm_corr_int_avg[m,,]<-apply(swm_corr_int[m,,],2,function(x){out<-mov_avg(x,wdow);return(out)})
}
occurrence_emp_corr_int_avg<-apply(swm_corr_int_avg,c(1,2),function(x){out<-noexceed_fun_trivar(x,obs_avg_est);return(out)})
emp_exc_corr_int_avg<-sum(occurrence_emp_corr_int_avg)/(length(yrs)*n)

#independent
swm_ind<-syn_flow_ind[,,sel_idx]
swm_ind_avg<-array(NA,dim(swm_corr))
for(m in 1:n){
  swm_ind_avg[m,,]<-apply(swm_ind[m,,],2,function(x){out<-mov_sum(x,wdow);return(out)})
}
occurrence_emp_ind_avg<-apply(swm_ind_avg,c(1,2),function(x){out<-noexceed_fun_trivar(x,obs_avg_est);return(out)})
emp_exc_ind_avg<-sum(occurrence_emp_ind_avg)/(length(yrs)*n)

saveRDS(c(emp_exc_obs_avg,emp_exc_sim_avg,emp_exc_corr_avg,emp_exc_corr_int_avg,emp_exc_ind_avg),paste('out/trivar-loflow_',wdow,'d-avg_',sel_sites[1],'-',sel_sites[2],'-',sel_sites[3],'_',sged_mv2,'_ar',ar,'_7Q10.rds',sep=''))


rm(list=ls());gc()

###############################################END###############################################