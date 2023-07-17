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

sged_mv1<-'prinorm-ls_pen-cbias'
sged_mv2<-'prinorm-ls_pen-cbias'
sged_ind<-'prinorm-ls_ind_ols-cbias'
ar<-3
sites<-14
n<-1000
#dsgn<-100
yrs<-88:113
plt_gev<-F
sel_sites<-c('ORO','YRS','FOL')
sel_idx<-c()
for(i in 1:length(sel_sites)){sel_idx[i]<-which(locs_full==sel_sites[i])}

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
  dsn<-qevd(dsgn_prob,loc=fit$results$par[1],scale=fit$results$par[2],shape=fit$results$par[3],type='GEV')
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

exceed_fun_trivar<-function(x,thresh){
    out<-0
    if(x[1]>thresh[1]&x[2]>thresh[2]&x[3]>thresh[3]){out<-1}
    return(out)
  }

#dsgn_events<-c(10,15,20,25,30,50,100)
dsgn_events<-100

for(d in 1:length(dsgn_events)){
  dsgn<-dsgn_events[d]
  #calculate design event on daily flows
  dsgn_prob<-1 - (1/dsgn)

  obs_annmax<-array(NA,c(length(yrs),length(sel_sites)))
  sim_annmax<-array(NA,c(length(yrs),length(sel_sites)))

  for(i in 1:length(yrs)){
    idx<-which(ix2$year==yrs[i])
    obs_annmax[i,]<-apply(obs[idx,],2,max)
    sim_annmax[i,]<-apply(sim[idx,],2,max)
  }

  obs_param<-apply(obs_annmax,2,gev_param)
  sim_param<-apply(sim_annmax,2,gev_param)
  
  obs_est<-apply(obs_annmax,2,gev_fun)
  sim_est<-apply(sim_annmax,2,gev_fun)
  
  print(obs_est)
  
  if(dsgn==100){
    obsum<-apply(obs,1,sum)
    obs_est<-obs[which.max(obsum),]
    print(ix[which.max(obsum)])
  }
  print(obs_est)
  
  obs_dsgn_est<-c()
  for(i in 1:length(sel_sites)){
    obs_dsgn_est[i]<-1/(1-pevd(obs_est[i],loc=obs_param[1,i],scale=obs_param[2,i],shape = obs_param[3,i],type='GEV'))
  }
  
  print(obs_dsgn_est)
    
  if(plt_gev==T){
    for(i in 1:length(sel_sites)){
      par(mfrow=c(1,2))
      hist(obs_annmax[,i],freq = F)
      lines(1:max(obs_annmax[,i]),devd(1:max(obs_annmax[,i]),loc=obs_param[1,i],scale=obs_param[2,i],shape = obs_param[3,i],type='GEV'),col='red')

      hist(sim_annmax[,i],freq = F)
      lines(1:max(sim_annmax[,i]),devd(1:max(sim_annmax[,i]),loc=sim_param[1,i],scale=sim_param[2,i],shape = sim_param[3,i],type='GEV'),col='red')
    }
  }

  #calculate design event on X-day flows
  wdow<-3

  obs_sum_annmax<-array(NA,c(length(yrs),length(sel_sites)))
  obs_sum<-apply(obs,2,function(x){out<-mov_sum(x,wdow);return(out)})

  sim_sum_annmax<-array(NA,c(length(yrs),length(sel_sites)))
  sim_sum<-apply(sim,2,function(x){out<-mov_sum(x,wdow);return(out)})

  for(i in 1:length(yrs)){
    idx<-which(ix2$year==yrs[i])
    obs_sum_annmax[i,]<-apply(obs_sum[idx,],2,max)
    sim_sum_annmax[i,]<-apply(sim_sum[idx,],2,max)
  }

  obs_sum_param<-apply(obs_sum_annmax,2,gev_param)
  sim_sum_param<-apply(sim_sum_annmax,2,gev_param)

  obs_sum_est<-apply(obs_sum_annmax,2,gev_fun)
  sim_sum_est<-apply(sim_sum_annmax,2,gev_fun)
  
  print(obs_sum_est)
  
  if(dsgn==100){
    obss_sum<-apply(obs_sum,1,sum)
    obs_sum_est<-obs_sum[which.max(obss_sum),]
    print(ix[which.max(obss_sum)])
  }
  print(obs_sum_est)
  
  obs_sum_dsgn_est<-c()
  for(i in 1:length(sel_sites)){
    obs_sum_dsgn_est[i]<-1/(1-pevd(obs_sum_est[i],loc=obs_sum_param[1,i],scale=obs_sum_param[2,i],shape = obs_sum_param[3,i],type='GEV'))
  }
  
  print(obs_sum_dsgn_est)

  if(plt_gev==T){
    for(i in 1:length(sel_sites)){
      par(mfrow=c(1,2))
      hist(obs_sum_annmax[,i],freq = F)
      lines(1:max(obs_sum_annmax[,i]),devd(1:max(obs_sum_annmax[,i]),loc=obs_sum_param[1,i],scale=obs_sum_param[2,i],shape = obs_sum_param[3,i],type='GEV'),col='red')
  
      hist(sim_sum_annmax[,i],freq = F)
      lines(1:max(sim_sum_annmax[,i]),devd(1:max(sim_sum_annmax[,i]),loc=sim_sum_param[1,i],scale=sim_sum_param[2,i],shape = sim_sum_param[3,i],type='GEV'),col='red')
    }
  }

  #Convert the sim thresholds to u's for each basin 
  u_simref_event<-c()
  for(i in 1:length(sel_sites)){
    u_simref_event[i]<-pgev(obs_est[i], loc=sim_param[1,i], scale=sim_param[2,i], sim_param[3,i], lower.tail = TRUE)
  }

  u_simref_sum_event<-c()
  for(i in 1:length(sel_sites)){
    u_simref_sum_event[i]<-pgev(obs_sum_est[i], loc=sim_sum_param[1,i], scale=sim_sum_param[2,i], sim_sum_param[3,i], lower.tail = TRUE)
  }

  u_obsref_event<-c()
  for(i in 1:length(sel_sites)){
    u_obsref_event[i]<-pgev(sim_est[i], loc=obs_param[1,i], scale=obs_param[2,i], obs_param[3,i], lower.tail = TRUE)
  }

  u_obsref_sum_event<-c()
  for(i in 1:length(sel_sites)){
    u_obsref_sum_event[i]<-pgev(sim_sum_est[i], loc=obs_sum_param[1,i], scale=obs_sum_param[2,i], obs_sum_param[3,i], lower.tail = TRUE)
  }

  #Now takes these u's and push them through a standard normal to get Z scores 
  z_obs_event<-qnorm(rep(dsgn_prob,length(sel_sites)))
  z_obs_sum_event<-qnorm(rep(dsgn_prob,length(sel_sites)))

  z_simref_event<-qnorm(u_simref_event)
  z_simref_sum_event<-qnorm(u_simref_sum_event)

  z_obsref_event<-qnorm(u_obsref_event)
  z_obsref_sum_event<-qnorm(u_obsref_sum_event)

  #Now we want to fit the copula on our historical flows. First create u's
  #obs
  u_obs<-array(NA,dim(obs_annmax))
  for(i in 1:length(sel_sites)){
    u_obs[,i]<-pgev(obs_annmax[,i], loc=obs_param[1,i], scale=obs_param[2,i], obs_param[3,i], lower.tail = TRUE)
  }

  u_obs_sum<-array(NA,dim(obs_sum_annmax))
  for(i in 1:length(sel_sites)){
    u_obs_sum[,i]<-pgev(obs_sum_annmax[,i], loc=obs_sum_param[1,i], scale=obs_sum_param[2,i], obs_sum_param[3,i], lower.tail = TRUE)
  }

  #sim
  u_sim<-array(NA,dim(sim_annmax))
  for(i in 1:length(sel_sites)){
    u_sim[,i]<-pgev(sim_annmax[,i], loc=sim_param[1,i], scale=sim_param[2,i], sim_param[3,i], lower.tail = TRUE)
  }

  u_sim_sum<-array(NA,dim(sim_sum_annmax))
  for(i in 1:length(sel_sites)){
    u_sim_sum[,i]<-pgev(sim_sum_annmax[,i], loc=sim_sum_param[1,i], scale=sim_sum_param[2,i], sim_sum_param[3,i], lower.tail = TRUE)
  }

  #Then create Z-scores
  #obs
  #daily
  z_obs<-array(NA,dim(obs_annmax))
  for(i in 1:length(sel_sites)){
    z_obs[,i]<-qnorm(u_obs[,i])
  }

  #summed
  z_obs_sum<-array(NA,dim(obs_sum_annmax))
  for(i in 1:length(sel_sites)){
    z_obs_sum[,i]<-qnorm(u_obs_sum[,i])
  }

  #sim
  #daily
  z_sim<-array(NA,dim(sim_annmax))
  for(i in 1:length(sel_sites)){
    z_sim[,i]<-qnorm(u_sim[,i])
  }

  #summed
  z_sim_sum<-array(NA,dim(sim_sum_annmax))
  for(i in 1:length(sel_sites)){
    z_sim_sum[,i]<-qnorm(u_sim_sum[,i])
  }

  #Now calculate the historical probability of exceeding the threshold in all basins
  exceedance_obs_annmax=1-pnorm(z_obs_event[1])-pnorm(z_obs_event[2])-pnorm(z_obs_event[3])+pmnorm(z_obs_event[1:2],mean=rep(0,2),varcov = cor(z_obs[,1:2]))+pmnorm(z_obs_event[c(1,3)],mean=rep(0,2),
    varcov = cor(z_obs[,c(1,3)]))+pmnorm(z_obs_event[2:3],mean=rep(0,2),varcov = cor(z_obs[,2:3]))-pmnorm(z_obs_event,mean=rep(0,3),varcov = cor(z_obs))
  exceedance_obs_sum_annmax=1-pnorm(z_obs_sum_event[1])-pnorm(z_obs_sum_event[2])-pnorm(z_obs_sum_event[3])+pmnorm(z_obs_sum_event[1:2],mean=rep(0,2),varcov = cor(z_obs_sum[,1:2]))+
    pmnorm(z_obs_sum_event[c(1,3)],mean=rep(0,2),varcov = cor(z_obs_sum[,c(1,3)]))+pmnorm(z_obs_sum_event[2:3],mean=rep(0,2),varcov = cor(z_obs_sum[,2:3]))-pmnorm(z_obs_sum_event,mean=rep(0,3),
    varcov = cor(z_obs_sum))

  exceedance_simref_annmax=1-pnorm(z_obs_event[1])-pnorm(z_obs_event[2])-pnorm(z_obs_event[3])+pmnorm(z_obs_event[1:2],mean=rep(0,2),varcov = cor(z_sim[,1:2]))+pmnorm(z_obs_event[c(1,3)],
    mean=rep(0,2),varcov = cor(z_sim[,c(1,3)]))+pmnorm(z_obs_event[2:3],mean=rep(0,2),varcov = cor(z_sim[,2:3]))-pmnorm(z_obs_event,mean=rep(0,3),varcov = cor(z_sim))
  exceedance_simref_sum_annmax=1-pnorm(z_obs_sum_event[1])-pnorm(z_obs_sum_event[2])-pnorm(z_obs_sum_event[3])+pmnorm(z_obs_sum_event[1:2],mean=rep(0,2),varcov = cor(z_sim_sum[,1:2]))+
    pmnorm(z_obs_sum_event[c(1,3)],mean=rep(0,2),varcov = cor(z_sim_sum[,c(1,3)]))+pmnorm(z_obs_sum_event[2:3],mean=rep(0,2),varcov = cor(z_sim_sum[,2:3]))-
    pmnorm(z_obs_sum_event,mean=rep(0,3),varcov = cor(z_sim_sum))

  exceedance_obsref_annmax=1-pnorm(z_obsref_event[1])-pnorm(z_obsref_event[2])-pnorm(z_obsref_event[3])+pmnorm(z_obsref_event[1:2],mean=rep(0,2),varcov = cor(z_obs[,1:2]))+
    pmnorm(z_obsref_event[c(1,3)],mean=rep(0,2),varcov = cor(z_obs[,c(1,3)]))+pmnorm(z_obsref_event[2:3],mean=rep(0,2),varcov = cor(z_obs[,2:3]))-pmnorm(z_obsref_event,mean=rep(0,3),varcov = cor(z_obs))
  exceedance_obsref_sum_annmax=1-pnorm(z_obsref_sum_event[1])-pnorm(z_obsref_sum_event[2])-pnorm(z_obsref_sum_event[3])+pmnorm(z_obsref_sum_event[1:2],mean=rep(0,2),varcov = cor(z_obs_sum[,1:2]))+
    pmnorm(z_obsref_sum_event[c(1,3)],mean=rep(0,2),varcov = cor(z_obs_sum[,c(1,3)]))+pmnorm(z_obsref_sum_event[2:3],mean=rep(0,2),varcov = cor(z_obs_sum[,2:3]))-
    pmnorm(z_obsref_sum_event,mean=rep(0,3),varcov = cor(z_obs_sum))
  
  occurrence_emp_sim<-apply(sim,1,function(x){out<-exceed_fun_trivar(x,obs_est);return(out)})
  emp_exc_sim<-sum(occurrence_emp_sim)/(length(yrs))
  
  occurrence_emp_sim_sum<-apply(sim_sum,1,function(x){out<-exceed_fun_trivar(x,obs_est);return(out)})
  emp_exc_sim_sum<-sum(occurrence_emp_sim_sum)/(length(yrs))

  1/exceedance_obs_annmax
  1/exceedance_obs_sum_annmax

  1/exceedance_simref_annmax
  1/exceedance_simref_sum_annmax

  1/exceedance_obsref_annmax
  1/exceedance_obsref_sum_annmax

  saveRDS(c(exceedance_obs_annmax,exceedance_obs_sum_annmax,exceedance_simref_annmax,exceedance_simref_sum_annmax,exceedance_obsref_annmax,exceedance_obsref_sum_annmax,emp_exc_sim,emp_exc_sim_sum),
    paste('out/trivar-exceedance_obs-sim_',sel_sites[1],'-',sel_sites[2],'-',sel_sites[3],'_',dsgn,'-yr-event.rds',sep=''))

#------------------------------------------------------------------------------------
  #calculate empirical occurrences in SWM samples
  #daily

  #corr
  swm_corr<-syn_flow_corr[,,sel_idx]
  occurrence_emp_corr<-apply(swm_corr,c(1,2),function(x){out<-exceed_fun_trivar(x,obs_est);return(out)})
  emp_exc_corr<-sum(occurrence_emp_corr)/(length(yrs)*n)

  #intermittent
  swm_corr_int<-syn_flow_corr_int[,,sel_idx]
  occurrence_emp_corr_int<-apply(swm_corr_int,c(1,2),function(x){out<-exceed_fun_trivar(x,obs_est);return(out)})
  emp_exc_corr_int<-sum(occurrence_emp_corr_int)/(length(yrs)*n)

  #independent
  swm_ind<-syn_flow_ind[,,sel_idx]
  occurrence_emp_ind<-apply(swm_ind,c(1,2),function(x){out<-exceed_fun_trivar(x,obs_est);return(out)})
  emp_exc_ind<-sum(occurrence_emp_ind)/(length(yrs)*n)

  saveRDS(c(emp_exc_corr,emp_exc_corr_int,emp_exc_ind),paste('out/trivar-exceedance_',sel_sites[1],'-',sel_sites[2],'-',sel_sites[3],'_',sged_mv2,'_ar',ar,'_',dsgn,'-yr-event.rds',sep=''))

#----------------------------------------------------------------------------
  #summed
  #corr
  swm_corr_sum<-array(NA,dim(swm_corr))
  for(m in 1:n){
    swm_corr_sum[m,,]<-apply(swm_corr[m,,],2,function(x){out<-mov_sum(x,wdow);return(out)})
  }
  occurrence_emp_corr_sum<-apply(swm_corr_sum,c(1,2),function(x){out<-exceed_fun_trivar(x,obs_sum_est);return(out)})
  emp_exc_corr_sum<-sum(occurrence_emp_corr_sum)/(length(yrs)*n)

  #intermittent
  swm_corr_int_sum<-array(NA,dim(swm_corr_int))
  for(m in 1:n){
    swm_corr_int_sum[m,,]<-apply(swm_corr_int[m,,],2,function(x){out<-mov_sum(x,wdow);return(out)})
  }
  occurrence_emp_corr_int_sum<-apply(swm_corr_int_sum,c(1,2),function(x){out<-exceed_fun_trivar(x,obs_sum_est);return(out)})
  emp_exc_corr_int_sum<-sum(occurrence_emp_corr_int_sum)/(length(yrs)*n)

  #independent
  swm_ind_sum<-array(NA,dim(swm_ind))
  for(m in 1:n){
    swm_ind_sum[m,,]<-apply(swm_ind[m,,],2,function(x){out<-mov_sum(x,wdow);return(out)})
  }
  occurrence_emp_ind_sum<-apply(swm_ind_sum,c(1,2),function(x){out<-exceed_fun_trivar(x,obs_sum_est);return(out)})
  emp_exc_ind_sum<-sum(occurrence_emp_ind_sum)/(length(yrs)*n)

  saveRDS(c(emp_exc_corr_sum,emp_exc_corr_int_sum,emp_exc_ind_sum),
    paste('out/trivar-exceedance_',wdow,'d-sum_',sel_sites[1],'-',sel_sites[2],'-',sel_sites[3],'_',sged_mv2,'_ar',ar,'_',dsgn,'-yr-event.rds',sep=''))

  print('daily trivar return period')
  print(c(1/exceedance_obs_annmax,'obs'))
  print(c(1/exceedance_obsref_annmax,'sim'))
  print(c(1/emp_exc_corr,'mv'))
  print(c(1/emp_exc_corr_int,'mv-log'))
  print(c(1/emp_exc_ind,'ind'))

  print('sum trivar return period')
  print(c(1/exceedance_obs_sum_annmax,'obs'))
  print(c(1/exceedance_obsref_sum_annmax,'sim'))
  print(c(1/emp_exc_corr_sum,'mv'))
  print(c(1/emp_exc_corr_int_sum,'mv-log'))
  print(c(1/emp_exc_ind_sum,'ind'))
  }

rm(list=ls());gc()

###############################################END###############################################