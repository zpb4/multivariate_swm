#date
setwd('z:/mv_swm/')
library(extRemes)
library(Hmisc)

ix<-seq(as.Date('1988-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)
cmb<-2:4
dsgn<-100
yrs<-88:113

arr_full<-readRDS('data/arr_full.rds')

obs<-arr_full[1,,]
sim<-arr_full[2,,]

obs_comb<-apply(obs[,cmb],1,sum)
sim_comb<-apply(sim[,cmb],1,sum)

obs_annmax<-array(NA,c(length(yrs),14))
obs_comb_annmax<-c()

sim_annmax<-array(NA,c(length(yrs),14))
sim_comb_annmax<-c()

for(i in 1:length(yrs)){
  idx<-which(ix2$year==yrs[i])

  obs_annmax[i,]<-apply(obs[idx,],2,max)
  obs_comb_annmax[i]<-max(obs_comb[idx])
  
  sim_annmax[i,]<-apply(sim[idx,],2,max)
  sim_comb_annmax[i]<-max(sim_comb[idx])
}

dsgn_prob<-1 - (1/dsgn)

gev_fun<-function(x){
  fit<-fevd(x,type='GEV')
  dsn<-qevd(dsgn_prob,loc=fit$results$par[1],scale=fit$results$par[2],shape=fit$results$par[3],type='GEV')
  return(dsn)}

sim_gev_fit<-fevd(sim_comb_annmax,type='GEV')
sim_comb_est<-qevd(dsgn_prob,loc=sim_gev_fit$results$par[1],scale =sim_gev_fit$results$par[2],shape = sim_gev_fit$results$par[3],type='GEV')
sim_est<-apply(sim_annmax,2,gev_fun)

obs_gev_fit<-fevd(obs_comb_annmax,type='GEV')
obs_comb_est<-qevd(dsgn_prob,loc=obs_gev_fit$results$par[1],scale =obs_gev_fit$results$par[2],shape = obs_gev_fit$results$par[3],type='GEV')
obs_est<-apply(obs_annmax,2,gev_fun)

saveRDS(obs_annmax,'data/obs_annmax.rds')
saveRDS(sim_annmax,'data/sim_annmax.rds')

saveRDS(obs_comb_annmax,paste('data/obs_comb_annmax_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))
saveRDS(sim_comb_annmax,paste('data/sim_comb_annmax_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))

saveRDS(obs_comb_est,paste('data/obs_design-est_',dsgn,'-yr-event_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))
saveRDS(sim_comb_est,paste('data/sim_design-est_',dsgn,'-yr-event_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))

saveRDS(obs_est,paste('data/obs_design-est_',dsgn,'-yr-event_all-sites.rds',sep=''))
saveRDS(sim_est,paste('data/sim_design-est_',dsgn,'-yr-event_all-sites.rds',sep=''))

rm(list=ls());gc()

###############################################END########################################################