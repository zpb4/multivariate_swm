setwd('z:/mv_swm/')

library(scales)
library(tidyverse)
library(data.table)
library(gridExtra)
library(ks)

arr_full<-readRDS('data/arr_full.rds')
obs<-arr_full[1,,]
sim<-arr_full[2,,]
emp_err<-arr_full[3,,]

ar<-3
arfit_mv<-'pen-cbias' #penalized or ols version of VAR model and residuals
arfit_ind<-'ols-cbias'
sites <- 14
n<-1000
percentile<-0.9
fill_dens<-F

#date/time indices
ix<-seq(as.Date('1988-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

locs_full<-c('SHA', 'ORO','YRS','FOL','MKM','NHG','NML','TLG','MRC','MIL','PNF','TRM','SCC','ISB')

syn_flow_corr<-readRDS(paste('out/syn_flow_prinorm-ls_',arfit_mv,'_ar',ar,'.rds',sep=''))
#syn_flow_ind<-readRDS(paste('out/syn_flow_prinorm-ls_ind_',arfit_ind,'_ar',ar,'.rds',sep=''))
#syn_flow_corr_int<-readRDS(paste('out/syn_flow_prinorm-ls_',arfit,'_ar',ar,'_intermittent.rds',sep=''))

#-------------------------------------------

pcntile_fun<-function(x,pcntile){
  two_tail<-(1-pcntile)/2
  upr_idx<-round((1-two_tail)*length(x))
  lwr_idx<-round((two_tail)*length(x))
  srt_dat<-sort(x)
  upr<-srt_dat[upr_idx]
  lwr<-srt_dat[lwr_idx]
  return(c(lwr,upr))
}

bnds_corr<-apply(syn_flow_corr,c(2,3),function(x){out<-pcntile_fun(x,percentile);return(out)})
bnds_ind<-apply(syn_flow_ind,c(2,3),function(x){out<-pcntile_fun(x,percentile);return(out)})
bnds_corr_int<-apply(syn_flow_corr_int,c(2,3),function(x){out<-pcntile_fun(x,percentile);return(out)})

sum_obs<-apply(obs[,2:4],1,sum)
max_idx<-which.max(sum_obs)
min_idx<-which.min(sum_obs)

#plot ORO, YRS, FOL timeseries plot
max_days<-(max_idx-15):(max_idx+15)
my.gray<-alpha('gray',alpha=0.4)
samp<-sample(1:n,1)
#samp<-863
print(ix[max_idx])

oro_dt<-data.table(x=-15:15,obs=obs[max_days,2]/1000,sim=sim[max_days,2]/1000,samp=syn_flow_corr[samp,max_days,2]/1000,
                   ymax=bnds_corr[1,max_days,2]/1000,ymin=bnds_corr[2,max_days,2]/1000)
yrs_dt<-data.table(x=-15:15,obs=obs[max_days,3]/1000,sim=sim[max_days,3]/1000,samp=syn_flow_corr[samp,max_days,3]/1000,
                   ymax=bnds_corr[1,max_days,3]/1000,ymin=bnds_corr[2,max_days,3]/1000)
fol_dt<-data.table(x=-15:15,obs=obs[max_days,4]/1000,sim=sim[max_days,4]/1000,samp=syn_flow_corr[samp,max_days,4]/1000,
                   ymax=bnds_corr[1,max_days,4]/1000,ymin=bnds_corr[2,max_days,4]/1000)


oro_ggp<-ggplot(data=oro_dt)+theme_minimal()+
  geom_ribbon(mapping=aes(x=x,ymax=ymax,ymin=ymin),fill='gray70',alpha=0.5)+
  geom_line(mapping=aes(x=x,y=obs,color='obs'),linewidth=0.75)+
  geom_line(mapping=aes(x=x,y=sim,color='sim'),linewidth=0.75,alpha=0.5)+
  geom_line(mapping=aes(x=x,y=samp,color='samp'),linewidth=0.5,alpha=0.5)+
  geom_label(x=-9,y=0.6*max(oro_dt),label='ORO',size=5)+
  annotate('text',x=-12.5,y=0.9*max(oro_dt),label='a)',size=6)+
  labs(x='',y='flow (kcfs)')+
  coord_cartesian(expand=F)+
  scale_color_manual(breaks=c('obs','sim','samp'),values=c('black','dodgerblue3','seagreen4'))+
  theme(legend.position=c(.8,.75),
        legend.title=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=12),
        legend.text = element_text(size=9))

yrs_ggp<-ggplot(data=yrs_dt)+theme_minimal()+
  geom_ribbon(mapping=aes(x=x,ymax=ymax,ymin=ymin),fill='gray70',alpha=0.5)+
  geom_line(mapping=aes(x=x,y=obs,color='obs'),linewidth=0.75)+
  geom_line(mapping=aes(x=x,y=sim,color='sim'),linewidth=0.75,alpha=0.5)+
  geom_line(mapping=aes(x=x,y=samp,color='samp'),linewidth=0.5,alpha=0.5)+
  geom_label(x=-9,y=0.6*max(yrs_dt),label='YRS',size=5)+
  labs(x='',y='flow (kcfs)')+
  coord_cartesian(expand=F)+
  scale_color_manual(breaks=c('obs','sim','samp'),values=c('black','dodgerblue3','seagreen4'))+
  theme(legend.position='none',
        legend.title=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=12),
        legend.text = element_text(size=9))

fol_ggp<-ggplot(data=fol_dt)+theme_minimal()+
  geom_ribbon(mapping=aes(x=x,ymax=ymax,ymin=ymin),fill='gray70',alpha=0.5)+
  geom_line(mapping=aes(x=x,y=obs,color='obs'),linewidth=0.75)+
  geom_line(mapping=aes(x=x,y=sim,color='sim'),linewidth=0.75,alpha=0.5)+
  geom_line(mapping=aes(x=x,y=samp,color='samp'),linewidth=0.5,alpha=0.5)+
  geom_label(x=-9,y=0.6*max(fol_dt),label='FOL',size=5)+
  labs(x='days from max flow',y='flow (kcfs)')+
  coord_cartesian(expand=F)+
  scale_color_manual(breaks=c('obs','sim','samp'),values=c('black','dodgerblue3','seagreen4'))+
  theme(legend.position='none',
        legend.title=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=12),
        legend.text = element_text(size=9))

#mat<-matrix(1:3,nrow=3)
#gplot<-grid.arrange(oro_ggp,yrs_ggp,fol_ggp,layout_matrix=mat)
#gplot
#ggsave(paste('h:/mv_swm/paper/figs/ts_',arfit_mv,'_ar',ar,'.png',sep=''),gplot,dpi=320,width=3,height=9,unit='in')

#-------------------------------------------
#design event plots
#10 year
cmb<-2:4
dsgn<-10

comb_out_sged_mv<-readRDS(paste('out/syn_design-est_sged_mv_prinorm-ls_',arfit_mv,'_ar',ar,'_',dsgn,'-yr-event_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))/1e3
comb_out_sged_ind<-readRDS(paste('out/syn_design-est_sged_ind_prinorm-ls_ind_',arfit_ind,'_ar',ar,'_',dsgn,'-yr-event_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))/1e3

obs_comb_est<-readRDS(paste('h:/mv_swm/data/obs_design-est_',dsgn,'-yr-event_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))/1e3
sim_comb_est<-readRDS(paste('h:/mv_swm/data/sim_design-est_',dsgn,'-yr-event_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))/1e3

smth_mv<-kde(comb_out_sged_mv,positive = T,gridsize = 10000)
smth_ind<-kde(comb_out_sged_ind,positive = T,gridsize = 10000)
max_y<-max(smth_ind$estimate,smth_mv$estimate)

dsgn_dat<-data.table(ind=c(comb_out_sged_ind),mv=c(comb_out_sged_mv),obs=obs_comb_est,sim=sim_comb_est)

if(fill_dens==F){
densplot_10<-ggplot(data=dsgn_dat)+theme_minimal()+
  geom_density(aes(x=ind),color='tan2',adjust=1.5,linewidth=1)+
  geom_hline(aes(yintercept=-1,color='ind'))+
  geom_density(aes(x=mv),color='seagreen4',adjust=1.5,linewidth=1)+
  geom_hline(aes(yintercept=-1,color='mv'))+
  geom_point(aes(x=obs,y=0,color='obs'),pch=17,size=3)+
  geom_point(aes(x=sim,y=0,color='sim'),pch=17,size=3)+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
  scale_color_manual(breaks=c('obs','sim','ind','mv'),values=c('black','dodgerblue3','tan2','seagreen4'),
                     guide=guide_legend(override.aes = list(shape=c(17,17,20,20),size=c(3,3,0,0),
                                                            linewidth=c(0,0,1,1),linetype=c(0,0,1,1))))+
  coord_cartesian(xlim=c(0,2*obs_comb_est),ylim=c(0,max_y),expand=T)+
  labs(x='flow (kcfs)',y='density')+
  annotate('text',x=0.05*2*obs_comb_est,y=0.975*max_y,label='b)',size=6)+
  theme(legend.position=c(0.75,.7),
        legend.title=element_blank(),
        #axis.text.x=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=12),
        legend.text = element_text(size=9),
        legend.spacing.y = unit(.1,'lines'))
}

#ggsave('h:/mv_swm/paper/figs/test.png')

if(fill_dens==T){
densplot_10<-ggplot(data=dsgn_dat)+theme_minimal()+
  geom_density(aes(x=ind),color='tan2',fill='yellow2',alpha=0.1,adjust=1.5,linewidth=1)+
  geom_hline(aes(yintercept=-1,color='ind'))+
  geom_density(aes(x=mv),color='seagreen4',fill='seagreen3',alpha=0.1,adjust=1.5,linewidth=1)+
  geom_hline(aes(yintercept=-1,color='mv'))+
  geom_point(aes(x=obs,y=0,color='obs'),pch=17,size=3)+
  geom_point(aes(x=sim,y=0,color='sim'),pch=17,size=3)+
  scale_color_manual(breaks=c('obs','sim','ind','mv'),values=c('black','dodgerblue3','tan2','seagreen4'),
                     guide=guide_legend(override.aes = list(shape=c(17,17,20,20),size=c(3,3,0,0),
                                                            linewidth=c(0,0,1,1),linetype=c(0,0,1,1))))+
  coord_cartesian(xlim=c(0,2*obs_comb_est),ylim=c(0,max_y),expand=T)+
  labs(x='flow (kcfs)',y='density')+
  annotate('text',x=0.05*2*obs_comb_est,y=0.975*max_y,label='b)',size=6)+
  theme(legend.position=c(0.75,.7),
        legend.title=element_blank(),
        #axis.text.x=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=12),
        legend.text = element_text(size=9),
        legend.spacing.y = unit(.1,'lines'))
}

#ggsave('h:/mv_swm/paper/figs/test.png')



#-------------------------------------------
#design event plots
#100 year
cmb<-2:4
dsgn<-100

comb_out_sged_mv<-readRDS(paste('out/syn_design-est_sged_mv_prinorm-ls_',arfit_mv,'_ar',ar,'_',dsgn,'-yr-event_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))/1e3
comb_out_sged_ind<-readRDS(paste('out/syn_design-est_sged_ind_prinorm-ls_ind_',arfit_ind,'_ar',ar,'_',dsgn,'-yr-event_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))/1e3

obs_comb_est<-readRDS(paste('h:/mv_swm/data/obs_design-est_',dsgn,'-yr-event_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))/1e3
sim_comb_est<-readRDS(paste('h:/mv_swm/data/sim_design-est_',dsgn,'-yr-event_sites_',cmb[1],'-',tail(cmb,1),'.rds',sep=''))/1e3

smth_mv<-kde(comb_out_sged_mv,positive = T,gridsize = 10000)
smth_ind<-kde(comb_out_sged_ind,positive = T,gridsize = 10000)
max_y<-max(smth_ind$estimate,smth_mv$estimate)

dsgn_dat<-data.table(ind=c(comb_out_sged_ind),mv=c(comb_out_sged_mv),obs=obs_comb_est,sim=sim_comb_est)

if(fill_dens==F){
  densplot_100<-ggplot(data=dsgn_dat)+theme_minimal()+
    geom_density(aes(x=ind),color='tan2',adjust=1.5,linewidth=1)+
    geom_hline(aes(yintercept=-1,color='ind'))+
    geom_density(aes(x=mv),color='seagreen4',adjust=1.5,linewidth=1)+
    geom_hline(aes(yintercept=-1,color='mv'))+
    geom_point(aes(x=obs,y=0,color='obs'),pch=17,size=3)+
    geom_point(aes(x=sim,y=0,color='sim'),pch=17,size=3)+
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    scale_color_manual(breaks=c('obs','sim','ind','mv'),values=c('black','dodgerblue3','tan2','seagreen4'),
                       guide=guide_legend(override.aes = list(shape=c(17,17,20,20),size=c(3,3,0,0),
                                                              linewidth=c(0,0,1,1),linetype=c(0,0,1,1))))+
    coord_cartesian(xlim=c(0,2.5*obs_comb_est),ylim=c(0,max_y),expand=T)+
    labs(x='flow (kcfs)',y='density')+
    annotate('text',x=0.05*2*obs_comb_est,y=0.975*max_y,label='c)',size=6)+
    theme(legend.position='none',
          legend.title=element_blank(),
          #axis.text.x=element_blank(),
          axis.text=element_text(size=9),
          axis.title=element_text(size=12),
          legend.text = element_text(size=9),
          legend.spacing.y = unit(.1,'lines'))
}

#densplot_100
#ggsave('h:/mv_swm/paper/figs/test.png')

if(fill_dens==T){
  densplot_100<-ggplot(data=dsgn_dat)+theme_minimal()+
    geom_density(aes(x=ind),color='tan2',fill='yellow2',alpha=0.1,adjust=1.5,linewidth=1)+
    geom_hline(aes(yintercept=-1,color='ind'))+
    geom_density(aes(x=mv),color='seagreen4',fill='seagreen3',alpha=0.1,adjust=1.5,linewidth=1)+
    geom_hline(aes(yintercept=-1,color='mv'))+
    geom_point(aes(x=obs,y=0,color='obs'),pch=17,size=3)+
    geom_point(aes(x=sim,y=0,color='sim'),pch=17,size=3)+
    scale_color_manual(breaks=c('obs','sim','ind','mv'),values=c('black','dodgerblue3','tan2','seagreen4'),
                       guide=guide_legend(override.aes = list(shape=c(17,17,20,20),size=c(3,3,0,0),
                                                              linewidth=c(0,0,1,1),linetype=c(0,0,1,1))))+
    coord_cartesian(xlim=c(0,2.5*obs_comb_est),ylim=c(0,max_y),expand=T)+
    labs(x='flow (kcfs)',y='density')+
    annotate('text',x=0.05*2*obs_comb_est,y=0.975*max_y,label='c)',size=6)+
    theme(legend.position=c(0.75,.7),
          legend.title=element_blank(),
          #axis.text.x=element_blank(),
          axis.text=element_text(size=9),
          axis.title=element_text(size=12),
          legend.text = element_text(size=9),
          legend.spacing.y = unit(.1,'lines'))
}

#ggsave('h:/mv_swm/paper/figs/test.png')

#-------------------------------------------------------
dsgn<-100
wdow<-3
sged_mv<-paste('prinorm-ls_',arfit_mv,sep='')
sged_ind<-paste('prinorm-ls_ind_',arfit_ind,sep='')
sel_sites<-c('ORO','YRS','FOL')

#Plot individual basins - Annual Maxima
#c(exceedance_obs_annmax,exceedance_obs_sum_annmax,exceedance_simref_annmax,exceedance_simref_sum_annmax,exceedance_obsref_annmax,exceedance_obsref_sum_annmax)
#c(emp_exc_corr,emp_exc_corr_int,emp_exc_ind)
#c(emp_exc_corr_sum,emp_exc_corr_int_sum,emp_exc_ind_sum)
obs_sim_trivar<-readRDS(paste('out/trivar-exceedance_obs-sim_',sel_sites[1],'-',sel_sites[2],'-',sel_sites[3],'_',dsgn,'-yr-event.rds',sep=''))
swm_trivar<-readRDS(paste('out/trivar-exceedance_',sel_sites[1],'-',sel_sites[2],'-',sel_sites[3],'_',sged_mv,'_ar',ar,'_',dsgn,'-yr-event.rds',sep=''))
swm_sum_trivar<-readRDS(paste('out/trivar-exceedance_',wdow,'d-sum_',sel_sites[1],'-',sel_sites[2],'-',sel_sites[3],'_',sged_mv,'_ar',ar,'_',dsgn,'-yr-event.rds',sep=''))

#exceedance prob plots
#daily flow
#bplot_vec<-c(swm_trivar[3],swm_trivar[2],swm_trivar[1])
#obs_est<-1/obs_sim_trivar[1]
#sim_est<-1/obs_sim_trivar[3]
bplot_vec<-c(swm_sum_trivar[3],swm_sum_trivar[2],swm_sum_trivar[1])
obs_est<-1/obs_sim_trivar[2]
sim_est<-1/obs_sim_trivar[4]
print(paste('obs',obs_est))
print(paste('sim',sim_est))

#return period plots
rtns<-1/bplot_vec

max_y<-max(rtns)*1.1
min_y<-min(rtns)*0.9

#trivar_hi_dat<-data.table(ind=rtns[1],mv_alog=rtns[2],mv_trunc=rtns[3])
trivar_hi_dat<-data.table(dat=rtns)

trivar_hi_gplot<-ggplot(data=trivar_hi_dat)+theme_minimal()+
  geom_col(aes(x=1:3,y=dat,fill=c('ind','mv-alog','mv-trunc')),width=0.5)+
  scale_fill_manual(breaks=c('ind','mv-alog','mv-trunc'),values = c('tan2','seagreen4','gray75'))+
  scale_x_continuous(breaks=1:3,labels=c('ind','mv-alog','mv-trunc'))+
  labs(x='',y='return period (years)')+
  annotate('text',x=0.5,y=0.95*max_y,label='d)',size=6)+
  theme(legend.position='none',
        legend.title=element_blank(),
        #axis.text.x=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=12),
        legend.text = element_text(size=9))


#-------------------------------------------------------
wdow<-7
sel_sites<-c('ORO','YRS','FOL')

#Plot individual basins - Annual Maxima
#c(emp_exc_obs_avg,emp_exc_sim_avg,emp_exc_corr_avg,emp_exc_corr_int_avg,emp_exc_ind_avg)
chinook_thresh<-readRDS(paste('out/trivar-loflow_',wdow,'d-avg_',sel_sites[1],'-',sel_sites[2],'-',sel_sites[3],'_',sged_mv,'_ar',ar,'_chinook-thresh.rds',sep=''))
obs_7q10<-readRDS(paste('out/trivar-loflow_',wdow,'d-avg_',sel_sites[1],'-',sel_sites[2],'-',sel_sites[3],'_',sged_mv,'_ar',ar,'_7Q10.rds',sep=''))

#exceedance prob plots
#daily flow
loflow_vec<-chinook_thresh[c(1,2,5,4,3)]

obs_est<-loflow_vec[1]
sim_est<-loflow_vec[2]

print(paste('obs',obs_est))
print(paste('sim',sim_est))

bplot_vec<-loflow_vec[3:5]

max_y<-max(bplot_vec)*1.1
min_y<-min(bplot_vec)*0.9


#trivar_hi_dat<-data.table(ind=rtns[1],mv_alog=rtns[2],mv_trunc=rtns[3])
trivar_lo_dat<-data.table(dat=bplot_vec)

trivar_lo_gplot<-ggplot(data=trivar_lo_dat)+theme_minimal()+
  geom_col(aes(x=1:3,y=dat,fill=c('ind','mv-alog','mv-trunc')),width=0.5)+
  scale_fill_manual(breaks=c('ind','mv-alog','mv-trunc'),values = c('tan2','seagreen4','gray75'))+
  scale_x_continuous(breaks=1:3,labels=c('ind','mv-alog','mv-trunc'))+
  labs(x='',y='7d low flow frequency (per year)')+
  annotate('text',x=0.5,y=0.95*max_y,label='e)',size=6)+
  theme(legend.position='none',
        legend.title=element_blank(),
        #axis.text.x=element_blank(),
        axis.text=element_text(size=9),
        axis.title=element_text(size=12),
        legend.text = element_text(size=9))
#trivar_lo_gplot
#--------------------------------------------------

mat<-rbind(c(1,1,1,2,2,2,2,3,3,3,3),c(1,1,1,2,2,2,2,3,3,3,3),c(4,4,4,2,2,2,2,3,3,3,3),
           c(4,4,4,5,5,5,5,6,6,6,6),c(7,7,7,5,5,5,5,6,6,6,6),c(7,7,7,5,5,5,5,6,6,6,6))
#plot
pl1<-oro_ggp
pl2<-densplot_10
pl3<-trivar_hi_gplot
pl4<-yrs_ggp
pl5<-densplot_100
pl6<-trivar_lo_gplot
pl7<-fol_ggp

grid.arrange(pl1,pl2,pl3,
             pl4,pl5,pl6,
             pl7,layout_matrix=mat)

#save for manuscript
gplot<-grid.arrange(pl1,pl2,pl3,
                    pl4,pl5,pl6,
                    pl7,layout_matrix=mat)

ggsave(paste('h:/mv_swm/paper/figs/',arfit_mv,'-ar',ar,'_ts-dsgn-trivar_fig4_ggplot.png',sep=''),gplot,dpi=320,width=8,height=5,unit='in')


###############################################END#######################################
