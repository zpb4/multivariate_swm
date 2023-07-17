setwd('z:/mv_swm/')

library(scales)
library(data.table)
library(gridExtra)
library(fGarch)
library(ggcorrplot)
library(markovchain)

arr_full<-readRDS('data/arr_full.rds')
obs<-arr_full[1,,]
sim<-arr_full[2,,]
emp_err<-arr_full[3,,]

ar<-3
arfit<-'pen-cbias' #penalized or ols version of VAR model and residuals
sites <- 14
n<-1000

#date/time indices
ix<-seq(as.Date('1988-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

locs_full<-c('SHA', 'ORO','YRS','FOL','MKM','NHG','NML','TLG','MRC','MIL','PNF','TRM','SCC','ISB')

syn_flow_corr<-readRDS(paste('out/syn_flow_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))
syn_resid_corr<-readRDS(paste('out/syn_resid_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))
syn_empcop_corr<-readRDS(paste('out/syn_empcop_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))
at_lst<-readRDS(paste('fit/at_lst_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))
gl_par_arr<-readRDS(paste('fit/gl_par_arr_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))
#-------------------------------------------
#plot histogram and fitted distr
mths<-c(12,6)
mth_lab<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
#site_sel<-c('SHA','FOL','MRC','ISB')
site_sel<-c('FOL','NHG')
gg_dist_out<-vector('list',4)
idx_mat<-matrix(1:4,ncol=2)
pl_lab_mat<-matrix(c('a)','','',''),ncol=2)
xlab_mat<-matrix(c('','','','residual'),ncol=2)

for(i in 1:length(mths)){
  seas<-which(ix2$mon==(mths[i]-1))
  for(k in 1:length(site_sel)){
    site<-which(locs_full==site_sel[k])
    #mat<-at_lst[[mths[i]]]
    Mu<-format(round(gl_par_arr[mths[i],site,1],2),nsmall=2)
    Sig<-format(round(gl_par_arr[mths[i],site,2],2),nsmall=2)
    Beta<-format(round((2/gl_par_arr[mths[i],site,3]-1),2),nsmall=2)
    Xi<-format(round(gl_par_arr[mths[i],site,4],2),nsmall=2)
    
    site<-k
    mat<-at_lst[[mths[i]]]
    inp<-mat[,k]
    
    dist_dat<-data.table(emp=inp)
    x<-seq(-10,10,.01)
    y<-dsged(seq(-10,10,.01),mean=gl_par_arr[mths[i],site,1],sd=gl_par_arr[mths[i],site,2],nu=gl_par_arr[mths[i],site,3],xi=gl_par_arr[mths[i],site,4])
    fit_dat<-data.table(x=x,y=y)
    
    g<-idx_mat[i,k]
    l<-pl_lab_mat[i,k]
    xlab<-xlab_mat[i,k]
    
    gg_dist_out[[g]]<-ggplot(data=dist_dat)+theme_minimal()+
      geom_histogram(aes(x=emp,y=after_stat(density)),binwidth = .3,color='gray30',fill='gray',linewidth=0.1)+
      coord_cartesian(xlim=c(-5,5),ylim=c(0,0.8),expand=F)+
      geom_line(data=fit_dat,aes(x=x,y=y),col='tomato3',linewidth=.5)+
      annotate('text',x=3,y=0.75,label=bquote(~mu==.(Mu)),size=2.5)+
      annotate('text',x=3,y=0.65,label=bquote(~sigma==.(Sig)),size=2.5)+
      annotate('text',x=3,y=0.55,label=bquote(~beta==.(Beta)),size=2.5)+
      annotate('text',x=3,y=0.45,label=bquote(~xi==.(Xi)),size=2.5)+
      annotate('text',x=-4,y=0.75,label=l,size=6)+
      geom_label(x=-3,y=0.5,label=site_sel[k],size=3.5)+
      labs(x=xlab,y='density')+
      annotate('text',x=-3,y=0.3,label=paste('italic(',mth_lab[mths[i]],')',sep=''),size=3,parse=T)+
      theme(
          axis.text=element_text(size=6),
          axis.title=element_text(size=8))
  }
}

#-------------------------------------------------------
#spearman correlation plot
syn_resid_corr<-readRDS(paste('out/syn_resid_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))
cbias_resids<-readRDS('fit/rresids-cbias_full.rds')

site_sel<-c('ORO','YRS','FOL','NHG','MRC','SCC')
site_idx<-c()
for(i in 1:length(site_sel)){site_idx[i]<-which(locs_full==site_sel[i])}

emp_err<-array(NA,dim(sim))
for(i in 1:12){
  seas<-which(ix2$mo==(i-1))
  emp_err[seas,]<-cbias_resids[[i]]
}

corr_emp<-cor(emp_err[,site_idx],method = 'spearman')
colnames(corr_emp)<-site_sel
rownames(corr_emp)<-site_sel
corr_pred<-cor(syn_resid_corr[sample(1:n,1),,site_idx],method = 'spearman')
corr_emp_gg<-corr_emp

corr_emp_gg[which(lower.tri(corr_emp_gg)==T)]<-corr_pred[which(lower.tri(corr_emp_gg)==T)]
corr_new<-corr_emp_gg[6:1,6:1]
corr_new[corr_new==1]<-NA

spear_corr<-ggcorrplot(corr_new[length(site_sel):1,],show.diag = T,lab_size=3,tl.cex=8,
                  lab=T,colors=c('white','white','dodgerblue3'),ggtheme=ggplot2::theme_void,show.legend = F)+
                  annotate('text',x=1,y=6,label='d)',size=6)

#-----------------------------------------------
#binary pearson plots
syn_flow_corr_int<-readRDS(paste('out/syn_flow_prinorm-ls_',arfit,'_ar',ar,'_autolog_obsim-min.rds',sep=''))

corr_obs<-cor(arr_full[5,,site_idx])
colnames(corr_obs)<-site_sel
rownames(corr_obs)<-site_sel
syn_bin_corr<-syn_flow_corr_int[sample(1:n,1),,site_idx]
syn_bin_corr[syn_bin_corr>0]<-1
corr_pred<-cor(syn_bin_corr)
corr_obs_gg<-corr_obs

corr_obs_gg[which(lower.tri(corr_obs)==T)]<-corr_pred[which(lower.tri(corr_obs_gg)==T)]
corr_new<-corr_obs_gg[length(site_sel):1,length(site_sel):1]
corr_new[corr_new==1]<-NA

bin_corr<-ggcorrplot(corr_new[length(site_sel):1,],show.diag = T,lab_size=3,tl.cex=8,
                  lab=T,colors=c('white','white','dodgerblue3'),ggtheme=ggplot2::theme_void,show.legend = F)+
                  annotate('text',x=1,y=6,label='e)',size=6)


#-------------------------------------------------------
#barplots
sflow_corr_trans<-readRDS(paste('out/mc-trans-probs_syn_flow_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))
sflow_corr_int_trans<-readRDS(paste('out/mc-trans-probs_syn_flow_prinorm-ls_',arfit,'_ar',ar,'_autolog_obsim-min.rds',sep=''))
obs_mc<-apply(arr_full[5,,],2,function(x){out<-markovchainFit(x)$estimate@transitionMatrix; return(as.vector(out))})

trans_probs<-c('0-0','1-0','0-1','1-1')

trans_idx<-1

bplot_vec<-array(NA,c(3,2*length(site_sel)))

for(i in 1:length(site_sel)){
  corr_int_samps<-sflow_corr_int_trans[,site_idx[i],trans_idx]
  corr_samps<-sflow_corr_trans[,site_idx[i],trans_idx]
  bplot_vec[,(i*2-1)]<-c(median(corr_int_samps),min(corr_int_samps),max(corr_int_samps))
  bplot_vec[,(i*2)]<-c(median(corr_samps),min(corr_samps),max(corr_samps))
}

xpos<-rep(c(-.75,-0.35),length(site_sel))+rep(1:length(site_sel),each=2)
xobs<-rep(-0.55,length(site_sel))+1:length(site_sel)

bplot_dat<-data.table(x=xpos,xob=xobs,obs=obs_mc[trans_idx,site_idx],med=bplot_vec[1,],min=bplot_vec[2,],max=bplot_vec[3,])

bplot_trans<-ggplot(data=bplot_dat)+theme_minimal()+
  geom_col(aes(x=xpos,y=med,fill=rep(c('mv-alog','mv-trunc'),6)),width=.35)+
  coord_cartesian(ylim=c(0.3,1))+
  geom_errorbar(aes(x=xpos,ymax=max,ymin=min),width=.15)+
  geom_point(aes(x=xob,y=obs,color='obs'),pch=18,size=3)+
  scale_fill_manual(breaks=c('mv-alog','mv-trunc'),values=c('seagreen4','gray75'))+
  scale_color_manual(breaks='obs',values='black')+
  scale_x_continuous(breaks=xobs,labels=site_sel)+
  labs(x='',y='0-0 transition probability')+
  annotate('text',x=0,y=1,label='b)',size=6)+
  theme(legend.position=c(.25,.8),
        legend.title=element_blank(),
        #axis.text.x=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        legend.text = element_text(size=10),
        legend.spacing.y = unit(.1,'lines'))

#bplot_trans
#-------------------------------------------------------
#Non-zero flow rate boxplot
obs_rate<-apply(arr_full[5,,],2,function(x){length(which(x>0))/length(x)})
rand_rates<-apply(syn_flow_corr,c(1,3),function(x){length(which(x>0))/length(x)})
int_rates<-apply(syn_flow_corr_int,c(1,3),function(x){length(which(x>0))/length(x)})

bplot_vec<-array(NA,c(3,2*length(site_sel)))

for(i in 1:length(site_sel)){
  corr_int_samps<-int_rates[,site_idx[i]]
  corr_samps<-rand_rates[,site_idx[i]]
  bplot_vec[,(i*2-1)]<-c(median(corr_int_samps),min(corr_int_samps),max(corr_int_samps))
  bplot_vec[,(i*2)]<-c(median(corr_samps),min(corr_samps),max(corr_samps))
}

xpos<-rep(c(-.75,-0.35),length(site_sel))+rep(1:length(site_sel),each=2)
xobs<-rep(-0.55,length(site_sel))+1:length(site_sel)

bplot_dat<-data.table(x=xpos,xob=xobs,obs=obs_rate[site_idx],med=bplot_vec[1,],min=bplot_vec[2,],max=bplot_vec[3,])

bplot_nzero<-ggplot(data=bplot_dat)+theme_minimal()+
  geom_col(aes(x=xpos,y=med),fill=rep(c('seagreen4','gray75'),6),width=.35)+
  coord_cartesian(ylim=c(0.5,1))+
  geom_errorbar(aes(x=xpos,ymax=max,ymin=min),width=.15)+
  geom_point(aes(x=xob,y=obs),color='black',pch=18,size=3)+
  labs(x='',y='non-zero flow fraction')+
  annotate('text',x=0,y=1,label='c)',size=6)+
  scale_x_continuous(breaks=xobs,labels=site_sel)+
  theme(#axis.text.x=element_blank(),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10))
#bplot_nzero
#---------------------------------------------------------------

laymat<-rbind(c(1,1,2,2,2,2,3,3,3,3),c(1,1,2,2,2,2,3,3,3,3),c(4,4,2,2,2,2,3,3,3,3),
           c(4,4,2,2,2,2,3,3,3,3),c(5,5,6,6,6,6,7,7,7,7),c(5,5,6,6,6,6,7,7,7,7),
           c(8,8,6,6,6,6,7,7,7,7),c(8,8,6,6,6,6,7,7,7,7))
#plot

ggp1<-gg_dist_out[[1]]
ggp2<-bplot_trans
ggp3<-bplot_nzero
ggp4<-gg_dist_out[[2]]
ggp5<-gg_dist_out[[3]]
ggp6<-spear_corr
ggp7<-bin_corr
ggp8<-gg_dist_out[[4]]

#grid.arrange(ggp1,ggp2,ggp3,ggp4,
             #ggp5,ggp6,ggp7,ggp8,
             #layout_matrix=laymat)
#save for manuscript
gplot<-grid.arrange(ggp1,ggp2,ggp3,ggp4,
                    ggp5,ggp6,ggp7,ggp8,
                    layout_matrix=laymat)

ggsave(paste('h:/mv_swm/paper/figs/',arfit,'-ar',ar,'_fit-dist_intermitt_corr_fig3_ggplot.png',sep=''),gplot,dpi=320,width=8,height=6,unit='in')
########################################END#################################################