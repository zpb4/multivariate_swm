setwd('z:/mv_swm/')

library(scales)
library(fGarch)
library(corrplot)
library(ggplot2)
library(ggcorrplot)

arr_full<-readRDS('data/arr_full.rds')
obs<-arr_full[1,,]
sim<-arr_full[2,,]
emp_err<-arr_full[3,,]

ar<-3
arfit<-'pen-cbias' #penalized or ols version of VAR model and residuals
arfit_at<-'pen-cbias'
sites <- 14
n<-1000

#date/time indices
ix<-seq(as.Date('1988-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

locs_full<-c('SHA', 'ORO','YRS','FOL','MKM','NHG','NML','TLG','MRC','MIL','PNF','TRM','SCC','ISB')

at_lst<-readRDS(paste('fit/at_lst_prinorm-ls_',arfit_at,'_ar',ar,'.rds',sep=''))
cbias_resids<-readRDS('fit/rresids-cbias_full.rds')

#-------------------------------------------

#annual empirical residual correlations
emp_resid<-array(NA,dim(sim))
for(i in 1:12){
  seas<-which(ix2$mo==(i-1))
  emp_resid[seas,]<-at_lst[[i]]
}

#annual empirical residual correlations
emp_err<-array(NA,dim(sim))
for(i in 1:12){
  seas<-which(ix2$mo==(i-1))
  emp_err[seas,]<-cbias_resids[[i]]
}


corr_emp<-cor(emp_err,method = 'spearman')
colnames(corr_emp)<-locs_full
rownames(corr_emp)<-locs_full
#corrplot(corr_emp,addCoefasPercent = F,addCoef.col = 'black',tl.col='black',method = 'shade',diag = F)
#mtext('Spearman residual correlation',3,line=4.5,font=2,cex=2)

corr_new<-corr_emp[sites:1,sites:1]
corr_new[corr_new==1]<-NA

gplot<-ggcorrplot(corr_new[sites:1,],show.diag = T,lab=T,colors=c('white','white','dodgerblue3'),
                  show.legend = F,ggtheme=ggplot2::theme_void)
ggsave('h:/mv_swm/paper/figs/all-sites-corrplot_fig1_ggplot.png',gplot,dpi=320,width=8,height=8,units = 'in')
#--------------------------------------------------
