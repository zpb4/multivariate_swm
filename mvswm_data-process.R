#Initial data processing script
setwd('h:/mv_swm')
source('mm-cfs_conversion.R')

#seasonal inference stuff
ix<-seq(as.Date('1988-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

#1) Read in observed and simulated inflows

#order- SHA, BND, ORO, YRS, FOL, MKM, NHG, NML, TLG, MRC, MIL, PNF, TRM, SCC, ISB
##exclude BND which has less days than the rest, use 'full' timespan
locs_full<-c('SHA', 'ORO','YRS','FOL','MKM','NHG','NML','TLG','MRC','MIL','PNF','TRM','SCC','ISB')
mos<-c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
area_mi2_full<-c(area_mi2_SHA,area_mi2_ORO,area_mi2_YRS,area_mi2_FOL,area_mi2_MKM,area_mi2_NHG,area_mi2_NML,area_mi2_TLG,area_mi2_MRC,
                 area_mi2_MIL,area_mi2_PNF,area_mi2_TRM,area_mi2_SCC,area_mi2_ISB)

arr_full<-array(0,c(5,length(ix),length(locs_full)))
for(i in 1:3){
  colnames(arr_full[i,,])<-locs_full
  rownames(arr_full[i,,])<-as.character(ix)
}

for(i in 1:length(locs_full)){
  ws<-locs_full[i]
  dat<-read.table(paste('data/SAC-SMA_CA/simflow_sacramento_',ws,'.txt',sep=''))
  obs<-dat$V5;obs[obs<0]<-0
  obs_log<-obs;obs_log[obs_log>0]<-1
  sim<-dat$V4
  
  #determine where obs are zero and non-zero
  zero_idx<-which(obs==0)
  nzero_idx<-1:length(sim)
  nzero_idx<-nzero_idx[-c(zero_idx)]
  
  arr_full[1,,i]<-obs*area_mi2_full[i]*mm_to_cfs
  arr_full[2,,i]<-sim*area_mi2_full[i]*mm_to_cfs
  arr_full[3,,i]<-sim-obs*area_mi2_full[i]*mm_to_cfs #differenced errors
  arr_full[4,nzero_idx,i]<-log(sim[nzero_idx]/obs[nzero_idx]) #only calculate log ratio errors for non-zero elements
  arr_full[5,,i]<-obs_log #create matrix of zeros and ones for logistic regression models
}

saveRDS(arr_full,'data/arr_full.rds')

rm(list=ls());gc()

###############################################END################################################################