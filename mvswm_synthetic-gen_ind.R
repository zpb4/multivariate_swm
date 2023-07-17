#Generate Independent AR SGED model
#setwd('z:/mv_swm')

#Required Packages
library(fGarch)
library(doParallel)
library(foreach)
library(abind)

#change for number of simulations desired
n <- 1000

#other parameters
ar<-3
arfit<-'pen-cbias' #penalized or ols version of VAR model and residuals
sites <- 14
incl_mean<-F

#date/time indices
ix<-seq(as.Date('1988-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

#Load required data
syn_resid<-array(NA,c(length(ix2),sites))
syn_flow<-array(NA,c(length(ix2),sites))
syn_flow_out<-array(NA,c(n,length(ix2),sites))
syn_resid_out<-array(NA,c(n,length(ix2),sites))
gl_par_arr<-readRDS(paste('fit/gl_par_arr_prinorm-ls_ind_',arfit,'_ar',ar,'.rds',sep=''))
#gl_par_arr<-readRDS(paste('fit/gl_par_arr_prinorm-ls_',arfit,'.rds',sep=''))
ar_coefs<-readRDS(paste('fit/ar_coefs_prinorm-ls_ind_',arfit,'_ar',ar,'.rds',sep=''))
norm_fit<-readRDS('fit/norm_fit_prinorm-ls-cbias_full.rds')
sd_arr<-readRDS('fit/sd_arr_prinorm-ls-cbias_full.rds')

#1) parallelization setup
parallel::detectCores()

n.cores <- parallel::detectCores()-1
my.cluster<-parallel::makeCluster(n.cores,type='PSOCK')
print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

#2) Define observations and simulations for synthetic generation
arr_full<-readRDS('data/arr_full.rds')
obs<-arr_full[1,,]
#sim<-arr_full[2,,]
sim_cbias<-readRDS('fit/sim_cbias.rds')

#3) Set up date/time indices for sequential (by month) VAR sampling
yr_idx<-ix2$year
yr_idx_lst<-vector('list',12)

for(i in 1:12){
  seas<-which(ix2$mon==(i-1))
  yr_idx_lst[[i]]<-yr_idx[seas]
}

yr_seq<-c(rep(88,3),rep(89:112,each=12),rep(113,9))
mo_seq<-c(10:12,rep(1:12,length(89:112)),1:9)

#3a. define results binding function
abnd<-function(x,y){abind(x,y,along=3)}

arsim<-function(n,ar,coef,start,innov){
  out_vec<-c(start,rep(0,n))
  for(i in (ar+1):length(out_vec)){
    out_vec[i]<-coef%*%out_vec[(i-1):(i-ar)]+innov[i-ar]
  }
  return(out_vec[(ar+1):(n+ar)])
}

#4) Synthetic SWM generation

#track start time
print(paste('start',Sys.time()))

#parallel for loop
out <- foreach(m = 1:n,.combine='abnd',.packages = 'fGarch') %dopar% {
  
  #appending matrix to maintain AR continuity
  app_mat<-array(0,c(ar,sites))
  
  for(i in 1:length(yr_seq)){
    seas<-which(ix2$mon==(mo_seq[i]-1) & ix2$year==yr_seq[i])
    sm<-sim_cbias[seas,]
    syn_resid_mat<-array(NA,c(length(seas),sites))
    ar_mat<-array(NA,dim(app_mat))
    
    for(j in 1:sites){
      ar_c<-ar_coefs[[mo_seq[i]]][[j]]
      #fitted SGED to generate new random residuals
      ats<-rsged(length(seas),mean=gl_par_arr[mo_seq[i],j,1],sd=gl_par_arr[mo_seq[i],j,2],nu=gl_par_arr[mo_seq[i],j,3],xi=gl_par_arr[mo_seq[i],j,4])
      #simulate from AR model individually by site
      if(incl_mean==F){
        ar_out<-arima.sim(n=length(seas),list(ar=ar_c$coef),innov = ats,n.start=ar,start.innov=app_mat[,j])}
        #ar_out<-arsim(length(seas),ar,ar_c$coef,app_mat[,j],ats)}
      
      if(incl_mean==T){
        ar_out<-arima.sim(n=length(seas),list(ar=ar_c$coef[1:ar],mean=ar_c$coef[ar+1]),innov = ats,n.start=ar,start.innov=app_mat[,j])}
      #de-normalize errors via fitted heteroscedastic model
      norm_val<-predict(norm_fit[[mo_seq[i]]][[j]],sm[,j])
      norm_val[norm_val<=0]<-sd_arr[mo_seq[i],j]
      ar_mat[,j]<-tail(ar_out,ar)
      syn_resid_mat[,j]<-ar_out * norm_val
    }
    
    syn_resid[seas,]<-syn_resid_mat
    syn_flow[seas,]<-sim_cbias[seas,] - syn_resid_mat
    app_mat<-ar_mat
  }
  return(cbind(syn_flow,syn_resid))
}

#5). Configure and save results
for(i in 1:n){
  syn_flow_out[i,,]<-out[,1:sites,i]
  syn_resid_out[i,,]<-out[,(sites+1):(sites*2),i]
}

#compare maximum synthetic and observed flows
print(max(obs))
print(max(syn_flow_out))

#ensure all synthetic values > 0
syn_flow_out[syn_flow_out<0]<-0

saveRDS(syn_flow_out,paste('out/syn_flow_prinorm-ls_ind_',arfit,'_ar',ar,'.rds',sep=''))
saveRDS(syn_resid_out,paste('out/syn_resid_prinorm-ls_ind_',arfit,'_ar',ar,'.rds',sep=''))
#saveRDS(syn_flow_out,paste('out/syn_flow_prinorm-ls_ind_',arfit,'_var-resid.rds',sep=''))
#saveRDS(syn_resid_out,paste('out/syn_resid_prinorm-ls_ind_',arfit,'_var-resid.rds',sep=''))

print(paste('end',Sys.time()))

rm(list=ls());gc()

###################################END######################################
