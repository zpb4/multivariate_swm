#Generate VAR correlated LRM model with bootstrapping
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
arfit<-'ols-cbias' #penalized or ols version of VAR model and residuals
sites <- 14
incl_mean<-F
bvar<-F

#date/time indices
ix<-seq(as.Date('1988-10-01'),as.Date('2013-09-30'),'day')
ix2<-as.POSIXlt(ix)

#Load required data
syn_resid<-array(NA,c(length(ix2),sites))
syn_flow<-array(NA,c(length(ix2),sites))
syn_empcop<-array(NA,c(length(ix2),sites))
syn_flow_out<-array(NA,c(n,length(ix2),sites))
syn_resid_out<-array(NA,c(n,length(ix2),sites))
syn_empcop_out<-array(NA,c(n,length(ix2),sites))
at_lst<-readRDS(paste('fit/at_lst_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))
gl_par_arr<-readRDS(paste('fit/gl_par_arr_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))
var_coefs<-readRDS(paste('fit/var_coefs_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))
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

#3) PCA to use for kNN sample
sim_pc_fit<-prcomp(sim_cbias)

for(i in 1:14){
  print(paste(sum(sim_pc_fit$sdev[1:i])/sum(sim_pc_fit$sdev)*100))
}

#6 PCs explain 80% of variance
sim_pc<-sim_pc_fit$x[,1:6]
samp_pc<-sim_pc_fit$x[,1:6] #use same sampling matrix historical sampling

#4) Set up date/time indices for sequential (by month) VAR sampling
yr_idx<-ix2$year
yr_idx_lst<-vector('list',12)

for(i in 1:12){
  seas<-which(ix2$mon==(i-1))
  yr_idx_lst[[i]]<-yr_idx[seas]
}

yr_seq<-c(rep(88,3),rep(89:112,each=12),rep(113,9))
mo_seq<-c(10:12,rep(1:12,length(89:112)),1:9)

#4a. define results binding function
abnd<-function(x,y){abind(x,y,along=3)}


#5) Synthetic SWM generation
#define knn parameters
knn<-round(sqrt(length(seas)))
tot<-sum(rep(1,knn) / 1:knn)
wts<-rep(1,knn) / 1:knn / rep(tot,knn)

#track start time
print(paste('start',Sys.time()))

#parallel for loop
out <- foreach(m = 1:n,.combine='abnd',.inorder=F,.packages = 'fGarch') %dopar% {
  
  #5a. kNN sampling of simulation indices
  knn_lst<-vector('list',12)
  
  for(i in 1:12){
    knn_vec<-c()
    seas<-which(ix2$mon==(i-1))
    n_sim<-sim_pc[seas,]
    smp_sim<-samp_pc[seas,]
    for(j in 1:length(n_sim[,1])){
      y<-sqrt(apply((n_sim - smp_sim)^2,1,sum)) #find NEP closest by Euclidean distance
      x<-sort(y)
      x<-x[1:knn] #find kNN
      s<-sample(x,1,prob=wts)
      id<-which(y==s)
      if(length(id)>1) {id <- sample(id,1)} #resample the sample for any duplicated values
      knn_vec[j]<-id
    }
    knn_lst[[i]]<-knn_vec
  }

  #5b. Create synthetic empirical copula matching rank structure of empirical data
  syn_cop<-vector('list',12)
  
  for(i in 1:12){
    seas<-which(ix2$mon==(i-1))
    mat<-at_lst[[i]]
    ecop<-apply(mat,2,function(x){rank(x,ties.method = 'random')})
    syn_ecop_mat<-array(NA,c(dim(mat)[1],sites))
    
    #empirical copula is populated with random generation from fitted SGED
    for(j in 1:sites){
      syn_at<-rsged(dim(mat)[1],mean=gl_par_arr[i,j,1],sd=gl_par_arr[i,j,2],nu=gl_par_arr[i,j,3],xi=gl_par_arr[i,j,4])
      
      #rank order randomly generated values to match empirical data
      r_syn_at<-rank(syn_at,ties.method = 'random')
      for(k in 1:length(r_syn_at)){
        syn_ecop_mat[k,j]<-syn_at[which(r_syn_at==ecop[k,j])]
      }
    }
    syn_cop[[i]]<-syn_ecop_mat
  }

  #5c. Use indices from 5a to sample from 5b
  syn_cop_knn<-vector('list',12)
  for(i in 1:12){
    id<-knn_lst[[i]]
    syn_cop_knn[[i]]<-syn_cop[[i]][id,]
  }
  
  #syn_empcop[[m]]<-syn_cop_knn #can be used for correlation analysis
    
  #starting matrix to append to VAR simulation
  app_mat<-syn_cop_knn[[mo_seq[1]]][1:ar,]
    
  #5d. Primary simulation script
  for(i in 1:length(yr_seq)){ #loops in sequence through months in timespan to maintain VAR/AR continuity
    mat2<-syn_cop_knn[[mo_seq[i]]][which(yr_idx_lst[[mo_seq[i]]]==yr_seq[i]),]
    coeff<-var_coefs[mo_seq[i],,]
    seas<-which(ix2$mon==(mo_seq[i]-1) & ix2$year==yr_seq[i])
    syn_empcop[seas,]<-mat2
    sm<-sim_cbias[seas,]
    syn_resid_mat<-matrix(0,ncol=sites,nrow=(dim(mat2)[1]+ar))
    syn_var_mat<-matrix(0,ncol=sites,nrow=(dim(mat2)[1]+ar))
    syn_var_mat[1:ar,]<-app_mat
    
    #Simulate from VAR model using synthetically generated residuals
    for(k in (ar+1):(dim(mat2)[1]+ar)){
      if(incl_mean==T & bvar==F){
        bvar_res<-c(t(syn_var_mat[(k-1):(k-ar),]),1) %*% t(coeff) + mat2[(k-ar),]}
      if(incl_mean==T & bvar==T){
        bvar_res<-c(1,t(syn_var_mat[(k-1):(k-ar),])) %*% t(coeff) + mat2[(k-ar),]}
      if(incl_mean==F & bvar==T){
        bvar_res<-c(t(syn_var_mat[(k-1):(k-ar),])) %*% t(coeff[,2:(sites*ar+1)]) + mat2[(k-ar),]}
      if(incl_mean==F & bvar==F){
        bvar_res<-c(t(syn_var_mat[(k-1):(k-ar),])) %*% t(coeff) + mat2[(k-ar),]}
        
      s<-sm[(k-ar),]
      var_res<-bvar_res[1,]
      syn_var_mat[k,]<-var_res #VAR matrix maintained separately from final residual matrix for continuity
      norm_val<-c()
      for(j in 1:sites){
        norm_val[j]<-predict(norm_fit[[mo_seq[i]]][[j]],s[j])
      }
      res<-var_res*norm_val
      sub_zero_idx<-which((s-res)<0)
      res[sub_zero_idx]<-s[sub_zero_idx]
      syn_resid_mat[k,]<-res
    }
  syn_resid[seas,]<-syn_resid_mat[(ar+1):k,]
  syn_flow[seas,]<-sim_cbias[seas,] - syn_resid_mat[(ar+1):k,]
  app_mat<-syn_var_mat[(k-ar+1):k,]
  }
  return(cbind(syn_flow,syn_resid,syn_empcop))
  print(paste(m,Sys.time()))
}

print(paste('end',Sys.time()))

#6) Reconfigure and save output
for(i in 1:n){
  syn_flow_out[i,,]<-out[,1:sites,i]
  syn_resid_out[i,,]<-out[,(sites+1):(sites*2),i]
  syn_empcop_out[i,,]<-out[,(sites*2+1):(sites*3),i]
}

#compare max synthetic flows to maximum observations
print(max(obs))
print(max(syn_flow_out))
print(min(syn_flow_out))

#ensure no synthetic flows < 0
syn_flow_out[syn_flow_out<0]<-0

saveRDS(syn_flow_out,paste('out/syn_flow_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))
saveRDS(syn_resid_out,paste('out/syn_resid_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))
saveRDS(syn_empcop_out,paste('out/syn_empcop_prinorm-ls_',arfit,'_ar',ar,'.rds',sep=''))

rm(list=ls());gc()

###################################END######################################
