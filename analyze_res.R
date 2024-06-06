
source("../analyze_function.R")
source("../gbcmetric.R")

dt= readRDS("../mcmc_chain.rds")

# generate w_est
chain_w = dt$mcmc_13$chain_w
p = dim(dt$true_X)[1]
L = (dim(chain_w)[2])/p
n = dim(dt$true_X)[2]
start = 1 
end = dim((chain_w))[1]
interval_w = apply(chain_w[start:end,],2,sample_quantile)
df_w = sum(interval_w[1,]*interval_w[2,]>0)
w_est= apply(chain_w[start:end,],2,mean)
zero_index_w = which(interval_w[1,]*interval_w[2,]<0,arr.ind = TRUE)
w_est[zero_index_w] = 0
w_est=matrix(w_est,nrow=p,ncol=L,byrow=TRUE)


####
true_L = 3
truth = Bi_obj(W=true_W,Z=t(true_Z))
S2= list()
for (i in 1:L) {
  tmp= list() 
  tmp[['r']] = which(truth@GenesMembership[,i]==TRUE) 
  tmp[['c']] = 1 
  S2[[i]] = tmp
}

fpr = 0
fnr =0
for (w_ind in 1:length(w_list)) {
  w_est = w_list[[w_ind]]
  L = dim(w_est)[2]
  est   = Bi_obj(W=w_est,Z=t(z_est))
  S1= list()
  for (i in 1:L) {
    tmp= list() 
    tmp[['r']] = which(est@GenesMembership[,i]==TRUE) 
    tmp[['c']] = 1 
    S1[[i]] = tmp
  }
  re_metric = gbcmetric(S1,S2,p,n,clus_dir=0)
  fpr = fpr + re_metric$FP_CE/p
  fnr = fnr + re_metric$FN_CE/(p*L-p)
}

fpr/100
fnr/100




