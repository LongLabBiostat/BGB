#source("E:/Post_doc/aim_2a/aim_2a/function/aim_2a.R")
source("../BGB_mcmcfunction.R")


###### 
d_ind = 1
dt_olp_ind = 1 # 

dt_olp_mt = matrix(0,nrow=8,ncol=2)
dt_olp_mt[,1]=rep(seq(1,4),2)
dt_olp_mt[,2]=c(1,1,1,1,2,2,2,2)

dt= dt_olp_mt[dt_olp_ind,1]
olp = dt_olp_mt[dt_olp_ind,2]

## data type
suffix_dt = list()
suffix_dt[[1]] = '_gau'
suffix_dt[[2]] = '_ber'
suffix_dt[[3]] = '_bi'
suffix_dt[[4]] = '_mix'

# overlap pattern
suffix_olp = list()
suffix_olp[[1]] = 'non_overlap'
suffix_olp[[2]] = 'overlap2'


suffix_graph = list()
suffix_graph[[1]] = "g_1"
suffix_graph[[2]] = "g_2"
suffix_graph[[3]] = "g_3"
suffix_graph[[4]] = "g_4"
suffix_graph[[5]] = "g_5"

for (mm in (2*d_ind-1):(2*d_ind)) {
for (gg in c(1,2,5)) {
    arg = mm
    bic = rep(Inf,6)
    set.seed(19402+arg)
    n = 100
    T = 1500
    start = 400
    end = T
    p1 = 30
    p2 = 30
    p3 = 30
    p = p1+p2+p3 
    data_dim = c(p1,p2,p3)
    H=3
    ind_s = c(1,p1+1,p2+p1+1,p1+p2+p3)
    ind_e = c(p1,p1+p2,p1+p2+p3)
    
    graph_ind = gg
    ## Data Generating
    save_data_mask = readRDS("../save_data.rds")
    data_1 = save_data_mask[[mm]]
    W = data_1$W
    Z = t(data_1$Z)
    X = data_1$X
    L = dim(W)[2]
    #grid_L  = c(L-1,L,L+1)
    grid_L  = c(L-2,L-1,L,L+1,L+2)
    #grid_L  = c(L-3,L-2,L-1,L,L+1,L+2,L+3,L+4,L+5,L+6,L+7)
    #grid_L  = c(L)
    
    # X: p x n matrix 
    mu = W %*% Z
    if(dt==1){
      trials = matrix(0,nrow = p, ncol = n)  
    }else{
      trials = data_1$trials
    }
    
    data_tbs = list('true_W'= W,'true_Z'=Z,'true_X'=X)
    #####################################################
    
    Sigma = rep(0.01,p)  # variance for each entry of m
    
    Q = diag(4,nrow = p,ncol = p) # proposal density
    
    ######## working graph####################################
    GG = readRDS("../graph_data.rds")
    graph = GG[[gg]]
    data_tbs[['G']] = graph 
    ########################################################
    
    w_ini_l = list()
    z_ini_l = list()
    L_seed = 1113
    for (l in 1:length(grid_L)) {
      set.seed(l%%3+L_seed)
      LL= grid_L[l]
      w_ini_l[[l]]=matrix(rnorm(LL*p,0,1),nrow = p,ncol = LL)
    }
    
    for (l in 1:length(grid_L)) {
      set.seed(l%%3+L_seed)
      LL= grid_L[l]
      z_ini_l[[l]] = matrix(rnorm(LL*n,0,1),nrow=LL,ncol=n)
    }
    
    ## initialization 
    
    grid_nu_1 = c(-1,0,1)   # mean 
    #  grid_nu_1 = c(1)   # mean 
    
    grid_nu_2 = c(0.5) # variance 
    #  grid_nu_2 = c(0.5) # variance 
    
    eta = 15
    eps = 0.1
    
    grid_nu_4 = c(1,2,8)
    nu_3 = 1
    
    for (l in 1:length(grid_L)) {
      L = grid_L[l]
      rho_ini = matrix(1,nrow=p,ncol=n )
      
      tau_temp = matrix(1,nrow=p,ncol=L )
      tau_ini  = tau_temp
      
      tauz_ini = matrix(1,nrow=L,ncol=n )
      
      # W and Z are initialized with the true value
      #w_ini = W 
      #z_ini = Z 
      
      # with wrong value 
      w_ini = w_ini_l[[l]] 
      z_ini = z_ini_l[[l]] 
      
      ## omega should be positive definite
      omega_ini = diag(1,nrow=p,ncol=p)  # compatible with graph
      
    
      ind = which(graph!=0,arr.ind = TRUE)  # ind of the nonzero elements
      
      if(dim(ind)[1]!=0){
        for (i in 1:dim(ind)[1]) {
          if(ind[i,1]>ind[i,2]){
            r_ind = ind[i,1]
            c_ind = ind[i,2]
            omega_ini[r_ind,c_ind] = 0.05
          }
        }
      }
      
      
      omega_temp = as.matrix(forceSymmetric(omega_ini,uplo = 'L'))
      
      inv_omega_temp = solve(omega_temp) 
      
      for (j in 1:length(grid_nu_1)) {
        for (i in 1:length(grid_nu_2)) {
          for (k in 1:length(grid_nu_4)) {
            nu_1 = grid_nu_1[j]
            nu_2 = grid_nu_2[i]
            nu_4 = grid_nu_4[k]
            alpha_ini = matrix(nu_1,nrow=p,ncol=L )
            
            xiz_ini = matrix(nu_4/nu_3,nrow=L,ncol=n )
            
            phi_ini = matrix(1,nrow=H,ncol=L)
            mcmc_box = empty_chain(n=n,p=p,T=T,L=L,rho_ini=rho_ini,phi_ini=phi_ini,w_ini=w_ini,z_ini=z_ini,alpha_ini=alpha_ini,tauz_ini=tauz_ini,xiz_ini=xiz_ini)
            
            chain_rho = mcmc_box$chain_rho
            chain_w = mcmc_box$chain_w
            chain_z = mcmc_box$chain_z
            chain_tauz = mcmc_box$chain_tauz
            chain_xiz = mcmc_box$chain_xiz
            chain_alpha = mcmc_box$chain_alpha
            chain_m = mcmc_box$chain_m
            chain_tau = mcmc_box$chain_tau
            chain_phi = mcmc_box$chain_phi
            
            
            w_temp   = matrix(as.numeric(chain_w[1,]),nrow=p ,ncol=L,byrow=T)
            z_temp   = matrix(as.numeric(chain_z[1,]),nrow=L ,ncol=n,byrow=T)
            tauz_temp =matrix(as.numeric(chain_tauz[1,]),nrow=L ,ncol=n,byrow=T)
            xiz_temp  =matrix(as.numeric(chain_xiz[1,]),nrow=L ,ncol=n,byrow=T)
            rho_temp   = matrix(as.numeric(chain_rho[1,]),nrow=p ,ncol=n,byrow=T)
            alpha_temp = matrix(as.numeric(chain_alpha[1,]),nrow=p ,ncol=L,byrow=T)
            
            m_temp = rep(0,p)
            trial_result = BFGA_MCMC(dt=dt-1,T=T,L=L,p=p,n=n,X=X,trials=trials,nu_1=nu_1,nu_2=nu_2,nu_3=nu_3,nu_4=nu_4,Sigma=Sigma,Q=Q,eta=eta,eps=eps,rho_temp=rho_temp,tau_temp=tau_temp,omega_temp=omega_temp,inv_omega_temp=inv_omega_temp,w_temp=w_temp,z_temp=z_temp,alpha_temp=alpha_temp,m_temp=m_temp,tauz_temp=tauz_temp,xiz_temp=xiz_temp,start=start,end=end,mu=mu)       
            
            
            trial_name = c('mcmc_result_1','mcmc_result_2','mcmc_result_3','mcmc_result_4','mcmc_result_5','mcmc_result_6')
            bic_list = c(trial_result$bic_11,trial_result$bic_12,trial_result$dic_11,trial_result$dic_12,trial_result$waic_11,trial_result$waic_12)
            for (ii in 1:6) {
              if ( bic_list[ii]<= bic[ii]){
                bic[ii] = bic_list[ii]
                data_tbs[[trial_name[ii]]] = trial_result
              }
              
            }
            
            #data_tbs[[trial_name]] = trial_result
            #tun_arg = tun_arg+1 
          }
        }
      }
    }
    saveRDS(data_tbs,paste0("/project/LongGroup/qiyiwen/aim_2a/low/",suffix_olp[[olp]],suffix_dt[[dt]],"/",suffix_graph[[graph_ind]],"_",arg,".rds"))  
  }
}
