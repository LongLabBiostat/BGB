library(Matrix)

## toy example with ber 

H=3
n=90
L=3

p1=p2=10
p3=10
p = p1+p2+p3
data_dim = c(p1,p2,p3)
ind_s = c(1,p1+1,p2+p1+1,p1+p2+p3)
ind_e = c(p1,p1+p2,p1+p2+p3)

# non-overlapping 
save_data = list()

W = matrix(0,p,L)
for (l in 1:L) {
  nz = 10
  sp = nz*(l-1)+1
  W[sp:(sp+nz-1),l] =  rnorm(nz,2,0.5)
}

heatmap(W, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

Z = matrix(0,nrow=n,ncol=L)
for (l in 1:L) {
  nz = 30
  sp = nz*(l-1)+1
  Z[sp:(sp+nz-1),l] =  rnorm(nz,2,0.5)
}
mu = W %*% t(Z)
heatmap(mu, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))


save_data = list()
for (mm in 1:10) {
  X = matrix(0,nrow=p,ncol=n)
  for (j in 1:p) {
    for (i in 1:n) {
      X[j,i] = rbinom(1,1,1/(1+exp(-mu[j,i])))
    }
  }
  data = list("Z" = Z, "W" = W,"X" = X,"trials"= matrix(1,nrow=p,ncol=n))
  save_data[[mm]] = data
}

saveRDS(save_data,paste0('save_data.rds'))

p =30
H=3
GG= list()
GG[[1]] = matrix(0, nrow = p,ncol = p)
pathway_list=list()
pathway_list[[1]] = rep(5,2)
pathway_list[[2]] = rep(5,2)
pathway_list[[3]] = rep(5,2)

#GG[[2]] = working_graph(x=2,pathway_list=pathway_list,H=H,data_dim=data_dim,ind_s=ind_s,ind_e=ind_e)

heatmap(GG[[2]] , Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

graph = GG[[2]]
for (h in 1:H) {
  graph_tmp= graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]]
  pathway_tmp = pathway_list[[h]]
  node = c(1)
  for (xx in 1:(length(pathway_tmp)-1)) {
    node = c(node,1+sum(pathway_tmp[1:xx]))
  }
  
  for (nn in 1:length(node)) {
    r_ind = node[nn]
    c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn]),1)
    
    # across pathway 
    if(sum(pathway_tmp[1:nn])+1 < sum(pathway_tmp)){
      graph_tmp[node[nn]:sum(pathway_tmp[1:nn]),(sum(pathway_tmp[1:nn])+1):sum(pathway_tmp)] = rbinom(pathway_tmp[nn]*(sum(pathway_tmp)-sum(pathway_tmp[1:nn])),1,0.1) 
    }
  }
  graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
}
graph = as.matrix(forceSymmetric(graph,uplo="U"))
diag(graph)=0

GG[[3]] = graph


sum(GG[[2]])
sum(GG[[3]])
saveRDS(GG,'toyber_graph_data.rds')

#########################################
# High-dimentional

H=3
n=90
L=6

p1=p2=p3=100
p = p1+p2+p3
data_dim = c(p1,p2,p3)
ind_s = c(1,p1+1,p2+p1+1,p1+p2+p3)
ind_e = c(p1,p1+p2,p1+p2+p3)

# non-overlapping 
save_data = list()

W = matrix(0,p,L)
for (l in 1:L) {
  nz = 50
  sp = nz*(l-1)+1
  W[sp:(sp+nz-1),l] =  rnorm(nz,2,0.5)*sample(c(-1,1),1,size = nz)
}

heatmap(W, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

Z = matrix(0,nrow=n,ncol=L)
for (l in 1:L) {
  nz = 15
  sp = nz*(l-1)+1
  Z[sp:(sp+nz-1),l] =  rnorm(nz,2,0.5)*sample(c(-1,1),1,size = nz)
}
mu = W %*% t(Z)
heatmap(mu, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

# signal to noise ratio
sqrt((norm(mu,type="F"))^2/(n*p))


# overlap=15
op = 15
v = 0.1
W = matrix(0,p,L)

for (l in 1:L) {
  if(l != L){
#    nz = 50 # for gaussian
    nz = 63 # for binomial 
    sp = (nz-op)*(l-1)+1 
    W[sp:(sp+nz-1),l] =  rnorm(nz,2,v)*sample(c(-1,1),1,size = nz)
  }else{
    sp = (nz-op)*(l-1)+1 
    W[sp:p,l] =  rnorm(p-sp+1,2,v)*sample(c(-1,1),1,size = p-sp+1)
  }
}

heatmap(abs(W), Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

# if Z no overlap 

Z = matrix(0,nrow=n,ncol=L)
for (l in 1:L) {
  nz = 15
  sp = nz*(l-1)+1
  Z[sp:(sp+nz-1),l] =  rnorm(nz,2,v)*sample(c(-1,1),1,size = nz)
}

heatmap(Z, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

mu = W %*% t(Z)
heatmap(mu, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

sqrt((norm(mu,type="F"))^2/(n*p))


##################


# Gaussian: non_overlap
for (mm in 1:100) {
  # signal to noise ratio around 1~1.3
  X = mu + matrix(rnorm(p*n,0,1.7),nrow=p,ncol=n)
  data = list("Z" = Z, "W" = W,"X" = X)
  save_data[[mm]] = data
}
heatmap(X, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

saveRDS(save_data,paste0('save_data.rds'))


# Gaussian: overlap=15
for (mm in 1:100) {
  # signal to noise ratio around 1~1.3
  X = mu + matrix(rnorm(p*n,0,1.8),nrow=p,ncol=n)
  data = list("Z" = Z, "W" = W,"X" = X)
  save_data[[mm]] = data
}
heatmap(X, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

saveRDS(save_data,paste0('save_data.rds'))


# # Ber 
# W = matrix(0,p,L)
# for (l in 1:L) {
#   nz = 100
#   sp = nz*(l-1)+1
#   W[sp:(sp+nz-1),l] =  rnorm(nz,2,0.5)*sample(c(-1,1),1,size = nz)
# }
# 
# heatmap(W, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))
# 
# Z = matrix(0,nrow=n,ncol=L)
# for (l in 1:L) {
#   nz = 30
#   sp = nz*(l-1)+1
#   Z[sp:(sp+nz-1),l] = rnorm(nz,2,0.5)*sample(c(-1,1),1,size = nz)
# }
# mu = W %*% t(Z)
# heatmap(mu, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))
# save_data = list()
# for (mm in 1:1) {
#   X = matrix(0,nrow=p,ncol=n)
#   for (j in 1:p) {
#     for (i in 1:n) {
#       X[j,i] = rbinom(1,1,1/(1+exp(-mu[j,i])))
#     }
#   }
#   data = list("Z" = Z, "W" = W,"X" = X,"trials"= matrix(1,nrow=p,ncol=n))
#   save_data[[mm]] = data
# }
# 
# 
# 
# 
# saveRDS(save_data,paste0('save_data.rds'))


## Binomial 

H=3
n=90
L=6

p1=p2=p3=100
p = p1+p2+p3
data_dim = c(p1,p2,p3)
ind_s = c(1,p1+1,p2+p1+1,p1+p2+p3)
ind_e = c(p1,p1+p2,p1+p2+p3)

# non-overlapping 
#save_data = list()
save_data = readRDS("E:/Post_doc/aim_2a/aim_2a/simulation_data/high/non_overlap_bi/save_data.rds")

W = matrix(0,p,L)
for (l in 1:L) {
  nz = 50
  sp = nz*(l-1)+1
  W[sp:(sp+nz-1),l] =  rnorm(nz,2,0.1)*sample(c(-1,1),1,size = nz)
}

#heatmap(W, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

Z = matrix(0,nrow=n,ncol=L)
for (l in 1:L) {
  nz = 15
  sp = nz*(l-1)+1
  Z[sp:(sp+nz-1),l] =  rnorm(nz,2,0.1)*sample(c(-1,1),1,size = nz)
}
#mu = W %*% t(Z)
#heatmap(mu, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))


## how to generate trial 
trials = matrix(0,nrow=p,ncol=n)

for (j in 1:p) {
  set.seed(1+j)
  trials[j,]= sample(seq(5,10,1),size=1)
}
save_data = list()
for (mm in 61:80) {
  # signal to noise ratio around 1~1.3
  X = matrix(0,nrow=p,ncol=n)
  for (j in 1:p) {
    for (i in 1:n) {
      X[j,i] = rbinom(1,trials[j,i],1/(1+exp(-mu[j,i])))
    }
  }
  data = list("Z" = Z, "W" = W,"X" = X,"trials"= trials)
  save_data[[mm]] = data
}

heatmap(X, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

saveRDS(save_data,paste0('save_data.rds'))

haha = readRDS("E:/Post_doc/aim_2a/aim_2a/simulation_data/high/non_overlap_bi/save_data.rds")
heatmap(haha[[1]]$X, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

class(save_data[[61]]$Z)

# overlapping

op = 15
v = 0.1
W = matrix(0,p,L)

for (l in 1:L) {
  if(l != L){
    #    nz = 50 # for gaussian
    nz = 63 # for binomial 
    sp = (nz-op)*(l-1)+1 
    W[sp:(sp+nz-1),l] =  rnorm(nz,2,v)*sample(c(-1,1),1,size = nz)
  }else{
    sp = (nz-op)*(l-1)+1 
    W[sp:p,l] =  rnorm(p-sp+1,2,v)*sample(c(-1,1),1,size = p-sp+1)
  }
}

#heatmap(abs(W), Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

# if Z no overlap 

Z = matrix(0,nrow=n,ncol=L)
for (l in 1:L) {
  nz = 15
  sp = nz*(l-1)+1
  Z[sp:(sp+nz-1),l] =  rnorm(nz,2,v)*sample(c(-1,1),1,size = nz)
}

#heatmap(Z, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

mu = W %*% t(Z)
heatmap(mu, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

trials = matrix(0,nrow=p,ncol=n)
for (j in 1:p) {
  set.seed(1+j)
  trials[j,]= sample(seq(10,15,1),size=1)
}
#save_data = list()
save_data = readRDS("E:/Post_doc/aim_2a/aim_2a/simulation_data/high/overlap2_bi/save_data.rds")
for (mm in 81:100) {
  # signal to noise ratio around 1~1.3
  X = matrix(0,nrow=p,ncol=n)
  for (j in 1:p) {
    for (i in 1:n) {
      X[j,i] = rbinom(1,trials[j,i],1/(1+exp(-mu[j,i])))
    }
  }
  data = list("Z" = Z, "W" = W,"X" = X,"trials"= trials)
  save_data[[mm]] = data
}

heatmap(X, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

saveRDS(save_data,paste0('save_data.rds'))


# nz=c(w_nz,z_nz)
# op: overlap no 
data_generate = function(nz_wz,op,v,p,n,L,t_min,t_max,replicate_no){
  op = op
  v = v
  W = matrix(0,p,L)
  for (l in 1:L) {
    if(l != L){
      #    nz = 50 # for gaussian
      nz = nz_wz[1] # for binomial 
      sp = (nz-op)*(l-1)+1 
      W[sp:(sp+nz-1),l] =  rnorm(nz,2,v)*sample(c(-1,1),1,size = nz)
    }else{
      sp = (nz-op)*(l-1)+1 
      W[sp:p,l] =  rnorm(p-sp+1,2,v)*sample(c(-1,1),1,size = p-sp+1)
    }
  }
  
  Z = matrix(0,nrow=n,ncol=L)
  for (l in 1:L) {
    nz = nz_wz[2]
    sp = nz*(l-1)+1
    Z[sp:(sp+nz-1),l] =  rnorm(nz,2,v)*sample(c(-1,1),1,size = nz)
  }
  
  mu = W %*% t(Z)
  
  trials = matrix(0,nrow=p,ncol=n)
  for (j in 1:p) {
    trials[j,]= sample(seq(t_min,t_max,1),size=1)
  }
  
  save_data = list()
  for (mm in 1:replicate_no) {
    # signal to noise ratio around 1~1.3
    X = matrix(0,nrow=p,ncol=n)
    for (j in 1:p) {
      for (i in 1:n) {
        X[j,i] = rbinom(1,trials[j,i],1/(1+exp(-mu[j,i])))
      }
    }
    data = list("Z" = Z, "W" = W,"X" = X,"trials"= trials)
    save_data[[mm]] = data
  }
  return(save_data)
}

toy_data = data_generate(nz_wz=c(63,15),op=15,v=0.1,p=300,n=90,L=6,t_min=10,t_max=30,replicate_no=100)

heatmap(toy_data[[2]]$X, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))
saveRDS(toy_data,'save_data.rds')

## graph generation

haha = readRDS("save_data.rds")

min(haha[[1]]$trials)
max(haha[[1]]$trials)

##############################################################

# Low-dimensional

# non-overlap 

H=3
n=100
L = 6
p1=p2=p3=30
p = p1+p2+p3
data_dim = c(p1,p2,p3)
ind_s = c(1,p1+1,p2+p1+1,p1+p2+p3)
ind_e = c(p1,p1+p2,p1+p2+p3)

# non-overlapping

W = matrix(0,p,L)
for (l in 1:L) {
  nz = 15
  sp = nz*(l-1)+1
  W[sp:(sp+nz-1),l] =  rnorm(nz,2,0.1)*sample(c(-1,1),1,size = nz)
}

heatmap(W, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

Z = matrix(0,nrow=n,ncol=L)
for (l in 1:L) {
  nz = 15
  sp = nz*(l-1)+1
  Z[sp:(sp+nz-1),l] =  rnorm(nz,2,0.1)*sample(c(-1,1),1,size = nz)
}

heatmap(Z, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))


mu = W %*% t(Z)
heatmap(mu, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))


## Gaussian
#save_data = list()
save_data = readRDS("E:/Post_doc/aim_2a/aim_2a/simulation_data/non_overlap/save_data.rds")
for (mm in 21:40) {
  # signal to noise ratio around 1~1.3
  X = mu + matrix(rnorm(p*n,0,1.5),nrow=p,ncol=n)
  #X = mu + matrix(rnorm(p*n,0,0.7),nrow=p,ncol=n)
  data = list("Z" = Z, "W" = W,"X" = X)
  save_data[[mm]] = data
}
saveRDS(save_data,paste0('save_data.rds'))


## binomial 




# bernoulli


save_data = list()
for (mm in 1:100) {
  # signal to noise ratio around 1~1.3
  X = matrix(0,nrow=p,ncol=n)
  for (j in 1:p) {
    for (i in 1:n) {
      # to make sure we have the biclustering 
      mu_itm = sample(c(mu[j,i],abs(mu[j,i])),1,prob = c(0.5,0.5))
      X[j,i] = rbinom(1,1,1/(1+exp(-mu_itm)))
    }
  }
  data = list("Z" = Z, "W" = W,"X" = X,"trials"= matrix(1,nrow=p,ncol=n))
  save_data[[mm]] = data
}
saveRDS(save_data,paste0('save_data.rds'))

heatmap(X, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))



## mix



## Binomial 
trials = matrix(0,nrow=p,ncol=n)
for (j in 1:p) {
  set.seed(1+j)
  trials[j,]= sample(seq(5,10,1),size=1)
}

for (mm in 1:100) {
  # signal to noise ratio around 1~1.3
  X = matrix(0,nrow=p,ncol=n)
  for (j in 1:p) {
    for (i in 1:n) {
      X[j,i] = rbinom(1,trials[j,i],1/(1+exp(-mu[j,i])))
    }
  }
  data = list("Z" = Z, "W" = W,"X" = X,"trials"= trials)
  save_data[[mm]] = data
}
saveRDS(save_data,paste0('save_data.rds'))

heatmap(X, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))



# overlapping: case I

save_data = list()

#overlap_size = 15
W = matrix(0,p,L)
for (l in 1:L) {
  if(l != L){
  nz = 30
  sp = 15*(l-1)+1 
  W[sp:(sp+nz-1),l] =  rnorm(nz,1,0.5)*sample(c(-1,1),1,size = nz)
  }else{
    sp = 15*(l-1)+1 
    W[sp:90,l] =  rnorm(15,1,0.5)*sample(c(-1,1),1,size = 15)
  }
}

heatmap(abs(W), Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

Z = matrix(0,nrow=n,ncol=L)
for (l in 1:L) {
  if(l != L){
  nz = 30
  sp = 15*(l-1)+1
  Z[sp:(sp+nz-1),l] =  rnorm(nz,1,0.5)*sample(c(-1,1),1,size = nz)
  }else{
    sp = 15*(l-1)+1 
    Z[sp:100,l] =  rnorm(25,1,0.5)*sample(c(-1,1),1,size = 25)
  }
}
heatmap(abs(Z), Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))


mu = W %*% t(Z)
heatmap(mu, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))


# # graph
# pathway_list=list()
# pathway_list[[1]] = c(15,15)
# pathway_list[[2]] = c(15,15)
# pathway_list[[3]] = c(15,15)
# graph = working_graph(x=2,pathway_list=pathway_list,H=H,data_dim=data_dim,ind_s=ind_s,ind_e=ind_e)
# heatmap(graph, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))


## Gaussian
# for (mm in 1:100) {
#   # signal to noise ratio around 1~1.3
#   X = mu + matrix(rnorm(p*n,0,0.7),nrow=p,ncol=n)
#   data = list("Z" = Z, "W" = W,"X" = X, "graph"=graph)
#   save_data[[mm]] = data
# }
# 
# saveRDS(save_data,paste0('save_data.rds'))

## Bernoulli 

# for (mm in 1:100) {
#   # signal to noise ratio around 1~1.3
#   X = matrix(0,nrow=p,ncol=n)
#   for (j in 1:p) {
#     for (i in 1:n) {
#       X[j,i] = rbinom(1,1,1/(1+exp(-mu[j,i])))
#     }
#   }
#   data = list("Z" = Z, "W" = W,"X" = X,"trials"= matrix(1,nrow=p,ncol=n))
#   save_data[[mm]] = data
# }
# saveRDS(save_data,paste0('save_data.rds'))

## Binomial 
trials = matrix(0,nrow=p,ncol=n)
for (j in 1:p) {
  set.seed(1+j)
  trials[j,]= sample(seq(5,10,1),size=1)
}

for (mm in 1:100) {
  # signal to noise ratio around 1~1.3
  X = matrix(0,nrow=p,ncol=n)
  for (j in 1:p) {
    for (i in 1:n) {
      X[j,i] = rbinom(1,trials[j,i],1/(1+exp(-mu[j,i])))
    }
  }
  data = list("Z" = Z, "W" = W,"X" = X,"trials"= trials)
  save_data[[mm]] = data
}
saveRDS(save_data,paste0('save_data.rds'))





# overlap case II
H=3
n=100
L = 6
p1=p2=p3=30
p = p1+p2+p3
data_dim = c(p1,p2,p3)
ind_s = c(1,p1+1,p2+p1+1,p1+p2+p3)
ind_e = c(p1,p1+p2,p1+p2+p3)



save_data = list()

#overlap_size = 5
W = matrix(0,p,L)
for (l in 1:L) {
  if(l != L){
    nz = 20
    sp = 15*(l-1)+1 
    W[sp:(sp+nz-1),l] =  rnorm(nz,2,0.1)*sample(c(-1,1),1,size = 20)
  }else{
    sp = 15*(l-1)+1 
    W[sp:90,l] =  rnorm(15,2,0.1)*sample(c(-1,1),1,size = 15)
  }
}

heatmap(abs(W), Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

# if Z no overlap 

Z = matrix(0,nrow=n,ncol=L)
for (l in 1:L) {
  nz = 15
  sp = nz*(l-1)+1
  Z[sp:(sp+nz-1),l] =  rnorm(nz,2,0.1)*sample(c(-1,1),1,size = nz)
}

# if Z has overlap
# Z = matrix(0,nrow=n,ncol=L)
# for (l in 1:L) {
#   if(l != L){
#     nz = 20
#     sp = 15*(l-1)+1
#     Z[sp:(sp+nz-1),l] =  rnorm(nz,2,0.5)*sample(c(-1,1),1,size = nz)
#   }else{
#     sp = 15*(l-1)+1 
#     Z[sp:100,l] =  rnorm(25,2,0.5)*sample(c(-1,1),1,size = 25)
#   }
# }
heatmap(abs(Z), Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))


mu = W %*% t(Z)
heatmap(mu, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))


# # graph
# pathway_list=list()
# pathway_list[[1]] = c(15,15)
# pathway_list[[2]] = c(15,15)
# pathway_list[[3]] = c(15,15)
# graph = working_graph(x=2,pathway_list=pathway_list,H=H,data_dim=data_dim,ind_s=ind_s,ind_e=ind_e)
# heatmap(graph, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))


# Gaussian:03/15
save_data = readRDS("E:/Post_doc/aim_2a/aim_2a/simulation_data/overlap2/save_data.rds")

for (mm in 21:40) {
  # signal to noise ratio around 1~1.3
  X = mu + matrix(rnorm(p*n,0,1.5),nrow=p,ncol=n)
  data = list("Z" = Z, "W" = W,"X" = X)
  save_data[[mm]] = data
}

saveRDS(save_data,paste0('save_data.rds'))

## Ber



save_data=list()
for (mm in 1:100) {
  # signal to noise ratio around 1~1.3
  X = matrix(0,nrow=p,ncol=n)
  for (j in 1:p) {
    for (i in 1:n) {
      mu_itm = sample(c(mu[j,i],abs(mu[j,i])),1,prob = c(0.5,0.5))
      X[j,i] = rbinom(1,1,1/(1+exp(-mu_itm)))
    }
  }
  data = list("Z" = Z, "W" = W,"X" = X,"trials"= matrix(1,nrow=p,ncol=n))
  save_data[[mm]] = data
}
saveRDS(save_data,paste0('save_data.rds'))

heatmap(X, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))


## Binomial 
trials = matrix(0,nrow=p,ncol=n)
for (j in 1:p) {
  set.seed(1+j)
  trials[j,]= sample(seq(5,10,1),size=1)
}

for (mm in 1:100) {
  # signal to noise ratio around 1~1.3
  X = matrix(0,nrow=p,ncol=n)
  for (j in 1:p) {
    for (i in 1:n) {
      X[j,i] = rbinom(1,trials[j,i],1/(1+exp(-mu[j,i])))
    }
  }
  data = list("Z" = Z, "W" = W,"X" = X,"trials"= trials)
  save_data[[mm]] = data
}
saveRDS(save_data,paste0('save_data.rds'))

heatmap(X, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))




## plot the data matrix
dt = readRDS("E:/Post_doc/aim_2a/aim_2a/simulation_data/overlap/save_data.rds")
dt_1 = dt[[1]]
heatmap(dt_1$W %*% t(dt_1$Z) , Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))


##

working_graph <- function(x,ol,pathway_list,H,data_dim,ind_s,ind_e) {
  graph = matrix(0,p,p)
  if(x==1){
    graph = matrix(0,p,p)
    return(graph)  
  }else if(x==2){
    for (h in 1:H) {
      graph_tmp= matrix(0,nrow=data_dim[h],ncol=data_dim[h])
      pathway_tmp = pathway_list[[h]]
      #node = c(1,1+pathway_tmp[1],1+sum(pathway_tmp[1:2]))
      node = c(1)
      for (xx in 1:(length(pathway_tmp)-1)) {
        node = c(node,1+sum(pathway_tmp[1:xx])-xx*ol)
      }
      for (nn in 1:length(node)) {
        r_ind = node[nn]
        c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn])-(nn-1)*ol,1)

        graph_tmp[r_ind,c_ind]=1
        
        graph_tmp[(node[nn]+1):(sum(pathway_tmp[1:nn])-(nn-1)*ol),(node[nn]+1):(sum(pathway_tmp[1:nn])-(nn-1)*ol)] = rbinom((pathway_tmp[nn]-1)^2,1,0.05)
      }
      graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
    }
    graph = as.matrix(forceSymmetric(graph,uplo="U"))
    diag(graph)=0
    return(graph)
  }else if(x==3){
    for (h in 1:H) {
      graph_tmp= matrix(0,nrow=data_dim[h],ncol=data_dim[h])
      pathway_tmp = pathway_list[[h]]
      #node = c(1,1+pathway_tmp[1],1+sum(pathway_tmp[1:2]))
      node = c(1)
      for (xx in 1:(length(pathway_tmp)-1)) {
        node = c(node,1+sum(pathway_tmp[1:xx]))
      }
      
      for (nn in 1:length(node)) {
        r_ind = node[nn]
        c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn]),1)
        graph_tmp[r_ind,c_ind]=1
        # within pathway 
        graph_tmp[(node[nn]+1):sum(pathway_tmp[1:nn]),(node[nn]+1):sum(pathway_tmp[1:nn])] = rbinom((pathway_tmp[nn]-1)^2,1,0.3)
        # across pathway 
        if(sum(pathway_tmp[1:nn])+1 < sum(pathway_tmp)){
          graph_tmp[node[nn]:sum(pathway_tmp[1:nn]),(sum(pathway_tmp[1:nn])+1):sum(pathway_tmp)] = rbinom(pathway_tmp[nn]*(sum(pathway_tmp)-sum(pathway_tmp[1:nn])),1,0.05) 
        }
      }
      graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
    }
    graph = as.matrix(forceSymmetric(graph,uplo="U"))
    diag(graph)=0
    return(graph)
  }
}


## generate graph 
## high 
# g1 = null 
# g2(pathway size =50) = informative
# g3(pathway size =50) = add noise edge 
# g4(pathway size =25) = informative
# g5(pathway size =25) = add noise edge 
p1=p2=p3=100
H=3
p = p1+p2+p3
data_dim = c(p1,p2,p3)
ind_s = c(1,p1+1,p2+p1+1,p1+p2+p3)
ind_e = c(p1,p1+p2,p1+p2+p3)
GG= list()
GG[[1]] = matrix(0, nrow = p,ncol = p)
pathway_list=list()

# pathway_list[[1]] = rep(50,2)
# pathway_list[[2]] = rep(50,2)
# pathway_list[[3]] = rep(50,2)

# pathway_list[[1]] = rep(25,4)
# pathway_list[[2]] = rep(25,4)
# pathway_list[[3]] = rep(25,4)







# overlap
GG= list()
GG[[1]] = matrix(0, nrow = p,ncol = p)

p1=300
H=1
L=6
p = p1
data_dim = c(p1)
ind_s = c(1,p1)
ind_e = c(p1)
ol=15
nz=50

pathway_list=list()
pathway_list[[1]] = c(rep(nz,5),p-nz*(L-1)+ol*(L-1)) 

GG[[2]] = working_graph(x=2,ol=15,pathway_list=pathway_list,H=H,data_dim=data_dim,ind_s=ind_s,ind_e=ind_e)


graph = GG[[2]]
prob_across = 0.01
for (h in 1:H) {
  graph_tmp= matrix(0,nrow=data_dim[h],ncol=data_dim[h])
  pathway_tmp = pathway_list[[h]]
  #node = c(1,1+pathway_tmp[1],1+sum(pathway_tmp[1:2]))
  node = c(1)
  for (xx in 1:(length(pathway_tmp)-1)) {
    node = c(node,1+sum(pathway_tmp[1:xx])-xx*ol)
  }
  
  for (nn in 1:length(node)) {
    r_ind = node[nn]
    c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn])-(nn-1)*ol,1)
    # within pathway 
    
    # across pathway 
    if(sum(pathway_tmp[1:nn])+1-(nn-1)*ol < sum(pathway_tmp)-(length(node)-1)*ol){
      graph_tmp[node[nn]:(sum(pathway_tmp[1:nn])-(nn)*ol),(sum(pathway_tmp[1:nn])-(nn-1)*ol+1):(sum(pathway_tmp)-(length(node)-1)*ol)] = rbinom((pathway_tmp[nn]-ol)*(sum(pathway_tmp)-(length(node)-1)*ol-(sum(pathway_tmp[1:nn])-(nn-1)*ol)),1,prob_across) 
      }
  }
  graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
}
graph = as.matrix(forceSymmetric(graph,uplo="U"))
diag(graph)=0

GG[[3]] = graph

## replace the within-pathway edge with across pathway 

whole_graph = matrix(0,300,300)
#ss = c(1,101,201)
#ee = c(100,200,300)
ss = c(1)
ee = c(300)
for (j in 1:1) {
  tmp_graph = haha[[8]]
  graph = tmp_graph[ss[j]:ee[j],ss[j]:ee[j]]
  replace_prop = 0.9
  
  ind_within = which(graph!=0,arr.ind = TRUE)
  ind_within = ind_within[ind_within[,1]>ind_within[,2],]
  rep_edge_within= sample(seq(1,dim(ind_within)[1]),floor((sum(graph)/2)*replace_prop),replace=FALSE)
  ind_within = ind_within[rep_edge_within,]

  #start with the full connect graph 
  p1=300
  H=1
  L=6
  p = p1
  data_dim = c(p1)
  ind_s = c(1,p1)
  ind_e = c(p1)
  ol=15
  nz=63
  pathway_list=list()
  pathway_list[[1]] = c(rep(nz,5),p-nz*(L-1)+ol*(L-1)) 
  graph_full_connect = working_graph(x=2,ol=15,pathway_list=pathway_list,H=H,data_dim=data_dim,ind_s=ind_s,ind_e=ind_e)
  ind_across = which(graph_full_connect==0,arr.ind = TRUE)
  ind_across = ind_across[ind_across[,1]>ind_across[,2],]
  rep_edge_across= sample(sample(seq(1,dim(ind_across)[1]),floor((sum(graph)/2)*replace_prop),replace=FALSE)) 
  ind_across = ind_across[rep_edge_across,]
  
  for(i in 1:(floor((sum(graph)/2)*replace_prop))) {
    graph[ind_within[i,1],ind_within[i,2]] = 0 
    graph[ind_across[i,1],ind_across[i,2]] = 1 
  }
  graph = as.matrix(forceSymmetric(graph,uplo="L"))
  diag(graph)=0
  whole_graph[ss[j]:ee[j],ss[j]:ee[j]] = graph 
  
}
 


## visual 
image(
  haha[[14]],
  axes = FALSE,
  col = colorRampPalette(c("white", "black"))(30),
  useRaster = TRUE,
  main = "With rasterisation"
)



## add more graph to the exist data
haha = readRDS('graph_data_high.rds')
#haha[[6]] = GG[[2]]
#haha[[7]] = GG[[2]]+ GG[[3]]
# for binomial overlap 
haha[[8]] = GG[[2]]          # no noisy 
haha[[9]] = GG[[2]]+ GG[[3]] # add noisy 
haha[[10]] = graph           # replace info with noisy
# for gaussian overlap
haha[[11]] = graph 
# for non-overlap 
haha[[12]] = whole_graph 

haha[[13]] = whole_graph 
## 

haha[[14]] = whole_graph 
haha[[15]] = whole_graph 
haha[[16]] = whole_graph 

##

haha[[17]] = whole_graph 
haha[[18]] = whole_graph 
haha[[19]] = whole_graph 

## 
haha[[20]] = whole_graph 
haha[[21]] = whole_graph 
haha[[22]] = whole_graph 

saveRDS(haha,'graph_data_high.rds')



length(haha)









#############################################################################3
## low 
p =90 

p1=p2=p3=30
p = p1+p2+p3
data_dim = c(p1,p2,p3)
ind_s = c(1,p1+1,p2+p1+1,p1+p2+p3)
ind_e = c(p1,p1+p2,p1+p2+p3)


GG= list()
GG[[1]] = matrix(0, nrow = p,ncol = p)
GG[[2]] = dt[[1]]$graph

pathway_list=list()
pathway_list[[1]] = c(15,15)
pathway_list[[2]] = c(15,15)
pathway_list[[3]] = c(15,15)
GG[[3]] = working_graph(x=3,pathway_list=pathway_list,H=H,data_dim=data_dim,ind_s=ind_s,ind_e=ind_e)

heatmap(GG[[3]], Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

saveRDS(GG,'graph_data.rds')
GG=readRDS('graph_data.rds')

#####################################

H=3

graph = GG[[2]]
for (h in 1:H) {
  graph_tmp= graph[(1+30*(h-1)):(30*h),(1+30*(h-1)):(30*h)]
  pathway_tmp = pathway_list[[h]]
  node = c(1)
  for (xx in 1:(length(pathway_tmp)-1)) {
    node = c(node,1+sum(pathway_tmp[1:xx]))
  }

  for (nn in 1:length(node)) {
    r_ind = node[nn]
    c_ind = seq(node[nn]+1,sum(pathway_tmp[1:nn]),1)
  
    # across pathway 
    if(sum(pathway_tmp[1:nn])+1 < sum(pathway_tmp)){
      graph_tmp[node[nn]:sum(pathway_tmp[1:nn]),(sum(pathway_tmp[1:nn])+1):sum(pathway_tmp)] = rbinom(pathway_tmp[nn]*(sum(pathway_tmp)-sum(pathway_tmp[1:nn])),1,0.1) 
    }
  }
  graph[ind_s[h]:ind_e[h],ind_s[h]:ind_e[h]] = graph_tmp
}
graph = as.matrix(forceSymmetric(graph,uplo="U"))
diag(graph)=0


image(
  graph_tmp,
  axes = FALSE,
  col = colorRampPalette(c("white", "black"))(30),
  useRaster = TRUE,
  main = "With rasterisation"
)




heatmap(graph, Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

identical(graph[16:30, 16:30], GG[[2]][16:30, 16:30]) 

sum(graph)-sum(GG[[4]])




## Based on GG[[2]], replace within-pathway edges with across pathway 
GG[[3]] = graph
##Based on GG[[2]], add across-pathway edges with p=0.05 
GG[[4]] = graph
##Based on GG[[2]], add across-pathway edges with p=0.1 
GG[[5]] = graph

saveRDS(GG,'graph_data.rds')



trial_name_2 = paste0('mcmc_result_2')
bic_2 = Inf
if (trial_result$bic_2 <= bic_2){
  bic_2 = trial_result$bic_2
  data_tbs[[trial_name_2]] = trial_result
}



heatmap(GG[[4]], Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))
heatmap(GG[[2]], Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))
heatmap(GG[[3]], Colv = NA, Rowv = NA,scale=NULL,col = cm.colors(256))

sum(GG[[4]])

sum(GG[[3]])

sum(GG[[2]])



####

w_list=list()
z_list=list()
x_list= list()
for (i in 1:100) {
  w_base_1=matrix(0,90,9)
  w_base_2=matrix(0,90,9)
  w_base_3=matrix(0,90,9)
  
  for (l in 1:9) {
    nz = 10
    sp = nz*(l-1)+1
    a = 1-rbinom(nz,1,0.5)
    w_base_1[sp:(sp+nz-1),l] =  a*rnorm(nz,1.5,0.1)+ (1-a)*rnorm(nz,7.5,0.1)  
  }
  
  for (l in 1:9) {
    nz = 10
    sp = nz*(l-1)+1
    a = 1-rbinom(nz,1,0.5)
    w_base_2[sp:(sp+nz-1),l] =  a*rnorm(nz,1.5,0.1)+ (1-a)*rnorm(nz,7.5,0.1)  
  }
  
  for (l in 1:9) {
    nz = 10
    sp = nz*(l-1)+1
    a = 1-rbinom(nz,1,0.5)
    w_base_3[sp:(sp+nz-1),l] =  a*rnorm(nz,1.5,0.1)+ (1-a)*rnorm(nz,7.5,0.1)  
  }
  
  true_W = rbind(w_base_1,w_base_2,w_base_3)
  
  true_Z=matrix(rnorm(200*9,0,1),9,200)
  X = true_W %*% true_Z +rnorm(200*270,0,1)
  w_list[[i]]= true_W
  z_list[[i]]= true_Z
  x_list[[i]]= X
}



