install.packages("superbiclust")
library(superbiclust)
library(fabia)

Bi_obj = function(W,Z){
  data_ini = fabia(matrix(rnorm(500,0,1),nrow=5,ncol=100),p=3,alpha=0.01,center=0,norm=0,scale=0.0)
  bi_ini= BiclustSet(data_ini)
  # change W Z and to boolean 
  p = dim(W)[1] 
  L = dim(W)[2]
  n = dim(Z)[1]
  boolean_W = matrix(as.logical(W),nrow = p, ncol=L)
  boolean_Z = matrix(as.logical(Z),nrow = n, ncol=L)
  bi_ini@GenesMembership =  boolean_W
  bi_ini@ColumnMembership = t(boolean_Z)
  bi_ini@Number =  L
  return(bi_ini)
}

##  construct object for gbcmetric 
metric_obj = function(bi_obj){
S1= list()
for (i in 1:bi_obj@Number) { 
  tmp = list()
  tmp[["r"]]=which(bi_obj@GenesMembership[,i]==TRUE)
  tmp[["c"]]=which(bi_obj@ColumnMembership[i,]==TRUE)
  S1[[i]]= tmp  
}
return(S1)
}

########## Sensitivity and Specificity #################

sensep = function(bitrue,bisslb){
  # gene-wise
  p = dim(bi_true@GenesMembership)[1]
  n = dim(bi_true@ColumnMembership)[1]
  L_S = dim(bi_true@GenesMembership)[2]
  L_E = dim(bi_sslb@GenesMembership)[2]
  search_col = seq(1,L_E,1)
  tp = c() 
  fp = c()
  if(L_S<=L_E){
    for (l_t in 1:L_S) {
      match_factor = bi_true@GenesMembership[,l_t] 
      tf = c() 
      for (l_e in search_col) {
        tf = c(tf,sum(match_factor*bi_sslb@GenesMembership[,l_e]))
      }
      
      if(length(which(tf==max(tf))==1)){
        tp = c(tp, max(tf))
        remove_col = which(tf==max(tf))
        fp = c(fp,sum(bi_sslb@GenesMembership[,search_col[remove_col]])-max(tf))
        search_col = search_col[-remove_col]
      }else{
        print("Problem")
        break
      }
    }
    
    if(length(search_col)>0){
      fp = c(fp,sum(bi_sslb@GenesMembership[,search_col]))   
    }
    
    sen=sum(tp)/sum(bi_true@GenesMembership)
    sep=(p*L_S-sum(bi_true@GenesMembership)-sum(fp))/(p*L_S-sum(bi_true@GenesMembership))
    return(c(sen,sep))
  }else{
    return(c(0,0))
  }
}



