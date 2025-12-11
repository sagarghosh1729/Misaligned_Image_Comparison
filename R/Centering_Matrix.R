#--------------- 1. Centering matrix transformation --------------#

Centering_matrix=function(k){
  C=matrix(0,nrow=k,ncol = k)
  for(i in 1:k){
    for(j in 1:k){
      if(i==j){
        C[i,j]=(k-1)/k
      }
      else{
        C[i,j]=(-1)/k
      }
    }
  }
  return(C)
}
