rm(list=ls())
#install.packages("gmp")
library(gmp)

get_expectations=function(max_coverage,nd){
  coverage=seq(0,max_coverage,by=1)
  nc=2*nd
  mat=matrix(rep(0,nc*length(coverage)),nrow=nc)
  fac_nc=factorialZ(nc)
  for(co in 1:length(coverage)){
    print(co)
    for(j in 1:nc){
      if(coverage[co]<j){
        mat[j,co]=0
      }else{
        num=fac_nc*Stirling2(coverage[co],j)#sn(coverage[co],j)
        den=as.bigz(factorialZ(nc-j)*pow.bigz(nc,coverage[co]))
        mat[j,co]=as.numeric(num/den)
      }
    }
  }
  
  expectations=rep(0,max_coverage)
  for(col in 1:ncol(mat)){
    expectations[col]=sum(mat[,col]*seq(1,nc,by=1))
  }
  return(expectations)
}

pool_sizes=c(41,41,41,41,41,9,31,31,31,31,31,31,32,34,23,34,31,31)

expectations_mat=matrix(rep(0,1001*length(pool_sizes)),nrow=length(pool_sizes))

for(row in 1:length(pool_sizes)){
  print(row)
  expectations_mat[row,]=get_expectations(1000,pool_sizes[row])
}

write.table(expectations_mat,file="uniquedraws_expectations.txt",col.names = F,row.names = F,quote=F)



