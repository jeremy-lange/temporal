rm(list=ls())
WF=function(N,freq_,s){
  
  num_with_gene=round(N*freq_)
  
  exact=F
  for(i in 1:121){
    
    if(i==1){
      year_75=rbinom(1,size=13,prob=num_with_gene/N)
      if(year_75!=0){
      #if(year_75!=0 && year_75!=1){
        break
      }
    }
    if(i==31){
      year_77=rbinom(1,size=12,prob=num_with_gene/N)
      if(year_77!=1){
      #if(year_77!=1 && year_77!=0 && year_77!=2){
        break
      }
    }
    if(i==46){
      year_78=rbinom(1,size=14,prob=num_with_gene/N)
      if(year_78!=4){
      #if(year_78!=3 && year_78!=4 && year_78!=5){
        break
      }
    }
    if(i==61){
      year_79=rbinom(1,size=23,prob=num_with_gene/N)
      if(year_79!=5){
      #if(year_79!=4 && year_79!=5 && year_79!=6){
        break
      }
    }
    if(i==76){
      year_80=rbinom(1,size=13,prob=num_with_gene/N)
      if(year_80!=5){
      #if(year_80!=4 && year_80!=5 && year_80!=6){
        break
      }
    }
    if(i==120){
      year_83=rbinom(1,size=10,prob=num_with_gene/N)
      if(year_83!=5){
      #if(year_83!=4 && year_83!=5 && year_83!=6){
        break
      }
      exact=T
    }
    
    wA=1+s;wa=1
    p=num_with_gene/N
    # calculate sampling probilities
    prob=wA*p/(wA*p+wa*(1-p))  
    
    # binomial sampling to the next generation
    num_with_gene=rbinom(1,N,prob)
    
  }
  return(exact)
  
}

num_sims=1000000

successes=list()
count=1
for(j in 1:num_sims){
  if(j %% 50000 == 0){
    print(j)
  }
  
  s=runif(1,0,0.2)
  p=runif(1,0,0.2)
  if(WF(19000,s,p)){
    print("found one!")
    successes[[count]]=c(p,s)
    count=count+1
  }
}

mat=matrix(unlist(successes),ncol=2,byrow=T)

write.table(mat,"out.txt",sep=" ",quote=F,row.names=F,col.names = c("initial_freq","selection_strength"))


