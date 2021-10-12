rm(list=ls())

ne=9500
num_targets=5000

sim_old=function(p,sizes,ne){#,ne){
  mult=1
  
  year1=rbinom(1,sizes[1],p)
  
  for(i in seq(1,30)){
    p=rbinom(1,2*mult*ne,p)/(2*mult*ne)
  }
  year2=rbinom(1,sizes[2],p)
  for(i in seq(1,15)){
    p=rbinom(1,2*mult*ne,p)/(2*mult*ne)
  }
  year3=rbinom(1,sizes[3],p)
  for(i in seq(1,15)){
    p=rbinom(1,2*mult*ne,p)/(2*mult*ne)
  }
  year4=rbinom(1,sizes[4],p)
  for(i in seq(1,15)){
    p=rbinom(1,2*mult*ne,p)/(2*mult*ne)
  }
  year5=rbinom(1,sizes[5],p)
  for(i in seq(1,45)){
    p=rbinom(1,2*mult*ne,p)/(2*mult*ne)
  }
  year6=rbinom(1,sizes[6],p)
  return(c(year1,year2,year3,year4,year5,year6,p))
}
run_sim=function(ne,fall_medians,spring_medians,fall_medians_mat,spring_medians_mat,ss){
  
  mult=1#for autosomes
  
  p=runif(1)
  
  ######################################################
  sizes=c(2,2,2,2,2,2)
  ######################################################
  years=sim_old(p,sizes,ne)
  p=years[7]
  p_keep=p
  years=years[1:6]
  
  for(i in seq(1,31*15)){
    p=rbinom(1,2*mult*ne,p)/(2*mult*ne)
  }
  
  fall=rbinom(c(1,1,1,1,1,1),2*c(41,41,41,41,41,9),p)
  
  rand_cov=sample(nrow(fall_medians_mat))[1]
  resample_fall=rbinom(c(1,1,1,1,1,1),as.numeric(fall_medians_mat[rand_cov,]),fall/(2*c(41,41,41,41,41,9)))
  for(i in seq(1,7)){
    p=rbinom(1,2*mult*ne,p)/(2*mult*ne)
  }
  
  spring=rbinom(c(1,1,1,1,1,1,1,1,1,1,1,1),2*c(31,31,31,31,31,31,31,34,23,34,31,31),p)
  
  resample_spring=rbinom(c(1,1,1,1,1,1,1,1,1,1,1,1),as.numeric(spring_medians_mat[rand_cov,]),spring/(2*c(31,31,31,31,31,31,31,34,23,34,31,31)))
  
  f_old=rbinom(1,ss,p_keep)/ss
  
  f_fall=sum(
    (resample_fall/fall_medians_mat[rand_cov,])*
      (2*c(41,41,41,41,41,9)/
         sum(2*c(41,41,41,41,41,9))
      ))
  
  f_spring=sum(
    (resample_spring/spring_medians_mat[rand_cov,])*
      (2*c(31,31,31,31,31,31,31,34,23,34,31,31)/
         sum(2*c(31,31,31,31,31,31,31,34,23,34,31,31))
      ))
  
  
  f_new=sum((resample_fall/fall_medians_mat[rand_cov,])*(2*c(41,41,41,41,41,9)),(resample_spring/spring_medians_mat[rand_cov,])*(2*c(31,31,31,31,31,31,31,34,23,34,31,31)))/(sum(2*c(41,41,41,41,41,9),2*c(31,31,31,31,31,31,31,34,23,34,31,31)))
  
  f=sum(c(years,
          (resample_fall/as.numeric(fall_medians_mat[rand_cov,]))*(2*c(41,41,41,41,41,9)),
          (resample_spring/as.numeric(spring_medians_mat[rand_cov,]))*(2*c(31,31,31,31,31,31,31,34,23,34,31,31))))/
    (sum(c(sizes,2*c(41,41,41,41,41,9),2*c(31,31,31,31,31,31,31,34,23,34,31,31))))
  
  if(f>0.5){
    fold=T
  }else{
    fold=F
  }
  
  if(fold){
    f_old=1-f_old
    f_fall=1-f_fall
    f_spring=1-f_spring
  }
  
  
  
  return(c(f_old,f_fall,f_spring))
}

#setwd("/Users/jeremylange/Documents/Projects/providence/manuscript/scripts/ne_analysis/")
setwd("/path/to/seasonal_median txt files")
fall_medians_mat=read.table("fall_medians_mat.txt")
spring_medians_mat=read.table("spring_medians_mat.txt")
fall_medians=c(106,  33,  38,  49,  40, 102)
spring_medians=c(20, 15, 17, 16, 21, 14, 14, 16, 15, 17, 17, 17)

#setwd("../PBS/")
setwd("/path/to/SNP PBS files")
chr=c("2L","2R","3L","3R")
emp_diffs=c()
emp_ss=c()
for(c in 1:length(chr)){
  print(c)
  
  file=paste("SNPoutput_Chr",chr[c],".txt",sep="")
  
  if(c==1){
    df=readLines(file)
    df=df[-1]
    df=matrix(unlist(strsplit(df,split="\t")),nrow=length(df),byrow=T)
  }else{
    tmp=readLines(file)
    tmp=tmp[-1]
    tmp=matrix(unlist(strsplit(tmp,split="\t")),nrow=length(tmp),byrow=T)
    df=rbind(df,tmp)
  }
}
#looking at high sample size SNP at >5% frequency
df=df[as.numeric(df[,8])>60 & as.numeric(df[,26])>150 & as.numeric(df[,27])>150 & as.numeric(df[,4])>=0.05,]

emp_diffs=((as.numeric(df[,5])+as.numeric(df[,6]))/2-as.numeric(df[,3]))
print(paste("empirical standard dev is ",sd(emp_diffs),sep=""))#0.07834901
emp_ss=as.numeric(df[,13])+as.numeric(df[,14])+as.numeric(df[,15])+as.numeric(df[,16])

minimum=seq(0.05,0.97,by=0.02);maximum=seq(0.07,1,by=0.02)
sfs=rep(0,length(minimum))
for(i in 1:length(sfs)){
  sfs[i]=sum(as.numeric(df[,4])>=minimum[i] & as.numeric(df[,4])<maximum[i])/nrow(df)
}

targets=round(sfs*num_targets)
output_list=list()
num_snps=1
while(sum(targets)>100){
  fs=run_sim(ne,fall_medians,spring_medians,fall_medians_mat,spring_medians_mat,sample(emp_ss,1))
  if(fs[1]<0.05 | fs[2]<0.05 | fs[3]<0.05 ){
    next
  }
  freq_class=max(which(fs[1]>=minimum))
  if(length(freq_class)!=1){
    break
  }
  if(targets[freq_class]>0){
    output_list[[num_snps]]=fs
    num_snps=num_snps+1
    targets[freq_class]=targets[freq_class]-1
  }
  if(num_snps %% 1000 == 0){
    print(num_snps)
    print(targets)
  }
}

output_mat=matrix(unlist(output_list),ncol=3,byrow=T)
hist(output_mat[,1])
(output_mat[,2]+output_mat[,3])/2-output_mat[,1]
hist((output_mat[,2]+output_mat[,3])/2-output_mat[,1])
print(sd((output_mat[,2]+output_mat[,3])/2-output_mat[,1]))