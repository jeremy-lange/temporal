rm(list=ls())

run_sim=function(df,inv_pos,pool_sizes){
  num_pools=(ncol(df)-1)/2
  num_snps=nrow(df)
  sim_df=matrix(rep(0,(num_snps*num_pools*2)),nrow=num_snps)
  for(row in 1:nrow(df)){
    for(pool in 1:num_pools){
      coverage=df[row,pool*2]+df[row,(pool*2+1)]
      if(coverage<10){
        coverage=10
      }
      inv_reads=rbinom(1,coverage,inv_pos[pool]/pool_sizes[pool])
      noninv_reads=coverage-inv_reads
      sim_df[row,(pool*2)-1]=inv_reads
      sim_df[row,(pool*2)]=noninv_reads
    }
  }
  sim_df=cbind(rep(0,num_snps),sim_df)
  return(sim_df)
}
get_inv_freq=function(df){
  num_pools=(ncol(df)-1)/2
  mean_freqs=rep(0,num_pools)
  for(i in seq(2,ncol(df)-1,by=2)){
    eligible=(df[,i]+df[,i+1])>=10
    mean_freqs[i/2] =mean(df[eligible,i]/(df[eligible,i]+df[eligible,i+1]))
  }
  if(num_pools==6){
    pool_sizes=c(41,41,41,41,41,9) 
  }else if(num_pools==12){
    pool_sizes=c(31,31,31,31,31,31,32,34,23,34,31,31)
  }else{
    pool_sizes=1
  }
  
  inv_freq=0
  for(i in 1:length(mean_freqs)){
    inv_freq=inv_freq+mean_freqs[i]*pool_sizes[i]/(sum(pool_sizes))
  }
  return(c(inv_freq,num_pools))
}

inversions=c("(2L)t",	"(2R)Ns",	"(3L)P",	"(3R)C",	"(3R)K",	"(3R)Mo",	"(3R)P")

fall_freq=c(0.1563, 0.06823, 0.005944, 0.01972, 0.01778,	0.05493,	0.03951)
spring_freq=c(0.1132,	0.02930,	0.001913,	0.02323,	0.002553,	0.08902,	0.08833)
diff=abs(fall_freq-spring_freq)

setwd("/path/to/*pool_freqs.txt")

master_func=function(inv,season){
  df_fall=read.table(paste("In",inv,"_fall_poolfreqs.txt",sep=""))
  df_spring=read.table(paste("In",inv,"_spring_poolfreqs.txt",sep=""))
  
  if(season=="fall"){
    vec=get_inv_freq(df_fall)
  }else{
    vec=get_inv_freq(df_spring)
  }
  
  vec=get_inv_freq(df_fall)
  p1=vec[1]
  vec=get_inv_freq(df_spring)
  p2=vec[1] 
  p=median(c(p1,p2))
  
  num_pools_fall=get_inv_freq(df_fall)[2]
  num_pools_spring=get_inv_freq(df_spring)[2]
  
  num_snps=nrow(df_fall)
  pool_sizes_fall=c(41,41,41,41,41,9) 
  pool_sizes_spring=c(31,31,31,31,31,31,32,34,23,34,31,31)
  pool_sizes_fall=2*pool_sizes_fall
  pool_sizes_spring=2*pool_sizes_spring
  inv_pos_fall=rbinom(num_pools_fall,pool_sizes_fall,p)
  inv_pos_spring=rbinom(num_pools_spring,pool_sizes_spring,p)
  sim_df_fall=run_sim(df_fall,inv_pos_fall,pool_sizes_fall)
  sim_df_spring=run_sim(df_spring,inv_pos_spring,pool_sizes_spring)
  
  localities_fall=data.frame(cbind(sim_df_fall[,2]/(sim_df_fall[,2]+sim_df_fall[,3]),
                                   sim_df_fall[,4]/(sim_df_fall[,4]+sim_df_fall[,5]),
                                   sim_df_fall[,6]/(sim_df_fall[,6]+sim_df_fall[,7]),
                                   sim_df_fall[,8]/(sim_df_fall[,8]+sim_df_fall[,9]),
                                   sim_df_fall[,10]/(sim_df_fall[,10]+sim_df_fall[,11]),
                                   sim_df_fall[,12]/(sim_df_fall[,12]+sim_df_fall[,13])
  ))
  
  localities_spring=data.frame(cbind((sim_df_spring[,2]+sim_df_spring[,4])/(sim_df_spring[,2]+sim_df_spring[,3]+sim_df_spring[,4]+sim_df_spring[,5]),
                                     (sim_df_spring[,6]+sim_df_spring[,8])/(sim_df_spring[,6]+sim_df_spring[,7]+sim_df_spring[,8]+sim_df_spring[,9]),
                                     (sim_df_spring[,10]+sim_df_spring[,12])/(sim_df_spring[,10]+sim_df_spring[,11]+sim_df_spring[,12]+sim_df_spring[,13]),
                                     (sim_df_spring[,14]+sim_df_spring[,16]+sim_df_spring[,18]+sim_df_spring[,20])/(sim_df_spring[,14]+sim_df_spring[,15]+sim_df_spring[,16]+sim_df_spring[,17]+sim_df_spring[,18]+sim_df_spring[,19]+sim_df_spring[,20]+sim_df_spring[,21]),
                                     (sim_df_spring[,22]+sim_df_spring[,24])/(sim_df_spring[,22]+sim_df_spring[,23]+sim_df_spring[,24]+sim_df_spring[,25])
  ))
  
  fall_means=c(mean(localities_fall[,1],na.rm=T),mean(localities_fall[,2],na.rm=T),mean(localities_fall[,3],na.rm=T),mean(localities_fall[,4],na.rm=T),mean(localities_fall[,5],na.rm=T),mean(localities_fall[,6],na.rm=T))
  spring_means=c(mean(localities_spring[,1],na.rm=T),mean(localities_spring[,2],na.rm=T),mean(localities_spring[,3],na.rm=T),mean(localities_spring[,4],na.rm=T),mean(localities_spring[,5],na.rm=T))
  diff_fall=max(fall_means)-min(fall_means)
  diff_spring=max(spring_means)-min(spring_means)
  
  sim_fall_invfreq=get_inv_freq(sim_df_fall)
  sim_spring_invfreq=get_inv_freq(sim_df_spring)
  return(c(sim_fall_invfreq[1],sim_spring_invfreq[1],diff_fall,diff_spring))
}

num_reps=10000
sim_diff=rep(0,num_reps)
temporal_diff=rep(0,num_reps)
inseason_differences_fall=rep(0,num_reps)
inseason_differences_spring=rep(0,num_reps)
inv="3RPayne"
season="fall"

if(inv=="2Lt" && season=="fall"){
  max_diff_fall=0.0839289
  max_diff_spring=0.10407192
  inv_index=1
}else if(inv=="2RNs" && season=="fall"){
  max_diff_fall=0.04460213
  max_diff_spring=0.104793358
  inv_index=2
}else if(inv=="3LP" && season=="fall"){
  max_diff_fall=0.015409009
  max_diff_spring=0.003498559
  inv_index=3
}else if(inv=="3RC" && season=="fall"){
  max_diff_fall=0.030966607
  max_diff_spring=0.067740783
  inv_index=4
}else if(inv=="3RK" && season=="fall"){
  max_diff_fall=0.030584075
  max_diff_spring=0.015277778
  inv_index=5
}else if(inv=="3RMo" && season=="fall"){
  max_diff_fall=0.05091697
  max_diff_spring=0.09041257
  inv_index=6
}else if(inv=="3RPayne" && season=="fall"){
  max_diff_fall=0.05391314
  max_diff_spring=0.181480789
  inv_index=7
}

for(i in 1:num_reps){
  if(i %% 1000 == 0){
    print(i)
  }
  ret_vec=master_func(inv,season)
  
  ##################################################
  sim_diff[i]=ret_vec[1]-ret_vec[2]
  inseason_differences_fall[i]=ret_vec[3]
  inseason_differences_spring[i]=ret_vec[4]
  ##################################################
  
}

inseason_pvalue_fall=sum(inseason_differences_fall>max_diff_fall)/num_reps
inseason_pvalue_spring=sum(inseason_differences_spring>max_diff_spring)/num_reps
print(paste("fall inseason p value ",inseason_pvalue_fall))
print(paste("spring inseason p value ",inseason_pvalue_spring))

if(season=="fall"){
  pvalue_fall=sum(abs(sim_diff)>diff[inv_index])/num_reps
  print(paste("seasonal p value ",pvalue_fall))
}else{
  pvalue_spring=sum(abs(sim_diff)>diff[inv_index])/num_reps
  print(paste("seasonal p value ",pvalue_spring))
}
