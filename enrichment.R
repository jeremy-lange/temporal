rm(list=ls())

run_removal=function(removal_type){
  find_region=function(removal_type,genome_wide,stat,sig_flank,flank_length){
    if(stat == 'snp'){
      col=5
    }else if(stat=='window'){
      col=4
    }else{
      print("provide 'window' or 'snp' for which statistic to use")
    }
    
    if(removal_type=='deterministic'){
      ordered=order(genome_wide[,col])
      find_window=T
      count_=0
      while(find_window){
        p_value=genome_wide[ordered[which(genome_wide[ordered,6]>=1)[1+count_]],col]
        min_=p_value-1/nrow(genome_wide)/2;max_=p_value+1/nrow(genome_wide)/2
        num_in_range=sum(genome_wide[,col]>min_ & genome_wide[,col]<max_)
        exp_=(max_-min_)*nrow(genome_wide)
        if(num_in_range>exp_){
          find_window=F
        }else{
          count_=count_+1
        }
      }
      row=which(as.numeric(genome_wide[,col]) == p_value)
    }else if(removal_type=='random'){
      row=sample(which(genome_wide[,col]<0.01 & genome_wide[,6]>=1),1)
    }else{
      print("provide removal type, random or deterministic")
    }
    
    if(length(row)>1){
      row=row[1]
    }
    genome_wide[row,]
    right_flank=1
    count=0
    while(count<(flank_length-1)){
      if(genome_wide[row,1] != genome_wide[row+right_flank,1] || 
         (row+right_flank)>nrow(genome_wide)){
        if((row+right_flank)>nrow(genome_wide)){
          right_flank=0
        }
        break
      }
      seq1=genome_wide[row+right_flank-1,7]
      seq2=genome_wide[row+right_flank,7]
      if(seq2-seq1!=1){
        #print("found right boundary")
        break
      }
      
      if(as.numeric(genome_wide[row+right_flank,col])>sig_flank){  
        count=count+1
      }else{
        count=0
      }
      right_flank=right_flank+1
    }
    
    left_flank=-1
    count=0
    while(count<(flank_length-1)){
      
      if(genome_wide[row,1] != genome_wide[row+left_flank,1] || 
         (row+left_flank)<1){
        if((row+left_flank)<1){
          left_flank=0
        }
        break
      }
      seq1=genome_wide[row+right_flank-1,7]
      seq2=genome_wide[row+right_flank,7]
      if(seq2-seq1!=1){
        #print("found left boundary")
        break
      }
      
      if(as.numeric(genome_wide[row+left_flank,col])>sig_flank){
        count=count+1
      }else{
        count=0
      }
      left_flank=left_flank-1
    }
    return(c(row,left_flank,right_flank))
  }
  
  get_enrichment=function(lower,upper,p_values){
    enrichments=c()
    for(i in 1:length(lower)){
      enrichments=c(enrichments,sum(p_values>lower[i] & p_values<upper[i])/((upper[i]-lower[i])*length(p_values)))
    }
    return(enrichments)
  }
  
  threepool=rbind(read.table("windowoutput_Chr2L_3pools.txt",header=T),
                  read.table("windowoutput_Chr2R_3pools.txt",header=T),
                  read.table("windowoutput_Chr3L_3pools.txt",header=T),
                  read.table("windowoutput_Chr3R_3pools.txt",header=T),
                  read.table("windowoutput_ChrX_3pools.txt",header=T))
  
  genome_wide=read.table("auto_pvalues.txt")
  
  genome_wide=cbind(genome_wide,seq(1,nrow(genome_wide)))
  genome_wide_keep=genome_wide
  
  flank_length_=10
  sleep=F
  stat='snp'
  #stat='window'
  if(stat == 'snp'){
    col=5
    threepool_col=6
  }else if(stat=='window'){
    col=4
    threepool_col=5
  }
  lower=c(0,0.05,0.1,0.15,0.2)
  upper=c(0.05,0.1,0.15,0.2,1)
  #for visualization
  #lower=c(0,     0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5)
  #upper=c(0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1)
  
  
  midpoints=(lower+upper)/2
  
  e=get_enrichment(lower,upper,genome_wide[genome_wide[,6]>=1,col])
  num_regions=0
  total_removed=0
  
  get_figure=function(save){
    par(mfrow=c(1,1))
    if(save){
      setwd("/Users/jeremylange/Documents/Projects/providence/manuscript/figures")
      png("snp_enrichment.png",
          width = 8*300,        
          height = 5*300,
          res = 300,            
          pointsize = 12)  
    }
    plot(e~log(midpoints),xaxt='n',xlab="Midpoint of P-value Bin",ylab="Enrichment",main="SNP P-value Enrichment",lwd=2)
    axis(1, at=log(midpoints), labels=midpoints)
    if(save){
      dev.off()
    }
  }
  get_figure(FALSE)
  
  while((e[1]>e[2])&(e[1]>1 | e[2]>1)){
    
    ret=find_region(removal_type,genome_wide,stat,0.1,flank_length_)#0.1 is reported in manuscript
    if(ret[3]<flank_length_){
      ret[3]=ret[3]+flank_length_-1
    }
    if(abs(ret[2])<flank_length_){
      ret[2]=ret[2]-flank_length_+1
    }
    
    q=sum(threepool[which(genome_wide[ret[1],1]==threepool[,1] &
                            genome_wide[ret[1],2]==threepool[,2] &
                            genome_wide[ret[1],3]==threepool[,3]),threepool_col]>threepool[,threepool_col])/
      nrow(threepool)
    if(q>0.95){
      num_regions=num_regions+1
      total_removed=total_removed+genome_wide[(ret[1]+ret[3]-flank_length_),3]-genome_wide[(ret[1]+ret[2]+flank_length_),2]+1
    }
    
    genome_wide[((ret[1]+ret[2]+flank_length_):(ret[1]+ret[3]-flank_length_)),]
    genome_wide=genome_wide[-((ret[1]+ret[2]+flank_length_):(ret[1]+ret[3]-flank_length_)),]
    e=get_enrichment(lower,upper,genome_wide[genome_wide[,6]>=1,col])
  }
  return(c(num_regions,total_removed/sum(genome_wide_keep[,3]-genome_wide_keep[,2])))
}


run_removal('deterministic')#20 regions, 4.833%
