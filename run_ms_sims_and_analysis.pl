#!/usr/bin/perl -w

my $reps=10000;
my $pop_size=9500;
my $locus_size;
my $mult=0.44;
my $id=1;
my $recbin=2;
my $min_freq=0.05;
my $site_ss=-1;
my @SiteSSs=();

use strict;
use List::Util 'shuffle';

my @NChoose2 = (0, 0);
my $j = 0;
for (my $i = 1; $i <= 124; $i++){
  $j += $i;
  push @NChoose2, $j;
}

$cmd="tar -xzf recrates.tgz";
system($cmd);
my @recrate_AoA=();
my $rec_file="recrate".$recbin."_lengths.txt";
open I, "<$rec_file" or die "cannot find $rec_file";
while(<I>){
   chomp;
   my @line=split /\s+/, $_;
   push @recrate_AoA, [@line];
}
close I;

my $num_row=@recrate_AoA;
$num_row--;
my @rows=(0 .. $num_row);
@rows=shuffle @rows;
my @shuffled_recrate_AoA=();
for(my $i=1;$i<@rows;$i++){
   push @shuffled_recrate_AoA, [@{$recrate_AoA[$rows[$i]]}];
}
my $recrate_count=0;


#################################################################
system("tar -xzf RI_samplesize_coverage_autosomes.tgz");
system("gunzip RI_samplesize_coverage_autosomes_bin".$recbin.".txt.gz");
my @cov_AoA=();
my $cov_input="RI_samplesize_coverage_autosomes_bin".$recbin.".txt";#  seasonal_effss_and_cov_withinds.txt";
open I, "<$cov_input" or die "cannot find $cov_input\n";
while(<I>){
   chomp;
   my @line=split /\s+/, $_;
   push @cov_AoA, [@line];
}
close I;

$num_row=@cov_AoA;
$num_row--;
@rows=(0 .. $num_row);
my @shuffled_cov_AoA=();
for(my $i=1;$i<@rows;$i++){
   push @shuffled_cov_AoA, [@{$cov_AoA[$rows[$i]]}];
}
my $cov_count=int(rand(@shuffled_cov_AoA));

my $cmd_line;
my $get_sim=0;
my $num_sites;
my @positions;
my @AoA=();
my @freqs=();
my @pass_arr;
our $MaxFST=0.95;
our $Num12Sum;
our $Den12Sum;
our $Num13Sum;
our $Den13Sum;
our $Num23Sum;
our $Den23Sum;
my @SiteFSTAoA=();
my $MaxPBS=0;
my @SitePBSs=();
my @window_stats_full=();
my @chr;
my $sites_analyzed=0;

my $rec_rate;
for(my $rep=0;$rep<$reps;$rep++){

$locus_size=$shuffled_recrate_AoA[$recrate_count][3];
$rec_rate=$shuffled_recrate_AoA[$recrate_count][4];
$recrate_count++;

if($recrate_count>=(@shuffled_recrate_AoA-1)){
   $recrate_count=0;
}

run_model((1,$pop_size,$locus_size,$rec_rate,$mult));

$get_sim=0;
@AoA=();
my $num_pops=3;
my @pop_starts=(14,0,0);
my @pop_stops=(14,0,0);
my $rep_count=0;
my $input="msms_output.txt";
my $num_snps=0;
my $total_samples;
my $PiSum;
my $window_comps;
open I, "$input" or die "cannot find $input !\n";
while(<I>){
   chomp;
   my @line=split /\s+/, $_;

   if(@line==0){
      next;
   }

   if($line[0] eq './ms'){
      $total_samples=$line[1];
      for(my $i=0;$i<@line;$i++){
         if($line[$i] eq '-I'){
            for(my $j=($i+9);$j<($i+15);$j++){
               $pop_stops[0]+=$line[$j];
            }
            $pop_stops[0]--;
            $pop_starts[1]=($pop_stops[0]+1);
            $pop_stops[1]=$pop_starts[1]+$line[$i+15]-1;
            $pop_starts[2]=$pop_stops[1]+1;
            $pop_stops[2]=$pop_starts[2]+$line[$i+16]-1;
         }
      }
   }

   if($line[0] eq '//'){
      $cov_count=int(rand(@shuffled_cov_AoA));
      $PiSum=0;
      $window_comps=0;
      $rep_count++;
      $get_sim=1;
      $cov_count=int(rand(@shuffled_cov_AoA));
      next;
   }
   if($get_sim==1){
      if($line[0] eq 'segsites:'){
         $num_sites=$line[1];
         next;
      }
      if($line[0] eq 'positions:'){
         @positions=@line;
         shift @positions;
         next;
      }
      @chr=split //, $line[0];
      push @AoA, [@chr];
      if(@AoA == $total_samples){

         for(my $col=0;$col<@{$AoA[0]};$col++){

            my @cov=@{$shuffled_cov_AoA[$cov_count]};
            $cov_count++;
            if($cov_count>=(@cov_AoA-1)){
               $cov_count=0;
            }

            $cov[1]=int($cov[21])+int($cov[22])+int($cov[23])+int($cov[24])+int($cov[25])+int($cov[26]);
            $cov[2]=int($cov[27])+int($cov[28])+int($cov[29])+int($cov[30])+int($cov[31])+int($cov[32])+int($cov[33])+int($cov[34])+int($cov[35])+int($cov[36])+int($cov[37])+int($cov[38]);
            my @fall_pool_sizes=(int($cov[21]),int($cov[22]),int($cov[23]),int($cov[24]),int($cov[25]),int($cov[26]));
            my @spring_pool_sizes=(int($cov[27]),int($cov[28]),int($cov[29]),int($cov[30]),int($cov[32]),int($cov[32]),int($cov[33]),int($cov[34]),int($cov[35]),int($cov[36]),int($cov[37]),int($cov[38]));

            my @tmp_pop1=();
            my @pop2=();
            my @pop3=();
            for(my $row=$pop_starts[0];$row<=$pop_stops[0];$row++){
               push @tmp_pop1, $AoA[$row][$col];
            }

            my @tmp1=splice @tmp_pop1, 0, 20;
            my @tmp2=splice @tmp_pop1, 0, 20;
            my @tmp3=splice @tmp_pop1, 0, 22;
            my @tmp4=splice @tmp_pop1, 0, 28;
            my @tmp5=splice @tmp_pop1, 0, 20;
            my @tmp6=splice @tmp_pop1, 0, 18;
            my @pop1=();
            for(my $i=0;$i<$cov[39];$i++){
               push @pop1, $tmp1[$i];
            }
            for(my $i=0;$i<$cov[40];$i++){
               push @pop1, $tmp2[$i];
            }
            for(my $i=0;$i<$cov[41];$i++){
               push @pop1, $tmp3[$i];
            }
            for(my $i=0;$i<$cov[42];$i++){
               push @pop1, $tmp4[$i];
            }
            for(my $i=0;$i<$cov[43];$i++){
               push @pop1, $tmp5[$i];
            }
            for(my $i=0;$i<$cov[44];$i++){
               push @pop1, $tmp6[$i];
            }

            for(my $row=$pop_starts[1];$row<=$pop_stops[1];$row++){
               push @pop2, $AoA[$row][$col];
            }
            for(my $row=$pop_starts[2];$row<=$pop_stops[2];$row++){
               push @pop3, $AoA[$row][$col];
            }

            my @pop2_keep=@pop2;
            my @pop3_keep=@pop3;

            my $summed_cov_fall=$cov[3]+$cov[4]+$cov[5]+$cov[6]+$cov[7]+$cov[8];
            my $summed_ss_fall=$cov[21]+$cov[22]+$cov[23]+$cov[24]+$cov[25]+$cov[26];

            my @fall_pool1=splice @pop2, 0, $fall_pool_sizes[0];
            my @fall_pool2=splice @pop2, 0, $fall_pool_sizes[1];
            my @fall_pool3=splice @pop2, 0, $fall_pool_sizes[2];
            my @fall_pool4=splice @pop2, 0, $fall_pool_sizes[3];
            my @fall_pool5=splice @pop2, 0, $fall_pool_sizes[4];
            my @fall_pool6=splice @pop2, 0, $fall_pool_sizes[5];

            @fall_pool1=read_sample([@fall_pool1],[$cov[3]]);
            @fall_pool2=read_sample([@fall_pool2],[$cov[4]]);
            @fall_pool3=read_sample([@fall_pool3],[$cov[5]]);
            @fall_pool4=read_sample([@fall_pool4],[$cov[6]]);
            @fall_pool5=read_sample([@fall_pool5],[$cov[7]]);
            @fall_pool6=read_sample([@fall_pool6],[$cov[8]]);

            #################################
            my $freq1=get_freq(@fall_pool1);
            my $freq2=get_freq(@fall_pool2);
            my $freq3=get_freq(@fall_pool3);
            my $freq4=get_freq(@fall_pool4);
            my $freq5=get_freq(@fall_pool5);
            my $freq6=get_freq(@fall_pool6);
            $freq1*=($cov[21]/$summed_ss_fall);
            $freq2*=($cov[22]/$summed_ss_fall);
            $freq3*=($cov[23]/$summed_ss_fall);
            $freq4*=($cov[24]/$summed_ss_fall);
            $freq5*=($cov[25]/$summed_ss_fall);
            $freq6*=($cov[26]/$summed_ss_fall);
            my $weighted_fall_freq=$freq1+$freq2+$freq3+$freq4+$freq5+$freq6;
            ###############################################

            my $summed_cov_spring=$cov[9]+$cov[10]+$cov[11]+$cov[12]+$cov[13]+$cov[14]+$cov[15]+$cov[16]+$cov[17]+$cov[18]+$cov[19]+$cov[20];
            my $summed_ss_spring=$cov[27]+$cov[28]+$cov[29]+$cov[30]+$cov[31]+$cov[32]+$cov[33]+$cov[34]+$cov[35]+$cov[36]+$cov[37]+$cov[38];

            my @spring_pool1=splice @pop3, 0, $spring_pool_sizes[0];
            my @spring_pool2=splice @pop3, 0, $spring_pool_sizes[1];
            my @spring_pool3=splice @pop3, 0, $spring_pool_sizes[2];
            my @spring_pool4=splice @pop3, 0, $spring_pool_sizes[3];
            my @spring_pool5=splice @pop3, 0, $spring_pool_sizes[4];
            my @spring_pool6=splice @pop3, 0, $spring_pool_sizes[5];
            my @spring_pool7=splice @pop3, 0, $spring_pool_sizes[6];
            my @spring_pool8=splice @pop3, 0, $spring_pool_sizes[7];
            my @spring_pool9=splice @pop3, 0, $spring_pool_sizes[8];
            my @spring_pool10=splice @pop3, 0, $spring_pool_sizes[9];
            my @spring_pool11=splice @pop3, 0, $spring_pool_sizes[10];
            my @spring_pool12=splice @pop3, 0, $spring_pool_sizes[11];

            @spring_pool1=read_sample([@spring_pool1],[$cov[9]]);
            @spring_pool2=read_sample([@spring_pool2],[$cov[10]]);
            @spring_pool3=read_sample([@spring_pool3],[$cov[11]]);
            @spring_pool4=read_sample([@spring_pool4],[$cov[12]]);
            @spring_pool5=read_sample([@spring_pool5],[$cov[13]]);
            @spring_pool6=read_sample([@spring_pool6],[$cov[14]]);
            @spring_pool7=read_sample([@spring_pool7],[$cov[15]]);
            @spring_pool8=read_sample([@spring_pool8],[$cov[16]]);
            @spring_pool9=read_sample([@spring_pool9],[$cov[17]]);
            @spring_pool10=read_sample([@spring_pool10],[$cov[18]]);
            @spring_pool11=read_sample([@spring_pool11],[$cov[19]]);
            @spring_pool12=read_sample([@spring_pool12],[$cov[20]]);

            $freq1=get_freq(@spring_pool1);
            $freq2=get_freq(@spring_pool2);
            $freq3=get_freq(@spring_pool3);
            $freq4=get_freq(@spring_pool4);
            $freq5=get_freq(@spring_pool5);
            $freq6=get_freq(@spring_pool6);
            my $freq7=get_freq(@spring_pool7);
            my $freq8=get_freq(@spring_pool8);
            my $freq9=get_freq(@spring_pool9);
            my $freq10=get_freq(@spring_pool10);
            my $freq11=get_freq(@spring_pool11);
            my $freq12=get_freq(@spring_pool12);

            $freq1*=($cov[27]/$summed_ss_spring);
            $freq2*=($cov[28]/$summed_ss_spring);
            $freq3*=($cov[29]/$summed_ss_spring);
            $freq4*=($cov[30]/$summed_ss_spring);
            $freq5*=($cov[31]/$summed_ss_spring);
            $freq6*=($cov[32]/$summed_ss_spring);
            $freq7*=($cov[33]/$summed_ss_spring);
            $freq8*=($cov[34]/$summed_ss_spring);
            $freq9*=($cov[35]/$summed_ss_spring);
            $freq10*=($cov[36]/$summed_ss_spring);
            $freq11*=($cov[37]/$summed_ss_spring);
            $freq12*=($cov[38]/$summed_ss_spring);

            my $weighted_spring_freq=$freq1+$freq2+$freq3+$freq4+$freq5+$freq6+$freq7+$freq8+$freq9+$freq10+$freq11+$freq12;

            my $derived1=0;
            for(my $i=0;$i<@pop1;$i++){
               $derived1+=$pop1[$i];
            }
            $freq1=$derived1/@pop1;
            my $summed_ss=$cov[0]+$cov[1]+$cov[2];
            my $freq=$freq1*$cov[0]/$summed_ss+$weighted_fall_freq*$cov[1]/$summed_ss+$weighted_spring_freq*$cov[2]/$summed_ss;
            my $fold=0;
            if($freq>0.5){
               $fold=1;
            }

            if($fold==1){
               $freq1=1-$freq1;
               $weighted_fall_freq=1-$weighted_fall_freq;
               $weighted_spring_freq=1-$weighted_spring_freq;
            }
            $freq2=$weighted_fall_freq;
            $freq3=$weighted_spring_freq;

            next if ((($freq1 + $freq2 + $freq3) / 3) < $min_freq);

            if($cov[0] < 50 || $cov[1] < 100 || $cov[2] < 100){
               next;
            }

            $PiSum+=(@pop1-$derived1)*$derived1;
            $window_comps+=$NChoose2[@pop1];

            $num_snps++;
            my @FSTs=();
            @pass_arr=($freq1,$freq2,$cov[0],$cov[1],1,2);
            my $FST=fst(@pass_arr);
            push @FSTs, $FST;
            @pass_arr=($freq1,$freq3,$cov[0],$cov[2],1,3);
            $FST=fst(@pass_arr);
            push @FSTs, $FST;

            @pass_arr=($freq2,$freq3,$cov[1],$cov[2],2,3);
            $FST=fst(@pass_arr);
            push @FSTs, $FST;

            push @SiteFSTAoA, [ @FSTs ];

            my $SitePBS = ((-1 * log (1 - $FSTs[0])) + (-1 * log (1 - $FSTs[1])) - (-1 * log (1 - $FSTs[2]))) / 2;
            if ($SitePBS > $MaxPBS && $cov[0]>50 && $cov[1] > 50 && $cov[2] > 50){
               $MaxPBS = $SitePBS;
               $site_ss=$cov[0];
            }
            push @SitePBSs, $SitePBS;
         }
         my @window_stats=window_stats();

         my $num=@{$AoA[0]};
         for(my $i=$num;$i<$locus_size;$i++){
            $cov_count++;
            if($cov_count>=(@cov_AoA-1)){
               $cov_count=0;
            }

            $window_comps+=$NChoose2[$shuffled_cov_AoA[$cov_count][0]];
         }

         $PiSum=$PiSum/($window_comps);
         push @window_stats, $MaxPBS;
         push @window_stats, $num_snps;
         push @window_stats, $PiSum;
         push @window_stats, $locus_size;
         push @window_stats, $site_ss;
         push @window_stats_full, [@window_stats];
         $Num12Sum=0;
         $Den12Sum=0;
         $Num13Sum=0;
         $Den13Sum=0;
         $Num23Sum=0;
         $Den23Sum=0;
         $MaxPBS=0;
         $get_sim=0;
         $num_snps=0;
         @AoA=();
      }
   }
}

}

my $output="prov_model_windowstats_bin".$recbin."_ne".$pop_size."_thetamult".$mult."_id".$id.".txt";
open O, ">$output";
for(my $i=0;$i<@window_stats_full;$i++){
   print O "@{$window_stats_full[$i]}\n";
}
close O;

system("rm recrate*");
system("rm RI_seasonal*");
system("rm seedms");
system("rm msms_output.txt");
system("rm RI_samplesize_coverage_autosomes_bin*.txt*");

sub get_freq{
   my @arr=@_;
   my $count=0;
   my $freq=-1;
   if(@arr==0){
      $freq=0;
   }

   if($freq==-1){
      for(my $i=0;$i<@arr;$i++){
         if(!defined $arr[$i]){
print "@arr\n";
print "i is $i\n";
my $num=@arr;
print "num is $num\n";
die "$freq\n";
            }
         $count+=$arr[$i];
      }

      $freq=$count/@arr;
   }
   return($freq);
}

sub run_model{

   my @input=@_;
   my $reps=$input[0];
   my $RI_popsize=$input[1];

   my $mu=5.21e-9;
   my $r=$input[3]*1e-8;
   my $sites=$input[2];
   my $mult=$input[4];
   my $ne=1000443.73;

   my $f;
   if($r==0){
      $f=4*$ne*6.25e-8;
   }else{
      $f=6.25e-8/($r);
   }

   my $lamda=518;

   my $theta=4*$ne*$mu*$sites*$mult;
   my $rho=4*$ne*$r*$sites;
   my $RI_scaled=$RI_popsize/($ne*1);

   #chr2r, with GC
   my $full_cmd="./ms 1312 ".$reps." -t ".$theta." -r ".$rho." ".$sites." -c ".$f." ".$lamda." -I 15 2 2 2 2 2 2 2 20 20 22 28 20 18 428 742 0 -en 0 1 2.329934943 -en 0 2 0.32654905 -en 0 3 0.462832587 -en 0 4 0.128821753 -en 0 5 0.273434259 -en 0 6 2.179843568 -en 0 7 0.450255983 -en 0 8 100 -en 0 9 100 -en 0 10 100 -en 0 11 100 -en 0 12 100 -en 0 13 100 -em 0 1 2 27.28290087 -em 0 2 1 27.28290087 -em 0 1 3 0.010337785 -em 0 3 1 0.010337785 -em 0 1 6 42.73962844 -em 0 6 1 42.73962844 -em 0 1 7 44.89191105 -em 0 7 1 44.89191105 -em 0 2 3 19.09486921 -em 0 3 2 19.09486921 -em 0 2 4 6.302255259 -em 0 4 2 6.302255259 -em 0 2 5 15.04827441 -em 0 5 2 15.04827441 -em 0 3 4 1.707928563 -em 0 4 3 1.707928563 -em 0 3 5 0.011128136 -em 0 5 3 0.011128136 -em 0 4 5 1.097629355 -em 0 5 4 1.097629355 -em 0 6 7 34.38499088 -em 0 7 6 34.38499088 -en 0 15 ".$RI_scaled." -ej 1.74922E-06 14 15 -ej 0.000117948 13 15 -ej 0.000129193 12 15 -ej 0.000132941 11 15 -ej 0.000136689 10 15 -ej 0.000140438 9 15 -ej 0.000147934 8 15 -es 0.000453549 15 0.17 -ej 0.000453549 16 4 -ej 0.000453549 15 3 -ej 0.00888381 6 1 -en 0.00888381 1 0.542234459 -em 0.00888381 1 7 44.89191105 -em 0.00888381 7 1 44.89191105 -em 0.00888381 1 2 45.20362531 -em 0.00888381 2 1 45.20362531 -em 0.00888381 1 3 1.021294138 -em 0.00888381 3 1 1.021294138 -ej 0.034268027 5 3 -en 0.034268027 3 0.255253002 -em 0.034268027 1 3 5.220992483 -em 0.034268027 3 1 5.220992483 -em 0.034268027 3 2 0.166922035 -em 0.034268027 3 2 0.166922035 -ej 0.061497279 4 3 -en 0.061497279 3 0.531134005 -em 0.061497279 2 3 0.17887182 -em 0.061497279 3 2 0.17887182 -ej 0.092526061 3 2 -en 0.092526061 2 0.356825546 -em 0.092526061 2 1 0.008717061 -em 0.092526061 1 2 0.008717061 -ej 0.095780424 2 1 -en 0.095780424 1 0.282788248 -em 0.095780424 1 7 1.986871724 -em 0.095780424 7 1 1.986871724 -ej 0.132010805 7 1 -en 0.132010805 1 1.011163346 -en 0.17881885 1 1 > full_autosome_model_withRI.txt";

   system($full_cmd);
   system("mv full_autosome_model_withRI.txt msms_output.txt");
}

sub sample_inds{
   my @input=@_;
   my @pop=@{$input[0]};
   my $ss=$input[1][0];
   my @rows=(0 .. 127);
   @rows=shuffle @rows;
   my @return_arr=();
   for(my $i=0;$i<$ss;$i++){
      push @return_arr, $pop[$rows[$i]];
   }
   return(@return_arr);
}

sub read_sample{
   my @input=@_;
   my @pop=@{$input[0]};
   my $cov=$input[1][0];
   my @resampled=();
   if(@pop==0){
      return(@resampled);
   }
   for(my $i=0;$i<$cov;$i++){
      my $random_number = int(rand(@pop));
      push @resampled, $pop[$random_number];
   }

   return(@resampled);
}

sub get_minor_freq{
   my @input=@_;
   my @pop1=@{$input[0]};
   my @pop2=@{$input[1]};
   my @pop3=@{$input[2]};
   my $derived1=0;
   my $derived2=0;
   my $derived3=0;
   my $fold=0;
   for(my $i=0;$i<@pop1;$i++){
      $derived1+=$pop1[$i];
   }
   for(my $i=0;$i<@pop2;$i++){
      $derived2+=$pop2[$i];
   }
   for(my $i=0;$i<@pop3;$i++){
      $derived3+=$pop3[$i];
   }
   my $ss=@pop1+@pop2+@pop3;
   my $freq=($derived1+$derived2+$derived3)/$ss;
   if($freq>0.5){
      $fold=1;
      $derived1=@pop1-$derived1;
      $derived2=@pop2-$derived2;
      $derived3=@pop3-$derived3;
   }
   return(($derived1,$derived2,$derived3,$fold));
}

sub get_coverage{
   my @AoA=@_;
   my $num_row=@AoA;
   $num_row--;
   my @rows=(0 .. $num_row);
   @rows=shuffle @rows;
   my @return_arr=();
   for(my $i=1;$i<@{$AoA[0]};$i++){
      push @return_arr, $AoA[$rows[0]][$i];
   }
   return(@return_arr);
}

sub fst{
   my @passed_arr=@_;
   my $MinorFreq1=$passed_arr[0];
   my $MinorFreq2=$passed_arr[1];
   my $SampleSize1=$passed_arr[2];
   my $SampleSize2=$passed_arr[3];
   my $pop1=$passed_arr[4];
   my $pop2=$passed_arr[5];

   #Reynolds calculations
   my $FST;
   my $denom = $MinorFreq1 + $MinorFreq2 - 2 * $MinorFreq1 * $MinorFreq2;
   my $SharedNum = $SampleSize1 * (2 * $MinorFreq1 - 2 * $MinorFreq1 ** 2) + $SampleSize2 * (2 * $MinorFreq2 - 2 * $MinorFreq2 ** 2);
   my $NumA = ($MinorFreq1 - $MinorFreq2) ** 2;
   my $FracNum = ($SampleSize1 + $SampleSize2) * $SharedNum;
   my $FracDen = 4 * $SampleSize1 * $SampleSize2 * ($SampleSize1 + $SampleSize2 - 1);
   my $frac = $FracNum / $FracDen;
   my $WholeNum = $NumA - $frac;
   my $DenFracNum = (4 * $SampleSize1 * $SampleSize2 - $SampleSize1 - $SampleSize2) * $SharedNum;
   my $DenFrac = $DenFracNum / $FracDen;
   my $WholeDen = $NumA + $DenFrac;
   if ($WholeDen != 0){
      $FST = $WholeNum / $WholeDen;
   }
   else{
      $FST = 0;
   }
   if ($FST == 1){
      $FST = $MaxFST;
   }

   if($pop1==1 && $pop2==2){
      $Num12Sum += $WholeNum;
      $Den12Sum += $WholeDen;
   }elsif($pop1==1 && $pop2==3){
      $Num13Sum += $WholeNum;
      $Den13Sum += $WholeDen;
   }elsif($pop1==2 && $pop2==3){
      $Num23Sum += $WholeNum;
      $Den23Sum += $WholeDen;
   }else{
      die "weird pops!\n";
   }

   return($FST);
}

sub window_stats{
  my @passed_arr=@_;

  #calculate window stats
  my @pis;

   my @FSTs = ();
   my $FST = $Num12Sum / $Den12Sum;
   if ($FST >= 1){
      $FST = $MaxFST;
   }
   push @FSTs, $FST;
      $FST = $Num13Sum / $Den13Sum;
   if ($FST >= 1){
      $FST = $MaxFST;
   }
   push @FSTs, $FST;
   $FST = $Num23Sum / $Den23Sum;
   if ($FST >= 1){
      $FST = $MaxFST;
   }
   push @FSTs, $FST;
   my $WinPBS = ((-1 * log (1 - $FSTs[0])) + (-1 * log (1 - $FSTs[1])) - (-1 * log (1 - $FSTs[2]))) / 2;
   push @FSTs, $WinPBS;
   return(@FSTs);
}
