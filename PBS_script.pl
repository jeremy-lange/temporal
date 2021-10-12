use strict;
no strict 'refs';

my $chr=$ARGV[0];

our @expected_unique_arr=get_expectations();
our @mean_coverage_fall=get_coverage("fall_meancoverage.txt");
our @mean_coverage_spring=get_coverage("spring_meancoverage.txt");

my %chr_hash = ("2L"  => 0,"2R"  => 1,"3L"  => 2,"3R"  => 3,"X"  => 4);

my $old_dir="/path/to/fas1k/";
my $het_regions="/path/to/het_regions/";
my $new_dir1="/path/to/sync1/";
my $new_dir2="/path/to/sync2/";
my $sync1="Chr".$chr.".sync";
my $sync2="Chr".$chr.".sync";
my $cwd=`pwd`;
chomp $cwd;

my @NChoose2 = (0, 0);
my $j = 0;
for (my $i = 1; $i <= 124; $i++){
  $j += $i;
  push @NChoose2, $j;
}

my $inv_freq1;
my $noninv_freq1;
my @inv_freq2;
my @noninv_freq2;
our @inv_lines=();
if($chr eq "2L"){
   @inv_freq2=(0.12998815,1-0.12998815);
   @inv_lines=([("RI75-1","RI77-51","RI77-55","RI78-13","RI79-13","RI79-16","RI79-19","RI83-1","RI83-3","RI83-7")],[("RI75-2","RI75-4","RI75-5","RI75-6","RI75-7","RI75-10","RI75-11","RI75-12","RI75-13","RI77-50","RI77-52","RI77-53","RI77-54","RI77-56","RI77-58","RI77-60","RI77-61","RI78-1A","RI78-1B","RI78-5","RI78-6","RI78-7","RI78-8","RI78-9","RI78-11","RI78-12","RI78-14","RI79-1","RI79-3","RI79-7","RI79-9","RI79-10","RI79-11","RI79-12","RI79-14","RI79-15","RI79-17","RI80-2","RI80-4","RI80-5","RI80-6","RI80-7","RI80-8","RI80-9","RI80-11","RI80-13","RI80-15","RI83-2","RI83-4","RI83-5","RI83-6","RI83-8","RI83-9")]);
}elsif($chr eq "2R"){
   @inv_freq2=(0.048706065,1-0.048706065);
   @inv_lines=([("RI75-2","RI75-10","RI75-13","RI77-51","RI77-60","RI78-1A","RI78-12","RI79-3","RI79-9","RI79-17","RI79-19","RI80-4","RI80-7","RI79-12")],[("RI75-1","RI75-4","RI75-5","RI75-6","RI75-7","RI75-11","RI75-12","RI77-50","RI77-52","RI77-53","RI77-54","RI77-55","RI77-56","RI77-58","RI77-61","RI78-5","RI78-6","RI78-7","RI78-8","RI78-9","RI78-11","RI78-13","RI78-14","RI79-1","RI79-7","RI79-10","RI79-11","RI79-13","RI79-14","RI79-15","RI79-16","RI80-2","RI80-5","RI80-6","RI80-8","RI80-9","RI80-11","RI80-13","RI80-15","RI83-1","RI83-2","RI83-3","RI83-4","RI83-5","RI83-6","RI83-7","RI83-8","RI83-9","RI78-1B")]);
}elsif($chr eq "3L"){
   @inv_freq2=(0.003997937,1-0.003997937);
   @inv_lines=([("RI77-50","RI79-13")],[("RI75-1","RI75-2","RI75-4","RI75-5","RI75-6","RI75-7","RI75-10","RI75-11","RI75-12","RI75-13","RI77-51","RI77-52","RI77-53","RI77-54","RI77-55","RI77-56","RI77-58","RI77-60","RI77-61","RI78-1A","RI78-5","RI78-6","RI78-7","RI78-8","RI78-9","RI78-11","RI78-12","RI78-13","RI78-14","RI79-1","RI79-3","RI79-7","RI79-9","RI79-10","RI79-11","RI79-12","RI79-14","RI79-15","RI79-16","RI79-17","RI79-19","RI79-20","RI80-2","RI80-4","RI80-5","RI80-6","RI80-7","RI80-8","RI80-9","RI80-11","RI80-13","RI80-15","RI83-1","RI83-2","RI83-3","RI83-4","RI83-5","RI83-6","RI83-7","RI83-8","RI83-9","RI78-1B")]);
}elsif($chr eq "3R"){
   #3RC, 3RP, 3RMo, 3RK, non_inv
   @inv_freq2=(0.021113585,0.063356715,0.071327415,0.010795037,(1-(0.021113585+0.010795037+0.071327415+0.063356715)));
   @inv_lines=([("RI75-6","RI77-52","RI77-61","RI78-1A","RI80-4","RI80-8","RI80-11","RI80-13")],[("RI75-4","RI75-5","RI77-50","RI77-55","RI78-1B","RI78-9","RI79-9","RI79-13","RI79-16","RI79-17","RI79-19","RI80-7","RI83-1","RI83-2","RI83-4","RI83-6","RI83-7")],[("RI75-2","RI77-53","RI77-58","RI78-8","RI78-14","RI79-20","RI80-2","RI80-6","RI80-9")],[("RI77-56","RI79-14","RI79-15")],[("RI77-54", "RI80-5","RI75-1","RI75-7","RI75-10","RI75-11","RI75-12","RI75-13","RI77-51","RI77-60","RI78-5","RI78-6","RI78-7","RI78-11","RI78-12","RI78-13","RI79-1","RI79-3","RI79-7","RI79-10","RI79-11","RI79-12","RI80-15","RI83-8","RI83-9")]);
}elsif($chr="X"){
   @inv_freq2=(1);
   @inv_lines=([("RI79-18","RI79-20","RI75-1","RI77-51","RI77-55","RI78-13","RI79-13","RI79-16","RI79-19","RI83-1","RI83-3","RI83-7","RI75-2","RI75-4","RI75-5","RI75-6","RI75-7","RI75-10","RI75-11","RI75-12","RI75-13","RI77-50","RI77-52","RI77-53","RI77-54","RI77-56","RI77-58","RI77-60","RI77-61","RI78-1A","RI78-1B","RI78-5","RI78-6","RI78-7","RI78-8","RI78-9","RI78-11","RI78-12","RI78-14","RI79-1","RI79-3","RI79-7","RI79-9","RI79-10","RI79-11","RI79-12","RI79-14","RI79-15","RI79-17","RI80-2","RI80-4","RI80-5","RI80-6","RI80-7","RI80-8","RI80-9","RI80-11","RI80-13","RI80-15","RI83-2","RI83-4","RI83-5","RI83-6","RI83-8","RI83-9")]);
}else{
   die "need a value chromosome input, $chr is incorrect!\n";
}
our $num_inversions=@inv_lines;

my $end_win;
if($chr eq "2L"){
   $end_win=23011543;
}elsif($chr eq "2R"){
   $end_win=21146707;
}elsif($chr eq "3L"){
   $end_win=24543556;
}elsif($chr eq "3R"){
   $end_win=27905052;
}elsif($chr eq "X"){
   $end_win=22422826;
}
my $OutputFile="windowoutput_Chr".$chr.".txt";
my $SNPFile="SNPoutput_Chr".$chr.".txt";
our $SampleMin=10;
my $FreqThresh = 0.05;
our $MaxFST = 0.95;
my @pops=("old","new1","new2");

chdir($old_dir);
my @annotation=glob("Chr".$chr."*annotations*");
my @genomes=glob("*Chr".$chr."*.fas1k");

my @InHandles=();
my @supp=('c','d','e','f');
for(my $c=0;$c<@supp;$c++){
  for(my $d=0;$d<@supp;$d++){
    for(my $e=0;$e<@supp;$e++){
      for(my $f=0;$f<@supp;$f++){
	push @InHandles, $supp[$c].$supp[$d].$supp[$e].$supp[$f];
      }
    }
  }
}

open $InHandles[0], "<$annotation[0]" or die "cannot find annotation file!\n";
my $file;
our @lines=();
for (my $f = 0; $f < @genomes; $f++){
   $file = $genomes[$f];
   my @arr=split /r/, $file;
   push @lines, $arr[0];
   open $InHandles[$f+1], "<$file" or print "Can't open $file\n";
}

chdir($new_dir1);
open NA, "<$sync1" or die "cannot find $sync1\n";

chdir($new_dir2);
open NB, "<$sync2" or die "cannot find $sync2\n";

chdir($het_regions);
my @het_region_files=glob("*Chr".$chr."*_hetRegions.txt");
our @het_regions=();
for(my $i=0;$i<@het_region_files;$i++){
  my @tmp=();
  open I, "<$het_region_files[$i]" or die "cannot fine $het_region_files[$i]\n";
  while(<I>){
    chomp;
    my @line=split /\s+/, $_;
    push @tmp, [@line];
  }
  close I;
  push @het_regions, [@tmp];
}
die if @het_regions != @genomes;

chdir($old_dir);
my $file1=$InHandles[0];
our $kb=-1;
our $pos;
our $SitesAnalyzed = 0;
our $Dxy12Sum = 0;
our $Dxy13Sum = 0;
our $Dxy23Sum = 0;
our $Pi1Sum = 0;
our $window_comps=0;
our $Pi2Sum = 0;
our $Pi3Sum = 0;
our $Num12Sum=0;
our $Den12Sum=0;
our $Num13Sum=0;
our $Den13Sum=0;
our $Num23Sum=0;
our $Den23Sum=0;
my @pass_arr;
my @keep;
our @positions=();
our @SiteWins=();
our @SiteFSTAoA=();
our @NucCountAoA=();
my $SitePBS=();
our $MaxPBS = 0;
our $max_snp_location;
our @max_snp_locations=();
our $window=0;
my @SitePBSs=();
our @WinSiteCounts = ();
our @PiAoA=();
our @DxyAoA=();
our @MaxPBSs=();
our @WinPBSs=();
our @WinZs=();
our @FSTAoA=();
our @ss1=();
our @ss2=();
our @ss3=();

our @unweighted_minors=();
our @weighted_minors=();
our @unweighted_ss=();
our @weighted_ss=();
our @minors2=();
our @minors3=();

my @WinStarts = ();
my @WinStops = ();

my $WindowFile="/path/to/window_file/admix_windows_200_Chr".$chr.".txt";
open W, "<$WindowFile" or die "Can't open window file $WindowFile\n";
while (<W>){
   chomp;
   last if m/^$/;
   my @line = split;
   push @WinStarts, $line[0];
   push @WinStops, $line[1];
}
close W;

while(<$file1>){
   chdir($old_dir);
   my @InputAoA=();
   chomp;
   $kb++;
if($kb>100){
last;
}
   last if m/^$/;
   my @line = split(//, $_);
   push @InputAoA, [ @line ];
   for (my $f = 0; $f < @genomes; $f++){
      $file = $InHandles[$f+1];
      $_ = (<$file>);
      chomp;
      @line = split(//, $_);
      push @InputAoA, [ @line ];
   }
   #get new samples
   my @new1_AoA=get_initial();
   chdir($new_dir1);
   $file='NA';
   $_=(<$file>);
   chomp;
   @line=split(/\s+/, $_);
   $pos=$line[1];
   if($pos-$kb*1000 > 1000){
      #dont need to do anything
   }else{
      @{$new1_AoA[($pos % 1000)-1]}=@line;
      while($pos<($kb*1000+1000)){
         $_=(<$file>);
         chomp;
         @line=split(/\s+/, $_);
         $pos=$line[1];
         last if $pos > ($kb*1000+1000);
         my $index=($pos % 1000)-1;
         @{$new1_AoA[$index]}=@line;
      }
   }

   my @new2_AoA=get_initial();
   chdir($new_dir2);
   $file='NB';
   $_=(<$file>);
   chomp;
   @line=split(/\s+/, $_);
   $pos=$line[1];
   if($pos-$kb*1000 > 1000){
      #dont need to do anything
   }else{
      @{$new2_AoA[($pos % 1000)-1]}=@line;
      while($pos<($kb*1000+1000)){
         $_=(<$file>);
         chomp;
         @line=split(/\s+/, $_);
         $pos=$line[1];
         last if $pos > ($kb*1000+1000);
         my $index=($pos % 1000)-1;
         @{$new2_AoA[$index]}=@line;
      }
   }

   #Check whether we've got our full data - if so, transfer the other part of this kb to NextAoA
   if ( ($WinStops[$window]) <= (($kb*1000)+999) ){

      my $SplicePos = $WinStops[$window] % 1000;
      my $SpliceNum = @{$InputAoA[0]} - $SplicePos;
      my @NextAoA1 = ();
      for (my $i = 0; $i < @InputAoA; $i++){
        @line = splice( @{$InputAoA[$i]}, $SplicePos, $SpliceNum );
        push @NextAoA1, [ @line ];
      }

      my @NextAoA2 = ();
      for(my $i=$SplicePos;$i<@new1_AoA;$i++){
         push @NextAoA2, [@{$new1_AoA[$i]}];
      }

      my @NextAoA3 = ();
      for(my $i=$SplicePos;$i<@new2_AoA;$i++){
         push @NextAoA3, [@{$new2_AoA[$i]}];
      }

      @pass_arr=([@InputAoA],[@new1_AoA],[@new2_AoA],[('end')]);
      master_routine(@pass_arr);

      @pass_arr=($SitesAnalyzed,$MaxPBS);
      window_stats(@pass_arr);

      print "running window $window\n";

      last if($WinStops[$window]==$end_win);
      @pass_arr=([@NextAoA1],[@NextAoA2],[@NextAoA3],[('start')]);
      master_routine(@pass_arr);
   }else{
      @pass_arr=([@InputAoA],[@new1_AoA],[@new2_AoA],[('middle')]);
      master_routine(@pass_arr);
   }

   last if($WinStops[$window]==$end_win);
}
chdir($old_dir);
for (my $h = 0; $h < @genomes; $h++){
   close $InHandles[$h];
}

chdir($cwd);
#Output
open O, ">$OutputFile" or die;
print O "Chr\tWinStart\tWinStop\tSitesAnalyzed\tWinPBSs\tMaxPBSs\tpi_$pops[0]\tpi_$pops[1]\tpi_$pops[2]\tDxy_$pops[0]_$pops[1]\tDxy_$pops[0]_$pops[2]\tDxy_$pops[1]_$pops[2]\tFST_$pops[0]_$pops[1]\tFST_$pops[0]_$pops[2]\tFST_$pops[1]_$pops[2]\tmax_snp_location\n";
for (my $i = 0; $i < @WinPBSs; $i++){
  print O "$chr\t$WinStarts[$i]\t$WinStops[$i]\t$WinSiteCounts[$i]\t$WinPBSs[$i]\t$MaxPBSs[$i]\t$PiAoA[$i][0]\t$PiAoA[$i][1]\t$PiAoA[$i][2]\t$DxyAoA[$i][0]\t$DxyAoA[$i][1]\t$DxyAoA[$i][2]\t$FSTAoA[$i][0]\t$FSTAoA[$i][1]\t$FSTAoA[$i][2]\t$max_snp_locations[$i]\n";
}
close O;

open S, ">$SNPFile" or die;
print S "Chr\tSNP_pos\t UnweightedMinor1\t WeightedMinor1\t Minor2\t Minor3\t UnweightedSS\t WeightedSS\tPBS\tFST_$pops[0]_$pops[1]\tFST_$pops[0]_$pops[2]\tFST_$pops[1]_$pops[2]\t$pops[0]_As\t$pops[0]_Cs\t$pops[0]_Gs\t$pops[0]_Ts\t$pops[1]_As\t$pops[1]_Cs\t$pops[1]_Gs\t$pops[1]_Ts\t$pops[2]_As\t$pops[2]_Cs\t$pops[2]_Gs\t$pops[2]_Ts\t ss1\t ss2\t ss3\n";
for (my $i = 0; $i < @SitePBSs; $i++){
  print S "$chr\t$positions[$i]\t$unweighted_minors[$i]\t$weighted_minors[$i]\t$minors2[$i]\t$minors3[$i]\t$unweighted_ss[$i]\t$weighted_ss[$i]\t$SitePBSs[$i]\t$SiteFSTAoA[$i][0]\t$SiteFSTAoA[$i][1]\t$SiteFSTAoA[$i][2]\t$NucCountAoA[$i][0]\t$NucCountAoA[$i][1]\t$NucCountAoA[$i][2]\t$NucCountAoA[$i][3]\t$NucCountAoA[$i][5]\t$NucCountAoA[$i][6]\t$NucCountAoA[$i][7]\t$NucCountAoA[$i][8]\t$NucCountAoA[$i][9]\t$NucCountAoA[$i][10]\t$NucCountAoA[$i][11]\t$NucCountAoA[$i][12]\t$ss1[$i]\t$ss2[$i]\t$ss3[$i]\n";
}
close S;

sub get_initial{
   my @initial_arr=();
   for(my $i=0;$i<1000;$i++){
      push @initial_arr, [(0,0)];
   }
   return(@initial_arr);
}

sub master_routine{

   my @passed_arr=@_;
   my @InputAoA=@{$passed_arr[0]};
   my @new1_AoA=@{$passed_arr[1]};
   my @new2_AoA=@{$passed_arr[2]};
   my $class=$passed_arr[3][0];

   my @Pop2Nucs;
   my @Pop3Nucs;
   my $MajorAllele;
   my $MinorAllele;
   my $SampleSize1;
   my $SampleSize2;
   my $SampleSize3;
   my $MajorFreq1;
   my $MajorFreq2;
   my $MajorFreq3;
   my $MinorFreq1;
   my $MinorFreq2;
   my $MinorFreq3;
   my @FSTs;
   my $FST;
   my $effective_chromo_pool1;
   my $effective_chromo_pool2;
   my @nucs=('A','C','G','T');

   for(my $s=0;$s<@{$InputAoA[0]};$s++){
     if($class eq 'middle' || $class eq 'end'){
       $pos = ($kb * 1000) + $s + 1;
     }elsif($class eq 'start'){
       $pos = ($kb * 1000) + (1000 - @{$InputAoA[0]}) + $s + 1;
     }

      my @AllNucCounts=();
      my @chr_only=@InputAoA[1..$#InputAoA];
      my @pass_arr=(@chr_only,[$s],[$pos]);
      my @return_arr=get_ind_freq(@pass_arr);

      my @snp_counts=@{$return_arr[1]};
      for(my $row=0;$row<@{$return_arr[0]};$row++){
         $InputAoA[$row+1][$s]=$return_arr[0][$row][$s];
      }
      my @inv_freqs=@{$return_arr[2]};

      my @Pop1Nucs=(0,0,0,0,0);
      for(my $i=0;$i<@snp_counts;$i++){
         for(my $j=0;$j<@{$snp_counts[$i]};$j++){
            $Pop1Nucs[$j]+=$snp_counts[$i][$j];
         }
      }
      $Pop1Nucs[0]=int($Pop1Nucs[0]);
      $Pop1Nucs[1]=int($Pop1Nucs[1]);
      $Pop1Nucs[2]=int($Pop1Nucs[2]);
      $Pop1Nucs[3]=int($Pop1Nucs[3]);
      $Pop1Nucs[4]=int($Pop1Nucs[4]);

      my $SampleSize1=$Pop1Nucs[0]+$Pop1Nucs[1]+$Pop1Nucs[2]+$Pop1Nucs[3];#unweighted sample size
      my $unweighted_ss=$SampleSize1;     
 

      next if ($SampleSize1 < $SampleMin);
      
      @Pop2Nucs=get_pool_freq(@{$new1_AoA[$s]});
      @Pop3Nucs=get_pool_freq(@{$new2_AoA[$s]});

      #get major, minor allele, and sample sizes
      @pass_arr=([@Pop1Nucs],[@Pop2Nucs],[@Pop3Nucs]);
      @return_arr=get_major_allele(@pass_arr);
      $MajorAllele=$return_arr[0][0];
      $MinorAllele=$return_arr[0][1];
      $SampleSize1=$return_arr[0][2];
      $SampleSize2=$return_arr[0][3];
      $SampleSize3=$return_arr[0][4];
      my @coverage2=@{$return_arr[1]};
      my @coverage3=@{$return_arr[2]};
      my @freqs_pool2=@{$return_arr[3]};
      my @freqs_pool3=@{$return_arr[4]};

      next if ($SampleSize2 < $SampleMin);
      next if ($SampleSize3 < $SampleMin);

      my $site_diffs=$Pop1Nucs[0]*$Pop1Nucs[1]+$Pop1Nucs[0]*$Pop1Nucs[2]+$Pop1Nucs[0]*$Pop1Nucs[3]+$Pop1Nucs[1]*$Pop1Nucs[2]+$Pop1Nucs[1]*$Pop1Nucs[3]+$Pop1Nucs[2]*$Pop1Nucs[3];
      my $site_comps=$NChoose2[$SampleSize1];
      $Pi1Sum+=$site_diffs;
      $window_comps+=$site_comps;

      my $freq=$Pop1Nucs[$MinorAllele]/($Pop1Nucs[$MajorAllele]+$Pop1Nucs[$MinorAllele]);#unweighted minor freq in old sample

      my $weighted_freq=0;
      my @minor_freqs=();
      my @weights=();
      my $summed_weights=0;
      for(my $i=0;$i<$num_inversions;$i++){
         if($inv_freqs[$i] == 0){
            push @weights, 1;
            $summed_weights++;
         }else{
            push @weights, $inv_freq2[$i]/$inv_freqs[$i];
            $summed_weights+=$inv_freq2[$i]/$inv_freqs[$i];
         }
      }

      for(my $i=0;$i<$num_inversions;$i++){
         if(($snp_counts[$i][$MinorAllele]+$snp_counts[$i][$MajorAllele])!=0){
            $weighted_freq+=($snp_counts[$i][$MinorAllele]/($snp_counts[$i][$MinorAllele]+$snp_counts[$i][$MajorAllele]))*($weights[$i]/$summed_weights);
         }
      }

      my $n_eff_num=0;
      my $n_eff_den=0;
      for(my $i=0;$i<@snp_counts;$i++){
         for(my $j=0;$j<($snp_counts[$i][$MinorAllele]+$snp_counts[$i][$MajorAllele]);$j++){
            $n_eff_num+=$weights[$i];
            $n_eff_den+=$weights[$i]**2;
         }
      }
      my $n_eff=($n_eff_num**2)/$n_eff_den;
      next if ($n_eff < $SampleMin);#for weighted analysis

      $MajorFreq1=1-$weighted_freq;
      $MajorFreq2=$freqs_pool2[$MajorAllele]/($freqs_pool2[$MajorAllele]+$freqs_pool2[$MinorAllele]);
      $MajorFreq3=$freqs_pool3[$MajorAllele]/($freqs_pool3[$MajorAllele]+$freqs_pool3[$MinorAllele]);
      $MinorFreq1=1-$MajorFreq1;

      $MinorFreq2=1-$MajorFreq2;
      $MinorFreq3=1-$MajorFreq3;

      $Pi1Sum += 2 * $MajorFreq1 * $MinorFreq1;
      $Pi2Sum += 2 * $MajorFreq2 * $MinorFreq2;
      $Pi3Sum += 2 * $MajorFreq3 * $MinorFreq3;
      $Dxy12Sum += $MajorFreq1 * $MinorFreq2 + $MajorFreq2 * $MinorFreq1;
      $Dxy13Sum += $MajorFreq1 * $MinorFreq3 + $MajorFreq3 * $MinorFreq1;
      $Dxy23Sum += $MajorFreq3 * $MinorFreq2 + $MajorFreq2 * $MinorFreq3;

      #Fst
      my @FSTs = ();
      @pass_arr=($MinorFreq1,$MinorFreq2,$n_eff,$SampleSize2,1,2);#fst with weighted freq in pop1
      $FST=fst(@pass_arr);
      push @FSTs, $FST;

      @pass_arr=($MinorFreq1,$MinorFreq3,$n_eff,$SampleSize3,1,3);#fst with weighted freq in pop1
      $FST=fst(@pass_arr);
      push @FSTs, $FST;

      @pass_arr=($MinorFreq2,$MinorFreq3,$SampleSize2,$SampleSize3,2,3);
      $FST=fst(@pass_arr);
      push @FSTs, $FST;
      #if site passes frequency threshold, record info for site output
      next if ((($MinorFreq1 + $MinorFreq2 + $MinorFreq3) / 3) < $FreqThresh);

      #####################################################
      push @AllNucCounts, @Pop1Nucs;
      push @AllNucCounts, @freqs_pool2;
      push @AllNucCounts, @freqs_pool3;
      $SitesAnalyzed++;
      push @unweighted_minors, $freq;
      push @unweighted_ss, $SampleSize1;
      push @weighted_ss, $n_eff;
      push @weighted_minors, $MinorFreq1;
      push @minors2, $MinorFreq2;
      push @minors3, $MinorFreq3;
      #####################################################

      push @positions, $pos;
      push @SiteWins, $window;
      push @SiteFSTAoA, [ @FSTs ];
      push @NucCountAoA, [ @AllNucCounts ];
      $SitePBS = ((-1 * log (1 - $FSTs[0])) + (-1 * log (1 - $FSTs[1])) - (-1 * log (1 - $FSTs[2]))) / 2;
      if ($SitePBS > $MaxPBS & $n_eff>=50 & abs($MinorFreq1 - $freq)<0.1 & $SampleSize2 >= 50 & $SampleSize3 >= 50){
         $MaxPBS = $SitePBS;
         $max_snp_location=$pos;
      }
      push @SitePBSs, $SitePBS;
      push @ss1, $SampleSize1;
      push @ss2, $SampleSize2;
      push @ss3, $SampleSize3;

      if($pos != $new1_AoA[$s][1] || $pos != $new2_AoA[$s][1] || $new2_AoA[$s][1] != $new1_AoA[$s][1]){
         print "$pos\n";
         print "@{$new1_AoA[$s]}\n";
         print "@{$new2_AoA[$s]}\n";
         print '($kb * 1000) + (1000 - @{$InputAoA[0]}) + $s + 1';
         print "\n";
         my $num=@{$InputAoA[0]};
         print "s $s kb $kb len $num\n\n";
         die "mismatch of positions!\n";
      }
   }

}

sub get_ind_freq{
   my @passed_arr=@_;
   my $s=$passed_arr[-2][0];
   my $pos=$passed_arr[-1][0];
   pop @passed_arr;
   pop @passed_arr;
   my @InputAoA=@passed_arr;
   my $As=0;
   my $Cs=0;
   my $Gs=0;
   my $Ts=0;
   my $Ns=0;
   my $num_chr=0;
   my @inv_freq=();
   for(my $i=0;$i<$num_inversions;$i++){
      push @inv_freq, 0;
   }
   my @snp_counts=();
   for(my $i=0;$i<$num_inversions;$i++){
      push @snp_counts, [(0,0,0,0,0)];
   }

   for(my $i = 0; $i < @InputAoA; $i++){

     my $region;
     for(my $h=0;$h<@{$het_regions[$i]};$h++){
       if($pos >= $het_regions[$i][$h][1] && $pos < $het_regions[$i][$h][2]){
         $region=$het_regions[$i][$h][0];
         last;
       }
     }
     my $counts;
     my $r;
     if($region eq 'I'){
       $counts=1;
       $r=rand(1);
     }else{
       $counts=2;
     }

     my $inv_kary=4;
     for(my $inv=0;$inv<$num_inversions;$inv++){
        if($lines[$i] ~~ @{$inv_lines[$inv]}){
           $inv_kary=$inv;

           if($InputAoA[$i][$s] ne 'N'){
              $inv_freq[$inv_kary]++;
              $inv_freq[4]++;
              $num_chr+=2;
           }
           last;
        }  
     }        

     if($InputAoA[$i][$s] eq 'A'){
        $snp_counts[$inv_kary][0]+=$counts;
        $As+=$counts;
     }elsif ($InputAoA[$i][$s] eq 'C'){
        $snp_counts[$inv_kary][1]+=$counts;
        $Cs+=$counts;
     }elsif ($InputAoA[$i][$s] eq 'G'){
        $snp_counts[$inv_kary][2]+=$counts;
        $Gs+=$counts;
     }elsif ($InputAoA[$i][$s] eq 'T'){
        $snp_counts[$inv_kary][3]+=$counts;
        $Ts+=$counts;
     }elsif($InputAoA[$i][$s] eq 'N'){
       $snp_counts[$inv_kary][4]+=$counts;
       $Ns+=$counts;
     }elsif($InputAoA[$i][$s] eq 'M'){
       if($region eq 'I'){
         if($r<=0.5){
           $InputAoA[$i][$s]='A';
         }else{
           $InputAoA[$i][$s]='C';
         }
         $As+=0.5;
         $snp_counts[$inv_kary][0]+=0.5;
         $Cs+=0.5;
         $snp_counts[$inv_kary][1]+=0.5;
       }else{
         $snp_counts[$inv_kary][0]++;
         $snp_counts[$inv_kary][1]++;
         $As++;
         $Cs++;
       }
     }elsif($InputAoA[$i][$s] eq 'R'){
       if($region eq 'I'){
         if($r<=0.5){
           $InputAoA[$i][$s]='A';
         }else{
           $InputAoA[$i][$s]='G';
         }
         $As+=0.5;
         $snp_counts[$inv_kary][0]+=0.5;
         $Gs+=0.5;
         $snp_counts[$inv_kary][2]+=0.5;
       }else{
         $As++;
         $Gs++;
         $snp_counts[$inv_kary][0]++;
         $snp_counts[$inv_kary][2]++;
       }
     }elsif($InputAoA[$i][$s] eq 'W'){
       if($region eq 'I'){
         if($r<=0.5){
           $InputAoA[$i][$s]='A';
         }else{
           $InputAoA[$i][$s]='T';
         }
         $As+=0.5;
         $snp_counts[$inv_kary][0]+=0.5;
         $Ts+=0.5;
         $snp_counts[$inv_kary][3]+=0.5;
       }else{
         $As++;
         $Ts++;
         $snp_counts[$inv_kary][0]++;
         $snp_counts[$inv_kary][3]++;
       }
     }elsif($InputAoA[$i][$s] eq 'S'){
       if($region eq 'I'){
         if($r<=0.5){
           $InputAoA[$i][$s]='C';
         }else{
           $InputAoA[$i][$s]='G';
         }
         $Cs+=0.5;
         $snp_counts[$inv_kary][1]+=0.5;
         $Gs+=0.5;
         $snp_counts[$inv_kary][2]+=0.5;
       }else{
         $Cs++;
         $Gs++;
         $snp_counts[$inv_kary][1]++;
         $snp_counts[$inv_kary][2]++;
       }
     }elsif($InputAoA[$i][$s] eq 'Y'){
       if($region eq 'I'){
         if($r<=0.5){
           $InputAoA[$i][$s]='C';
         }else{
           $InputAoA[$i][$s]='T';
         }
         $Cs+=0.5;
         $snp_counts[$inv_kary][1]+=0.5;
         $Ts+=0.5;
         $snp_counts[$inv_kary][3]+=0.5;
       }else{
         $Cs++;
         $Ts++;
         $snp_counts[$inv_kary][1]++;
         $snp_counts[$inv_kary][3]++;
       }
     }elsif($InputAoA[$i][$s] eq 'K'){
       if($region eq 'I'){
         if($r<=0.5){
           $InputAoA[$i][$s]='G';
         }else{
           $InputAoA[$i][$s]='T';
         }
         $Gs+=0.5;
         $snp_counts[$inv_kary][2]+=0.5;
         $Ts+=0.5;
         $snp_counts[$inv_kary][3]+=0.5;
       }else{
         $Ts++;
         $Gs++;
         $snp_counts[$inv_kary][2]++;
         $snp_counts[$inv_kary][3]++;
       }
     }
   }

   if($num_chr==0){
      $num_chr=1;
   }
   for(my $i=0;$i<@inv_freq;$i++){
      $inv_freq[$i]=$inv_freq[$i]/$num_chr;
   }
   my @return_arr=([@InputAoA],[@snp_counts],[@inv_freq]);

   return(@return_arr);
}

sub get_pool_freq{
    #our script is A:C:G:T:N
    #popoolation sync is A:T:C:G
    #fall data will have length 9
    #spring data will have length 15

    my @pops=@_;

    my $num_pools=@pops-3;
    my @Pop2Nucs=(0,0,0,0);
    my @Pop2Nucs_AoA=();
    for(my $i=0;$i<$num_pools;$i++){
       push @Pop2Nucs_AoA, [@Pop2Nucs];
    }

    for(my $i=3;$i<@pops;$i++){
      my @arr=split /:/, $pops[$i];
      $Pop2Nucs_AoA[$i-3][0]=$arr[0];
      $Pop2Nucs_AoA[$i-3][1]=$arr[2];
      $Pop2Nucs_AoA[$i-3][2]=$arr[3];
      $Pop2Nucs_AoA[$i-3][3]=$arr[1];
   }
   return(@Pop2Nucs_AoA);
}

sub get_major_allele{
   my @Pop1Nucs=@{$_[0]};
   my @Pop2Nucs=@{$_[1]};#fall
   my @Pop3Nucs=@{$_[2]};#spring
   my $SampleSize1=$Pop1Nucs[0]+$Pop1Nucs[1]+$Pop1Nucs[2]+$Pop1Nucs[3];
   my @SampleSize2=(0,0,0,0,0,0);
   my @coverage2=(0,0,0,0,0,0);
   my $summed_pool2=0;
   for(my $i=0;$i<@Pop2Nucs;$i++){
     my $coverage=0;
     foreach(@{$Pop2Nucs[$i]}){
       $coverage+=$_;
     }
     my $mean_cov=$mean_coverage_fall[$i][$chr_hash{$chr}];

     if($coverage>=0.33*$mean_cov && $coverage<=3*$mean_cov){
       $coverage2[$i]=$coverage;
       $SampleSize2[$i]=$expected_unique_arr[$i][$coverage];
     }
     $summed_pool2+=$SampleSize2[$i];
   }

   my @SampleSize3=(0,0,0,0,0,0,0,0,0,0,0,0);
   my @coverage3=(0,0,0,0,0,0,0,0,0,0,0,0);
   my $summed_pool3=0;
   for(my $i=0;$i<@Pop3Nucs;$i++){
     my $coverage=0;
     foreach(@{$Pop3Nucs[$i]}){
       $coverage+=$_;
     }

     my $mean_cov=$mean_coverage_spring[$i][$chr_hash{$chr}];
     if($coverage>=0.5*$mean_cov && $coverage<=2*$mean_cov){
       $coverage3[$i]=$coverage;
       $SampleSize3[$i]=$expected_unique_arr[$i+6][$coverage];
     }
     $summed_pool3+=$SampleSize3[$i];
   }

   #in pool, get frequency at read level, multiply by samplesize[i]/sum(samplesize)
   my $AFreqSum_pool2=0;
   my $CFreqSum_pool2=0;
   my $GFreqSum_pool2=0;
   my $TFreqSum_pool2=0;
   for(my $i=0;$i<@Pop2Nucs;$i++){
      if($coverage2[$i]>0){
         $AFreqSum_pool2+=($Pop2Nucs[$i][0]/$coverage2[$i])*($SampleSize2[$i]/$summed_pool2);
         $CFreqSum_pool2+=($Pop2Nucs[$i][1]/$coverage2[$i])*($SampleSize2[$i]/$summed_pool2);
         $GFreqSum_pool2+=($Pop2Nucs[$i][2]/$coverage2[$i])*($SampleSize2[$i]/$summed_pool2);
         $TFreqSum_pool2+=($Pop2Nucs[$i][3]/$coverage2[$i])*($SampleSize2[$i]/$summed_pool2);
      }
   }
   my @freqs_pool2=($AFreqSum_pool2,$CFreqSum_pool2,$GFreqSum_pool2,$TFreqSum_pool2);

   my $AFreqSum_pool3=0;
   my $CFreqSum_pool3=0;
   my $GFreqSum_pool3=0;
   my $TFreqSum_pool3=0;
   for(my $i=0;$i<@Pop3Nucs;$i++){
      if($coverage3[$i] > 0){
         $AFreqSum_pool3+=($Pop3Nucs[$i][0]/$coverage3[$i])*($SampleSize3[$i]/$summed_pool3);
         $CFreqSum_pool3+=($Pop3Nucs[$i][1]/$coverage3[$i])*($SampleSize3[$i]/$summed_pool3);
         $GFreqSum_pool3+=($Pop3Nucs[$i][2]/$coverage3[$i])*($SampleSize3[$i]/$summed_pool3);
         $TFreqSum_pool3+=($Pop3Nucs[$i][3]/$coverage3[$i])*($SampleSize3[$i]/$summed_pool3);
      }
   }
   my @freqs_pool3=($AFreqSum_pool3,$CFreqSum_pool3,$GFreqSum_pool3,$TFreqSum_pool3);

   #Define major and minor alleles (mask third/fourth alleles)
   my $AFreqSum = ($Pop1Nucs[0] / $SampleSize1) + $AFreqSum_pool2 + $AFreqSum_pool3;
   my $CFreqSum = ($Pop1Nucs[1] / $SampleSize1) + $CFreqSum_pool2 + $CFreqSum_pool3;
   my $GFreqSum = ($Pop1Nucs[2] / $SampleSize1) + $GFreqSum_pool2 + $GFreqSum_pool3;
   my $TFreqSum = ($Pop1Nucs[3] / $SampleSize1) + $TFreqSum_pool2 + $TFreqSum_pool3;
   my @FreqSums = ($AFreqSum, $CFreqSum, $GFreqSum, $TFreqSum);

   my @SortedValues = @FreqSums;
   @SortedValues = sort { $b <=> $a } @SortedValues;
   my $MajorAllele = -1;
   my $MinorAllele;
   for(my $j = 0; $j < @FreqSums; $j++){
      if (($FreqSums[$j] == $SortedValues[0]) && ($MajorAllele < -0.5)){
         $MajorAllele = $j;
         next;
      }
      if ($FreqSums[$j] == $SortedValues[1]){
         $MinorAllele = $j;
         next;
      }
   }
   my $SampleSize1 = $Pop1Nucs[$MajorAllele] + $Pop1Nucs[$MinorAllele];

   my $SampleSize2=$summed_pool2;
   if($SampleSize2!=0){
   }else{
     $SampleSize2=-99999;
   }
   my $SampleSize3=$summed_pool3;
   if($SampleSize3 != 0){
   }else{
     $SampleSize3=-99999;
   }
   my @return_arr=([$MajorAllele,$MinorAllele,$SampleSize1,$SampleSize2,$SampleSize3],[@coverage2],[@coverage3],[@freqs_pool2],[@freqs_pool3]);
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
  push @WinSiteCounts, $SitesAnalyzed;
   if ($SitesAnalyzed == 0){
      @pis = (-999, -999, -999);
      push @PiAoA, [ @pis ];
      push @DxyAoA, [ @pis ];
   }
   else{
      $Pi1Sum = $Pi1Sum / $window_comps;
      $Pi2Sum = $Pi2Sum / $SitesAnalyzed;
      $Pi3Sum = $Pi3Sum / $SitesAnalyzed;
      @pis = ($Pi1Sum, $Pi2Sum, $Pi3Sum);
      push @PiAoA, [ @pis ];
      $Dxy12Sum = $Dxy12Sum / $SitesAnalyzed;
      $Dxy13Sum = $Dxy13Sum / $SitesAnalyzed;
      $Dxy23Sum = $Dxy23Sum / $SitesAnalyzed;
      my @Dxys = ($Dxy12Sum, $Dxy13Sum, $Dxy23Sum);
      push @DxyAoA, [ @Dxys ];
   }
   $SitesAnalyzed = 0;
   push @MaxPBSs, $MaxPBS;
   push @max_snp_locations, $max_snp_location;
   $MaxPBS = 0;
   if (($Den12Sum == 0) || ($Den13Sum == 0) || ($Den23Sum == 0)){
      push @WinPBSs, -999;
      push @WinZs, -999;
      my @FSTs = (-999, -999, -999);
      push @FSTAoA, [ @FSTs ];
   }
   else{
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
      push @FSTAoA, [ @FSTs ];
      my $WinPBS = ((-1 * log (1 - $FSTs[0])) + (-1 * log (1 - $FSTs[1])) - (-1 * log (1 - $FSTs[2]))) / 2;
      push @WinPBSs, $WinPBS;
      my $Z = -1 * log (1 - $FSTs[2]);
      push @WinZs, $Z;
   }
   $Pi1Sum = 0;
   $window_comps=0;
   $Pi2Sum = 0;
   $Pi3Sum = 0;
   $Dxy12Sum = 0;
   $Dxy13Sum = 0;
   $Dxy23Sum = 0;
   $Num12Sum = 0;
   $Num23Sum = 0;
   $Num13Sum = 0;
   $Den12Sum = 0;
   $Den23Sum = 0;
   $Den13Sum = 0;
   $window++;
}

sub get_expectations{
  my $input="uniquedraws_expectations_updatedspringss.txt";
  open I, "<$input" or die "cannot find $input";
  my @AoA=();
  while(<I>){
    chomp;
    my @line=split /\s+/, $_;
    push @AoA, [@line];
  }
  close I;
  return(@AoA);
}

sub get_coverage{
  my $input=shift;
  open I, "<$input" or die "cannot find $input";
  my @AoA=();
  while(<I>){
    chomp;
    my @line=split /\s+/, $_;
    push @AoA, [@line];
  }
  close I;
  return(@AoA);
}

