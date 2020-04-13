#!/usr/bin/perl -w
use strict;

die("Usage: perl $0  <gtPool> <minMajorCount> <ishybrid.err1> <ishybrid.err2> ...\n
Used to filter genotype that not included in gtPool\n") if @ARGV <3;
my($gtPool,$minMajorCount,@files)=@ARGV;

#### get db 
my%all_GT=();
my%pos=();
open IN, "statMajorDistri.db.log " or die  "cannot open $!";
while(<IN>){
    chomp;
    next if /^#/;
##genotype       aid     num4Sample      majorCount      majorAbdList    stuterCount     stutterAbdList  stutterRatioList
    my($gt,$aid,$Num4sample,$majorCount,$majorAbdList,$stuterCount,$stutterAbdList,$stutterRatioList)=(split /\s+/,$_);
    next if $majorCount <$minMajorCount;
    $all_GT{$aid}{$gt}="";
}
close IN;

foreach my$inf(@files){
#    my$tmp_f="$inf\_bak";
#    print `mv $inf $tmp_f`;
print STDERR "Processing $inf \n";
    open IN,"$inf " or die "cannot open $inf $!";
    open OT, "> $inf.f " or die "cannot write $inf.f $!";
    while(<IN>){
        chomp;
        my($sam,$aid,$gt,$nul,$majorAbd,$stutterAbd,$major,$stutter)=(split /\s+/,$_);
        if($gt!~/hybrid/){
            print OT "$_\n";
            next;
        }
        
        if(exists $all_GT{$aid}{$stutter}){
            print OT "$_\n";
            next;
        }
        if(($stutterAbd/$majorAbd)>0.7){
            print OT "$_\n";
            next;
        }
        print OT "$sam\t$aid\tNA\t-\t$majorAbd\t0\t$major\t-\n";
        
=cut do not ignore indel and lianxu SNP; 
        my$diff=1;
        foreach my$gt_t(keys %{$all_GT{$aid}}){
            $diff=&compareGT($stutter,$gt_t);
            last if $diff==0;
        }

        if($diff==0){
            print OT "$_\n";
        }else{
            print OT "$sam\t$aid\tinbred\t-\t$majorAbd\t0\t$major\t-\n";
        }
=cut
    }
    close IN;
    close OT;
print STDERR "Finishing $inf \n";
}


sub compareGT{
    my($major_,$seGT_)=@_;
    %pos=();
    &getSNP($major_,"M");
    &getSNP($seGT_,"S");

    my$snpNum_=0;
    foreach my$snpPos(sort {$a<=>$b} keys %pos){
        my@tmp= keys %{$pos{$snpPos}};
        if(@tmp==1){
            $snpNum_ ++ if $pos{$snpPos}{$tmp[0]} ne "ignored";
            next;
        }

        next if $pos{$snpPos}{"M"} eq "ignored"; 
        next if $pos{$snpPos}{"S"} eq "ignored";
        $snpNum_ ++ if $pos{$snpPos}{"M"} ne $pos{$snpPos}{"S"};
    }
    return ($snpNum_);
}

sub getSNP{
    my($genotype,$order)=@_;
    $genotype=~s/^\d+,//;
    $genotype=~s/\|$//;
    
    my$genotype_=$genotype;
    my$cyc=0;
    my%posList=();
    while($genotype ne ""){
        $cyc++;
        if($genotype=~/^0/){
            $genotype=~s/^0//;
            $posList{$cyc}="indel";
            next;
        }

        if($genotype=~/^R/){
            $genotype=~s/^R//;
            next;
        }
        
        if($genotype=~/^[DI]:[ATCG]{1,},S:[ATCG]/){
            $genotype=~s/^([DI]:[ATCG]{1,},S:[ATCG])//;
            $posList{$cyc}="indel";
            next;
        }
        
        if($genotype=~/^[DI]:/){
            $genotype=~s/^([DI]:[ATCG]{1,})//;
            $posList{$cyc}="indel";
            next;
        }
        
        if($genotype=~/^S:/){
            $genotype=~s/^(S:[ATCG])//;
            $posList{$cyc}="snp";
        }
    }

    my$indel=0;
    foreach (sort {$a<=>$b} keys %posList){
        $indel=0;
        if( $posList{$_} eq "indel"){
            $pos{$_}{$order}="ignored";
            next;
        }else{
            $pos{$_}{$order}=$posList{$_};
        }
        next;

        
        for(my$i=1;$i<=1;$i++){
            my$tmp=$_+$i;
            if(exists $posList{$tmp} ){
                if($i<=5){
                    $indel ++ if $posList{$tmp} eq "snp";
                    $indel ++ if $posList{$tmp} eq "indel";
                }else{
                    $indel ++ if $posList{$tmp} eq "indel";
                }
                }
                }
        
        for(my$i=1;$i<=1;$i++){
            my$tmp=$_-$i;
            if(exists $posList{$tmp} ){
                if($i<=5){
                    $indel ++ if $posList{$tmp} eq "snp";
                    $indel ++ if $posList{$tmp} eq "indel";
                }else{
                    $indel ++ if $posList{$tmp} eq "indel";
                }
                }
                }
         
        $pos{$_}{$order}="snp"  if $indel ==0;
        $pos{$_}{$order}="ignored" if $indel >0;
#        delete $posList{$_} if $indel>0;
    }
}
