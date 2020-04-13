#!/usr/bin/perl -w
use strict;

die("Usage: perl $0 < abd.stat.xls1> <abd.stat.xls2> .....\
=================================================================\
#\
#Used to check the real genotype that have at least one SNP compared with major genotype for each amplicon site;\
#Authors: Fang Zhiwei (fangzhiwei1126\@126.com)\
#\
=================================================================\n
Revised2: To reduce the memary required, we clear all hash \%snp and \%gt foreach cycle at line 20 and 21; if you
          want to output Parent2F1.list, please commented line 21 (\%snp,\%gt)=();     \@20180510;

Revised: All SNP were ignored if two or more SNPs site were adjacent to each other; \@20171130;\n ") if @ARGV <1;


my(%snp,%gt)=();
my%pos=();
my($pre,$sam,$aid,$gtNum,$majorAbd,$major,$otherAbd,$otherGT)=();
foreach(@ARGV){
    (%snp,%gt)=();
    $pre=(split /\./,$_)[0];
    my$f=$_; 
    
    print STDERR "Processing $f at ".`date`;

    open IN , "$f " or die "cannot open $f $!";
    open OT,"> $f.check " or die "cannot write $f.check $!";
    while(<IN>){
        chomp; 
        ($sam,$aid,$gtNum,$majorAbd,$major,$otherAbd,$otherGT)=(split /\s+/,$_);
        $gt{$aid}{$pre}{"M"}="$major\t$majorAbd";
        if($gtNum==0){
            $gt{$aid}{$pre}{"s0"}="-\t-";
            $gt{$aid}{$pre}{"s1"}="-\t-"; 
            $snp{$aid}{$pre}{"s0"}="-";
            $snp{$aid}{$pre}{"s1"}="-";
            print OT "$_\n";
            next;
        }
        my@abd=(split /\,/,$otherAbd); 
        my@genotypes=(split /\;/,$otherGT);
        my(@abd_tmp,@genotypes_tmp)=();
        $gtNum=0;
        for(my$ii=0;$ii<=$#abd;$ii++){
            next if ($abd[$ii]/$majorAbd)<=0.02;
            $gtNum++;
            push @abd_tmp,$abd[$ii];
            push @genotypes_tmp,$genotypes[$ii];
        }
        if($gtNum==0){
            $gt{$aid}{$pre}{"s0"}="-\t-";
            $gt{$aid}{$pre}{"s1"}="-\t-";
            $snp{$aid}{$pre}{"s0"}="-";
            $snp{$aid}{$pre}{"s1"}="-";
            print OT "$_\n";
            next;
        }

        @abd=@abd_tmp;
        @genotypes=@genotypes_tmp;

        unshift @abd,$majorAbd;
        unshift @genotypes,$major;

        my$count=-1; 
        my($abdList,$gtList)=("","");
        for(my$i=1;$i<=$#genotypes;$i++){
            my$snpNum=-1;
            next if $genotypes[$i]=~/snp:N/;
            for(my$j=$i-1;$j>=0;$j--){
                next if $genotypes[$j]=~/snp:N/;
                %pos=();
                $snpNum=&compareGT($genotypes[$j],$genotypes[$i]);
                next if $snpNum!=0;
                last;
            }

            next if $snpNum==0;
            $count ++;
            $gt{$aid}{$pre}{"s$count"}="$genotypes[$i]\t$abd[$i]";
            $snp{$aid}{$pre}{"s$count"}=$snpNum;
            $abdList.="$abd[$i],";
            $gtList .="$genotypes[$i];";            
        }
        if($count ==-1){
            print OT "$sam\t$aid\t0\t$majorAbd\t$major\n";
        }else{
            $count ++;
            print OT"$sam\t$aid\t$count\t$majorAbd\t$major\t$abdList\t$gtList\n";
        }


        for(my$kk=1;$kk>=0;$kk--){  ### output two second genotypes;
            $gt{$aid}{$pre}{"s$kk"}="-\t-" unless exists $gt{$aid}{$pre}{"s$kk"};
            $snp{$aid}{$pre}{"s$kk"}="-"   unless exists $snp{$aid}{$pre}{"s$kk"};
        }
    }
    close IN ;
    close OT;
}

#open IN , "Parent2F1.list"; 
#
#while(<IN>){
#    chomp;
#    my($p1,$p2,$f1)=(split /\s+/,$_);
#    next if $f1=~/NMR5/;
#    foreach my$aid (keys %gt){
#        next unless exists $gt{$aid}{$f1}{"M"};
#        next unless exists $gt{$aid}{$p1}{"M"};
#        next unless exists $gt{$aid}{$p2}{"M"};
#        print "$aid\t",$gt{$aid}{$p1}{"M"}."\t",$gt{$aid}{$p1}{"s0"}."\t".$snp{$aid}{$p1}{"s0"}."\t$p1";
#        print     "\t",$gt{$aid}{$p2}{"M"}."\t",$gt{$aid}{$p2}{"s0"}."\t".$snp{$aid}{$p2}{"s0"}."\t$p2";
#        print     "\t",$gt{$aid}{$f1}{"M"}."\t",$gt{$aid}{$f1}{"s0"}."\t".$snp{$aid}{$f1}{"s0"};
#        print                              "\t",$gt{$aid}{$f1}{"s1"}."\t".$snp{$aid}{$f1}{"s1"}."\t$f1\n";;
#    }
#}

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
        
        if($genotype=~/^[DI]:[ATCGN]{1,},S:[ATCGN]/){
            $genotype=~s/^([DI]:[ATCGN]{1,},S:[ATCGN])//;
            $posList{$cyc}="indel";
            next;
        }
        
        if($genotype=~/^[DI]:/){
            $genotype=~s/^([DI]:[ATCGN]{1,})//;
            $posList{$cyc}="indel";
            next;
        }
        
        if($genotype=~/^S:/){
            $genotype=~s/^(S:[ATCGN])//;
            $posList{$cyc}="snp";
        }
    }

    my$indel=0;
    foreach (sort {$a<=>$b} keys %posList){
        $indel=0;
        if( $posList{$_} eq "indel"){
            $pos{$_}{$order}="ignored";
            next;
        }

        my$tmp=$_+1;
        if(exists $posList{$tmp} ){
            $indel ++ if $posList{$tmp} eq "indel";
            $indel ++ if $posList{$tmp} eq "snp";
            }

        $tmp=$_-1;
        if(exists $posList{$tmp} ){
            $indel ++ if $posList{$tmp} eq "indel";
            $indel ++ if $posList{$tmp} eq "snp";
            }
         
        $pos{$_}{$order}="snp"  if $indel ==0;
        $pos{$_}{$order}="ignored" if $indel >0;
#        delete $posList{$_} if $indel>0;
    }
}
