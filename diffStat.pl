#!/usr/bin/perl -w
use strict;

die("Usage: perl $0 <otPre> <cutoff> <hap1> <hap2>....

#==============================================================================
@
#Used to determine the real haplotype by compare each to the major haplotype
#Authors: Fang Zhiwei (fangzhiwei1126\@126.com
#
#==============================================================================
    otPre: prefix of output file'
    cutof: the SNP number compared with major genotype when a minor genotype determined to be real
    hap1: record the reads bunmber supported each haplotype; in this formated:
        ampliconID  readsNum1   haplotype1
        ampliconID  readsNum2   haplotype2

    Revised: add the parameter of SNPnumber; at 20161229;
    
    
    Revised: If a less abundant (i.e, abd1) genotype has no SNP with a higher (i.e. abd2) 
             abundant genotype in the UNCOVERED regions, then the less one will be deem 
             the same as the higher one and (abd1+abd2) will be the abundance of higher
             abundant genotype. Actually, we ignored indels in this region. at 20161223;
        \n") if @ARGV<2;
my($otPre,$cutoff,@hapF)=@ARGV;


my%SSR=();
open IN, " /result/pub/software/software/bin/SSR/pipeline//amplicon_SSR.list " or die "cannot open  ../amplicon_SSR.list $!";
while(<IN>){
    chomp;
    $_=(split /\s+/,$_)[0];
    $SSR{$_}="";
}
close IN;

my($sam,$aid,$lastAid,$abd,$type,$samID)=();
my%h=();

foreach my$inF(@hapF){
    %h=();
    $samID=(split /\./, (split /\//,$inF)[-1])[0];

    open IN, " $inF " or die  "cannto open $inF ";
    chomp(my$line=<IN>);
    ($sam,$lastAid,$abd,$type)=(split /\s+/,$line)[0,1,2,3];
######$h{$type}=$abd if ($abd>=5 and defined $type);
######$h{$type}=$abd if defined $type;#\@20180901;
    $h{$type}=$abd if ($abd>=2 and defined $type);#\@20180907, one reads could be ignored;

    while(<IN>){
        chomp;
        next if /^#*\s*\d*\s*$/;
        ($sam,$aid,$abd,$type)=(split /\s+/,$_);

####### next if $abd<5;### minor haplotype should be have at least 5 reads surpported; \@20180901;
        next if $abd<2;### minor haplotype should be have at least 2 reads surpported; \@20180907;
        
        
        if($aid eq $lastAid){
            $h{$type}=$abd;
            next;
        }
        &process();
    }
    &process();
    close IN;
}


################  sub program #######
sub process{

##  No minor haplotype exist;
        my($majorType,@allType)=(sort {$h{$b}<=>$h{$a}} keys %h);
        unless (defined $majorType){
            %h=();
            $lastAid=$aid;
            $h{$type}=$abd;
            return(1);
        }
        
###  Combine some genotype, see revise1 in usage;        
        my%diffSite=();
        my%samType=();
        unshift @allType,$majorType;
        foreach my$tmp_type(@allType){
            my@allSite=(split /\|/,$tmp_type);
            foreach my$Pos_tmp(@allSite){
                next if $Pos_tmp=~/,0/;
                next if $Pos_tmp=~/,Re*f*/;
                $Pos_tmp=(split /\,/,$Pos_tmp)[0];
                $diffSite{$Pos_tmp}="";
            }
        }

        my$ko2pos="";
        for(my$km=0;$km<$#allType;$km++){
            next if exists $samType{$allType[$km]};            
            my@site4km=(split /\|/,$allType[$km]);

            for(my$kn=$km+1;$kn<=$#allType;$kn++){
                next if exists $samType{$allType[$kn]};                
                my@site4kn=(split /\|/,$allType[$kn]);
                
                my$lastFlag=0;
                for(my$ko=0;$ko<=$#site4km;$ko++){
                    next if $site4km[$ko] eq $site4kn[$ko];
                    $ko2pos=(split /\,/,$site4km[$ko])[0]; ## this should be judged befor next two lines;
                    $lastFlag =1 if(exists $diffSite{$ko2pos});
                    last if(exists $diffSite{$ko2pos});
                    next if $site4km[$ko]=~/,0/;
                    next if $site4kn[$ko]=~/,0/;
                    $lastFlag =1;
                    last;
                }

                next if $lastFlag==1;
                $h{$allType[$km]} +=$h{$allType[$kn]};
                delete $h{$allType[$kn]};
                $samType{$allType[$kn]}="";
            }
        }


###   No minor haplotype exist;
        ($majorType,@allType)=(sort {$h{$b}<=>$h{$a}} keys %h);
        if(@allType ==0){
            my$majorType_=$majorType;
            $majorType_=~s/\|\d+,//g;
            $majorType_=~s/Ref/R/g;
            $majorType_=~s/delete/D/g;
            $majorType_=~s/insert/I/g;
            $majorType_=~s/snp/S/g;
            print "$samID\t$lastAid\t0\t$h{$majorType}\t$majorType_\n" unless exists $SSR{$lastAid};
            print "$samID\t$lastAid\t0\t$h{$majorType}\t$majorType_\tSSR\n" if exists $SSR{$lastAid};
            %h=();
            $lastAid=$aid;
            $h{$type}=$abd;
            return(1);
        }
       
##  record all genotype with more than two diff site with major genotype;
        my@majorSite=(split /\|/,$majorType);
        my%record=();

        foreach my$otherType (@allType){
            my@otherSite=(split /\|/,$otherType);
            my$n=0;

            for(my$i=0;$i<=$#otherSite;$i++){
                next if $majorSite[$i]=~/[ID]/;
                next if $majorSite[$i]=~/insert/;
                next if $majorSite[$i]=~/delete/;
                next if $otherSite[$i]=~/[ID]/;
                next if $otherSite[$i]=~/insert/;
                next if $otherSite[$i]=~/delete/;
                next if $majorSite[$i]=~/\d+,0/;
                next if $otherSite[$i]=~/\d+,0/;
                $n++ if $majorSite[$i] ne $otherSite[$i];
            }
            next if $n<=$cutoff;
            $record{$otherType}=$h{$otherType} ;
        }

##   out put all minor genotype have at least 5 reads surpported and >=2 diff sites with major type;
        my($ot,$otabd)=("","");
        my$m=0;
        
        foreach (sort {$record{$b}<=>$record{$a}} keys %record){
            $ot .="$_;";
            $otabd.="$record{$_},";
            $m++;
        }
        my$majorType__=$majorType;
        $majorType__=~s/\|\d+,//g;
        $majorType__=~s/Ref/R/g;
        $majorType__=~s/delete/D/g;
        $majorType__=~s/insert/I/g;
        $majorType__=~s/snp/S/g;

        $ot=~s/\|\d+,//g;
        $ot=~s/Ref/R/g;
        $ot=~s/delete/D/g;
        $ot=~s/insert/I/g;
        $ot=~s/snp/S/g;

        print "$samID\t$lastAid\t$m\t$h{$majorType}\t$majorType__\t$otabd\t$ot\n" unless exists $SSR{$lastAid};;
        print "$samID\t$lastAid\t$m\t$h{$majorType}\t$majorType__\t$otabd\t$ot\tSSR\n" if exists $SSR{$lastAid};;
        %record=();

        %h=();
        $lastAid=$aid;
        $h{$type}=$abd;
        return(1);
}

