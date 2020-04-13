#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my$usage="Usage: perl $0 -i1 <compareHap1> -i2  <compareHap2> -O <outPre> -p <posList> -m <abdance for major GT>\n
# ==============================================================================
#
#   Used to compare the genotypes of two strains, and get the site that different between two strains (only for haploid);
#    Authors: Fang Zhiwei (fangzhiwei1126\@126.com)
#
# ==============================================================================

    compareHap: was the genotyope derived from;
    outPre    : Prefix of output files;
    posList   : the pos of base on reference, as some region within amplicon were ignored when call SNP;
                so the fiftieth base must not be the fiftieth, but could be 56th or 58th or others; due 
                how many base have been ignored;
    majorAbd  : the least abundance of major GT that could be used to estimate the diff between two strains;

    
    Revised: All SNP were ignored if two or more SNPs site were adjacent to each other; \@20171130;

    compareHap1: generate as follows:
                perl ../diffStat.pl shan 3 singleSam/compar_all1613Sam_D62A.XLS | sort -k3gr| grep -v SSR   >D62A.abd.stat_compareHap.xls

    Revised: the GT of A sample's amplicon site should be delete when it found no difference with B sample's GTs in the same amplicon site,
             otherwise it could be used to compare with next GT of B sample's, in the case that no difference found again, A sample will be 
             in this amplicon site. this is wrong!!!. \@20180410;

    Revised: we compre the GT one position by another, and  all indel and uncovered region were set to be no difference for the Parents and F1; \@20171031;

    Revised: we add otPre to named the output files; F1 label was used to indicate the F1 generation; 
             The order for input files MUST be parent1 parent2 F1;   revised at 20170725;

    Attention: the amplicon site will be compared only if the marjor genotype has at least 20 reads supported, 
                a stutter genotype will be considered only if it has more 10% of the reads of marjor genotype.

                PLEASE CHECK IF THIS IS WHAT U WANT; AND this is just used for haploid,nor for multiploid.
    \n";


#die("Some difference were caused by the uncovered region which will be marked as 0, but other sample were Ref or others,this should be handle.\n");
my($infile1,$infile2,$otPre,$posList,$cutoff4abd)=();
my@sam=();

GetOptions(
    'i1=s' => \$infile1,
    'i2=s' => \$infile2,
    'O=s'  => \$otPre,
    'm=i'  => \$cutoff4abd,
);

die($usage) unless defined $infile1;
die($usage) unless defined $infile2;

$cutoff4abd=20 unless (defined $cutoff4abd);


my(%h,%type)=();
my%check=();
my@files=($infile1,$infile2);

foreach my$f(@files){
    my$pre=(split /\./,(split /\//,$f)[-1])[0];  
    die("Same sample names were found, please check it\n") if exists $check{$pre};
    $check{$pre}="";
    push @sam,$pre;
    
    open IN , "$f " or die  "cannot open $f "; 
    while(<IN>){
        chomp;  
        my($id,$aid,$numOfStutter,$abd,$major,$otherAbd,$otherType) =(split /\s+/,$_);
        next if $abd <$cutoff4abd;
        $major=~s/(\D)\,S:/$1/g;
        $type{$aid}{$major}="";

        $h{$aid}{$pre}{$major}=$abd;
        $type{$aid}{$major}=$abd;
        next if $numOfStutter ==0;
        my@stutterAbd=(split /\,/,$otherAbd); 
        my@stutter=(split /;/,$otherType);

        for(my$i=0;$i<=$#stutter;$i++){
#            next if ($stutterAbd[$i]/$abd)<0.1;
            $stutter[$i]=~s/(\D)\,S:/$1/g;

            $type{$aid}{$stutter[$i]}=""; 
            $h{$aid}{$pre}{$stutter[$i]}=$stutterAbd[$i];
        }
    }
    
    close IN; 
}

print STDERR "Finished read in the input files\n";

my$title=join"\t",@sam; 
my$ot=  "#aid\t$title\t$sam[0]_GT\t$sam[1]\_GT\n#we record all genotypes occured in all strains\
#for each amplicon site, and output the reads number support each genotype in each strains. The \
#reads number are in the same order of genotypes for each strains.\n";

open OT , "> compre_.$otPre.xls "  or die  "cannot write compre_parent2F1.NMR1.xls $!";
print OT "$ot\n";

my%common=();
my%commonGT=();
my%pos=();
my%delete=();
my($commonSiteNum,$diffSiteNum)=(0,0);
my$exist=0;
my$isDiff=0;
foreach my$aid (keys %type){
    my$flag=0;
    foreach my$sam(@sam){ 
        $flag =1 unless exists $h{$aid}{$sam};
        #$flag =1 if(scalar keys %{$h{$aid}{$sam}})>2; mutiple GT is ok;
    }    
    next if $flag ==1;
    
    $commonSiteNum++;
    my@all_gt=(sort {$h{$aid}{$sam[-1]}{$b}<=>$h{$aid}{$sam[-1]}{$a}} keys %{$h{$aid}{$sam[-1]}}); #output the abdunce from high to low based on F1;
    foreach my$gt(keys %{$type{$aid}}){
        push @all_gt,$gt unless exists $h{$aid}{$sam[-1]}{$gt};
    }

    $isDiff=0;
    foreach my$gt(@all_gt){
        next if $gt=~/^-/;#uncertain GT;
        next if exists $delete{$gt};
        $delete{$gt}="";
        
        $exist=0;
        foreach my$sam(@sam){
            if(exists $h{$aid}{$sam}{$gt}){
                $common{$aid}{$sam} .="$h{$aid}{$sam}{$gt},";
                $commonGT{$aid}{$sam} .="$gt;";
                delete $h{$aid}{$sam}{$gt};# add at 20180410;
                $exist++;
                next;
            }

            my$snpNum=0;
            $flag=0;
            foreach my$gt4each(keys %{$h{$aid}{$sam}}){
                %pos=();
                $snpNum=&compareGT($gt,$gt4each);
                next if $snpNum !=0;
                $common{$aid}{$sam} .="$h{$aid}{$sam}{$gt4each},";
                $commonGT{$aid}{$sam} .="$gt4each;";
                $delete{$gt4each}="";
                delete $h{$aid}{$sam}{$gt4each}; # add at 20180410;
                $flag=1;
                $exist++;
                last;
            }
            $common{$aid}{$sam} .="0,"   if $flag==0;
            $commonGT{$aid}{$sam} .="-;" if $flag==0;
        }
        $isDiff ++ if $exist ne @sam;
    }

    $ot ="$common{$aid}{$sam[0]}\t$common{$aid}{$sam[1]}";
    $ot.="\t$commonGT{$aid}{$sam[0]}\t$commonGT{$aid}{$sam[1]}";

    if($isDiff ==0){
        #print OT "$sam[0]\t$sam[1]\t$aid\t0\t$ot\n";
    }else{
        print OT  "$sam[0]\t$sam[1]\t$aid\t1\t$ot\n";
        $diffSiteNum ++;
    }

#    my@p1abd=(split /\,/,$common{$aid}{$sam[0]}); 
#    my@p2abd=(split /\,/,$common{$aid}{$sam[1]});
#
#    if(@p1abd ==1){
#        #print OT "$sam[0]\t$sam[1]\t$aid\t0\t$ot\n";
#        next;
#    }
#    
#    if(@p1abd ==2){
#        if(($p1abd[0]*$p2abd[0])==0){
#            print OT  "$sam[0]\t$sam[1]\t$aid\t1\t$ot\n";
#            $diffSiteNum ++;
#        }elsif(($p1abd[1]*$p2abd[1])==0){
#            print OT "$sam[0]\t$sam[1]\t$aid\t1\t$ot\n";
#            $diffSiteNum++;
#        }else{
#            #print OT "$sam[0]\t$sam[1]\t$aid\t0\t$ot\n";
#        }
#        next;
#    }
#
#    if(@p1abd >2){
#        print OT "$sam[0]\t$sam[1]\t$aid\t1\t$ot\n";
#        $diffSiteNum ++;
#    }
}
close OT;


print STDERR "$sam[0]\t$sam[1]\tCommonSiteNum\tdiffSiteNum\tdiffRatio\n";

if($commonSiteNum==0){
    print STDERR"$sam[0]\t$sam[1]\tnoCommonSite\t-\t-\n";
}else{
    print STDERR "$sam[0]\t$sam[1]\t$commonSiteNum\t$diffSiteNum\t";
    printf STDERR "%.3f\n",($diffSiteNum/$commonSiteNum);
}



#################################   subroutin   ##########################
sub compareGT{
    my($major_,$seGT_)=@_;
    $major_=~s/^-//;
    $seGT_ =~s/^-//;
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
        }
        
        for(my$i=1;$i<=1;$i++){
            my$tmp=$_+$i;
            if(exists $posList{$tmp} ){
                if($i<=1){
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
                if($i<=1){
                    $indel ++ if $posList{$tmp} eq "snp";
                    $indel ++ if $posList{$tmp} eq "indel";
                }else{
                    $indel ++ if $posList{$tmp} eq "indel";
                }
                }
                }
         
        $pos{$_}{$order}="snp"  if $indel ==0;
        $pos{$_}{$order}="ignored" if $indel >0;
        #delete $posList{$_} if $indel>0;
    }
}
