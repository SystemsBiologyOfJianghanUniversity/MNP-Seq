#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my$usage="Usage:  perl $0 -M <Minmum depth> -M2 <depth_cutoff> -S <uncoveredLen> -H1 <Hybrid type1> -H2 <Hybrid type1> -I1 <Inbred type1> -I2 <Inbred type2> -P <genotype Pools> -o <otPre> -r <ratio of hybrid site> -pl <parameter list> -3thGT <thirdRatio> -D <redo=1;not=0> -snp <yes=1,no=1>\n

# ==============================================================================
#
#    Used to identify if a amplicon site is hybrid or inbred based on SSR/SNP\n
#    Authors: Fang Zhiwei (fangzhiwei1126\@126.com)
#
# ==============================================================================
    Minmum depth: the mimum reads number that covered the amplicon site in each sample; default=20;
    depth_cutoff: Hybrid type2 and Inbred type2 will be used when the maojor abd  is more than this;
    uncoveredLen: the site will be determined to be NA if the ratio of third genotype to majorGT is more than uncoveredLen; default =5;
    Hybrid type1: Minmum depth> cutoff &&  stuttuer Distance>cutoff && coverDepth4major:coverDepth4minor > cutoff; default=0.7;
    Inbred type1: Minmum depth> cutoff &&  stuttuer Distance>cutoff && coverDepth4major:coverDepth4minor < cutoff; default=0.5;
    
    Hybrid type2: Minmum depth> cutoff &&  stuttuer Distance<cutoff && coverDepth4major:coverDepth4minor > cutoff && repeatNumber4MajorGT> repeatNumber4MinorGT; default=0.7;
    Inbred type2: Minmum depth> cutoff &&  stuttuer Distance<cutoff && coverDepth4major:coverDepth4minor < cutoff && repeatNumber4MajorGT> repeatNumber4MinorGT; default=0.6;
    
    parameterList: used for hybrid samples that more than 10% site were found to be hybrid in samples when genotyping  as a inred lines; it is in format: MinDept,uncoveredLen,hybrid1,hybrid2,hy                   brid3,inbred1,inbred2,inbred3;
    thirdRatio: the site will be determined to be NA if the ratio of third genotype to majorGT is more than thirdRatio;
    redo:       if it is required to redo genotyping when the ratio of hybrid sites more than 10%; 1=yes,0=no; default=1;
    snp: if all SNPs were ignored if two or more SNPs site were adjacent to each other; default=1: ignored;

    Revised5: redo was set to be 0=no redo genotyping; line 328 and 335 should be (>I1 and <H1) or(>I2 or <H2);\@20180424;


    Revised4: All SNP were ignored if two or more SNPs site were adjacent to each other; \@20171130;
   
    Revised3:  we changed the method used to count the variation: a snp will be ignored if there is a indel beside the SNP; if it is a SNP on one GT but indel on the other one when two GT compared, both variation ignored;
              Actually, an candidate genotype will be discarded ifs no SNP found when compared with major genotyoe or any other candidate genotype that with higher abundance than it; this should be included when a  at \@20171113;

    Revised2:
    ratio of hybrid: if one samples have more than xx% site to be hybrid, then this sample will be reanalized based on parameter set setted for hybride samples. \@20170829;

    The genotype will be determined to be NA  if a amplicon site has less than <Minmum depth> reads covered, or the ratio of coverDepth4minor to coverDepth4major GT is between cutof for inbred and cutoff for hybrid.
   
   Revised: we introduced the genotype pool derived from all samples sequenced in 34 bath; this is used to discard stutter genotype; \@20170825

    Attention: two or more second most abndant stutter maybe exists, both of which could have less repeat number compared with mayjor GT, or one of them have more repeat number compared major GT. In this program, we first sort all genotype based on abundance of each genotype, and selected the second most abundant genotype from all that have equal number of reads supported. 20170814. 

    Attention: we donot count the effect on stutture due to the GC content, both number and length of repeat unit, the SSR length and so on. Here, only the distance, direction of stutter and depth of reads covered were considered.
    \n";

my($MinDept,$depth_cutoff,$uncoveredLen,$hybrid1,$hybrid2,$hybrid3,$inbred1,$inbred2,$inbred3,$help,$h,$genotypePool,@allFiles,$otPre,$ratio,$parameterList,$thirdRatio,$redo,$ignore)=();
@allFiles=@ARGV;
my%pos=();

GetOptions(
    'M=i'  => \$MinDept,
    'M2=i' => \$depth_cutoff,
    'S=i'  => \$uncoveredLen,
    'H1=f' => \$hybrid1,
    'H2=f' => \$hybrid2,
    'I1=f' => \$inbred1,
    'I2=f' => \$inbred2,
    'pl=s' => \$parameterList,
    'P=s'  => \$genotypePool,
    'o=s'  => \$otPre,
    'r=f'  => \$ratio,
    '3thGT=f' => \$thirdRatio,
    'D=i'    => \$redo,
    'snp=i' => \$ignore,
    'help' => \$help,
);


if(defined $MinDept){
    print "\$MinDept is ok: $MinDept\n";
    shift @allFiles;
    shift @allFiles;
}else{$MinDept=20;}

if(defined $depth_cutoff){
    shift @allFiles;
    shift @allFiles;
}else{$depth_cutoff=50;}

if( defined $uncoveredLen){
    print "\$stutterDistance is ok: $uncoveredLen\n";
    shift @allFiles;
    shift @allFiles;
}else{$uncoveredLen=5;}

if(defined $hybrid1){
    print "\$hybrid1 is ok: $hybrid1\n";
    shift @allFiles;
    shift @allFiles;
}else{$hybrid1=0.3;}

if(defined $hybrid2){
    print "\$hybrid2 is ok: $hybrid2\n";
    shift @allFiles;
    shift @allFiles;
}else{$hybrid2=0.2;}

if(defined $inbred1){
    print "\$inbred1 is ok: $inbred1\n";
    shift @allFiles;
    shift @allFiles;
}else{$inbred1=0.05;}

if(defined $inbred2){
    print "\$inbred2 is ok: $inbred2\n";
    shift @allFiles;
    shift @allFiles;
}else{$inbred2=0.05;}

if( defined $thirdRatio){
    shift @allFiles;
    shift @allFiles;
}else{
    $thirdRatio=0.2;
}

if( defined $genotypePool){
    shift @allFiles;
    shift @allFiles;
}

if( defined $otPre){
    shift @allFiles;
    shift @allFiles;
} 

if(defined $ratio){
    shift @allFiles;
    shift @allFiles;
}else{
    $ratio=0.1;
}

if(defined $parameterList){
    shift @allFiles;
    shift @allFiles;
}

if(defined $redo){
    print STDERR "\$redo is ok: $redo\n";
    shift @allFiles;
    shift @allFiles;
}else{
    $redo=0;
    print STDERR "\$redo is default: $redo\n";
}

if (defined $ignore){
    shift @allFiles;
    shift @allFiles;
}else{
    $ignore=1;
}
die($usage) if(defined $help);
die($usage) if(@allFiles==0);
die($usage) unless defined $otPre;
print STDERR "No genotype pool supplied\n" unless defined $genotypePool;

#print STDERR "$infile\n@allFiles\n";

my%h=();
my%gtPool=();
my$cycFlag=0;
my%reAnalyz=();# record sample that have too much hybrid sites;

open OT , "> $otPre.xls " or die  "cannot write $otPre.xls $!";
open OT2, "> $otPre.err " or die  "cannot write $otPre.err $!";

my@type=qw/inbred0 inbred1 inbred2 inbred3 hybrid1 hybrid2 hybrid3 NA0 NA1 NA2 NA3/;
my$title=join "\t",@type; 
print OT "\t$title\n";

&getPools() if defined $genotypePool;
&getStat() if @allFiles >0;
&otStat();

exit (1) if $redo==0;
exit (1) if (scalar keys %reAnalyz) ==0;

print STDERR "ok\n";

$cycFlag=1;
%h=();

#my$date=`date`;
my$date=int(rand(1000000));
print`cp $otPre.err $otPre.err.$date`;
open OT2, " >$otPre.err " or die "cannot open $otPre.err $!";
open IN, " $otPre.err.$date " or die "cannot open $otPre.err.$date $!";
while(<IN>){
    chomp;
    if(/processing/i){print OT2 "$_\n";next;}
    my$sampleID=(split /\s+/,$_)[0];
    print OT2  "$_\n" unless exists $reAnalyz{$sampleID};
}
close IN;
print `rm -rf $otPre.err.$date`;

if(defined $parameterList){
    ($MinDept,$uncoveredLen,$hybrid1,$hybrid2,$hybrid3,$inbred1,$inbred2,$inbred3)=(split /\,/,$parameterList);
}else{
    ($MinDept,$uncoveredLen,$hybrid1,$hybrid2,$hybrid3,$inbred1,$inbred2,$inbred3)=(20,20,0.7,0.4,0.3,0.4,0.2,0.2);
}


&getStat()   if @allFiles >0;
&otStat();


############################################### v sub ############################################
sub otStat {
    my($hybridNum,$total,$otLine)=(0,0,"");
    foreach my$sam(keys %h){
        $otLine= "$sam";
        ($hybridNum,$total)=(0,0);
        foreach (@type){
            unless (exists $h{$sam}{$_}){
                $otLine .= "\t0";
                next;
            }
            $otLine .= "\t$h{$sam}{$_}";
            $total += $h{$sam}{$_} unless $_=~/NA/;;
            $hybridNum +=$h{$sam}{$_} if $_=~/hybrid/;
        }
        if($total ==0){
            print "Sample $sam has noamplicon site detected\n";
            next;
        }
        if($cycFlag==1){
            print OT "$otLine\n";
            next;
        }
        print OT "$otLine\n" if ($hybridNum/$total)<= $ratio;
        $reAnalyz{$sam}=""   if ($hybridNum/$total)  >$ratio;
    }
#    close OT;
    close OT2;
}


 sub getStat{
     foreach my$file (@allFiles){

        my%tmp=();
        my(%cov,%record)=(); # these two were kept for other use;
        my$pre=$file;
        $pre=(split /\./,(split /\//,$pre)[-1])[0];
        
        if($cycFlag ==1){ next unless exists $reAnalyz{$pre};} ## used for second cycle; 
        
        print OT2 "Processing $file \n";

        open IN , " $file " or  die "cannt open $file $!";
        open OTTEST, "> test_$pre.xls" or die  "cannot open test_$pre.xls $!";
        while(<IN>){
            chomp;
            my($pre_,$aid,$gtNum,$majorAbd,$majorGT,$otherAbd,$otherGT) =(split /\s+/,$_);
            
            my$majorGT_=$majorGT;
            my$partialFlag=0;
            while(1){
                last unless $majorGT_=~/\d+,\D*0+/;
                $majorGT_=~s/\d+,\D*(0+)//;
                if(length($1)>=$uncoveredLen){
                    print OT2 "$pre\t$aid\tNA0\tpartial\t-\t-\t-\t-\n";
                    $h{$pre}{"NA0"}++;
                    $partialFlag=1;
                    last;
                }
            }
            next if $partialFlag ==1;

            if($majorAbd<$MinDept){
                print  OT2 "$pre\t$aid\tNA0\tfewReads\t-\t-\t-\t-\n";
                $h{$pre}{"NA0"}++;
                next;
            }


            if($gtNum ==0){
                print OT2 "$pre\t$aid\tinbred0\t-\t$majorAbd\t0\t$majorGT\t-\n";
                print OTTEST "$pre\t$aid\t0\t$majorAbd\t$majorGT\n";
                $h{$pre}{"inbred0"}++ ;
                next;
            }
###### it is not necessary to check all the minor GT  as they have at least one SNP compared with majorGT;
#            my(@candiateAbd_,@candidateGT_)=("","");
#            my@candiateAbd=(split/\,/,$otherAbd);
#            my@candiateGT=(split /\;/,$otherGT);
#            for(my$cyc=0;$cyc<=$#candiateGT;$cyc++){
#                $partialFlag=0;
#                my$candGT=$candiateGT[$cyc];
#                while(1){
#                    last unless $candGT=~/\d+,\D*0+/;
#                    $candGT=~s/\d+,\D*(0+)//;
#                    if(length($1)>=$uncoveredLen){
#                        $partialFlag=1;
#                       last;
#                   }
#                }
#                next if $partialFlag==1;
#                push @candiateAbd_,$candiateAbd[$cyc];
#                push @candidateGT_,$candiateGT[$cyc];
#            }
###########################################################################################################


            my@stuGTs=(split/;/,$otherGT);
            my@stuAbd=(split /\,/,$otherAbd);

            if(($stuAbd[0]/$majorAbd)<$inbred1){
                print OT2 "$pre\t$aid\tinbred1\t-\t$majorAbd\t0\t$majorGT\t-\n";
                print OTTEST "$pre\t$aid\t0\t$majorAbd\t$majorGT\n";
                next;
            }

            my($certain,$unCertain)=("","");
            my($cAbd,$ucAbd)=("","");
            my($cNum,$ucNum)=(0,0);

            for(my$cc=0;$cc<=$#stuGTs;$cc++){
                next if ($stuAbd[$cc]/$majorAbd)<$inbred1;
                
                if( defined $genotypePool){
                    my$SNPnum=&checkDB($stuAbd[$cc],$aid);
                    next if $SNPnum!=0;
                }

                if(($stuAbd[$cc]/$majorAbd)>$hybrid1){
                    $cNum ++;
                    $cAbd.="$stuAbd[$cc],";
                    $certain .="$stuGTs[$cc];";
                }else{
                    $cNum ++;
                    $cAbd .="$stuAbd[$cc],";
                    $certain .="-$stuGTs[$cc];";
                    #$ucNum ++;
                    #$ucAbd .="$stuAbd[$cc],";
                    #$unCertain.="-$stuGTs[$cc];";
                }
            }
            if($cNum==0){
                 print OT2 "$pre\t$aid\tinbred1\t-\t$majorAbd\t0\t$majorGT\t-\n";
                 print OTTEST "$pre\t$aid\t0\t$majorAbd\t$majorGT\n";
                 next;
            }
            print OT2 "$pre\t$aid\thybrid\t-\t$majorAbd\t$cAbd\t$majorGT\t$certain\n";
            print OTTEST "$pre\t$aid\t$cNum\t$majorAbd\t$majorGT\t$cAbd\t$certain\n";

#            if(($ucNum+$cNum)==0){
#                print OT2 "$pre\t$aid\tinbred1\t-\t$majorAbd\t$stuAbd[0]\t$majorGT\t$stuGTs[0]\n";
#                next;
#            }
#
#           if(($ucNum*$cNum)!=0){
#                print OT2 "$pre\t$aid\thybrid\t-\t$majorAbd\t$cAbd\t$majorGT\t$certain\t$ucAbd\t$unCertain\n";
#                next;
#            }
#
#            if($ucNum!=0){
#                print OT2 "$pre\t$aid\thybrid\t-\t$majorAbd\t$cAbd\t$majorGT\t-\t$ucAbd\t$unCertain\n";
#                next;
#            }
#
#            if($cNum!=0){
#                print OT2 "$pre\t$aid\thybrid\t-\t$majorAbd\t$cAbd\t$majorGT\t$certain\t$ucAbd\t-\n";
#            }
        }
        close IN;
    }
}


sub getPools{
    open IN , " $genotypePool " or die  "cannot open $genotypePool $!";
    while(<IN>){
        chomp; 
        next unless /AMP/;
        my($aid,$gt)=(split /\s+/,$_)[0,3];
        $gtPool{$aid}{$gt}="";
    }
    close IN;
}

sub checkDB {
    my($queryGT,$qAid)=@_;
    my$qNum=-1;
    foreach (keys %{$gtPool{$qAid}}){
        $qNum=&compareGT($queryGT,$_);
        last if $qNum==0;
    }
    return($qNum);
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
    foreach ( sort {$a<=>$b} keys %posList){
        $indel=0;
        if( $posList{$_} eq "indel"){
            $pos{$_}{$order}="ignored";
            next;
        }
        
        for(my$i=1;$i<=1;$i++){
            my$tmp=$_+$i;
            if(exists $posList{$tmp} ){
                if($i<=1){
                    $indel ++ if (($posList{$tmp} eq "snp") and ($ignore==1)) ;
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
                    $indel ++ if (($posList{$tmp} eq "snp") and ($ignore==1));
                    $indel ++ if $posList{$tmp} eq "indel";
                }else{
                    $indel ++ if $posList{$tmp} eq "indel";
                }
            }
        }
        
        $pos{$_}{$order}="snp"  if $indel ==0;
        $pos{$_}{$order}="ignored" if $indel>0;
    }
}
