#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use threads;
use FindBin '$Bin';

my$usage="perl $0 -split <N> -pos <amplicon776.list> -fa <all776amplicon.fa> -p <776primer.fa> -o <otPre> -d <targetList> -t <gtPool> -db <dbFlag> <sam1.input> <sam2.input> <sam3.input> ....
#======================================================================

#  Used to genotype of amplicon-sequencing.
#  Auther: fangzhiwei1126\@126.com
#======================================================================

    split:        how many jobs U want to submit, max=30;
    targetList:   which samples should be compared with, it is in this formate:
                      sample1.abd.stat.xls.check
                      sample2.abd.stat.xls.check
                      ...
    ampliconList: chr   aid    start   end
    776primer.fa: it is in this formate:
                    >aid1_F
                    atcgatcgatcg
                    >aid1_R
                    atcgatcg
                    >adi2_F
                    atcgatcg
                    >aid2_R
                    atcgatcg
    targetList: sample that you want to compare with for all sample listed as sam1.input sam2.input et al; it is in formate as follows:
                /dir/samID1.abd.stat.xls
                /dir/samID2.abd.stat.xls
                .......
    dbFlag    : the directory of database files; If you want to compare the MNP fingerprints of test samples with other known varieties
                listed in a database, you should specified the directory of the database;
    gtPool    : a collection of true alleles for each amplicon that have been detected and validated in all the standard samples.

";

die("$usage") if @ARGV ==0;
my($split,$ampliconList,$ampliconFa,$ampliconPrimerFa,$otPre,$targetList,$dbFlag,$gtPool)=();

GetOptions(
    'split=i' => \$split,
    'pos=s'   => \$ampliconList,
    'fa=s'    => \$ampliconFa,
    'p=s'     => \$ampliconPrimerFa,
    'o=s'     => \$otPre,
    'd=s'     => \$targetList,
    'db=s'    => \$dbFlag,
    't=s'     => \$gtPool,
);

if(defined $split){
#    shift @ARGV;
#    shift @ARGV;
    $split =30 if $split>=30;
}else{
    $split=1;
}

if(defined $ampliconList){
#    shift @ARGV;
#    shift @ARGV;
}else{ die"No amplicon list found!\n";}

if(defined $ampliconFa){
#    shift @ARGV;
#    shift @ARGV;
}else{ die"No amplicon Fasta file found.\n";}

if(defined $ampliconPrimerFa){
#    shift @ARGV;
#    shift @ARGV;
}else{ die"No amplicon Primer found.\n";}

if(defined $otPre){
#    shift @ARGV;
#    shift @ARGV;
}else{
    my@second = localtime();
    $otPre=join "",@second;
}

my$i=1;
my$j=0;
my@samList=();
my@threadNum=();


open OT, ">x$i" or die "cannot write x$i $!";
my$maxJob=@ARGV/$split;
$maxJob=int($maxJob);

my%h=();
print `mkdir  singleSam`;
foreach(@ARGV){
    my$pre=(split /\./,(split /\//,$_)[-1])[0];
    push @samList,$pre;
    my$inf=$_;
    if($inf=~/gz$/){
        open IN, "gzip -cd $inf | " or  die ("cannot open $inf $!");
    }else{
        open IN, " $inf " or  die ("cannot open $inf $!");
    }
    open OT, "> singleSam/compar_all2336Sam_$pre.XLS " or die "cannot write commpar_all2336Sam_$pre.XLS $!";
    %h=();
    while(<IN>){
        chomp;
        next if /^#/;
        next if /^\s*$/;
        my($aid,$gt)=(split /\s+/,$_)[1,-1];
        next if $gt=~/N/;
        $h{$aid}{$gt} ++;
    }
    close IN;
    foreach my$aid(keys %h){
        foreach (sort {$h{$aid}{$b}<=>$h{$aid}{$a}} keys %{$h{$aid}}){
            next if $h{$aid}{$_} <2;
            print OT "$pre\t$aid\t$h{$aid}{$_}\t$_\n";
        }
    }
    close OT;
}

@threadNum=();
$i=1;$j=0;
open OT,">abd$i.sh" or die "cannot write abd$i.sh $!";
foreach my$pre(@samList){
    $j++;
    print OT  "perl $Bin/diffStat.pl shan 0 singleSam/compar_all2336Sam_$pre.XLS \| sort -k3gr\| grep -v SSR   >$pre.abd.stat.xls\n";#\@20190219;
    if($j>=$maxJob){
        $j=0;
        close OT; 
        $threadNum[$i]= threads->new(\&submit,"abd$i.sh");
        $i++;
        open OT,">abd$i.sh" or die "cannot write abd$i.sh $!";
    }
}
close OT;
$threadNum[$i] = threads->new(\&submit,"abd$i.sh");

&finish();
 
my@samList2=();  
if(defined $targetList){ 
    open IN , " $targetList" or die "cannot open $targetList $!";;
    while(<IN>){
        chomp; 
        next if /^#/;
        die "File named $_ cannot be found!\n" unless -f $_;
        my$tmp_fileName=(split /\./, (split /\//,$_)[-1])[0];#20180621;
        if (-f  "$tmp_fileName.abd.stat.xls"){
            print STDERR "File $tmp_fileName are crossed and renamed as $tmp_fileName\_target, please check!";
            print `ln -s  $_ "$tmp_fileName\_target.abd.stat.xls"`;
            push @samList2, "$tmp_fileName\_target";
            next;
        }

        print `ln -s  $_`;
        $_=(split /\./,(split /\//,$_)[-1])[0];
        push @samList2, $_;
    }
    close IN;
}

if(defined $dbFlag){
    foreach my$tmp_file(`ls $dbFlag/*.abd.stat.xls`){ #\@20190219;
        chomp($tmp_file);
        my$tmp_fileName=(split /\./, (split /\//,$tmp_file)[-1])[0];
        if (-f  "$tmp_fileName.abd.stat.xls"){
            print STDERR "File $tmp_fileName are crossed and renamed as $tmp_fileName\_db, please check!";
            print `ln -s  $tmp_file "$tmp_fileName\_db.abd.stat.xls"`;
            push @samList2,"$tmp_fileName\_db";
            next;
        }
        push @samList2,$tmp_fileName;
        print `ln -s $tmp_file`;

        print "$tmp_fileName\t$tmp_file\n";
    }
}


my$list=join".abd.stat.xls ",(@samList,@samList2);
print `perl $Bin/checkSNP.pl $list.abd.stat.xls`; 
$list=join".abd.stat.xls.check ",(@samList,@samList2);

die ("No gtPool named $gtPool found !\n") unless -f $gtPool;
print  `perl $Bin/isHybrid4MNP_peng.pl -o ishybrid_test -P $gtPool $list.abd.stat.xls.check`;
print  `perl $Bin/checkSNP2_afterGTping.pl $gtPool 4 ishybrid_test.err`;




%h=();
open IN , "  ishybrid_test.err.f " or die "cannot open ishybrid_test.err $!";
while(<IN>){
    chomp; 
    my($id,$aid,$type,$nul,$abd,$abd2,$gt,$gt2)=(split /\s+/,$_);
    $id=~s/compar_all2336Sam_//;
    if($type=~/inbred/){
        $h{$id}{$aid} .="0\t$abd\t$gt";
    }elsif($type=~/hybrid/){
        $h{$id}{$aid} .="1\t$abd\t$gt\t$abd2\t$gt2";
    }
}
close IN;

foreach my$id (keys %h){ 
    open OT, ">test_$id.xls" or die "canot write test_$id.xls $!"; 
    foreach (keys %{$h{$id}}){
        print OT "compar_all1613Sam_$id\t$_\t$h{$id}{$_}\n";
    }
    close OT;
} 

my$num=100; 
my$pair=($#samList*($#samList-1))/2;
$pair +=(scalar @samList)*(scalar @samList2);
if($pair>100){
    $num=int($pair/$split);
    $num ++;
}
@threadNum=();
my$m=0;my$n=1;

open OT,"> compare$n.sh " or die "cannot write compare$n.sh !";
print OT "\`rm -rf comparePart$n.list\`\n";
for(my$i=0;$i<=$#samList;$i++){ 
    for(my$j=$i+1;$j<=$#samList;$j++){
      
        print OT "perl $Bin/cauculteDustance4MNP.pl -i1  test_$samList[$i].xls -i2 test_$samList[$j].xls -O $samList[$i].vs.$samList[$j]  -m 20; \`cat compre_.$samList[$i].vs.$samList[$j].xls >> comparePart$n.list \`; \`rm -rf compre_.$samList[$i].vs.$samList[$j].xls\`;\n";   #\@20190219;
        $m++;
        next if $m<$num;
        close OT; 
        $threadNum[$n]=threads->new(\&submit,"compare$n.sh");
        $n++;
        open OT,"> compare$n.sh " or die "cannot write compare$n.sh !";
        print OT "\`rm -rf comparePart$n.list\`\n";
        $m=0;
    }
    next if @samList==0;
    for(my$j=0;$j<=$#samList2;$j++){

        print OT "perl  $Bin/cauculteDustance4MNP.pl -i1  test_$samList[$i].xls -i2 test_$samList2[$j].xls -O $samList[$i].vs.$samList2[$j]  -m 20; \`cat compre_.$samList[$i].vs.$samList2[$j].xls >> comparePart$n.list \`; \`rm -rf compre_.$samList[$i].vs.$samList2[$j].xls\`;\n";   #\@20190219;
        $m++;
        next if $m<$num;
        close OT;
        $threadNum[$n]=threads->new(\&submit,"compare$n.sh");
        $n++;
        open OT,"> compare$n.sh " or die "cannot write compare$n.sh !";
        print OT "\`rm -rf comparePart$n.list\`\n";
        $m=0;
    }
}
close OT;
$threadNum[$n]=threads->new(\&submit,"compare$n.sh");


&finish();

system("echo \"sample1\tsample2\tCommonSiteNum\tdiffSiteNum\tdiffRatio\" >$otPre\_compre.dist.list");
print `cat compare*.sh.err|grep -v Finished|grep -v CommonSiteNum >>$otPre\_compre.dist.list`;
system("echo \"sample1\tsample2\tampliconID\tsame=0/not=1\tabundance4sample1\tabundance4sample2\tgenotype4sample1\tgenotype4sample2\" > $otPre\_compre.amplicon.list");
print `cat comparePart*.list |grep -v "#" >$otPre\_compre.amplicon.list`;
print `rm -rf compare*.sh.err comparePart*.list`;

foreach my$tmpID(@samList2){
    print "rm -rf $tmpID.abd.stat.xls $tmpID.abd.stat.xls.check test_$tmpID.xls\n";
}
    

sub finish{
    while(scalar(threads->list())!=0){
        foreach my$mythread(threads->list(threads::all)){
            if($mythread ->is_joinable()){
                $mythread ->join();
            }else{
                sleep(3);
            }
        }
    }
}

sub submit{
    my$shell=shift @_;
    print `bash $shell 1>$shell.log 2>$shell.err `;

}
