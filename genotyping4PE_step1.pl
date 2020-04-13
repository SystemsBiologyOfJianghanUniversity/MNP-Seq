#!/usr/bin/perl -w 
use strict;
use Getopt::Long;
use FindBin '$Bin';
my$usage="perl $0 -1 <reads1.fq > -2 <reads2.fq> -r <ref> -o <otPre> -m <overlapLen> -e <errorRate> -a <ampliconList> -f <amplicon.fa> -p <primer.fa> -b <sam/bam>\n

#==============================================================================
#
#       Used to process amplicon sequencing reads and generate variant_stat.xls file.
#
#       Authors: Fang Zhiwei (fangzhiwei1126\@126.com)
#==============================================================================


    ref:          reference genome;
    otPre:        prefix for output files;
    overlapLen:   the minmum length of overlap region between a pair of reads;
    errorRate:    the error rate of overlap region between a pair of reads.
    ampliconList: tareget region list;
    amplicon.fa:  the sequence of amplicon, derived from reference genme sequece;
    primter.fa:   primer sequece used for multiplex PCR amplification;
    sam/bam:      a bam file that contained the sequence reads, and will be used as input.
                  Atttention you should not specified a bam file and fastq files simultaneously.
     
    \n";

my($fq1,$fq2,$ref,$otPre,$overlapLen,$errorRate,$ampliconList,$ampliconFa,$ampliconP,$bam)=();
GetOptions(
    '1=s' => \$fq1,
    '2=s' => \$fq2,
    'r=s' => \$ref,
    'o=s' => \$otPre,
    'm=i' => \$overlapLen,
    'e=f' => \$errorRate,
    'a=s' => \$ampliconList,
    'f=s' => \$ampliconFa,
    'p=s' => \$ampliconP,
    'b=s' => \$bam,
);
unless(defined $bam){
die("$usage") unless defined $fq1;
die("$usage") unless defined $fq2;
}
die("$usage") unless defined $ref;
die("$usage") unless defined $otPre;


$overlapLen=20 unless defined $overlapLen;
$errorRate=0.01 unless defined $errorRate;

$ampliconList="/result/pub/software/software/bin/SSR/pipeline/amplicon776.list" unless defined $ampliconList;
$ampliconFa="/result/pub/software/software/bin/SSR/pipeline/all776amplicon.fa" unless defined $ampliconFa;
$ampliconP= "/result/pub/software/software/bin/SSR/pipeline/776primer.fa" unless defined $ampliconP;


my$overlap="$Bin/overlap_pair_trim.new";
my$stat_var="$Bin/statVariant.pl";
my$combin="$Bin/combinVariant.pl";

my$bwt_index=$ref;
if($bwt_index=~/.fna$/){$bwt_index=~s/.fna//;}elsif($bwt_index=~/.fa$/){
$bwt_index=~s/.fa$//;}elsif($bwt_index=~/.fasta$/){$bwt_index=~s/.fasta$//;}
print "#######$bwt_index\n";
unless (-f "$bwt_index.4.bt2"){ die "No BOWTIE2 index found! \n";}
unless (-f "$bwt_index.sa"   ){ die "No BWA index found! \n";}

if(defined $bam){
    my%shan=();
    open IN, " samtools view $bam |" or die "cannot open $bam $!";
    open OT,"> $otPre.fa " or die "cannot open $otPre.fa $!";
    while(<IN>){
        chomp;
        my($id,$seq)=(split /\s+/,$_)[0,9];
        next if exists $shan{$id};
        $shan{$id}="";
        print OT ">$id\n$seq\n";
    }
    close IN;
    %shan=();
    print`bowtie2 -f --un $otPre.unmap.fa -x $bwt_index -U $otPre.fa -S $otPre.sam 1>$otPre\_bwt.log 2>$otPre\_bwt.err`;
    print`bwa mem $bwt_index $otPre.unmap.fa >> $otPre.sam 2>>$otPre\_bwt.err`;
    print`perl $stat_var $otPre.sam alleles_IonXpress_001.xls $otPre $ampliconList $ampliconFa $ampliconP 1`;
    print `gzip $otPre.variant_stat.xls $otPre.AllVarant2genotype.readsNum.xls $otPre.reads2haplo.list `;
    print `rm -rf $otPre.fa $otPre.sam $otPre.unmap.fa`;
    exit(0);

}

my%h=();
if(-f "$otPre.overlap.fa"){
    print STDERR "Warning: $otPre.overlap.fa exists and no overlap conducted!\n";
    open OT, "> $otPre.overlap.fil.fa " or die "cannot open 1 $!";
    open IN ,"  $otPre.overlap.fa " or die "cannot open 1 $!";
    $/=">"; <IN>;while(<IN>){chomp;$_=~s/^\@//;my$id=(split /\n+/,$_)[0];my($rid,$s,$mis)=(split /\s+/,$id)[0,2,3];$h{"\@$rid"}="";next if $mis !=0;
    unless($id=~/T$/){print OT ">$_" ;}else{ print OT ">$_" if $s ==0;}}close IN;close OT; $/="\n";
}elsif(-f "$otPre.overlap.fa.gz"){
    print STDERR "Warning: $otPre.overlap.fa.gz exists and no overlap conducted!\n";
    open OT, "> $otPre.overlap.fil.fa " or die "cannot open 1 $!";
    open IN ,"  gzip -cd $otPre.overlap.fa.gz | " or die "cannot open 1 $!";
    $/=">"; <IN>;while(<IN>){chomp;$_=~s/^\@//;my$id=(split /\n+/,$_)[0];my($rid,$s,$mis)=(split /\s+/,$id)[0,2,3];$h{"\@$rid"}="";next if $mis !=0;
    unless($id=~/T$/){print OT ">$_" ;}else{ print OT ">$_"if $s ==0;}}close IN; close OT; $/="\n";
}else{
    print` $overlap -a $fq1 -b $fq2 -c $fq1.left -d $fq2.left -o $otPre.overlap.fa -q $otPre.overlap.qv -m $overlapLen -e $errorRate -n 0.01 -s 2`;
    open IN ,"  $otPre.overlap.fa " or die "cannot open 1 $!";  
    open OT, "> $otPre.overlap.fil.fa " or die "cannot open 1 $!";
    open OT2, "> $otPre.overlap.fa2 " or die "cannot open 1 $!";
    #my%h=();
    $/=">"; <IN>;
    while(<IN>){
        chomp;
        $_=~s/^\@//;
        print OT2 ">$_";
        my$id=(split /\n+/,$_)[0];
        my($rid,$s,$mis)=(split /\s+/,$id)[0,2,3]; 
        $h{"\@$rid"}="";
        next if $mis !=0;
        unless($id=~/T$/){
            print OT ">$_" ;
        }else{ 
            print OT ">$_" if $s ==0; 
        }
    }
    close IN;
    close OT;
    $/="\n";
    print `rm -rf $otPre.overlap.fa && mv $otPre.overlap.fa2 $otPre.overlap.fa`;
}
#goto GO;

my($id,$rid,$rseq,$rqv)=();
open FQ1, " $fq1 " or die "cannot open $fq1 $!";
open FQ2, " $fq2 " or die "cannot open $fq2 $!";
open OT1, "> $fq1.left " or die " cannot write $fq1.left $!";
open OT2, "> $fq2.left " or die " cannot write $fq2.left $!";
while(<FQ1>){
    chomp; 
    $rid=(split /\s+/,$_)[0];
    $id=$_;
    $rseq=<FQ1>;
    <FQ1>;
    $rqv=<FQ1>;
    print OT1"$id\n$rseq+\n$rqv" unless (exists $h{$rid});

    chomp($id=<FQ2>);
    $rseq=<FQ2>;
     <FQ2>;
     $rqv=<FQ2>;
     print OT2"$id\n$rseq+\n$rqv" unless (exists $h{$rid});
}
close FQ1;
close FQ2;
close OT1;
close OT2;

GO:   
print`bowtie2 -f --un $otPre.unmap.fa -x $bwt_index -U $otPre.overlap.fil.fa -S $otPre.sam 1>$otPre\_bwt.log 2>$otPre\_bwt.err`;
print`bwa mem $bwt_index $otPre.unmap.fa >> $otPre.sam 2>>$otPre\_bwt.err`;
print`perl $stat_var $otPre.sam $Bin/alleles_IonXpress_001.xls $otPre\_bwt_overlap $ampliconList $ampliconFa $ampliconP 1`;

print `bwa mem $bwt_index $fq1.left $fq2.left > $otPre.sam_bwa`;
print `bowtie2 -q -x $bwt_index -1 $fq1.left -2 $fq2.left -S $otPre.sam_bwt`;

print `cat $otPre.sam_bwa $otPre.sam_bwt >$otPre.sam_bwta`;


open IN, " $otPre.sam_bwta " or die "cannot open $otPre.sam_bwta $!";
open OT, "> $otPre.sam_com " or die "cannot write $otPre.sam_com $!";
my$flag=0;
while(<IN>){
    chomp;
    if($_=~/^\@/){
        next if $flag!=0;
        $flag=1;
        print OT "$_\n";
        while(1){
            chomp(my$line=<IN>);
            if($line=~/^\@/){
                print  OT "$line\n";
                next;
            }
            $_=$line;
            last;
        }
    }
    ($rid,my$mapFlag,my$dis)=(split /\s+/,$_)[0,1,8];
    next if $mapFlag>=2000; 
    $dis=abs$dis;
    next  if ($dis==0 or $dis>=400);
    $h{$rid}++;
    next if $h{$rid} >2; 
    $_=(split /\s+/,$_,2)[-1];
    print OT "$rid\_$h{$rid}\t$_\n";
}
close IN;
close OT; 

print`perl $stat_var $otPre.sam_com alleles_IonXpress_001.xls $otPre\_bwta_fq $ampliconList $ampliconFa $ampliconP 1 1`;
print`perl $combin $otPre\_bwta_fq.variant_stat.xls $otPre\_bwta_fq 2>$otPre\_bwta_fq.err`;
print `cat $otPre\_bwt_overlap.variant_stat.xls $otPre\_bwta_fq_variant_stat_combin.xls >$otPre.variant_stat.xls`;

