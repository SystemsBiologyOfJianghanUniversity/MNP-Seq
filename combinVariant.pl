#!/usr/bin/perl -w
use strict;

die("Usage: perl $0 <variant_stat.xls> <otPre>\n") unless @ARGV ==2;
my($varF,$otPre)=@ARGV;

my(%line)=();
open IN, " $varF " or die "canot open $varF $!";
open OT, "> $otPre\_variant_stat_combin.xls " or die "2 $!";
while(<IN>){
    chomp;
    next if /^\#/;
    next if /^\s*$/;
    next if /^\@/;
    my($rID,$aid,$chr,$s,$e,$variant)=(split /\s+/,$_);
    die("reads ID is in wrong format\n") unless $rID=~/(\S+)\_([12])/;
    $line{$1}{$2}=$_;
}
close IN;

foreach my$rid(keys %line){
    unless((scalar keys %{$line{$rid}})==2){
        foreach my$k(keys %{$line{$rid}}){
            my($aid,$chr)=(split /\s+/,$line{$rid}{$k})[1,2];
            print STDERR "No paired reads found:\t$rid\t$aid\t$chr\n";
        }
        next;
    }
    my@reads1=(split /\s+/,$line{$rid}{1});
    my@reads2=(split /\s+/,$line{$rid}{2});
    
    if($reads1[1] ne $reads2[1]){
        print STDERR "Different aid found for paired reads: $rid\n";
        next;
    }elsif($reads1[2] ne $reads2[2]){        
        print STDERR "Different chr found for paired reads: $rid\n";
        next;
    }

    my($alignS,$alignE)=(sort {$a<=>$b} ($reads1[3],$reads1[4],$reads2[3],$reads2[4]))[0,-1];

    my@var1=(split /\|/,$reads1[-1]);
    my@var2=(split /\|/,$reads2[-1]);

    my$align="";
    my$flag=0;
    for(my$i=0;$i<=$#var1;$i++){
        if($var1[$i]=~/,0$/){
            $align.="$var2[$i]|";
        }elsif($var2[$i]=~/,0$/){
            $align.="$var1[$i]|";
        }else{
            $align.="$var1[$i]|";
            $flag=1 if($var1[$i] ne $var2[$i]);       
        }
    }

    if($flag==0){
        print OT "$rid\t$reads1[1]\t$reads1[2]\t$alignS\t$alignE\t$align\n";
        next;
    }
    print STDERR "#unconsistentSite:\t$rid\t$reads1[1]\t$reads1[2]\t$alignS\t$alignE\t$align\n";
}


        


