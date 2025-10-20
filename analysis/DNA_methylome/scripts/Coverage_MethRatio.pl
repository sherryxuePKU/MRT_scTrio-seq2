#!/usr/bin/perl -w
use strict;
use warnings;
my $file=$ARGV[0];
open OUT," >$ARGV[1]" or die $!;
my $type=$ARGV[3];

print OUT "Sample\t${type}_Total(1X)\t${type}_TotalUnique(1X)\t${type}_MeanCov(1X)\t${type}_MethRatio(1X)\t${type}_MethReadRatio(1X)\t${type}_Total$type(3X)\t${type}_TotalUnique(3X)\t${type}_MeanCov(3X)\t${type}_MethRatio(3X)\t${type}_MethReadRatio(3X)\t${type}_Total$type(5X)\t${type}_TotalUnique(5X)\t${type}_MeanCov(5X)\t${type}_MethRatio(5X)\t${type}_MethReadRatio(5X)\t${type}_Total$type(10X)\t${type}_TotalUnique(10X)\t${type}_MeanCov(10X)\t${type}_MethRatio(10X)\t${type}_MethReadRatio(10X)\n";


my @depth=(1,3,5,10);
my $length=@depth;
my (@sum,@num,$name,@meth,@meth_ratio);
$name=$ARGV[2];
my $j;

#$num[0]=$num[1]=$num[2]=$num[3]=0;
#open(FH,"gzip -dc $file |") or die "fail open the file";
open(FH, $file) or die "fail open the file";
while(<FH>){
	chomp;
	my $line=$_;
	#next if ($line =~/#/);
  next if ($line =~/chrBase/); ## added by xxh
	next if ($line =~/Lambda/); ## added by xxh
#	print "$line\n";
	my @words=split(/\s+/,$line);
	for($j=0;$j<$length;$j++){
		if($words[4]>=$depth[$j]){
			#print "@words\n";
			$sum[$j]+=$words[4];
			#print "total\t$words[4]\t";
			$num[$j]++;
			$meth[$j]+=$words[4]*$words[5];
			#print "meth\t$words[5]\t";
			$meth_ratio[$j]+=$words[5];
			#print "methRatio\t$words[7]\n";
		}
	}
}


my $ratio_1 = $sum[0]/$num[0];
my $Meth_Ratio_1 = $meth_ratio[0]/$num[0];
my $Meth_Read_1 = $meth[0]/$sum[0];

my ($ratio_3,$Meth_Ratio_3,$Meth_Read_3);

if($num[1]==0){
	$sum[1] = $ratio_3 = $Meth_Ratio_3 = $Meth_Read_3 = 0 ;
}else {
	$ratio_3 = $sum[1]/$num[1];
	$Meth_Ratio_3 = $meth_ratio[1]/$num[1];
	$Meth_Read_3 = $meth[1]/$sum[1];
}
#my $ratio_3 = $sum[1]/$num[1];
#my $Meth_Ratio_3 = $meth_ratio[1]/$num[1];
#my $Meth_Read_3 = $meth[1]/$num[1];

my ($ratio_5,$Meth_Ratio_5,$Meth_Read_5);
if($num[2]==0){
	$sum[2] = $ratio_5 = $Meth_Ratio_5 = $Meth_Read_5 =0;
}else{
	$ratio_5 =  $sum[2]/$num[2];
	$Meth_Ratio_5 = $meth_ratio[2]/$num[2];
	$Meth_Read_5 = $meth[2]/$sum[2];
}

#my $ratio_5 =  $sum[2]/$num[2];
#my $Meth_Ratio_5 = $meth_ratio[2]/$num[2];
#my $Meth_Read_5 = $meth[2]/$num[2];

my ($ratio_10, $Meth_Ratio_10, $Meth_Read_10);
if($num[3]==0){
	$sum[3] = $ratio_10 = $Meth_Ratio_10 = $Meth_Read_10 = 0;
}else{
	$ratio_10 = $sum[3]/$num[3];	
	$Meth_Ratio_10 = $meth_ratio[3]/$num[3];
	$Meth_Read_10 = $meth[3]/$sum[3];
}
#my $ratio_10 = $sum[3]/$num[3];
#my $Meth_Ratio_10 = $meth_ratio[3]/$num[3];
#my $Meth_Read_10 = $meth[3]/$num[3];

print OUT "$name\t$sum[0]\t$num[0]\t$ratio_1\t$Meth_Ratio_1\t$Meth_Read_1\t$sum[1]\t$num[1]\t$ratio_3\t$Meth_Ratio_3\t$Meth_Read_3\t$sum[2]\t$num[2]\t$ratio_5\t$Meth_Ratio_5\t$Meth_Read_5\t$sum[3]\t$num[3]\t$ratio_10\t$Meth_Ratio_10\t$Meth_Read_10\n";
