#!/usr/bin/perl -w
#use strict;

# methylkit format
# id	chromosome	position	strand	total_read	meth_ratio	unmeth_ratio
# chr1.10473      chr1    10473   R       1         0.00  100.00

my %promoter;
open IN,"$ARGV[0]" or die $!; # ../../../../annotation/hg19.Intragenic.xls
while(<IN>){
	chomp;
	my @a = split(/\t/,$_);
	my $s5 = $a[1];
	my $s3 = $a[2];
	if ($a[4] eq "-"){
		$s3 = $a[1];
		$s5 = $a[2];
	}
	$promoter{$a[0]}{$a[3]}{$s5} = $s3;
}
close IN;

#my $DEPTH = 3;
my $DEPTH = 1; ## added by xxh
my %CpG_meth;
my %CpG_num;
open IN,"$ARGV[1]" or die $!; # ../8cell.SingleCmet
while(<IN>){
	chomp;
	next if (/#/);
	my @a = split /\s+/,$_;
	#next if ($a[4] < $DEPTH);
  next if ($a[1] eq "chrM" or $a[0] =~ /random/);
	$CpG_meth{$a[1]}{$a[2]} = $a[4]*$a[5]; # added by xxh for methylkit
  $CpG_num{$a[1]}{$a[2]} = $a[4]
  #$CpG{$a[1]}{$a[2]} = $a[5]; #methylation level, modified in 20210201 
}
close IN;


my $bin_size = 300;
my $bin_num = 200;

foreach my $chr (sort keys %promoter){
	foreach my $gene (sort keys %{$promoter{$chr}}){
		foreach my $s5 (sort keys %{$promoter{$chr}{$gene}}){
			my $s3 = $promoter{$chr}{$gene}{$s5};
			my $length = abs($s5-$s3);
			my $bin = int($length/$bin_num);
			if ($s5 > $s3){
				print "$gene";
				for(my $pos=$s5+15000;$pos-($bin_size - 1)>=$s5;$pos-=$bin_size){
                                        my ($sum,$n) = (0,0);
                                        foreach my $i (($pos-($bin_size-1))..$pos){
                                                if (exists $CpG_meth{$chr}{$i}){
                                                        $sum += $CpG_meth{$chr}{$i};
                                                        #$n ++;
                                                        $n += $CpG_num{$chr}{$i};
                                                }
                                        }
                                        my $mean = ($n != 0) ? $sum/$n : "NA";
                                        print "\t$mean";
                                }
				my $pos = $s5;
                        	for(my $i=1;$i<=$bin_num;$i++){
                                	my ($sum,$n) = (0,0);
                                	foreach my $j (($pos-$bin)..$pos){
                                        	if (exists $CpG_meth{$chr}{$j}){
                                                	$sum += $CpG_meth{$chr}{$j};
                                                  #$n ++;
                                                  $n += $CpG_num{$chr}{$j};
                                        	}
					}
                                	my $mean = ($n != 0) ? $sum/$n : "NA";
                                	print "\t$mean";
					$pos -= $bin;
                        	}
				for(my $pos=$s3;$pos-($bin_size - 1)>=$s3-15000;$pos-=$bin_size){
                                        my ($sum,$n) = (0,0);
                                        foreach my $i (($pos-$bin_size)..$pos){
                                                if (exists $CpG_meth{$chr}{$i}){
                                                        $sum += $CpG_meth{$chr}{$i};
                                                        #$n ++;
                                                        $n += $CpG_num{$chr}{$i};
                                                }
                                        }
                                        my $mean = ($n != 0) ? $sum/$n : "NA";
                                        print "\t$mean";
                                }
                        	print "\n";
			}
			else{
				print "$gene";
				for(my $pos=$s5-15000;$pos+($bin_size - 1)<=$s5;$pos+=$bin_size){
                                        my ($sum,$n) = (0,0);
                                        foreach my $i ($pos..($pos+($bin_size - 1))){
                                                if (exists $CpG_meth{$chr}{$i}){
                                                        $sum += $CpG_meth{$chr}{$i};
                                                        #$n ++;
                                                        $n += $CpG_num{$chr}{$i};
                                                }
                                        }
                                        my $mean = ($n != 0) ? $sum/$n : "NA";
                                        print "\t$mean";
                                }
				my $pos = $s5;
                                for(my $i=1;$i<=$bin_num;$i++){
                                        my ($sum,$n) = (0,0);
                                        foreach my $j ($pos..($pos+$bin)){
                                                if (exists $CpG_meth{$chr}{$j}){
                                                	$sum += $CpG_meth{$chr}{$j};
                                                  #$n ++;
                                                  $n += $CpG_num{$chr}{$j};
                                                }
                                        }
                                        my $mean = ($n != 0) ? $sum/$n : "NA";
                                        print "\t$mean";
					$pos += $bin;
                                }
				for(my $pos=$s3;$pos+($bin_size - 1)<=$s3+15000;$pos+=$bin_size){
                                        my ($sum,$n) = (0,0);
                                        foreach my $i ($pos..($pos+($bin_size - 1))){
                                                if (exists $CpG_meth{$chr}{$i}){
                                                        $sum += $CpG_meth{$chr}{$i};
                                                        #$n ++;
                                                        $n += $CpG_num{$chr}{$i};
                                                }
                                        }
                                        my $mean = ($n != 0) ? $sum/$n : "NA";
                                        print "\t$mean";
                                }
				print "\n";
			}
		}
	}
}
