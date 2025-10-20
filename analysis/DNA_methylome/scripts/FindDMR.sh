#!/bin/sh

############################################
### PBAT Pipeline                        ###
### Author: Xiaohui Xue                  ###
### Last Modification: 2021-07-26        ###
############################################

dir=$1
sp=$2
prefix=$3

echo "workdir=$dir"
echo "sample=$sp"

### Software
bin_dir=~/tangfuchou_coe/xuexiaohui/software
conda_PBAT_dir=$bin_dir/miniconda3/envs/PBAT/bin
samtools_exe=$conda_PBAT_dir/samtools
methylDackel_exe=$conda_PBAT_dir/MethylDackel
perl_exe=$conda_PBAT_dir/perl
bedtools_exe=$conda_PBAT_dir/bedtools

### Database
db_dir=~/tangfuchou_coe/xuexiaohui/database
ref=$db_dir/hg38/Bismark/hg38.genome_lambda.fa
annodir=$db_dir/hg38/Annotation/sub_group
genebody_ref=$db_dir/hg38/Annotation/hg38.gencode.p5.allGene.bed
promoter_ref=$db_dir/hg38/Annotation/hg38.gencode.p5.allGene_promoter.bed
element='ALR Alu ERV1 ERVK ERVL-MaLR ERVL Exon HCP ICP Intergenic Intragenic Intron L1 L2 LCP LINE LTR MIR SINE SVA'

### Script
code_dir=~/tangfuchou_coe/xuexiaohui/script/pipeline/PBAT_xxh
tile_pl=$code_dir/Tile_Meth_Methylkit_v2.pl

### Merge group bams
outdir=$dir/07.dmr/${prefix}
mkdir -p $outdir
bamlist=$outdir/${prefix}_bamlist.txt
merged_bam=$outdir/${prefix}_merged.bam

cat $sp | awk '{print "'$dir'""/02.bam/"$1"/"$1".sort.rmdup.bam"}' > $bamlist

function do_01_MergeGroupBam(){
  $samtools_exe merge -@ 3 -b $bamlist $merged_bam
}

### Convert bam to methylkit
function do_02_Bam2Methylkit(){
  methylkit_prefix=$outdir/${prefix}_merged
  
  $methylDackel_exe extract \
  $ref -@ 1 --methylKit \
  --CHG --CHH \
  $merged_bam \
  -o $methylkit_prefix.tmp
  
  grep "Lambda" $methylkit_prefix.tmp* >\
  $methylkit_prefix.Lambda.methylkit
  
  for tp in CpG CHH CHG
  do
    cat $methylkit_prefix.tmp_${tp}.methylKit |\
    awk '{if(($1!="chrBase" && $1!~/^Lambda/) && ($6 <10 || $6 >90))print $0}' >\
    ${methylkit_prefix}.${tp}.methylkit
  	rm $methylkit_prefix.tmp_${tp}.methylKit
  done
}


### Annotation
function do_03_AnnoSite(){
  methylkit_prefix=$outdir/${prefix}_merged
  tmp=$outdir/${prefix}_merged_CpG.tmp.bed
	
  cat $methylkit_prefix.CpG.methylkit |\
  awk '{print $2"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6}' > $tmp
  
	for i in $element
	do
		output=$methylkit_prefix.${i}_CpG.bed
		$bedtools_exe intersect \
		-b $tmp -a $annodir/hg38.$i.xls -wa -wb |\
		$bedtools_exe groupby -i - -g 1-4 -c 11 -o mean > $output
	done
 
  output=$methylkit_prefix.gene_CpG.bed
  $bedtools_exe intersect \
  -b $tmp -a $genebody_ref -wa -wb |\
	$bedtools_exe groupby -i - -g 1-4 -c 11 -o mean > $output
 
  output=$methylkit_prefix.CGI_CpG.bed
  $bedtools_exe intersect \
  -b $tmp -a $annodir/hg38.CGI.xls -wa -wb |\
	$bedtools_exe groupby -i - -g 1-4 -c 10 -o mean > $output
 
  output=$methylkit_prefix.promoter_gene_CpG.bed
  $bedtools_exe intersect \
  -b $tmp -a $promoter_ref -wa -wb |\
	$bedtools_exe groupby -i - -g 1-4 -c 11 -o mean > $output
	rm -rf $tmp
}

