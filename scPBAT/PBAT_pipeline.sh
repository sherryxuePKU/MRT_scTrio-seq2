#!/bin/sh

#!/bin/sh

############################################
### PBAT Pipeline                        ###
### Author: Xiaohui Xue                  ###
### Last Modification: 2021-07-26        ###
############################################

wkdir=$1
sp=$2

echo "workdir=$wkdir"
echo "sample=$sp"

# Step_00 Prepare directory
indir=$wkdir/00.raw_data
trim_dir=$wkdir/01.trim
bam_dir=$wkdir/02.bam
methylkit_dir=$wkdir/03.methylkit
info_dir=$wkdir/04.info
tile_dir=$wkdir/05.tile
cnv_dir=$wkdir/06.cnv
#beddir=$wkdir/05.bedfile
echo "=================================="
echo "Prepare directory for PBAT Analysis"
echo "raw data : $indir"
echo "clean data : $trim_dir"
echo "bam data : $bam_dir"
echo "methylation data : $methylkit_dir"
echo "info data : $info_dir"
echo "tile data : $tile_dir"
echo "cnv data : $cnv_dir"

### Software
bin_dir=~/tangfuchou_coe/xuexiaohui/software
conda_PBAT_dir=$bin_dir/miniconda3/envs/PBAT/bin
conda_R4_dir=$bin_dir/miniconda3/envs/R4.0/bin
trim_galore=$conda_PBAT_dir/trim_galore
bismark_exe=$conda_PBAT_dir/bismark
samtools_exe=$conda_PBAT_dir/samtools
perl_exe=$conda_PBAT_dir/perl
methylDackel_exe=$conda_PBAT_dir/MethylDackel
bowtie_dir=$bin_dir/bowtie2/bowtie2-2.3.5
bedtools_exe=$conda_PBAT_dir/bedtools
freec_exe=$bin_dir/FREEC-11.5/src/freec
Rscript_exe=$conda_R4_dir/Rscript

### Script
code_dir=~/tangfuchou_coe/xuexiaohui/script/pipeline/PBAT_xxh
changeID_pl=$code_dir/ChangeReadID.pl
coverage_pl=$code_dir/Coverage_MethRatio.pl
tile_pl=$code_dir/Tile_Meth_Methylkit_v2.pl
basic_summary_pl=$code_dir/basic_sammary_PBAT_xxh_V2.pl
site_anno_r=$code_dir/SiteAnno_merge.R
AbD_pl=$code_dir/AbsoluteDistance.pl
merge_mtx_r=$code_dir/CpGProfile_merge.R

### Database
db_dir=~/tangfuchou_coe/xuexiaohui/database
bismark_db=$db_dir/hg38/Bismark
ref=$db_dir/hg38/Bismark/hg38.genome_lambda.fa
annodir=$db_dir/hg38/Annotation/sub_group
genebody_ref=$db_dir/hg38/Annotation/hg38.gencode.p5.allGene.bed
CGI_ref=$db_dir/hg38/Annotation/hg38.CGI.bed
element='ALR Alu CGI ERV1 ERVK ERVL-MaLR ERVL Exon HCP ICP Intergenic Intragenic Intron L1 L2 LCP LINE LTR MIR Promoter SINE SVA'
freec_db=$db_dir/hg38/FREEC/C2T_ChrFile
chr_len=$db_dir/hg38/FREEC/hg38.genome.len
mappability_ref=$db_dir/hg38/FREEC/out100m2_hg38.gem
subgroup=$db_dir/hg38/Annotation/subgroup.txt
imprint_M=$annodir/hg38.Maternal_Imprinted_Gene.xls
imprint_P=$annodir/hg38.Paternal_Imprinted_Gene.xls

### Pipeline
function do_01_Trim() {
	echo "=================================="
	echo "Trimming begin at: "`date +%Y-%m-%d,%H:%M:%S`
	mkdir -p $trim_dir/$sp
	r1=$indir/$sp/${sp}_R1.fastq.gz
	r2=$indir/$sp/${sp}_R2.fastq.gz
	outdir=$trim_dir/$sp
	
	$trim_galore  --quality 20 --stringency 3 --length 50 \
	--clip_R1 9 --clip_R2 9  --paired       \
	--trim1 --phred33 --gzip --cores  4     \
	--output_dir $outdir  $r1 $r2
	echo "Trimming end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_02_Align() {
	echo "=================================="
	echo "Align begin at: "`date +%Y-%m-%d,%H:%M:%S`
	mkdir -p $bam_dir/$sp
	trim_fq1=${sp}_R1_val_1.fq.gz
	trim_fq2=${sp}_R2_val_2.fq.gz
 
	trim_sam_dir=$trim_dir/$sp
	bam_sam_dir=$bam_dir/$sp

	bam=$bam_sam_dir/${sp}_R1_val_1_bismark_bt2_pe
	unmap_fq1=$bam_sam_dir/${trim_fq1}_unmapped_reads_1.fq.gz
	unmap_fq2=$bam_sam_dir/${trim_fq2}_unmapped_reads_2.fq.gz

	sam1_unmap=$bam_sam_dir/unmap1/${trim_fq1}_unmapped_reads_1_bismark_bt2
	sam2_unmap=$bam_sam_dir/unmap2/${trim_fq2}_unmapped_reads_2_bismark_bt2
  
	$bismark_exe      --fastq  --bowtie2  --non_directional   --unmapped                 \
		--phred33-quals --path_to_bowtie2 $bowtie_dir                            \
		--output_dir $bam_sam_dir --temp_dir $bam_sam_dir $bismark_db             \
		-1 $trim_sam_dir/$trim_fq1 -2 $trim_sam_dir/$trim_fq2                && \
		$samtools_exe view -ubS -t $ref $bam.bam         |\
		$samtools_exe sort -m 900M -T $sp -o $bam.sort.bam

	$bismark_exe      --fastq  --bowtie2  --non_directional   --unmapped                 \
		--phred33-quals --path_to_bowtie2 $bowtie_dir                            \
		--output_dir $bam_sam_dir/unmap1 --temp_dir $bam_sam_dir/unmap1         \
		$bismark_db   $unmap_fq1                                               && \
		$samtools_exe view -uSb -t $ref ${sam1_unmap}.bam  |\
		$samtools_exe sort -m 900M -T $sp -o ${sam1_unmap}.sort.bam

	$bismark_exe      --fastq  --bowtie2  --non_directional   --unmapped                 \
		--phred33-quals --path_to_bowtie2 $bowtie_dir                            \
		--output_dir $bam_sam_dir/unmap2 --temp_dir $bam_sam_dir/unmap2         \
		$bismark_db   $unmap_fq2                                               && \
		$samtools_exe view -uSb -t $ref ${sam2_unmap}.bam  |\
		$samtools_exe sort -m 900M -T $sp -o ${sam2_unmap}.sort.bam

	$perl_exe $changeID_pl $bam.sort.bam       $bam.sort.ReID.bam            && \
	$samtools_exe rmdup    $bam.sort.ReID.bam                                   \
						   $bam.sort.ReID.rmdup.bam                          && \

	$perl_exe $changeID_pl $sam1_unmap.sort.bam $sam1_unmap.sort.ReID.bam    && \
	$samtools_exe rmdup -s $sam1_unmap.sort.ReID.bam                            \
						   $sam1_unmap.sort.ReID.rmdup.bam                   && \

	$perl_exe $changeID_pl $sam2_unmap.sort.bam $sam2_unmap.sort.ReID.bam    && \
	$samtools_exe rmdup -s $sam2_unmap.sort.ReID.bam                            \
						   $sam2_unmap.sort.ReID.rmdup.bam                   && \

	$samtools_exe merge -@ 4 -f $bam_sam_dir/$sp.rmdup.bam                    \
						   $bam.sort.ReID.rmdup.bam                             \
						   $sam1_unmap.sort.ReID.rmdup.bam                      \
						   $sam2_unmap.sort.ReID.rmdup.bam                   && \
		
	$samtools_exe sort -@ 4 -m 900M -T $sp -o $bam_sam_dir/$sp.sort.rmdup.bam          \
		$bam_sam_dir/$sp.rmdup.bam                                   && \

	$samtools_exe index -@ 4 $bam_sam_dir/$sp.sort.rmdup.bam

	rm  $sam1_unmap.sort.bam      $sam2_unmap.sort.bam                          \
		$sam1_unmap.sort.ReID.bam $sam2_unmap.sort.ReID.bam                     \
		$sam1_unmap.sort.ReID.rmdup.bam $sam2_unmap.sort.ReID.rmdup.bam         \
		$bam_sam_dir/*reads*fq.gz   $bam_sam_dir/$sp.rmdup.bam             \
		$bam.sort.bam  $bam.sort.ReID.bam  $bam.sort.ReID.rmdup.bam       \
		$bam_sam_dir/unmap1/*_unmapped_reads.fq.gz      \
		$bam_sam_dir/unmap2/*_unmapped_reads.fq.gz
	echo "Align end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_03_methylkit() {
	echo "=================================="
	echo "Extract methylation site begin at: "`date +%Y-%m-%d,%H:%M:%S`
	mkdir -p $methylkit_dir/$sp

	$methylDackel_exe extract \
	$ref -@ 1 --methylKit --CHG --CHH \
	$bam_dir/$sp/$sp.sort.rmdup.bam \
	-o $methylkit_dir/$sp/$sp.tmp

	grep "Lambda" $methylkit_dir/$sp/$sp.tmp* >\
	$methylkit_dir/$sp/$sp.Lambda.methylkit
	for tp in CpG CHH CHG
	do
	grep -v "Lambda" $methylkit_dir/$sp/${sp}.tmp_${tp}.methylKit  |\
	grep -v "chrBase" > $methylkit_dir/$sp/${sp}_${tp}.methylkit
	rm $methylkit_dir/$sp/${sp}.tmp_${tp}.methylKit
	done

	cat $methylkit_dir/$sp/${sp}_* |\
	sort >$methylkit_dir/$sp/$sp.methylkit
	echo "Extract methylation site end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_04_StatMeth(){
	echo "=================================="
	echo "Stat methylation info begin at: "`date +%Y-%m-%d,%H:%M:%S`
	#mkdir -p $info_dir/$sp
	
	for tp in CpG CHH CHG
	do
		#mkdir -p $info_dir/${tp}_MethRatio
		$perl_exe $coverage_pl $methylkit_dir/$sp/${sp}_${tp}.methylkit \
		$methylkit_dir/$sp/${sp}_${tp}.cov_meth.txt $sp $tp
	done
	echo "Stat methylation info end at: "`date +%Y-%m-%d,%H:%M:%S`	
}

function do_05_MergeStatMeth(){
	echo "=================================="
	echo "Merge stat meth begin at:"`date +%Y-%m-%d,%H:%M:%S`
  sample=$sp
  batch=`echo $sample | awk -F "/" '{print $NF}' | awk -F "_sample.xls" '{print $1}'`
  #sample=$dir/${batch}_sample.xls
  indir=$methylkit_dir
  outdir=$info_dir/$batch
  #batch=`echo $dir | awk -F "/" '{print $NF}'`
  
  mkdir -p $outdir
  for tp in CpG CHH CHG
  do
    echo -e -n "Sample\t${tp}_Total(1X)\t${tp}_TotalUnique(1X)\t${tp}_MeanCov(1X)\t${tp}_MethRatio(1X)\t${tp}_MethReadRatio(1X)\t${tp}_Total$tp(3X)\t${tp}_TotalUnique(3X)\t${tp}_MeanCov(3X)\t${tp}_MethRatio(3X)\t${tp}_MethReadRatio(3X)\t${tp}_Total$tp(5X)\t${tp}_TotalUnique(5X)\t${tp}_MeanCov(5X)\t${tp}_MethRatio(5X)\t${tp}_MethReadRatio(5X)\t${tp}_Total$tp(10X)\t${tp}_TotalUnique(10X)\t${tp}_MeanCov(10X)\t${tp}_MethRatio(10X)\t${tp}_MethReadRatio(10X)\n" > $outdir/${tp}_title.txt
    cat $indir/*/*_${tp}.cov_meth.txt | grep -f $sample | sort -k1,1 >$outdir/${batch}_${tp}.stat_meth.tmp.txt
    cat $outdir/${tp}_title.txt $outdir/${batch}_${tp}.stat_meth.tmp.txt >$outdir/${batch}_${tp}.stat_meth.txt
    #rm -r $indir/${tp}_MethRatio
    rm $outdir/${tp}_title.txt $outdir/${batch}_${tp}.stat_meth.tmp.txt
  done
  
  join $outdir/${batch}_CHH.stat_meth.txt $outdir/${batch}_CHG.stat_meth.txt >$outdir/${batch}_nonCpG.stat_meth.txt
  join $outdir/${batch}_CpG.stat_meth.txt $outdir/${batch}_nonCpG.stat_meth.txt >$outdir/${batch}.stat_meth.txt
  rm $outdir/${batch}_CHH.stat_meth.txt $outdir/${batch}_CHG.stat_meth.txt
	echo "Merge stat meth end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_06_StatBasic(){
	echo "=================================="
	echo "Stat basic info begin at: "`date +%Y-%m-%d,%H:%M:%S`
	bam=$bam_dir/$sp/$sp.sort.rmdup.bam
	outfile=$bam_dir/$sp/$sp.stat_basic.txt
	
	echo -n -e "$sp\t" >$outfile
	## deduplicated reads number
	dedup_reads=`$samtools_exe flagstat $bam | grep "in total" | awk '{print $1}'`
	echo -e -n "$dedup_reads\t" >> $outfile
	## genomic coverage
	cov=`$samtools_exe depth $bam | awk '{if ($1!="Lambda"){g+=1;gsum+=$3};if($1=="Lambda"){l+=1;lsum+=$3}} END {lper=lsum/(gsum+lsum)*100;print g/2934860425*100"\t"l/48510*100"\t"gsum"\t"lsum"\t"lper}'`
	echo "$cov" >> $outfile
	echo "Stat basic info begin at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_07_SiteAnno(){
	echo "=================================="
	echo "Annotate CpG site begin at: "`date +%Y-%m-%d,%H:%M:%S`
	methylkit=$methylkit_dir/$sp/${sp}_CpG.methylkit
	tmp=$methylkit_dir/$sp/${sp}_CpG.tmp.bed
	
	cat $methylkit | awk '{print $2"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6}' > $tmp
	for i in $element
	do
		output=$methylkit_dir/$sp/${sp}.${i}_CpG.bed
		$bedtools_exe intersect \
		-a $tmp -b $annodir/hg38.$i.xls \
		-u > $output
	done
	rm -rf $tmp
	echo "Annotate CpG site end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_08_MergeSiteAnno(){
	echo "=================================="
	echo "Merge site annotation begin at: "`date +%Y-%m-%d,%H:%M:%S`
  batch=`echo $sp | awk -F "/" '{print $NF}' | awk -F "_sample.xls" '{print $1}'`
  $Rscript_exe $site_anno_r $dir $batch
	echo "Merge site annotation end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_09_ComputeMatrix(){
	echo "=================================="
	echo "ComputeMatrix begin at: "`date +%Y-%m-%d,%H:%M:%S`
	region=$1
	input=$methylkit_dir/$sp/${sp}_CpG.methylkit
	output=$methylkit_dir/$sp/${sp}_CpG.$region.absoluteDistance.txt
	
	case $region in
	"genebody") ref=$genebody_ref
	;;
	"CGI") ref=$CGI_ref
	;;
	esac
	
	$perl_exe $AbD_pl $ref $input > $output
	echo "ComputeMatrix end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_10_MergeMatrix() {
	echo "=================================="
	echo "Merge profile matrix begin at: "`date +%Y-%m-%d,%H:%M:%S`
  batch=`echo $sp | awk -F "/" '{print $NF}' | awk -F "_sample.xls" '{print $1}'`
  $Rscript_exe $merge_mtx_r $dir $batch genebody
  $Rscript_exe $merge_mtx_r $dir $batch CGI
	echo "Merge profile matrix end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_11_Tile_split(){
	echo "=================================="
	echo "Split tile begin at: "`date +%Y-%m-%d,%H:%M:%S`
	window=$1
	depth=$2
	mkdir -p $tile_dir/$sp
	
	$perl_exe $tile_pl $methylkit_dir/$sp/${sp}_CpG.methylkit \
	$window $depth CpG \
	$tile_dir/$sp/${sp}.${window}bp_${depth}X_CpG.txt
	echo "Split tile end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_12_CNV(){
	echo "=================================="
	echo "Calling CNV begin at: "`date +%Y-%m-%d,%H:%M:%S`
	window=$1
	mkdir -p $cnv_dir/$sp/${window}M
	
	bam=$bam_dir/$sp/$sp.sort.rmdup.bam
	config=$cnv_dir/$sp/${window}M/$sp.${window}M_config.tmp.txt
	echo "
	[general]

	chrLenFile = $chr_len
	ploidy = 2
	breakPointThreshold = .6
	chrFiles = $freec_db
	window = ${window}000000
	maxThreads = 2
	outputDir = $cnv_dir/$sp/${window}M
	BedGraphOutput = TRUE
	samtools = $samtools_exe
	minExpectedGC = 0.17
	maxExpectedGC = 0.27
	sex=XY
	#gemMappabilityFile = $mappability_ref
	[sample]

	mateFile = $bam
	inputFormat = BAM
	mateOrientation = 0
	" > $config
	
	$freec_exe -conf $config
	echo "Calling CNV end at: "`date +%Y-%m-%d,%H:%M:%S`
}

function do_13_MergeCNVRes() {
  echo "=================================="
	echo "Merge CNV result begin at:"`date +%Y-%m-%d,%H:%M:%S`
  sample=$sp
  batch=`echo $sample | awk -F "/" '{print $NF}'|awk -F "_sample.xls" '{print $1}'`
  indir=$cnv_dir
  outdir=$info_dir/$batch
  window=1
  #batch=`echo $dir | awk -F "/" '{print $NF}'`
  
  mkdir -p $outdir
  fsp=`head -n 1 $sample`
  cut -f 1,2 $indir/$fsp/${window}M/*_sample.cpn > $outdir/$batch.${window}M_Merge_read.txt
  cut -f 1,2 $indir/$fsp/${window}M/*_sample.cpn > $outdir/$batch.${window}M_Merge_GC.txt
  cut -f 1,2 $indir/$fsp/${window}M/*_sample.cpn > $outdir/$batch.${window}M_Merge_cpn.txt
  echo -e -n "Chr\tStart" > $outdir/$batch.${window}M_CNV_header.txt
  
  cat $sample | while read sp
  do
    echo -e -n "\t$sp" >> $outdir/$batch.${window}M_CNV_header.txt
  
    cut -f 4 $indir/$sp/${window}M/*bam_ratio.txt | grep -v "MedianRatio" >$indir/$sp/${window}M/$sp.MedianRatio.txt
    paste $outdir/$batch.${window}M_Merge_GC.txt $indir/$sp/${window}M/$sp.MedianRatio.txt >$outdir/$batch.${window}M_Merge_GC.tmp.txt
    mv $outdir/$batch.${window}M_Merge_GC.tmp.txt $outdir/$batch.${window}M_Merge_GC.txt
    rm $indir/$sp/${window}M/$sp.MedianRatio.txt
  
    cut -f 5 $indir/$sp/${window}M/*bam_ratio.txt | grep -v "CopyNumber" >$indir/$sp/${window}M/$sp.CopyNumber.txt
    paste $outdir/$batch.${window}M_Merge_cpn.txt $indir/$sp/${window}M/$sp.CopyNumber.txt >$outdir/$batch.${window}M_Merge_cpn.tmp.txt
    mv $outdir/$batch.${window}M_Merge_cpn.tmp.txt $outdir/$batch.${window}M_Merge_cpn.txt
    rm $indir/$sp/${window}M/$sp.CopyNumber.txt
  
    cut -f 3 $indir/$sp/${window}M/*bam_sample.cpn >$indir/$sp/${window}M/$sp.read.txt
    paste $outdir/$batch.${window}M_Merge_read.txt $indir/$sp/${window}M/$sp.read.txt >$outdir/$batch.${window}M_Merge_read.tmp.txt
    mv $outdir/$batch.${window}M_Merge_read.tmp.txt $outdir/$batch.${window}M_Merge_read.txt
    rm $indir/$sp/${window}M/$sp.read.txt
  
  done
  
  echo -e -n "\n" >> $outdir/$batch.${window}M_CNV_header.txt
  
  cat $outdir/$batch.${window}M_CNV_header.txt $outdir/$batch.${window}M_Merge_GC.txt >$outdir/$batch.${window}M_Merge_GC_Final.txt
  mv $outdir/$batch.${window}M_Merge_GC_Final.txt $outdir/$batch.${window}M_Merge_GC.txt
  
  cat $outdir/$batch.${window}M_CNV_header.txt $outdir/$batch.${window}M_Merge_cpn.txt >$outdir/$batch.${window}M_Merge_cpn_Final.txt
  mv $outdir/$batch.${window}M_Merge_cpn_Final.txt $outdir/$batch.${window}M_Merge_cpn.txt
  echo "Merge CNV result end at:"`date +%Y-%m-%d,%H:%M:%S`
}

function do_14_MergeStatBasic(){
	echo "=================================="
	echo "Merge stat basic begin at:"`date +%Y-%m-%d,%H:%M:%S`
	# mkdir -p $info_dir
	# prefix=$1
	# $perl_exe $basic_summary_pl $wkdir $info_dir/$prefix.stat_basic.txt
	
	## raw reads and raw bases
	trim_file_r1=$trim_dir/$sp/${sp}_R1.fastq.gz_trimming_report.txt
	total_read=`grep "Total reads processed:" $trim_file_r1 | awk -F ":" '{gsub(/^[ ]+|,/,"",$2);print $2*2}'`
	total_base=`grep "Total basepairs processed:" $trim_file_r1  | awk -F ":" '{gsub(/,| |bp/,"",$2);print $2}'`
	
	## clean reads
	pe_file=$bam_dir/$sp/${sp}_R1_val_1_bismark_bt2_PE_report.txt
	se_file_r1=$bam_dir/$sp/unmap1/${sp}_R1_val_1.fq.gz_unmapped_reads_1_bismark_bt2_SE_report.txt
	se_file_r2=$bam_dir/$sp/unmap2/${sp}_R2_val_2.fq.gz_unmapped_reads_2_bismark_bt2_SE_report.txt
	clean_read=`grep "Sequence pairs analysed in total:" $pe_file | awk -F ":" '{print $2*2}'`
		
	## mapped reads
	pe_map_read=`grep "Number of paired-end alignments with a unique best hit:" $pe_file | awk -F ":" '{gsub(/\t/,"",$2);print $2*2}'`
	se_r1_map_read=`grep "Number of alignments with a unique best hit from the different alignments:" $se_file_r1 | awk -F ":" '{gsub(/\t/,"",$2);print $2}'`
	se_r2_map_read=`grep "Number of alignments with a unique best hit from the different alignments:" $se_file_r2 | awk -F ":" '{gsub(/\t/,"",$2);print $2}'`
	map_read=`echo ${pe_map_read}+${se_r1_map_read}+${se_r2_map_read} | bc`
	
	## dedup reads
	cov_depth_file=$bam_dir/$sp/$sp.stat_basic.txt
	dedup_read=`cat $cov_depth_file | awk '{print $2}'`
	dedup_base=`cat $cov_depth_file | awk '{print $5+$6}'`
	genomic_cov=`cat $cov_depth_file | awk '{print $3}'`
	genomic_base=`cat $cov_depth_file | awk '{print $5}'`
	genomic_depth=`cat $cov_depth_file | awk '{print $5/2934860425}'`
	lambda_percent=`cat $cov_depth_file | awk '{print $7}'`
	
	## conversion ratio
	lambda_methylkit=$methylkit_dir/$sp/$sp.Lambda.methylkit
	conv_ratio=`cat $lambda_methylkit | awk '{sum+=$5;unmeth+=$5*$7}END{print unmeth/sum}'`
	
	outfile=$info_dir/basic_summary/$sp.basic_summary.txt
	mkdir -p $info_dir/basic_summary
	
	echo -e -n "Sample\tTotal_base\tTotal_read\tClean_read\tMapped_read\tDedup_read\tDedup_base\tGenomic_base\tGenomic_coverage\tGenomic_depth\tLambda_percent\tConversion_ratio\n" > $outfile
	echo -e -n "$sp\t$total_base\t$total_read\t$clean_read\t$map_read\t$dedup_read\t$dedup_base\t$genomic_base\t$genomic_cov\t$genomic_depth\t$lambda_percent\t$conv_ratio\n" >> $outfile

	echo "Merge stat basic end at:"`date +%Y-%m-%d,%H:%M:%S`
}
