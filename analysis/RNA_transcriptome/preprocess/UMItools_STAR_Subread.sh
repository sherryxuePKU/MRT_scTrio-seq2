############################################
### Tanglab Single-Cell RNA-Seq Pipeline ###
### Author: Zorro Dong                   ###
### Date: 2018-07-25                     ###
############################################

#USAGE: sh UMItools_STAR_Subread.sh sample

sample=$1

######Reference and scripts######
BASE_DIR=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen
DATABASE=${BASE_DIR}/database/refdata-cellranger-GRCh38-3.0.0
PIPELINE=${BASE_DIR}/script/GRCh38
genomeDir=${DATABASE}/star
ref=${DATABASE}/fasta/genome.fa
gtf=${DATABASE}/genes/genes.gtf
barcode=${PIPELINE}/barcode_96_8bp.txt
trim_script=${PIPELINE}/trim_TSO_polyA.pl
norm_script=${PIPELINE}/scRNA_Normalization.R

#######Tools#######
umi_tools=${BASE_DIR}/software/anaconda3/bin/umi_tools
STAR=${BASE_DIR}/software/STAR-2.7.2b/bin/Linux_x86_64/STAR
perl=/usr/bin/perl
seqtk=${BASE_DIR}/software/seqtk/seqtk
featureCounts=${BASE_DIR}/software/subread-2.0.0-Linux-x86_64/bin/featureCounts
samtools=/appsnew/bioapps/samtools-1.9/bin/samtools
Rscript=${BASE_DIR}/software/anaconda3/bin/Rscript

######Step1 Extract Barcode and UMI######
$umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
                  --stdin $sample/${sample}*2.f*q.gz \
                  --stdout $sample/$sample.R1.extracted.fq.gz \
                  --read2-stdout \
                  --read2-in $sample/${sample}*1.f*q.gz \
                  --filter-cell-barcode \
                  --whitelist=$barcode

######Step2 Trim TSO and PolyA######
$perl $trim_script $sample/$sample.R1.extracted.fq.gz $sample/$sample.R1.trim.fq.gz 0
$seqtk trimfq $sample/$sample.R1.trim.fq.gz | gzip - > $sample/$sample.R1.clean.fq.gz

######Step3 Mapping With STAR######
$STAR --runThreadN 4 \
     --genomeDir $genomeDir \
     --readFilesIn $sample/$sample.R1.clean.fq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 1 \
     --outFileNamePrefix $sample/$sample. \
     --outSAMtype BAM SortedByCoordinate

######Step4 Add feature using featureCounts of subread#######
$featureCounts -a $gtf -o $sample/gene_assigned -R BAM $sample/$sample.Aligned.sortedByCoord.out.bam -T 4

#######Step5 sort and index bam file######
$samtools sort -m 15000000000 $sample/$sample.Aligned.sortedByCoord.out.bam.featureCounts.bam -o $sample/$sample.assigned_sorted.bam
$samtools index $sample/$sample.assigned_sorted.bam && rm $sample/$sample.R1.extracted.fq.gz $sample/$sample.R1.trim.fq.gz $sample/$sample.R1.clean.fq.gz $sample/$sample.Aligned.sortedByCoord.out.bam.featureCounts.bam $sample/$sample.Aligned.sortedByCoord.out.bam

######Step6 get umi count table######
$umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I $sample/$sample.assigned_sorted.bam -S $sample/$sample.UMI_counts.tsv

#######Step7 Normalization umi count table######
$Rscript $norm_script $sample
