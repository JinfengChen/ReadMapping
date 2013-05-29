#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR
tabix="/rhome/cjinfeng/software/tools/tabix/tabix-0.2.6"
export PERL5LIB=$PERL5LIB:/opt/vcftools/0.1.8.1/lib/perl5/site_perl/


for bam in `ls MSU7_BWA/HEG4_2.3.MSU7_BWA.Chr*.bam`
do
perl step2_CallVariant.pl --ref ../input/MSU_r7.fa --bam $bam --samtools --project $bam
$tabix/bgzip $bam.flt.vcf
$tabix/tabix -p vcf $bam.flt.vcf.gz
done
/opt/vcftools/0.1.8.1/bin/vcf-concat `ls MSU7_BWA/*.vcf.gz` > MSU7_BWA/HEG4_2.3.MSU7_BWA.merge.vcf
echo "All Done"

