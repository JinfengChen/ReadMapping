#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR
tabix="/rhome/cjinfeng/software/tools/tabix/tabix-0.2.6"
export PERL5LIB=$PERL5LIB:/opt/vcftools/0.1.8.1/lib/perl5/site_perl/
#cp var.flt.vcf var.flt1.vcf
#cp var.flt.vcf var.flt2.vcf
#$tabix/bgzip var.flt1.vcf
#$tabix/bgzip var.flt2.vcf
#$tabix/tabix -p vcf var.flt1.vcf.gz
#$tabix/tabix -p vcf var.flt2.vcf.gz
/opt/vcftools/0.1.8.1/bin/vcf-concat `ls *.vcf.gz` > merge1.vcf


