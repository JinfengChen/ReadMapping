#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR
tabix="/rhome/cjinfeng/software/tools/tabix/tabix-0.2.6"
export PERL5LIB=$PERL5LIB:/opt/vcftools/0.1.8.1/lib/perl5/site_perl/

for file in `ls MSU7_BWA/*.vcf.gz`
do
  gunzip -d $file
done

echo "All Done"

