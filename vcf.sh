#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR

/opt/samtools/0.1.18/bin/samtools mpileup -uf ../input/MSU_r7.fa MSU7_BWA/HEG4_2.3.MSU7_BWA.Chr12.bam | /opt/samtools/0.1.18/bin/bcftools view -bvcg - > var.raw.bcf
/opt/samtools/0.1.18/bin/bcftools view var.raw.bcf | /opt/samtools/0.1.18/bin/vcfutils.pl varFilter -D100 > var.flt.vcf

echo "done"

