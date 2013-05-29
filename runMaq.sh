#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR

perl /rhome/cjinfeng/HEG4_cjinfeng/MappingReads/bin/step1_Mapping.pl -ref ../input/MSU_r7.fa -1 /rhome/cjinfeng/HEG4_cjinfeng/fastq/errorcorrection/soapec/HEG4_2_200bp/HEG4_2.3_p1.fq -2 /rhome/cjinfeng/HEG4_cjinfeng/fastq/errorcorrection/soapec/HEG4_2_200bp/HEG4_2.3_p2.fq -min 100 -max 500 -cpu 12 -tool maq -project ./HEG4_2.3.map

echo "Done"

