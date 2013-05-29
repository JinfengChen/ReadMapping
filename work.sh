echo "generate sh and run read mapping for MSU7"
perl runMapping.pl --ref ../input/MSU_r7.fa --project MSU7 > log 2> log2 &
perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --resource nodes=1:ppn=12,mem=4G,walltime=20:00:00 MSU7.sh > log 2> log2 &

echo "generate sh and run read mapping for HEG4_RAW"
perl runMapping.pl --ref ../input/HEG4_RAW.fa --project HEG4_RAW > log 2> log2 &
perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --resource nodes=1:ppn=12,mem=4G,walltime=20:00:00 HEG4_RAW.sh > log 2> log2 &

echo "HEG4_luluv2"
perl runMapping.pl --ref ../input/HEG4_luluv2.fa --project HEG4_luluv2 > log 2> log2 &
perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --resource nodes=1:ppn=12,mem=4G,walltime=20:00:00 HEG4_luluv2.sh > log 2> log2 & 

echo "MSU7_BWA for SNP calling"
perl runMapping.pl --ref ../input/MSU_r7.fa --project MSU7_BWA > log 2> log2 &
perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --resource nodes=1:ppn=30,mem=6g,walltime=10:00:00 MSU7_BWA.sh > log 2> log2 &
perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --resource nodes=1:ppn=1,mem=13g,walltime=100:00:00 MSU7_BWA2.sh > log 2> log2 &
perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --resource nodes=1:ppn=8,mem=13g,walltime=100:00:00 MSU7_BWA3.sh > log 2> log2 &
perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --resource nodes=1:ppn=1,mem=13g,walltime=100:00:00 MSU7_BWA4.sh > log 2> log2 &

perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --resource nodes=1:ppn=8,mem=10g,walltime=100:00:00 MSU7_BWA5.sh > log 2> log2 &
perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --resource nodes=1:ppn=30,mem=10g,walltime=100:00:00 MSU7_BWA.3_5k.sh > log 2> log2 &
