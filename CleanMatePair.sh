#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR

#/usr/local/bin/samtools view -h FC67_3.MSU7_BWA.bam | sort | perl -e '{while(<>){if($_=~/^@/){print $_}else{my @read1=split("\t",$_);my @read2=split("\t",<>);if($read1[5]=~/^\d+M$/ and $read2[5]=~/^\d+M$/){my $r1=join("\t",@read1);my $r2=join("\t",@read2);print $r1,$r2}}}}' | samtools view -Sb -o FC67_3.MSU7_BWA.clean.bam -

#prefix=FC67_3.MSU7_BWA
#prefix=FC70_5.MSU7_BWA
prefix=FC20_3.MSU7_BWA

##clean chemiric reads in mating pairs mapping bam file
/usr/local/bin/samtools view -H $prefix.bam > $prefix.header
/usr/local/bin/samtools view $prefix.bam | sort | perl -e '{while(<>){if($_=~/^@/){print $_}else{my @read1=split("\t",$_);my @read2=split("\t",<>);while($read1[0] ne $read2[0]){@read1=@read2;@read2=split("\t",<>)};if($read1[5]=~/^\d+M$/ and $read2[5]=~/^\d+M$/){my $r1=join("\t",@read1);my $r2=join("\t",@read2);print $r1,$r2}}}}' > $prefix.clean.sam 
cat $prefix.header $prefix.clean.sam > $prefix.clean.header.sam
cat $prefix.clean.header.sam | samtools view -Sb -o $prefix.clean.bam -
rm $prefix.header $prefix.clean.sam $prefix.clean.header.sam


