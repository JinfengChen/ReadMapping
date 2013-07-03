#!/usr/bin/perl
=header
The scripts is designed to call SNP/INDEL with SAMtools or GATK.
-ref:     reference sequence
-bam:     sorted bam alignment file of mapped reads after rmdup
-realign  Realigning the BAM file using GATK's RealignerTargetCreator and IndelRealigner (optional)
-recal    Recalibration (optional, not added)
-samtools call variant using samtools (optional)
-GATK     call variant using GATK (optional)
-dbSNP    vcf of dbSNP used to recalibrate base quality and variant quality
-dbINDEL  vcf of dbINDEL used to recalibate variant qulaity of indels
-snpvcf   SNP vcf of GATK call used to do hard filter or VQSR
-indelvcf INDEL vcf of GATK call used to do hard filter or VQSR
-project: project name that used for result file (optional)
=cut

use Getopt::Long;
my %opt;
GetOptions(\%opt,"ref:s","bam:s","mem:s","samtools","realign","recal","dbSNP:s","snpvcf:s","indelvcf:s","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 -ref ../input/all.con -bam out.bam\n";
   exit();
}
$opt{realign} = 1 if ($opt{recal});
$opt{project} ||= "HEG4_BWA";
$opt{mem} ||= 10;
my $bwa="/opt/tyler/bin/bwa";
my $SAMtool="/opt/samtools/0.1.18/bin/samtools";
my $BCFtool="/opt/samtools/0.1.18/bin/bcftools";
my $vcfutils="/opt/samtools/0.1.18/bin/vcfutils.pl";
my $GATK="/opt/GATK/2.4-3-g2a7af43/GenomeAnalysisTK.jar";
my $java="/opt/java/jdk1.6.0_38/bin/java";
my $picard="/opt/picard/1.81";

### reference
my $name=$1 if $opt{ref}=~/(.*)\.f*a$/;
unless (-e "$name.dict"){
   `$java -jar $picard/CreateSequenceDictionary.jar R=$opt{ref} O=$name.dict`;
}
unless (-e "$opt{ref}.fai"){
   `$SAMtool faidx $opt{ref}`;
}

### realign
if ($opt{realign}){
   print "Realigning the BAM file using GATK's RealignerTargetCreator and IndelRealigner\n";
   `$java -Xmx$opt{mem}g -jar $GATK -T RealignerTargetCreator -R $opt{ref} -I $opt{bam} -o $opt{project}.gatk.intervals > $opt{project}.gatk.interval.log 2> $opt{project}.gatk.interval.log2` unless (-e "$opt{project}.gatk.intervals");
   `$java -Xmx$opt{mem}g -jar $GATK -T IndelRealigner -R $opt{ref} -I $opt{bam} -targetIntervals $opt{project}.gatk.intervals -o $opt{project}.realign.bam > $opt{project}.gatk.realign.log 2> $opt{project}.gatk.realign.log2` unless (-e "$opt{project}.realign.bam"); 
}
### recalibration, assume that do realign before this (recalibrator need to be done for whole-genome, as designed in the method)
### the memory useage is linear increace with the read group number in bam
if ($opt{recal}){
   print "BaseRecalibration with BaseRecalibrator and PrintReads\n";
   `$java -Xmx$opt{mem}g -jar $GATK -T BaseRecalibrator -R $opt{ref} -I $opt{project}.realign.bam -knownSites $opt{dbSNP} -o $opt{project}.recal.grp > $opt{project}.gatk.recal.log 2> $opt{project}.gatk.recal.log2` unless (-e "$opt{project}.recal.grp");
   `$java -Xmx$opt{mem}g -jar $GATK -T PrintReads -R $opt{ref} -I $opt{project}.realign.bam -BQSR $opt{project}.recal.grp -o $opt{project}.realign.recal.bam > $opt{project}.gatk.print.log 2> $opt{project}.gatk.print.log2` unless (-e "$opt{project}.realign.recal.bam");
}

### call variants
if ($opt{samtools}){
   print "Call variant using SAMtools/BCFtools!\n";
   if ($opt{recal}){         ##do realign and do recalibrate
      print "Using realigned and recalibrated bam file!\n";
      `$SAMtool mpileup -uf $opt{ref} $opt{project}.realign.recal.bam | $BCFtool view -bvcg - > $opt{project}.realign.recal.raw.bcf 2> $opt{project}.realign.recal.raw.bcf.log2`;
      `$BCFtool view $opt{project}.realign.recal.raw.bcf | $vcfutils varFilter -D200 > $opt{project}.realign.recal.flt.vcf 2> $opt{project}.realign.recal.flt.vcf.log2`;
      print "Done!\n";
   }elsif($opt{realign}){    ##only do realign
      print "Using realigned bam file!\n";
      `$SAMtool mpileup -uf $opt{ref} $opt{project}.realign.bam | $BCFtool view -bvcg - > $opt{project}.realign.raw.bcf 2> $opt{project}.realign.raw.bcf.log2`;
      `$BCFtool view $opt{project}.realign.raw.bcf | $vcfutils varFilter -D200 > $opt{project}.realign.flt.vcf 2> $opt{project}.realign.flt.vcf.log2`;
      print "Done!\n";
   }else{                    ##use original bam
      print "Using original bam file!\n";
      `$SAMtool mpileup -uf $opt{ref} $opt{bam} | $BCFtool view -bvcg - > $opt{project}.raw.bcf 2> $opt{project}.raw.bcf.log2`;
      `$BCFtool view $opt{project}.raw.bcf | $vcfutils varFilter -D200 > $opt{project}.flt.vcf 2> $opt{project}.flt.vcf.log2`;
      print "Done!\n";
   }
}elsif($opt{GATK}){ ### GATK
   print "Call variant using GATK!\n";   
   if ($opt{recal}){         ##do realign and do recalibrate
      print "Using realigned and recalibrated bam file!\n";
      `$SAMtool index $opt{project}.realign.recal.bam` unless (-e "$opt{project}.realign.recal.bam.bai");
      ####SNPs
      `$java -Xmx2g -jar $GATK -T UnifiedGenotyper -R $opt{ref} -I $opt{project}.realign.recal.bam -o $opt{project}.realign.recal.gatk.snp.raw.vcf -stand_call_conf 30 -stand_emit_conf 10 -glm SNP > $opt{project}.gatk.realign.recal.snp.log 2> $opt{project}.gatk.realign.recal.snp.log2` unless (-e "$opt{project}.realign.recal.gatk.snp.raw.vcf");
      $opt{snpvcf} ||= "$opt{project}.realign.recal.gatk.snp.raw.vcf";
      ####INDELs
      `$java -Xmx2g -jar $GATK -T UnifiedGenotyper -R $opt{ref} -I $opt{project}.realign.recal.bam -o $opt{project}.realign.recal.gatk.indel.raw.vcf -stand_call_conf 30 -stand_emit_conf 10 -glm INDEL > $opt{project}.gatk.realign.recal.indel.log 2> $opt{project}.gatk.realign.recal.indel.log2` unless (-e "$opt{project}.realign.recal.gatk.indel.raw.vcf");
      $opt{indelvcf} ||= "$opt{project}.realign.recal.gatk.indel.raw.vcf";
      print "Done!\n";
   }elsif($opt{realign}){    ##only do realign
      print "Using realigned bam file!\n";
      `$SAMtool index $opt{project}.realign.bam` unless (-e "$opt{project}.realign.bam.bai");
      ####SNPs
      `$java -Xmx2g -jar $GATK -T UnifiedGenotyper -R $opt{ref} -I $opt{project}.realign.bam -o $opt{project}.realign.gatk.snp.raw.vcf -stand_call_conf 30 -stand_emit_conf 10 -glm SNP > $opt{project}.gatk.realign.snp.log 2> $opt{project}.gatk.realign.snp.log2` unless (-e "$opt{project}.realign.gatk.snp.raw.vcf");
      $opt{snpvcf} ||= "$opt{project}.realign.gatk.snp.raw.vcf";
      ####INDELs
      `$java -Xmx2g -jar $GATK -T UnifiedGenotyper -R $opt{ref} -I $opt{project}.realign.bam -o $opt{project}.realign.gatk.indel.raw.vcf -stand_call_conf 30 -stand_emit_conf 10 -glm INDEL > $opt{project}.gatk.realign.indel.log 2> $opt{project}.gatk.realign.indel.log2` unless (-e "$opt{project}.realign.gatk.indel.raw.vcf");
      $opt{indelvcf} ||= "$opt{project}.realign.gatk.indel.raw.vcf";
      print "Done!\n";
   }else{                    ##use original bam
      `$SAMtool index $opt{project}.bam` unless (-e "$opt{project}.bam.bai");
      ####SNPs
      `$java -Xmx2g -jar $GATK -T UnifiedGenotyper -R $opt{ref} -I $opt{bam} -o $opt{project}.gatk.snp.raw.vcf -stand_call_conf 30 -stand_emit_conf 10 -glm SNP > $opt{project}.gatk.snp.log 2> $opt{project}.gatk.snp.log2` unless (-e "$opt{project}.gatk.snp.raw.vcf");
      $opt{snpvcf} ||= "$opt{project}.gatk.snp.raw.vcf";
      ####INDELs
      `$java -Xmx2g -jar $GATK -T UnifiedGenotyper -R $opt{ref} -I $opt{bam} -o $opt{project}.gatk.indel.raw.vcf -stand_call_conf 30 -stand_emit_conf 10 -glm INDEL > $opt{project}.gatk.indel.log 2> $opt{project}.gatk.indel.log2` unless (-e "$opt{project}.gatk.indel.raw.vcf");
      $opt{indelvcf} ||= "$opt{project}.gatk.indel.raw.vcf";
      print "Done!\n";
   }
}
  
if ($opt{VQSR}){ ### VQSR apply to large sample size, small need to use hard filter
      print "Variant quality score recalibration!\n";
      ###INDEL, $opt{indelvcf}=*.indels.raw.vcf
      if ($opt{indelvcf}){
      my $temp=$1 if ($opt{vcf}=~/(.*)\.raw\.vcf$/);
      `$java -Xmx20g -jar $GATK -T VariantRecalibrator -R $opt{ref} -input $opt{indelvcf} --maxGaussians 4 --percentBadVariants 0.05 -resource:dbsnp,known=true,training=true,truth=true,prior=10.0 $opt{dbINDEL} -an QD -an FS -an MQ -an MQRankSum -an HaplotypeScore -an ReadPosRankSum --ts_filter_level 99.0 -mode INDEL -recalFile $opt{project}.VQSR.indel.recal -tranchesFile $opt{project}.VQSR.indel.tranches -rscriptFile $opt{project}.VQSR.indel.plots.R`;   
      `$java -Xmx20g -jar $GATK -T ApplyRecalibration -R $opt{ref} -input $opt{indelvcf} --ts_filter_level 99.0 -tranchesFile $opt{project}.VQSR.indel.tranches -recalFile $opt{project}.VQSR.indel.recal -mode INDEL -o $temp.VQSR.vcf`;
      }
      ###SNP, $opt{snpvcf}=*.snps.raw.vcf
      if ($opt{snpvcf}){
      my $temp=$1 if ($opt{vcf}=~/(.*)\.raw\.vcf$/);
      `$java -Xmx20g -jar $GATK -T VariantRecalibrator -R $opt{ref} -input $opt{snpvcf} --maxGaussians 4 --percentBadVariants 0.05 -resource:dbsnp,known=true,training=true,truth=true,prior=10.0 $opt{dbSNP} -an QD -an FS -an MQ -an MQRankSum -an HaplotypeScore -an ReadPosRankSum --ts_filter_level 99.0 -mode SNP -recalFile $opt{project}.VQSR.snp.recal -tranchesFile $opt{project}.VQSR.snp.tranches -rscriptFile $opt{project}.VQSR.snp.plots.R`; 
      `$java -Xmx20g -jar $GATK -T ApplyRecalibration -R $opt{ref} -input $opt{snpvcf} --ts_filter_level 99.0 -tranchesFile $opt{project}.VQSR.snp.tranches -recalFile $opt{project}.VQSR.snp.recal -mode SNP -o $temp.VQSR.vcf`;
      }
}
   
if ($opt{hardfilter}){
     print "hardfilter variant by GATK VariantFiltration!\n";
  if ($opt{indelvcf}){
     ###INDELs, $opt{indelvcf}=*.indels.raw.vcf
     my $indel=$1 if ($opt{indelvcf}=~/(.*)\.raw\.vcf$/);
     `$java -Xmx2g -jar $GATK -T VariantFiltration -R $opt{ref} --variant $opt{indelvcf} -o $indel.hardfilter.vcf \
  --filterExpression "QD < 2.0" \
  --filterName "QDFilter" \
  --filterExpression "ReadPosRankSum < -20.0" \
  --filterName "ReadPosFilter" \
  --filterExpression "FS > 200.0" \
  --filterName "FSFilter" \
  --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
  --filterName "HARD_TO_VALIDATE" \
  --filterExpression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" \
  --filterName "QualFilter"`;
     `$java -Xmx2g -jar $GATK -T SelectVariants -R $opt{ref} -variant $indel.hardfilter.vcf -o $indel.vcf --excludeFiltered`;
  }
  if ($opt{snpvcf}){
     ###SNPs, $opt{snpvcf}=*.snps.raw.vcf
     my $snp=$1 if ($opt{snpvcf}=~/(.*)\.raw\.vcf$/);
     `$java -Xmx2g -jar $GATK -T VariantFiltration -R $opt{ref} --variant $opt{snpvcf} -o $snp.hardfilter.vcf \
  --filterExpression "QD < 2.0" \
  --filterName "QDFilter" \
  --filterExpression "MQ < 40.0" \
  --filterName "MQFilter" \
  --filterExpression "FS > 60.0" \
  --filterName "FSFilter" \
  --filterExpression "HaplotypeScore > 13.0" \
  --filterName "HaplotypeScoreFilter" \
  --filterExpression "MQRankSum < -12.5" \
  --filterName "MQRankSumFilter" \
  --filterExpression "ReadPosRankSum < -8.0" \
  --filterName "ReadPosRankSumFilter" \
  --filterExpression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" \
  --filterName "StandardFilters" \
  --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
  --filterName "HARD_TO_VALIDATE"`;
     `$java -Xmx2g -jar $GATK -T SelectVariants -R $opt{ref} -variant $snp.hardfilter.vcf -o $snp.vcf --excludeFiltered`;
  }
} 










