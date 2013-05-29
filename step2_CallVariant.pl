#!/usr/bin/perl
=header
The scripts is designed to call SNP/INDEL with SAMtools or GATK.
-ref:     reference sequence
-bam:     sorted bam alignment file of mapped reads after rmdup
-realign  Realigning the BAM file using GATK's RealignerTargetCreator and IndelRealigner (optional)
-recal    Recalibration (optional, not added)
-samtools call variant using samtools (optional, default using GATK)
-project: project name that used for result file (optional)
=cut

use Getopt::Long;
my %opt;
GetOptions(\%opt,"ref:s","bam:s","samtools","realign","recal","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 -ref ../input/all.con -bam out.bam\n";
   exit();
}

$opt{project} ||= "HEG4_BWA";
my $bwa="/opt/tyler/bin/bwa";
my $SAMtool="/opt/samtools/0.1.18/bin/samtools";
my $BCFtool="/opt/samtools/0.1.18/bin/bcftools";
my $vcfutils="/opt/samtools/0.1.18/bin/vcfutils.pl";
my $GATK="/opt/GATK/2.4-3-g2a7af43/GenomeAnalysisTK.jar";
my $java="/opt/java/jdk1.6.0_38/bin/java";


### realign
if ($opt{realign}){
   print "Realigning the BAM file using GATK's RealignerTargetCreator and IndelRealigner\n";
   `$java -jar $GATK -T RealignerTargetCreator -R $opt{ref} -I $opt{bam} -o $opt{project}.gatk.intervals > $opt{project}.gatk.interval.log 2> $opt{project}.gatk.interval.log2`;
   `$java -jar $GATK -T IndelRealigner -R $opt{ref} -I $opt{bam} -targetIntervals $opt{project}.gatk.intervals -o $opt{project}.realigned.bam > $opt{project}.gatk.realign.log 2> $opt{project}.gatk.realign.log2`; 
}
### recalibration
if ($opt{recal}){
   print "Recalibration\n";
   "not added"
}

### call variants
if ($opt{samtools}){
   print "Call variant using SAMtools/BCFtools!\n";
   unless (-e "$opt{ref}.fai"){
      `samtools faidx $opt{ref}`;
   }
   unless ($opt{realign}){
      my $cmd1="$SAMtool mpileup -uf $opt{ref} $opt{bam} | $BCFtool view -bvcg - > $opt{project}.raw.bcf 2> $opt{project}.raw.bcf.log2";
      my $cmd2="$BCFtool view $opt{project}.raw.bcf | $vcfutils varFilter -D100 > $opt{project}.flt.vcf 2> $opt{project}.flt.vcf.log2";
      `$cmd1`;
      `$cmd2`;
      print "$cmd1\n$cmd2\n";
      print "Done!\n";
   }else{
      print "Using realigned bam file!\n";
      `$SAMtool mpileup -uf $opt{ref} $opt{project}.realigned.bam | $BCFtool view -bvcg - > $opt{project}.realigned.raw.bcf 2> $opt{project}.realigned.raw.bcf.log2`;
      `$BCFtool view $opt{project}.realigned.raw.bcf | $vcfutils varFilter -D100 > $opt{project}.realigned.flt.vcf 2> $opt{project}.realigned.flt.vcf.log2`;
      print "Done!\n";
   }
}else{ ### GATK
   print "Call variant using GATK!\n";
   unless ($opt{realign}){
      `$SAMtool index $opt{bam}` unless (-e "$opt{bam}.bai");
      `$java -jar $GATK -T UnifiedGenotyper -R $opt{ref} -I $opt{bam} -o $opt{project}.gatk.calls  -stand_call_conf 30 -stand_emit_conf 10 > $opt{project}.gatk.call.log 2> $opt{project}.gatk.call.log2`; 
      print "Done!\n";
   }else{
      print "Using realigned bam file!\n";
      `$SAMtool index $opt{project}.realigned.bam` unless (-e "$opt{project}.realigned.bam.bai");
      `$java -jar $GATK -T UnifiedGenotyper -R $opt{ref} -I $opt{project}.realigned.bam -o $opt{project}.realigned.gatk.calls  -stand_call_conf 30 -stand_emit_conf 10 > $opt{project}.gatk.realigned.call.log 2> $opt{project}.gatk.realigned.call.log2`;
      print "Done!\n";
   }
}

