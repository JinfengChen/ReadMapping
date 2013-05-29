#!/usr/bin/perl
=header
The scripts is designed to run bwa to map solexa sequencing read to reference genome.
--ref:   reference sequence
--1:     paired read with -2, SRR034638_1.fastq
--2:     paired read with -1, SRR034638_2.fastq
--tool:  mapping tools: bwa, maq, ssaha, soap
-project: project name that used for result file 
=cut

use Getopt::Long;
my %opt;
GetOptions(\%opt,"ref:s","1:s","2:s","tool:s","min:s","max:s","cpu:s","bam","verbose","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 -ref ../input/all.con -1 1.fq -2 2.fq -min 0 -max 500 -cpu 12 -tool soap\n";
   exit();
}

$opt{project} ||= "out";
$opt{tool} ||= "bwa";
$opt{cpu} ||=12;
$opt{min} ||= 0;
$opt{max} ||= 500; 

my $bwa="/opt/tyler/bin";
my $soap="/usr/local/bin";
my $ssaha="/home/jfchen/software/ssaha2_v2.5.5_x86_64/";
my $maq="/opt/tyler/bin/maq";
 
my $SAMtool="/usr/local/bin/samtools";
my $rmdup="/opt/picard/1.81/MarkDuplicates.jar";

if (exists $opt{1} and exists $opt{2}){
   if ($opt{tool}=~/bwa/){
      print "Run pair-end mapping by BWA!\n";
      unless (-e "$opt{ref}.sa"){
         `$bwa/bwa index $opt{ref} > $opt{project}.index.log 2> $opt{project}.index.log2`;
      }
      print "Align Read 1!\n";
      `$bwa/bwa aln -t $opt{cpu} $opt{ref} $opt{1} > $opt{1}.sai 2> $opt{1}.bwa.log2`;
      print "Align Read 2!\n";
      `$bwa/bwa aln -t $opt{cpu} $opt{ref} $opt{2} > $opt{2}.sai 2> $opt{2}.bwa.log2`;
      print "Pairing!\n";
      `$bwa/bwa sampe -a $opt{max} $opt{ref} $opt{1}.sai $opt{2}.sai $opt{1} $opt{2} > $opt{project}.sam 2> $opt{project}.sampe.log2`;
      print "SAM 2 BAM!\n";
      `$SAMtool view -bS -o $opt{project}.raw.bam $opt{project}.sam > $opt{project}.convert.log 2> $opt{project}.convert.log2`;
      print "Sort Bam!\n";
      `$SAMtool sort $opt{project}.raw.bam $opt{project}.sort > $opt{project}.sort.log 2> $opt{project}.sort.log2`;
      print "Remove duplicate!\n";
      `java -jar $rmdup ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT INPUT=$opt{project}.sort.bam OUTPUT=$opt{project}.bam METRICS_FILE=$opt{project}.dupli > $opt{project}.rmdup.log 2> $opt{project}.rmdup.log2`;
      unless ($opt{verbose}){
          `rm $opt{project}.sam $opt{project}.raw.bam $opt{project}.sort.bam`;
          `rm $opt{project}.*.log* $opt{project}.dupli $opt{1}.sai $opt{1}.bwa.log2 $opt{2}.sai $opt{2}.bwa.log2`;
      }
      print "Done!\n";
   }elsif($opt{tool}=~/soap/){ # soap
      print "Run pair-end mapping by soap!\n";
      unless (-e "$opt{ref}.index.sai"){
         `$soap/2bwt-builder $opt{ref} > $opt{project}.builder.log 2> $opt{project}.builder.log2`;
         `$SAMtool faidx $opt{ref}`; # generate $opt{ref}.fai used in samtools view
      }
      `$soap/soap -a $opt{1} -b $opt{2} -D $opt{ref}.index -o $opt{project}.soap.PE -2 $opt{project}.soap.SE -p $opt{cpu} -m $opt{min} -x $opt{max} > $opt{project}.soap.log 2> $opt{project}.soap.log2` if ($opt{max} < 2000);
      `$soap/soap -a $opt{1} -b $opt{2} -D $opt{ref}.index -o $opt{project}.soap.PE -2 $opt{project}.soap.SE -p $opt{cpu} -m $opt{min} -x $opt{max} -R > $opt{project}.soap.log 2> $opt{project}.soap.log2` if ($opt{max} >= 2000);
      `cat $opt{project}.soap.PE $opt{project}.soap.SE > $opt{project}.soap`;
      if ($opt{bam}){
      print "Convert SOAP to SAM\n";
      `perl /opt/tyler/bin/soap2sam.pl $opt{project}.soap.SE > $opt{project}.soap.SE.sam`;
      `perl /opt/tyler/bin/soap2sam.pl -p $opt{project}.soap.PE > $opt{project}.soap.PE.sam`;
      print "Convert SAM to BAM, sort and merge\n";
      `$SAMtool view -bS -t $opt{ref}.fai -o $opt{project}.raw.SE.bam $opt{project}.soap.SE.sam > $opt{project}.convert.SE.log 2> $opt{project}.convert.SE.log2`;
      `$SAMtool sort $opt{project}.raw.SE.bam $opt{project}.SE.sort > $opt{project}.SE.sort.log 2> $opt{project}.SE.sort.log2`;
      `$SAMtool view -bS -t $opt{ref}.fai -o $opt{project}.raw.PE.bam $opt{project}.soap.PE.sam > $opt{project}.convert.PE.log 2> $opt{project}.convert.PE.log2`;
      `$SAMtool sort $opt{project}.raw.PE.bam $opt{project}.PE.sort > $opt{project}.PE.sort.log 2> $opt{project}.PE.sort.log2`;
      `$SAMtool merge -f $opt{project}.sort.bam $opt{project}.SE.sort.bam $opt{project}.PE.sort.bam`;  
      print "Remove duplicate!\n";
      `java -jar $rmdup ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT INPUT=$opt{project}.sort.bam OUTPUT=$opt{project}.bam METRICS_FILE=$opt{project}.dupli > $opt{project}.rmdup.log 2> $opt{project}.rmdup.log2`;
      }
      unless ($opt{verbose}){
          `rm $opt{project}.*.sam $opt{project}.raw.*.bam $opt{project}.*.sort.bam $opt{project}.sort.bam`;
          `rm $opt{project}.*.log* $opt{project}.dupli`;
          
      }
      print "Done!\n";
   }elsif($opt{tool}=~/maq/){ # maq
      ## build reference index
      unless (-e "$opt{ref}.bfa"){
         `$maq fasta2bfa $opt{ref} $opt{ref}.bfa`;
      }
      `$maq fastq2bfq $opt{1} $opt{1}.bfq` unless (-e "$opt{1}.bfq");
      `$maq fastq2bfq $opt{2} $opt{2}.bfq` unless (-e "$opt{2}.bfq");
      `$maq match -a $opt{max} $opt{project}.Maq.map $opt{ref}.bfa $opt{1}.bfq $opt{2}.bfq`;
      `$maq mapview $opt{project}.Maq.map > $opt{project}.Maq.map.view`; 
   }
}


