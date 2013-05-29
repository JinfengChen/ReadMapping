#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"ref:s","project:s","help");


my $help=<<USAGE;
perl $0 --ref
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{project} ||= "map";

my $reflib=readlib("in_libs.HEG4_RAW.csv");
my $refgroup=readgroup("in_groups.HEG4_RAW.csv");
my $script="/rhome/cjinfeng/HEG4_cjinfeng/MappingReads/bin/step1_Mapping.pl";

open OUT, ">$opt{project}.sh" or die "$!";
foreach my $read (sort keys %$refgroup){
   #print "$read\n";
   my $fq1=$read;
   my $fq2=$read;
   $fq1=~s/\?/1/;
   $fq2=~s/\?/2/;
   my $head;
   if ($read=~ "HEG4.*\/(.*?)\.fq" ){
      $head=$1;
   }
   $head=~s/\_p*\?//;
   print "$head\n";
   #print "$fq1\n$fq2\n";
   #print "$refgroup->{$read}\t$reflib->{$refgroup->{$read}}->[0]\t$reflib->{$refgroup->{$read}}->[1]\n";
   my $min=$reflib->{$refgroup->{$read}}->[0]-$reflib->{$refgroup->{$read}}->[1];
   my $max=$reflib->{$refgroup->{$read}}->[0]+$reflib->{$refgroup->{$read}}->[1];
   my $cmd="perl $script -ref $opt{ref} -1 $fq1 -2 $fq2 -min $min -max $max -cpu 30 -tool bwa -project ./$head\.$opt{project}";
   print OUT "$cmd\n";
}
close OUT;

#####
#illuminaGAIIx_200_2.3,HEG4_RAW,HEG4,fragment,1,148,25,,,inward,0,0
#illuminaGAIIx_500,HEG4_RAW,HEG4,jumping,1,,,433,27,inward,0,0
sub readlib
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split(",",$_);
    $hash{$unit[0]}=$unit[3]=~/fragment/ ? [$unit[5],$unit[6]] : [$unit[7],$unit[8]];
}
close IN;
return \%hash;
}
 

#jump_500_P1,illuminaGAIIx_500,/rhome/cjinfeng/HEG4_cjinfeng/fastq/errorcorrection/soapec/HEG4_0_500bp/FC52_7_?.fq
sub readgroup
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split(",",$_);
    $hash{$unit[2]}=$unit[1];
}
close IN;
return \%hash;
}
 
