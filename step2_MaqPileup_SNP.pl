#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"ref:s","maq:s","project:s","help");


my $help=<<USAGE;
perl $0 --ref --maq 
Give a map result of maq, pileup and output a SNP file, which same with parents used in MPR.
SNPid refbase	SNP
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}


$opt{project} ||= "test";
my $maq="/opt/tyler/bin/maq";
unless (-e "$opt{ref}.bfa"){
   `$maq fasta2bfa $opt{ref} $opt{ref}.bfa`;
}
`$maq pileup -vP -q 40 $opt{ref}.bfa $opt{maq} | awk '\$4 > 0' > $opt{maq}.pileup` unless (-e "$opt{maq}.pileup");
parsepileup("$opt{maq}.pileup");

##chromosome05    94      T       5       @,.c.,  @H;G@H  @~~~~~  63,45,46,48,49,
sub parsepileup 
{
my ($file)=@_;
my $maxdepth=60; #max depth of one SNP to avoid repetitive sequence regions
my $minbaseq=20; #min base quality for at least one base of one allele in SNP site
my $minreads=4;  #min of reads support a allele in SNP site
my $minsumbq=80; #min of sum base quality for each allele in SNP site
my ($snpID);
open OUT, ">$opt{project}.SNPs.parents" or die "$!";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    next if ($unit[3] >= $maxdepth and $unit[3] == 0); ## possible repetitive regions or uncovered
    $unit[0]=~s/\D+//ig;
    $unit[1]=sprintf("%08d",$unit[1]);
    $snpID=$unit[0].$unit[1].$unit[2];
    my @base=split("",$unit[4]);
    my @qual=split("",$unit[5]);
    my %allele;
    for(my $i=1;$i<@base;$i++){
       if ($base[$i]=~/\,/ or $base[$i]=~/\./){
          my $qscore=ord($qual[$i])-33;
          push @{$allele{$unit[2]}}, $qscore;
       }else{
          $base[$i]=~tr/atcg/ATCG/;
          my $qscore=ord($qual[$i])-33;
          push @{$allele{$base[$i]}}, $qscore;
       }
    }
    if (keys %allele == 1){  ##only have one alleles in SNP site
       my $SNP=join("\t",keys %allele);
       #print OUT "$snpID\t$SNP\tcheck\n";
       my %maxbq; my %sumbq; my %maxreads;
       foreach my $a (keys %allele){
   
          $maxreads{$a}=@{$allele{$a}};
          $maxbq{$a}=max(\@{$allele{$a}});
          $sumbq{$a}=sum(\@{$allele{$a}});
          #print OUT "$a\t$maxreads{$a}\t$maxbq{$a}\t$sumbq{$a}\n";
       } 
       my $flag1; my $flag2; my $flag3;
       foreach my $v (values %maxreads){
          $flag1++ if $v >= $minreads; ## reads depth larger than $minreads==4 for each allele
       }
       foreach my $q (values %maxbq){ 
          $flag2++ if $q >= $minbaseq; ## at least have a base quality larger than $minbaseq==20 for each allele  
       }
       foreach my $q (values %sumbq){
          $flag3++ if $q >= $minsumbq; ## sum of base quality larger than $minsumbq==60 for each allele
       }
       #print OUT "$flag1\t$flag2\t$flag3\n";
       if ($flag1 == 1 and $flag2 == 1 and $flag3 ==1){ ## need to meet all three criteria for both allele
          print OUT "$snpID\t$unit[2]\t$SNP\n";
       }
    }
}
close IN;
close OUT;
}


sub max
{
my ($num)=@_;
my $max=0;
foreach  (@$num) {
        next if ($_ eq "NA");
        $max = $_ > $max ? $_ : $max;
}
return $max;
}


sub sum
{
my ($num)=@_;
my $loop=0;
my $total;
foreach (@$num) {
        next if ($_ eq "NA");
        $total+=$_;
        $loop++;
}
return 0 if ($loop == 0);
return $total;
}

