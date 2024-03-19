# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Filters the vcf for Polyphen, Sift, and gnomad AF

use List::Util qw( min max );
open(IN, $ARGV[0]);
while($line = <IN>){
  # Site cannot be functional in another gene
  if($line !~ /HIGH|MODERATE/){
    chomp($line);
    @colon = split ";", $line;
    $AF = 0;
    if($line =~ /gnomAD/){
      foreach $piece (@colon){
        if($piece =~ /gnomAD_/){
          @equal = split "=", $piece;
          $equal[1] =~ s/\s.+$//g;
          $AF = $equal[1];
        }
      }
    }
    if($AF < 0.001){
      print $line,  "\n";
    }
  }
}
close IN;
