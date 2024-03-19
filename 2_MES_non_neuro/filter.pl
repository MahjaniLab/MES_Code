# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Filters the vcf for Polyphen, Sift, and gnomad AF

use List::Util qw( min max );
open(IN, $ARGV[0]);
while($line = <IN>){
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
    undef($Poly);
    undef($Sift);
    if($line =~ /Polyphen2_HVAR_score/){
      foreach $piece (@colon){
        if($piece =~ /Polyphen2_HVAR_score/){
          $piece =~ s/dbNSFP_Polyphen2_HVAR_score=//g;
          @comma = split ",", $piece;
          $Poly = max @comma;
        }
      }
    }
   if($line =~ /SIFT4G/){
      foreach $piece (@colon){
        if($piece =~ /SIFT4G/){
          $piece =~ s/dbNSFP_SIFT4G_score=//g;
          @comma = split ",", $piece;
          $Sift = min @comma;
        }
      }
    }
    if($line =~ /HIGH/){$Sift = 0;}
    if($Poly >= 0.5 or $Sift <= 0.05){
      print $line,  "\n";
    }
  }
}
close IN;
