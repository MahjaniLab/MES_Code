# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Queries variants in the larger SPARK cohort

#!usr/bin/perl

$GENE1 = $ARGV[0];
$GENE2 = $ARGV[1];

open(IN, "variant_metrics.txt");
while($line = <IN>){ if($line !~ /Variant/){
  if($line =~ /\_$GENE1\t/){ 
    chomp($line);
    @split = split "\t", $line;
    @split_var = split "_", $split[0];
    $polyphen{$split_var[0]} = $split[1];
    $sift{$split_var[0]} = $split[2];
    $isptv{$split_var[0]} = $split[4];
    $af{$split_var[0]} =$split[3];
    $consequence{$split_var[0]} = $split[5];
    $gene{$split_var[0]} = $split_var[1];
    @colon = split ":", $split_var[0];
    $GENE1_CHR = $colon[0];

  }
  if($line =~ /\_$GENE2\t/){
    chomp($line);
    @split = split "\t", $line;
    @split_var = split "_", $split[0];
    $polyphen{$split_var[0]} = $split[1];
    $sift{$split_var[0]} = $split[2];
    $isptv{$split_var[0]} = $split[4];
    $af{$split_var[0]} =$split[3];
    $consequence{$split_var[0]} = $split[5];
    $gene{$split_var[0]} = $split_var[1];
    @colon = split ":", $split_var[0];
    $GENE2_CHR = $colon[0];
  }
}}
close IN;

# All individuals in SPARK with WES
open(INDIV, "individuals.txt");

while($INDIV = <INDIV>){ 
  chomp($INDIV);
  $in1 = "individual/" . $GENE1_CHR . "/" . $INDIV . ".var";
  open(IN, $in1);
  while($line = <IN>){
    # chr22,17841798,C,T,0/1
    chomp($line);
    @split = split ",", $line;
    $split[0] =~ s/chr//g;
    $variant = $split[0] . ":" . $split[1] . ":" . $split[2] . ":" . $split[3];
    if($isptv{$variant} =~ /./ and $gene{$variant} =~ /^$GENE1$/){
      if($af{$variant} !~ /./){ $af{$variant} = 0;}        
      if($polyphen{$variant} =~ /null/){ $polyphen{$variant} = 0;}
      if($sift{$variant} =~ /null/){ $sift{$variant} = 1;}

      if($af{$variant} < 0.001 and $consequence{$variant} !~ /syn/){
        if($isptv{$variant} =~ /true/ or $polyphen{$variant} >= 0.5 or $sift{$variant} <= 0.05){
          print $GENE1, "\t", $INDIV, "\t", $variant, "\t", $split[4], "\t", $isptv{$variant}, "\t", 
	    $af{$variant}, "\t", $consequence{$variant}, "\t", $polyphen{$variant}, "\t", $sift{$variant}, "\t", "\n";
        }
      }
    }
  }
  close IN;

  $in2 = "individual/" . $GENE2_CHR . "/" . $INDIV . ".var";
  open(IN, $in2);
  while($line = <IN>){
    # chr22,17841798,C,T,0/1
    chomp($line);
    @split = split ",", $line;
    $split[0] =~ s/chr//g;
    $variant = $split[0] . ":" . $split[1] . ":" . $split[2] . ":" . $split[3];
    if($isptv{$variant} =~ /./ and $gene{$variant} =~ /^$GENE2$/){
      if($af{$variant} !~ /./){ $af{$variant} = 0;}
      if($polyphen{$variant} =~ /null/){ $polyphen{$variant} = 0;}
      if($sift{$variant} =~ /null/){ $sift{$variant} = 1;}

      if($af{$variant} < 0.001 and $consequence{$variant} !~ /syn/){
        if($isptv{$variant} =~ /true/ or $polyphen{$variant} >= 0.5 or $sift{$variant} <= 0.05){
          print $GENE2, "\t", $INDIV, "\t", $variant, "\t", $split[4], "\t", $isptv{$variant}, "\t",
            $af{$variant}, "\t", $consequence{$variant}, "\t", $polyphen{$variant}, "\t", $sift{$variant}, "\t", "\n";
        }
      }
    }
  }
  close IN;
}
close INDIV;

