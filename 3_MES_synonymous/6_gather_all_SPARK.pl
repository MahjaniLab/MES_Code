# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Similar to the code for delterious variants but for synonymous variants.

#!usr/bin/perl

$GENE1 = $ARGV[0];
$GENE2 = $ARGV[1];

open(IN, "variant_metrics.txt");
while($line = <IN>){ if($line !~ /Variant/){
  if($line =~ /\_$GENE1\t/){ 
    chomp($line);
    @split = split "\t", $line;
    @split_var = split "_", $split[0];
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
    $af{$split_var[0]} =$split[3];
    $consequence{$split_var[0]} = $split[5];
    $gene{$split_var[0]} = $split_var[1];
    @colon = split ":", $split_var[0];
    $GENE2_CHR = $colon[0];
  }
}}
close IN;

## OVERWRITE IF VAR IS BOTH MIS AND SYN

open(IN, "variant_metrics.txt");
while($line = <IN>){ if($line !~ /Variant/){
  if($line =~ /\_$GENE1\t/ and $line !~ /syn/){
    chomp($line);
    @split = split "\t", $line;
    @split_var = split "_", $split[0];
    $af{$split_var[0]} =1;
    $consequence{$split_var[0]} = $split[5];
    $gene{$split_var[0]} = $split_var[1];
    @colon = split ":", $split_var[0];
    $GENE1_CHR = $colon[0];
  }
  if($line =~ /\_$GENE2\t/ and $line !~ /syn/){
    chomp($line);
    @split = split "\t", $line;
    @split_var = split "_", $split[0];
    $af{$split_var[0]} =1;
    $consequence{$split_var[0]} = $split[5];
    $gene{$split_var[0]} = $split_var[1];
    @colon = split ":", $split_var[0];
    $GENE2_CHR = $colon[0];
  }
}}
close IN;

open(INDIV, "individuals.txt");

while($INDIV = <INDIV>){ 
  chomp($INDIV);
  $in1 = "individual/" . $GENE1_CHR . "/" . $INDIV . ".var";
  open(IN, $in1);
  while($line = <IN>){
    chomp($line);
    @split = split ",", $line;
    $split[0] =~ s/chr//g;
    $variant = $split[0] . ":" . $split[1] . ":" . $split[2] . ":" . $split[3];
    if($gene{$variant} =~ /^$GENE1$/){
      if($af{$variant} !~ /./){ $af{$variant} = 0;}        
      if($af{$variant} < 0.001 and $consequence{$variant} =~ /syn/){
          print $GENE1, "\t", $INDIV, "\t", $variant, "\t", $split[4], "\t", $isptv{$variant}, "\t", 
	    $af{$variant}, "\t", $consequence{$variant}, "\t", $polyphen{$variant}, "\t", $sift{$variant}, "\t", "\n";
      }
    }
  }
  close IN;

  $in2 = "individual/" . $GENE2_CHR . "/" . $INDIV . ".var";
  open(IN, $in2);
  while($line = <IN>){
    chomp($line);
    @split = split ",", $line;
    $split[0] =~ s/chr//g;
    $variant = $split[0] . ":" . $split[1] . ":" . $split[2] . ":" . $split[3];
    if($gene{$variant} =~ /^$GENE2$/){
      if($af{$variant} !~ /./){ $af{$variant} = 0;}
      if($af{$variant} < 0.001 and $consequence{$variant} =~ /syn/){
          print $GENE2, "\t", $INDIV, "\t", $variant, "\t", $split[4], "\t", $isptv{$variant}, "\t",
            $af{$variant}, "\t", $consequence{$variant}, "\t", $polyphen{$variant}, "\t", $sift{$variant}, "\t", "\n";
      }
    }
  }
  close IN;
}
close INDIV;

