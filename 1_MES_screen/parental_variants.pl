# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Calculates rate of ultra-rare delterious variants in SPARK parents per gene

#!usr/bin/perl
open(IN, "variant_metrics.txt");
while($line = <IN>){ if($line !~ /Variant/){
  chomp($line);
  @split = split "\t", $line;
  @split_var = split "_", $split[0];
  $polyphen{$split_var[0]} = $split[1];
  $sift{$split_var[0]} = $split[2];
  $isptv{$split_var[0]} = $split[4];
  $af{$split_var[0]} =$split[3];
}}
close IN;

open(FAM, "families.txt");

$count = 0;
$grade = "FAIL";
open(OUT, ">parent_variants.txt");

while($fam = <FAM>){ if($fam !~ /family/){
  chomp($fam);
  @split_fam = split "\t", $fam;
  $in = "intermediates/" . $split_fam[0] . "_valid_probands.fixed.txt";

  open(IN, $in);
  while($line = <IN>){
    chomp($line);
    @split = split "\t", $line;
    
    if($split[6] =~ /0.1/){
      if($af{$split[0]} < 0.01 and $split[2] !~ /syn/ and $isptv{$split[0]} =~ /./){
  	if($polyphen{$split[0]} =~ /null/){ $polyphen{$split[0]} = 0;}
        if($sift{$split[0]} =~ /null/){ $sift{$split[0]} = 1;}
        if($isptv{$split[0]} =~ /true/ or $polyphen{$split[0]} >= 0.5 or $sift{$split[0]} <= 0.05){
          print OUT $split[0], "\t", $split[1], "\t", $split[4], "\n";
        }
      }
    }
    if($split[7] =~ /0.1/){
      if($af{$split[0]} < 0.01 and $split[2] !~ /syn/ and $isptv{$split[0]} =~ /./){
        if($polyphen{$split[0]} =~ /null/){ $polyphen{$split[0]} = 0;}
        if($sift{$split[0]} =~ /null/){ $sift{$split[0]} = 1;}
        if($isptv{$split[0]} =~ /true/ or $polyphen{$split[0]} >= 0.5 or $sift{$split[0]} <= 0.05){
          print OUT $split[0], "\t", $split[1], "\t", $split[5], "\n";
        }
      }
    }
  }
  close IN;
}}
close OUT;


open(IN, "parent_variants.txt");
open(OUT, ">parental_gene_mut_stats.txt");
while($line = <IN>){
  @split = split " ", $line;
  $gene{$split[1]}++;
  $samp{$split[2]}++;
}
close IN;

foreach $g (keys %gene){
  print OUT $g, "\t", $gene{$g}/scalar keys %samp, "\n";
}
close OUT;

