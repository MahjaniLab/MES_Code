# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Filters the vcf for the genes of interest

open(GENES, "gene_names_Snpeffcompat.txt");
#open(GENES, "missing.txt");
while($line = <GENES>){
  chomp($line);
  $genes{$line}++;
  $out = "counts/" . $line . ".counts";
  open(OUT, ">$out");
  close OUT;
}
close GENES;

# Using the annotation information stored in the single sample file
# Also filters AF based on just samples in this cohort
  open(SINGLE, "singlesamp.ann.dbNSFP.deleterious.vcf");
  while($line = <SINGLE>){
    @split = split /\|/, $line;
    for($a = 1; $a < scalar @split; $a++){
      if($genes{$split[$a]} =~ /./){
        if($split[$a - 1] =~ /HIGH|MODERATE/){
          @tab_split = split "\t", $line;
          @comma_split = split ";", $tab_split[7];
          for $com (@comma_split){
            if($com =~ /^AF\=/){
              $com =~ s/AF=//g;
              if($com < 0.001){
                $key = $tab_split[0] . "_" . $tab_split[1] . "_" . $tab_split[3] . "_" . $tab_split[4];
                $key_to_gene{$key} = $split[$a];
              }
            } 
          }
        }
      }
    }
  close SINGLE;
}

# Parses the filtered vcf of all samples and extracts counts 
open(IN, "filtered.input.vcf");
while($line = <IN>){
  if($line =~ /^\#CHROM/){
    @split = split "\t", $line;
    for($a = 9; $a < scalar @split; $a++){$IDS{$a} =  $split[$a];}
  }
  if($line !~ /^\#/){
    @split = split "\t", $line;
    $key = $split[0] . "_" . $split[1] . "_" . $split[3] . "_" . $split[4];
    if($key_to_gene{$key} =~ /./){
      $out = "counts/" . $key_to_gene{$key} . ".counts";
      open(OUT, ">>$out");
      for($a = 9; $a < scalar @split; $a++){
        if($split[$a] =~ /^0\/1/){
          print OUT $key, "\t", $IDS{$a}, "\t", $split[$a], "\n";
        }
      }
      close OUT;
    }
  }
}
