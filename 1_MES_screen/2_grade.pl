# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Finds A* inheritence patterns on a per-gene basis

#!usr/bin/perl

# A file with variant scores and stats from SPARK
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

# Families without the top ASD genes impacted
open(FAM, "edited_families.txt");

$count = 0;
$grade = "FAIL";
open(OUT, ">intermediates/1_valid_probands.passed.txt");
open(FAIL, ">intermediates/1_valid_probands.failed.txt");

while($fam = <FAM>){ if($fam !~ /family/){
  chomp($fam);
  @split_fam = split "\t", $fam;
  $in = "intermediates/" . $split_fam[0] . "_valid_probands.txt";

  open(IN, $in);
  while($line = <IN>){
    chomp($line);
    $count++;
    if($count == 1){
      $line =~ s/\t$//;
      chomp($line);
      print OUT $line, "\tgrade\n";
      print FAIL $line, "\tgrade\n";
    }
    @split = split "\t", $line;

    if($count > 1 and $line !~ /Consequence/){
      $grade = "FAIL";
      # All probands have the variant
      if($split[14] == $split[10]){ 
        
        # Allele is inherited:
        if(($split[6] =~ /0.0/ and $split[7] =~ /0.1/) or ($split[7] =~ /0.0/ and $split[6] =~ /0.1/)){
        if($af{$split[0]} !~ /./){ $af{$split[0]} = 0;}

          # Rare variant, not synonymous, and exists
          if($af{$split[0]} < 0.001 and $split[2] !~ /syn/ and $isptv{$split[0]} =~ /./){
            # Is a PTV or presumed deleterious  
            if($polyphen{$split[0]} =~ /null/){ $polyphen{$split[0]} = 0;}
            if($sift{$split[0]} =~ /null/){ $sift{$split[0]} = 1;}

            if($isptv{$split[0]} =~ /true/ or $polyphen{$split[0]} >= 0.5 or $sift{$split[0]} <= 0.05){
                $grade = "PASS";
                $gene_save{$split[1]}++;
                $pass++; 
                print OUT $line, "\t", $grade, "\n";
              }
          }
        }
      } 
      if($grade =~ /FAIL/){
         print FAIL $line, "\t", $grade, "\n";
      }
    }
  }
  close IN;
}}
close OUT;

open(OUT, ">intermediates/graded.statistics.txt");
for $gene (keys %gene_save){
  if($gene_save{$gene} >= 3){
    $examples++;
  }
}

print OUT "Number of passed examples: ", $pass, "\n";
print OUT "Number of genes with 3+ examples: ", $examples, "\n";

close OUT;

open(OUT, ">intermediates/graded.statistics.genes.txt");

for $gene (keys %gene_save){
  if($gene_save{$gene} >= 3){
    print OUT $gene, "\t", $gene_save{$gene}, "\n";
  }
}

close OUT;
