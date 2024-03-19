# Moderate effect size genes project
# Author: Madison Caballero
# Desription: Created merged_table.txt, which has co-occurance and normalization stats

#!usr/bin/perl

# Table with the number of families with B* given A* inheritence
open(OUT, ">merged_table.txt");
open(GENES, "intermediates/graded.statistics.genes.txt");
while($gene = <GENES>){
  $count++; $all++;
  @gene_split = split "\t", $gene;

  open(IN, "intermediates/$gene_split[0].table");
  while($line = <IN>){
    chomp($line);
    @split = split "\t", $line;
    
    $key = $split[0] . "_" . $split[2];
    $store{$key} = $line;
  }
  close IN;
}
close GENES;

foreach $key (keys %store){
  @split = split "_", $key;
  $alt_key = $split[1] . "_" . $split[0];
  if($split[1] !~ /^$split[0]$/){
    if($done{$alt_key} !~ /./ and $done{$key} !~ /done/ and $store{$key} =~ /./ and $store{$alt_key} =~ /./){
       $done{$alt_key} = "done"; $done{$key} = "done";
       print OUT $store{$key}, "\t", $store{$alt_key}, "\n";
    }
  }
}
close OUT;


# awk '($12 < 1) && ($6 < 1)' ~/MES_code/1_MES_screen/merged_table.txt > ~/MES_code/1_MES_screen/merged_table.candidates_only.txt

