# Moderate effect size genes project
# Author: Madison Caballero
# Desription: A set of steps to annotate and extract variants in the non-neuro cohorts.
# These steps were used for 1kGP, BioMe, and All of Us.
# Similar to that for deleterious variants

# Single sample extraction, keeps all loci
bcftools view --threads 4 -s [random_sample_id] input.vcf.gz > singlesamp.vcf

# For speed, I just extract one sample. This keeps all loci even if 0/0
bcftools view --threads 4 -s [random_sample_id] input.vcf.gz > singlesamp.vcf

# SNPeff annotate for function
java -Xmx8g -jar snpEff.jar -v GRCh38.86 -canon singlesamp.vcf > singlesamp.ann.vcf

# SNPsift to get dbNSFP information
java -jar SnpSift.jar dbnsfp -v -db dbNSFP4.1a.txt.gz \
-f gnomAD_exomes_AF,SIFT_score,SIFT4G_score,Polyphen2_HDIV_score,Polyphen2_HVAR_score,genename,Ensembl_geneid \
singlesamp.ann.vcf > singlesamp.ann.dbNSFP.vcf

# Get gene IDs that are SNPeff compatible
awk '{print "synonymous_variant|LOW|" $1}' genes.compat.txt > grab.txt

# Grab variants that are syn for those genes
# There is a later filtering step in case grep grabs partial matches
grep -f grab.txt singlesamp.ann.dbNSFP.vcf > singlesamp.ann.dbNSFP.functional.noheader.vcf
perl filter.pl singlesamp.ann.dbNSFP.functional.noheader.vcf > singlesamp.ann.dbNSFP.synonymous.noheader.vcf

# header check
grep -P "^" singlesamp.regen.1.ann.dbNSFP.vcf > header.txt
cat header.txt singlesamp.ann.dbNSFP.synonymous.noheader.vcf > \
	singlesamp.ann.dbNSFP.synonymous.vcf

# Grab regions
    bcftools view \
    -T singlesamp.ann.dbNSFP.synonymous.vcf \
    -O v --threads 8 \
    -o filtered.input.vcf \
    input.vcf.gz 

# Get counts
perl grab_counts.pl

