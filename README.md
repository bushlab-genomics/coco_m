
# Example of commands:
```
#directory to latest cocos_m:
path_to_cocos_m=/storage/cocos_m
#directory to vep:
path_to_vep=/storage/software/ensembl-tools-release-80
${path_to_vep}/scripts/variant_effect_predictor/variant_effect_predictor.pl  -i test.vcf --cache --dir_cache ${path_to_vep}/cache --offline --force_overwrite --dir_plugins path_to_cocos_m --plugin cocos_m -o test.out --stats_text
```
**when using vep version 80, there is a 'negative count' warning. It won't affect the function. This warning has been further removed since version 87.**

# Additional output:

**IMPACT : impact of variant** 
- high severity if there is a PTC predicted with canonical rules    
- mild severity if there is a PTC predicted with noncanonical rules

**PTC_POSITION :** the new premature codon position in cDNA with variant inserted

**TERM_TYPE :** termination type of sequence    
- PTC,ncanonical,escaping NMD    
- PTC,noncanonical,escaping NMD

**ALT_AA_SEQ :** altered amino acid sequence

# Rationale for COCOS:

1. record exon junction positions
   - get all original exon junction positions
   - insert variant
   - get all modified exon junction positions
2. record truncated cDNA with no 5'UTR but 3'UTR 
   - scan the new cDNA and find the new PTC position
3. determine the consequence:
   - if new PTC position doesn't exist: non-stop decayif new PTC position is 55nt early than last exon-exon junction: NMD
   - if new PTC position is early than last exon-exon junction but within 55 nt: escape NMD with canonical rules            
   - if there is very long exons (greater than approximately 407nt) inhibit NMD: escape NMD with noncanonical rules 
   - if new PTC position is within 150nt from the start codon: escape NMD with noncanonical rules 

