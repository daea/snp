# snp

Data Sources
- 'AT4G34000.1.vcf' files are exports of 1001 genomes database with only variant sites for that gene
-'gene aliases' are taken from TAIR and are for quick lookup / annotation of sequences that only show AGIs
-'snps_AT4G34000' contains formatted SNPs exported from ePlants molecule viewer
- I have both Araport11 and TAIR10 .gff files that I pull from, (they're big, not on github, but are downloadable from TAIR and Araport)

Scripts
-'handle_json' is a temporary solution to my small amount of JSON/UTF8 knowledge
-'read_eplant_snps' converts the raw json export (needs to be formatted slightly) into a dictionary format that can be queried easily
-'read_1001g', reads the .vcfs and builds our sequence data (reference sequences, protein coding, SNPs, allele frequencies, genotypes etc.)
-'pull_domains' eventual location where I will pull domains from ePLANT using cURL to compare the locations of SNPs

raw_data
-data that has been formatted is pulled out of raw data and placed in the working directory (Not all on github)

Modules of Note
- PyVcf for parsing of .vcf files
- re for regular expressions
- Biopython (I haven't explicitly used it for anything yet but very soon).
- BCBio.GFF reader (kinda hard to use, I pretty much did this manually)
