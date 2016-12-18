# Matt Cumming
# Control File (for a little more automation)

import os, sys, glob
import read_eplant_snps

snp_ePlant_files = glob.glob('*_raw.txt')
snp_ebi_files = glob.glob('*.vcf') 
snp_1001_files = glob.glob('*_1001.vcf')

def make_snp_classes():
	'code to create SNP classes with non-syn snps here'



'''
This is just to automate the pulling of SNPS from ePlant molecule viewer files.
'''


