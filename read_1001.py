# Matt Cumming 2016 #
# Provart Lab #

import csv, re
import vcf
from Bio.Seq import MutableSeq, Seq as MutableSeq, Seq
import pprint as pp
import textwrap

''' For a gene with a given .vcf, and a reference .gff, this will create.
- full length CDS
- a list of potential SNPs
- first site
- last site
- AGI ID
- protein coding region (exons spliced out based on gff annotation)
- dictionary of snps
''' 

def mutate(sequence, site, mutation):
	try:
		old = list(sequence)
	except TypeError:	
		print "Input sequence should be a string or list"
	old[site-1] = mutation
	return old	

def codon(sequence, frame=0):
	return textwrap.wrap(sequence[frame:], 3)	



class Snp:
	def __init__(self, vcffile, gff):
		self.AGI = vcffile[:-4]
		print self.AGI	
		self.read_vcf(vcffile)
		self.find_exons(self.ref_seq, gff)

	def read_vcf(self, vcffile):
		''' Gets the sequence from the vcf files and grabs the associated SNPs'''
		self.fullsites = []
		self.snpsites = [] 
		self.first_site = 0
		self.last_site = 0
		gene = vcf.Reader(open(vcffile, 'r'))
		
		for record in gene:
			if record.POS > self.last_site:	#	Find Last
				self.last_site = record.POS
			if self.first_site == 0:
				self.first_site = record.POS	#	Find First	
		
			self.fullsites.append(record.REF[0]) #	Make Reference Sequence
		#	if record.REF not in ['A','T','C', 'G']:
		#		print 'Site: {0}, Base: {1}'.format(record.POS, record.REF)
			if record.is_snp == True:
				self.snpsites.append(
					[record.POS, record.REF, record.ALT, record.aaf]
				) #	Make SNP Sequence
				 	
		self.ref_seq = str(''.join(self.fullsites))	
		#print "Length of self.ref_seq: {0}".format(len(self.ref_seq))
	
	def find_exons(self, DNA_sequence, gff):
		''' Pull annotation from GFF files. For a given AGI_ID '''
		concat = MutableSeq('')
		self.exons = {}
		
		with open(gff) as tsv:
			for line in csv.reader(tsv, delimiter="\t"): 
				if self.AGI in line[8] and line[2] == 'CDS':
					exon = re.search(r'CDS:(\d+)', line[8])
					try:
						self.exons[int(exon.group(1))] = MutableSeq(
							self.ref_seq[int(line[3]) - 
							self.first_site: int(line[4]) - 
							self.first_site + 1]
						)
					except IndexError:
						continue			

			for i in self.exons.keys():
		#		print self.exons[i],":", i ,"\n", len(self.exons[i])
				concat = concat + self.exons[i]
			self.protein_coding = concat
				
	def polymorphisms(self):
		self.polymorphs = {}
		codons = codon(str(self.protein_coding))		
		for i in self.snpsites:
			try:
				print ("S: {0}, B: {2}, C: {3}, CB: {4}, R: {5}".format( 
					i[0],
					self.first_site,
					i[0] - self.first_site + 1,
					((i[0] - self.first_site + 1) / 3) + 1,  
					(i[0] - self.first_site + 1) % 3,
					MutableSeq(codons[(i[0] - self.first_site + 1)/3]).translate()
					)
				)
			except IndexError:
				continue

			self.polymorphs = {
				i[0] : {
					"site" : i[0],
					"reference" : i[1],
					"alt" : i[2],
					"freq": i[3]
				}
			}				
		
				
	#Access Functions
	def get_DNAseq(self):
		return self.ref_seq[:]

	def get_protein(self):
		return self.protein_seq[:]

	def get_CDSstart(self):
		return self.first_site[:]

	def get_CDSend(self):
		return self.last_site[:]

	def get_snpsites(self):
		return self.snpsites[:]

	def get_protein_coding(self):
		return self.protein_coding[:]

	def get_protein_seq(self):
		return self.get_protein_coding().translate(cds = True)	


gff = 'Araport11_GFF3_genes_transposons.201606.gff.txt'
vcf_file = 'AT4G34000.1.vcf'
test = Snp(vcf_file, gff)
sites = test.polymorphisms()



#test_seq = 'ctttttcactttcttgcttgacaggcttatacgatggaactggaagcagaaattgcgcaactcaaagaattgaatgaagagttgcagaagaaacaagtgtgtctcgcttcttccctatcacaattaagaatctcgagattttcatattttcttgaggttgtattcactgaccaaatgtttcatgcaggttgaaatcatggaaaagcagaaaaatcaggtactgtcttgatttgaatatcctctatggttgttggctaggcttttaactctcactcataatgaattacacttttggacagtattctaagcttttgagtagaatagtgtaagctataccatgaagtgagacatatcatcacatttttgatttcccactctgcataaagtatttaagatttgtgaatatgttgcaatgccaatttggatatttcatgagactaatctgacgagcatggatttaacggaagttggctcatttgttgttgcagcttctggagcctctgcgccagccatggggaatgggatgcaaaaggcaatgcttgcgaaggacattgacgggtccctggtagagcttataatggcgtctaaggaacccaacaaagcgccgaagttatagaacaactcagaagatagaaagctagctttgtacgtagtttaggcaggttctgtgggtgattgtaaatcttgaagtgtggcggatttgacagagatagataaacacatatctgttctattttcctaaatcttttggttttatcttcctgatgtaatggatctttatcatttgtcttgaacatctttgtgacttaaccagagtgaatttatcttgtatcttgtctgcaattttttctttagctcatttgggcccaagtcaacacctttaaa'

#mutate(test_seq, 1, 'A')
#print codon(test_seq, 1)



################################################################################
#Eventual Control Functions

def make_genelist(filelist):
	'''Create a list of Snp objects, from a list of vcf files'''
	global gene_list
	gene_list = {}
	for file in filelist:
		gene = vcf.Reader(open(file, 'r'))
		gene_list = {file[:-4] : Snp(gene)}
						

