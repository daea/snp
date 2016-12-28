# Matt Cumming 2016 #
# Provart Lab #

import csv, re, glob, sys, os, vcf, textwrap, dill
from Bio.Seq import MutableSeq, Seq as MutableSeq, Seq
import pprint as pp
import process_vcfs as pvcf

''' For a gene with a given .vcf, and a reference .gff, this will create.
- full length CDS
- a list of potential SNPs
- first site
- last site
- AGI ID
- protein coding region (exons spliced out based on gff annotation)
- dictionary of snpsites
-dictionary containing affected codons for a particular site, and the amino acid translations (non-syn)
''' 


# Helper Functions

def mutate(sequence, site, mutation):
	try:
		old = list(sequence)
	except TypeError:	
		print "Input sequence should be a string or list"
	old[site-1] = mutation
	return old	



def codon(sequence, frame=0):
	return textwrap.wrap(sequence[frame:], 3)	




# Snp Class

class Snp:
	def __init__(self, vcffile, gff):
		self.AGI = os.path.basename(vcffile)[:-4]
		self.fullsites, self.snpsites, self.first_site, self.last_site, self.ref_seq = pvcf.load_pickle(self.AGI, os.path.dirname(vcffile) + '/')
		self.find_exons(self.ref_seq, gff)
		self.polymorphisms()

	def find_exons(self, DNA_sequence, gff):
		''' Pull annotation from GFF files. For a given AGI_ID '''
		concat = MutableSeq('')
		self.exons = {}
		
		with open(gff) as tsv:
			for line in csv.reader(tsv): 
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
				concat = concat + self.exons[i]
			self.protein_coding = concat
				
	def print_polymorphs(self):
		print "{0} has the following polymorphic sites: \n".format(self.AGI)
		for i in self.polymorphs.keys():
			print "Site: {0}, Reference: {1}, Polymorphism: {2}, Non Synonymous?: {3}, Residue: {4}, Mutant: {5}".format(
					self.polymorphs[i]["site"], 
					self.polymorphs[i]["reference"],
					self.polymorphs[i]["alt"], 
					self.polymorphs[i]["non_syn"],	
					self.polymorphs[i]["old_residue"], 
					self.polymorphs[i]["new_residue"]
			) 
	
	def polymorphisms(self):
		self.polymorphs = {}
		codons = codon(str(self.protein_coding))		
		for i in self.snpsites:
			try:
				site = i[0] - self.first_site + 1
				codon_site = site % 3,
				old_res = MutableSeq(codons[site/3]).translate()
				new_res = MutableSeq("".join(mutate(codons[site/3], site%3, i[2]))).translate() 
				self.polymorphs[i[0]] = {
						"site" : i[0],
						"reference" : i[1],
						"alt" : i[2],
						"freq": i[3],
						"old_residue" : old_res,
						"new_residue" : new_res,
						"non_syn" : old_res == new_res
					}
			except IndexError:
				continue

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

###############################################################################
	
gff = 'Araport11_GFF3_genes_transposons.201606.gff.txt'
directory = "./1001g_variants/"
file_list, pickled = pvcf.get_files(directory)

