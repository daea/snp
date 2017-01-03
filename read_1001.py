# Matt Cumming 2016 #
# Provart Lab #

import csv, re, glob, sys, os, vcf, textwrap, dill
from Bio.Seq import MutableSeq, Seq as MutableSeq, Seq
import process_vcfs as pvcf
import process_gff as pgff
from collections import namedtuple


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
		self.vcffile = vcffile
		self.gfffile = gff
		# Need to grab the AGI, currently using the filename
		self.AGI = os.path.basename(vcffile)[:-4]
		self.fullsites, self.snpsites, self.first_site, self.last_site, self.ref_seq = pvcf.load_pickle(self.AGI, os.path.dirname(vcffile) + '/')
		self.find_exons(self.ref_seq, gff)
		#self.polymorphisms()
		

	def find_exons(self, DNA_sequence, gff):
		''' Pull annotation from GFF files. For a given AGI_ID '''

		records = pgff.parseGFF3(self.gfffile)
		concat = MutableSeq('')
		self.exons = {}
		self.gff_records = []	
		for agi in records:
			try:
				if self.AGI == agi.attributes['Parent']: 
					agi.attributes['Parent']
					self.gff_records.append(agi)
			except KeyError:
				next
		
		for i in self.gff_records:
			exonnum = float(i.attributes['ID'].split(':')[2])
			try:
				self.exons[exonnum] = MutableSeq(
				self.ref_seq[int(i.start) - 
				self.first_site: int(i.end) - 
				self.first_site + 1]
				)
			except IndexError:
				print ''			
	
		print self.exons.keys()	
		print float(self.exons.keys()).sort()	
		#for i in self.exons.keys().sort():
		#	concat = concat + self.exons[i]
		#self.protein_coding = concat
		#print concat	



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
		try:
			return self.get_protein_coding().translate(cds = True)	
		except:
			print self.get_protein_coding()
###############################################################################
gff = 'Araport11_GFF3_genes_transposons.201606.gff.gz'
directory = "./1001g_variants/"
file_list, pickled = pvcf.get_files(directory)

Snp('./1001g_variants/AT1G45249.1.vcf' ,gff)
