# Matt Cmming 2016 #
# Provart Lab

import vcf, glob, dill, re, os
import pprint as pp


def get_files(directory):
	''' Given a directory, returns a list of .vcfs and .pkls'''
	vcffiles = []
	pickles = []
	
	for root, dirs, files in os.walk(directory, topdown=True):
		for name in files:
			filename = os.path.join(root, name)
			if filename[-4:] == ".pkl":
				pickles.append(filename)
			elif filename[-4:] == ".vcf":
				vcffiles.append(filename)
	return vcffiles, pickles




def read_vcf(vcffile):
	''' Gets the sequence from the vcf files and grabs the associated SNPs'''
	fullsites = []
	snpsites = []
	first_site = 0
	last_site = 0
	gene = vcf.Reader(open(vcffile, 'r'))
	for record in gene:

		# Find first and last sites in the vcf for index operations later.
		if record.POS > last_site:
			last_site = record.POS
		if first_site == 0:
			first_site = record.POS

		# Append all the sites to create a reference gene.
		''' record.REF[0] is used because indels are listed at the first
			site at which they occur, and the REFerence sequence contains
			the entire indel if a deletion is listed in the record.ALT. So
			we only want the first base for the reference sequence here
		''' 
		fullsites.append(record.REF[0]) #  
		
		# Create a list to be returned later containing all the info we need.
		if record.is_snp == True:
			alts = re.search(r'(\w)', str(record.ALT))
			snpsites.append(
				[record.POS,
				 record.REF,
				 alts.group(1),
				 record.aaf,
				 record.is_sv]
			)
			
		# Return ref_seq as a string.	
		ref_seq = str(''.join(fullsites))
	return fullsites, snpsites, first_site, last_site, ref_seq



	
def store_genes(filename):
	directory = os.path.dirname(filename)
	AGI = os.path.basename(filename)[:-4]
	destfile = directory + '/' + AGI + ".pkl"
	
	try:
		print "Processing {0}: ".format(AGI),
		with open(destfile, 'wb') as f:
			info = read_vcf(filename)
			dill.dump(info, f)
		print "{0} successfully saved.".format(destfile)
		return info

	except IOError:
		print "The file was not saved.".format(destfile)	
		return False



def	update_allvcfs(vcf_directory):		
	vcfs, pickles  = get_files(vcf_directory)
	# First we need to remove files that have already been processed from our to-do list	
	done = []
	[done.append(os.path.basename(i)[:-4]) for i in pickles]	
	for name in vcfs:
		todo = os.path.basename(name)[:-4] 
		if todo in done:
			continue		
		else:	
			return store_genes(name)




def load_pickle(AGI, directory):
	try:
		with open(directory + AGI + ".pkl") as f:
			info = dill.load(f)
		return info
	except:
		info = store_genes(directory)
		return info
	
		
		
		

