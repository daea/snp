# Matt Cumming 2015
# parsing the ePlant export of SNPs from the "Molecule Viewer"

import json			# for handling json data from ePlant
from pprint import pprint	# for nice printing, when checking ePlant nested data
import handle_json		# small function code to turn unicode strings into bytes 
import dpath.util		# allows us to traverse nested dictionaries as if they were filesystems
import os, sys, re
# You need to pull the ePlant data, and make sure the SNPs file is in the correct json format using a json validator

class SNP:
	""" Class will need to be on a per Transcription factor / Enzyme basis. Has to include SNPs by site, frequency, to / from residue """
	def __init__(self, id):
		self.id = id


# pull our json file
def load_entry(file):
	""" Returns the JSON file of SNPs as a nested dictionary. This just pulls the SNP file from whatever the text file is named, and converts it to strings so we can do stuff to it easier. dpath returns a dictionary at the level that you specify; here it's the level just under 'data' """

	with open(file) as file:
		cheese = byteify(json.loads(file.read()))
	# Pull just the accessions from the SNP data
	out_var = dpath.util.get(cheese, '/SNP Raw Data/data')
	return out_var



def list_accessions(dictionary):
	""" Return a list of all the accessions that are included in ePlant that have SNPs for the given sequence.This is mainly included as a double check for me."""
	mylist = []
	for i in dictionary.keys():
		mylist.append(i)
	return sorted(mylist)



def get_snps(infile):
	""" This is where we are going to have to build our dataframe containing our SNPs and other info. Eventually it will have to include all the other metadata that you get from the ePlant export """		
	snp_freq = {}
	snps = load_entry(infile)
	a = list_accessions(snps)
	for accession in a:
		for base in snps[accession].keys():
			#print base
			if base in snp_freq.keys():
			#	print snp_freq[base]
				if snps[accession][base] in snp_freq[base]: 
					snp_freq[base][snps[accession][base]] += 1
				else:
					snp_freq[base][snps[accession][base]] = 1
			else:
				snp_freq[base] = {snps[accession][base]: 1}
				
	return snp_freq	



def byteify(input):
	''' Converts the default unicode that you get from a json load to bytecode so that I can make it into strings easier.  You call this on json.load or json.loads(file containing json). It returns a json dictionary that you can traverse using dpath pretty easily.'''
	if isinstance(input, dict):
		return {byteify(key): byteify(value)
			for key, value in input.iteritems()}
	elif isinstance(input, list):
		return [byteify(element) for element in input]
	elif isinstance(input, unicode):
		return input.encode('utf-8')
	else:
		return input	
	


def nice_format(filename):
	''' This just strips all the info out of the ePlant SNP export that I don't need and grabs only the SNP Frequencies''' 
	with open (filename, 'r') as myfile: 
		data = "".join(line.rstrip() for line in myfile)	
	start = re.search(r'SNP Raw Data', data)
	finish = re.search(r'"retrieved"', data)
	f = open('out.txt', 'w')
	f.write( "{" + 
		 '"' + 
		 start.group(0) + 
		 '":' +
		 data[start.start()+13:finish.start()-1] + 
		 "}}")
	f.close()
	return ( "{" + data[start.start():finish.start()-1] + "}}")	
	
	
#nice_format('ePlant_snps_AT4G34000_raw.txt')




