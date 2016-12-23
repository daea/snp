# Matt Cumming 2016
# Format the raw snp export from eplant

import re, sys

def nice_format(filename):
	''' This just strips all the info out of the ePlant SNP export that I don't need and grabs only the SNP Frequencies''' 
	with open (filename, 'r') as myfile: 
		data = "".join(line.rstrip() for line in myfile)	
	start = re.search(r'SNP', data)
	finish = re.search(r'"retrieved"', data)
	return ("{" + data[start.start():finish.start()-1] + "}}")
