# Matt Cumming 2016
# Provart Lab 2016

'''
Simple script to create a dictionary where the AGI id is linked to all it's aliases in an easy to understand format. I just want to be able to read what genes I'm looking at. It reads from the gene aliases file you provide it.
Usage: agis = read_AGI(file_name); agis['AGI ID']['full_name']
'''
import re, sys

old =  r'[a-zA-Z0-9\. ?]\t'

def agi_alias(agi_id):
	''' Just finds all the instances of each AGI in the given line returns them as a list, splits on tabs.'''
	alias = re.findall(r'[a-zA-Z0-9\. ?] ?[a-zA-Z0-9\. ?]*', agi_id)	
	return [i.rstrip('\t') for i in alias]

def read_AGI(file):	
	''' Reads the TAIR outputted list of Gene aliases. Returns a dictionary where the keys are the AGI IDs.
	    Ex. agis['AT4G34000'] gives you the information for ABF3. You can return lists of the associated symbol,
	    and full_names, through agis['AGI ID']['symbol'], agis['AGI ID']['full_name']
	'''
	agis = {}
	with open(file) as f:
		for line in f:
			a = agi_alias(str(line.rstrip()))	
			if a[0] in agis:
				agis[a[0]]['symbol'].append(a[1])
				agis[a[0]]['full_name'].append(a[2])
			else:	
 				agis[a[0]] = {'locus_name': [a[0]],
					      'symbol': [a[1]],
					      'full_name': [a[2]]
					     }
	return agis		
			
agis = read_AGI(sys.argv[1])

