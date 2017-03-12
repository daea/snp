#!/usr/bin/env python

"""
A simple parser for the GFF3 format.

Test with transcripts.gff3 from
http://www.broadinstitute.org/annotation/gebo/help/gff3.html.

Format specification source:
http://www.sequenceontology.org/gff3.shtml

Version 1.0
"""
from __future__ import with_statement
from collections import namedtuple
import gzip   # Matt: For decompression of .gz if encountered
import urllib # Matt: For removal of %characters, replaces with single string equiv
from Bio.Seq import MutableSeq, Seq as MutableSeq, Seq 

__author__  = "Uli Koehler"   
# Matt: Shamelessly stolen from techoverflow.net/blog/2013/11/30/parsing-gff3-in-python/
__license__ = "Apache License v2.0"

#Initialized GeneInfo named tuple. Note: namedtuple is immutable
gffInfoFields = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", gffInfoFields) 

def parseGFFAttributes(attributeString):
    """Parse the GFF3 attribute column and return a dict"""
    if attributeString == ".": return {}
    ret = {}
    for attribute in attributeString.split(";"):
        key, value = attribute.split("=")
        ret[urllib.unquote(key)] = urllib.unquote(value)
    return ret

def parseGFF3(filename):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.
    
    Supports transparent gzip decompression.
    """
    #Parse with transparent decompression
    openFunc = gzip.open if filename.endswith(".gz") else open
    with openFunc(filename) as infile:
        for line in infile:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            #If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            #Normalize data
            normalizedInfo = { 
                "seqid": None if parts[0] == "." else urllib.unquote(parts[0]), 
                "source": None if parts[1] == "." else urllib.unquote(parts[1]),
                "type": None if parts[2] == "." else urllib.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand": None if parts[6] == "." else urllib.unquote(parts[6]),
                "phase": None if parts[7] == "." else urllib.unquote(parts[7]),
                "attributes": parseGFFAttributes(parts[8]) # defined above
            }
            #Alternatively, you can emit the dictionary here, if you need mutability:
            #    yield normalizedInfo
            yield GFFRecord(**normalizedInfo)
            
'''
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="The GFF3 input file (.gz allowed)")
    parser.add_argument("--print-records", action="store_true", help="Print all GeneInfo objects, not only")
    parser.add_argument("--filter-type", help="Ignore records not having the given type")
    args = parser.parse_args()
    #Execute the parser
    recordCount = 0
    for record in parseGFF3(args.file):
        #Apply filter, if any
        if args.filter_type and record.type != args.filter_type:
            continue
        #Print record if specified by the user
        if args.print_records: print record
        #Access attributes like this: my_strand = record.strand
        recordCount += 1
    print "Total records: %d" % recordCount
'''

def slice_ref_seq(sequence, start, end, phase, locus):
	a = MutableSeq(
		sequence[
			int(start) - int(locus) + int(phase) - 1:
			int(end) - int(locus) 		
			]
		)
	print int(start), int(phase), int(end),
	print int(start) - int(locus) + int(phase), int(end) - int(locus)
	return a


def _join_records(recordtype, splice_variant, sequence, locus='null'):
	''' Takes a list of gff records as input, joins them in the correct orientation.'''
	
	records = retrieve_gff_records(recordtype, splice_variant)
	print records
	feats = []
	concat = ''

	# First see how many records we recovered
	# If one, just return whatever it slices	
	if len(records) <= 1:
		# If a locus isn't specified just slice from the beginning 
	
		if locus == 'null':
			locus = int(records[0].start)
		if records[0].strand == '+':	
			feats.append(slice_ref_seq(sequence, records[0].start, records[0].end, record[0].phase, locus))	
		
		elif records[0].strand == '-':
			feats.append(slice_ref_seq(sequence, records[0].start, records[0].end, record[0].phase, locus).reverse_complement())	
			
	else:
		if locus == 'null':
			print "please specify a starting locus"
		else:	
			for feature in records:
				exonnum = float(feature.attributes['ID'].split(':')[2])
				if feature.strand == '-':
					try:
						feats.append(slice_ref_seq(sequence, feature.start, feature.end, feature.phase, locus).reverse_complement())
					except IndexError:
						feature.attributes['ID']
				
				elif feature.strand == '+':
					try:	
						feats.append(slice_ref_seq(sequence, feature.start, feature.end, feature.phase, locus))
						#print feature.start, feature.end, feature.phase, locus, int(feature.start) - int(locus), int(feature.end) - int(locus) + 1	
					except IndexError:
						print 'Cheesy McSneezy'
				else:
					print "There was an issue assigning the feature, {0}, to your list. Check ".format(feature)
		print feats

	try:
		for i in feats:
			concat = concat + i 
	except TypeError:
		print ''	

	print concat.translate()

def retrieve_gff_records(recordtype, splice_variant):
	gfflist = []

	for i in records:
		print i
		try:
			if i.type == recordtype and i.attributes['Parent'] == splice_variant:
				gfflist.append(i)
		
		except IndexError:
			print 'Error in your selection, check your AGI ID, and recordtype'
	if len(gfflist) == 0:
		print "No records match the variant or recordtype selected"
	
	return gfflist	


at4g34000 = 'TGAATCGATTTTTGTTGTCTATTACTGATTGGTTTTCTTGTTCAGATTCACTGATTCGAAGAGAATCATGATTTTTTTTTCCCGCTGAATAATAAGCATATGATTGGGTGTTTTGGAGATTTGTTTACTGATTAAAAGGAGATTCCTTTCCATTTTCACCATTTGCTCTGTTTGACTTCATTGTGCTTATATTTCATTTAGATCTTTTGTTTGGGTTTAGCTTTGGAACTGATAAAAATCTGATTTTGTCTCACGGCTTTGGATTTGGTTCTTAAATTTTGGTACTTTAAAACTGGATAAAGATCAGTGCTTTTTTAGATTCTTCGTTTGTTGATGAATTTATGGATGTATGTATAATTAAACCATAATCTCTCTGCTTGTTTGTTTTCTTATAGGTAAATATCCAGAAGCTTGATCCTCCTAGTTGTACGAAAGCTTGAGTAATGGGGTCTAGATTAAACTTCAAGAGCTTTGTTGATGGTGTGAGTGAGCAGCAGCCAACGGTGGGGACTAGTCTTCCATTGACTAGGCAGAACTCTGTGTTCTCGTTAACCTTTGATGAGTTTCAGAACTCATGGGGTGGTGGAATTGGGAAAGATTTTGGGTCTATGAACATGGATGAGCTCTTGAAGAACATTTGGACTGCAGAGGAAAGTCATTCAATGATGGGAAACAATACCAGTTACACCAACATCAGCAATGGTAATAGTGGAAACACTGTTATTAACGGCGGTGGTAACAACATTGGTGGGTTAGCTGTTGGTGTGGGAGGAGAAAGTGGTGGTTTTTTCACTGGTGGGAGTTTGCAGAGACAAGGTTCACTTACCTTGCCTCGGACGATTAGTCAGAAAAGGGTTGATGATGTCTGGAAGGAGCTGATGAAGGAGGATGACATTGGAAATGGTGTTGTTAATGGTGGGACAAGCGGAATTCCGCAGAGGCAACAAACGCTGGGAGAGATGACTTTGGAGGAGTTTTTGGTCAGGGCTGGTGTGGTTAGGGAAGAACCTCAACCGGTGGAGAGTGTAACTAACTTCAATGGCGGATTCTATGGATTTGGCAGTAATGGAGGTCTTGGGACAGCTAGTAATGGGTTTGTTGCAAACCAACCTCAAGATTTGTCAGGAAATGGAGTAGCGGTGAGACAGGATCTGCTGACTGCTCAAACTCAGCCACTACAGATGCAGCAGCCACAGATGGTGCAGCAGCCACAGATGGTGCAGCAGCCGCAACAACTGATACAGACGCAGGAGAGGCCTTTTCCCAAACAGACCACTATAGCATTTTCCAACACTGTTGATGTGGTTAACCGTTCTCAACCTGCAACACAGGTGACTAAAGAACCTAAGCTACTACATTTTGTAGGAAATACAATCATCTGAAAAATATTAGACGTAGCCTTCATGTTTTAAGATAAGCTGGTTTTGGATATGCATGTATGTTTCAGTTATCTGAGCATGACTTGTTTTTACTTTTTTCGCAGTGCCAGGAAGTGAAGCCTTCAATACTTGGAATTCATAACCATCCTATGAACAACAATCTACTGCAAGCTGTCGATTTTAAAACAGGAGTAACGGTTGCAGCAGTATCTCCTGGAAGCCAGATGTCACCTGATCTGACTCCAAAGAGCGCCCTGGATGCATCTTTGTCCCCTGTTCCTTACATGTTTGGGCGAGTGAGAAAAACAGGTGCAGTTCTGGAGAAAGTGATTGAGAGAAGGCAAAAAAGGATGATAAAGAATAGGGAATCAGCTGCAAGATCCCGCGCTCGCAAGCAAGTGAGTGTTTGTTTAAATTTTGGAGATTAAAGAAACCTTAAAACTGTGACCATGTTATTTACTTTTTCACTTTCTTGCTTGACAGGCTTATACGATGGAACTGGAAGCAGAAATTGCGCAACTCAAAGAATTGAATGAAGAGTTGCAGAAGAAACAAGTGTGTCTCGCTTCTTCCCTATCACAATTAAGAATCTCGAGATTTTCATATTTTCTTGAGGTTGTATTCACTGACCAAATGTTTCATGCAGGTTGAAATCATGGAAAAGCAGAAAAATCAGGTACTGTCTTGATTTGAATATCCTCTATGGTTGTTGGCTAGGCTTTTAACTCTCACTCATAATGAATTACACTTTTGGACAGTATTCTAAGCTTTTGAGTAGAATAGTGTAAGCTATACCATGAAGTGAGACATATCATCACATTTTTGATTTCCCACTCTGCATAAAGTATTTAAGATTTGTGAATATGTTGCAATGCCAATTTGGATATTTCATGAGACTAATCTGACGAGCATGGATTTAACGGAAGTTGGCTCATTTGTTGTTGCAGCTTCTGGAGCCTCTGCGCCAGCCATGGGGAATGGGATGCAAAAGGCAATGCTTGCGAAGGACATTGACGGGTCCCTGGTAGAGCTTATAATGGCGTCTAAGGAACCCAACAAAGCGCCGAAGTTATAGAACAACTCAGAAGATAGAAAGCTAGCTTTGTACGTAGTTTAGGCAGGTTCTGTGGGTGATTGTAAATCTTGAAGTGTGGCGGATTTGACAGAGATAGATAAACACATATCTGTTCTATTTTCCTAAATCTTTTGGTTTTATCTTCCTGATGTAATGGATCTTTATCATTTGTCTTGAACATCTTTGTGACTTAACCAGAGTGAATTTATCTTGTATCTTGTCTGCAATTTTTTCTTT'


records = parseGFF3('TAIR10_GFF3_genes_transposons.AIP.gff.gz')
_join_records('CDS', 'AT4G34000.1', at4g34000, int('16295558'))

