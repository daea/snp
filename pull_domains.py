# Matt Cumming
# Provart Lab

import subprocess, json, re
import handle_json
import dpath.util
import read_1001 as r1001
import pycurl, cStringIO


class Domains():
	
	url = "http://bar.utoronto.ca/eplant/cgi-bin/CDDannot.cgi"
	features = {}
	domains = {}	

	def __init__(self, seq):
		self.query_seq = seq
		rawcdd = pull_cdd_bar(seq, url)		
		decode_json(handle_json.json_loads_byteified(rawcdd))
				


	def decode_json(tata):
		#cheese = byteify(json.loads(tata))
		out_var = dpath.util.get(cheese, '/')
		return out_var



	def pull_cdd_bar(seq, url):
		sequence = 'FASTAseq=' + seq 
		
		buf = cStringIO.StringIO()
		c = pycurl.Curl()
		c.setopt(c.URL, url)
		c.setopt(c.POSTFIELDS, sequence)	
		c.setopt(c.WRITEFUNCTION, buf.write)
		c.perform()
		cheese = buf.getvalue()
		buf.close()
		return cheese


	def read_cdd(decoded):
		for i in decoded.keys():	
			#print i
			decoded_feat[i] = decoded[i].split(',')
			for k in decoded_feat[i]:
				residue = k[0]
				site = k[1:]	
				domains[self.query_seq] = {i:[site, residue]}	
				print domains	


abf3 = r1001.Snp('./1001g_variants/AT4G34000.1.vcf', r1001.gff)

abf3_feat = {}
domains = {}


