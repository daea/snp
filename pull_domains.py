# Matt Cumming
# Provart Lab

import subprocess, json, re
import handle_json
import dpath.util
import read_1001 as r1001
import pycurl, cStringIO
import process_vcfs as pvcf

class Annot():
	
	features = {}
	cdd = {}	

	def __init__(self, seq):
		
		self.pfamurl = 'http://bar.utoronto.ca/eplant/cgi-bin/PfamAnnot.cgi'	
		self.cddurl =  'http://bar.utoronto.ca/eplant/cgi-bin/CDDannot.cgi'
		self.query_seq = seq
		self.rawcdd = self.pull_bar(self.query_seq, self.cddurl)
		self.rawpfam = self.pull_bar(self.query_seq, self.pfamurl)		
		self.cdd_sites = self.read_cdd(handle_json.json_loads_byteified(self.rawcdd))
		#self.read_pfam(handle_json.json_loads_byteified(self.rawpfam))	


	def decode_json(self, tata):
		#cheese = byteify(json.loads(tata))
		out_var = dpath.util.get(cheese, '/')
		return out_var



	def pull_bar(self, seq, url):
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



	def read_cdd(self, decoded):
		decoded_feat = {}
		for i in decoded.keys():	
			decoded_feat[i] = decoded[i].split(',')
			for k in decoded_feat[i]:
				residue = k[0]
				site = k[1:]	
				self.cdd[self.query_seq] = {i:[site, residue]}	


for i in pvcf.get_files('./1001g_variants/'):
	for k in i:
		temp = r1001.Snp(k, r1001.gff)
		temp_annot = Annot(str(temp.get_protein_seq())) 
		print temp_annot.query_seq, temp_annot.cdd_sites, temp_annot.rawpfam



#abf3 = r1001.Snp('./1001g_variants/AT4G34000.1.vcf', r1001.gff)
#Annot(str(abf3.get_protein_seq()))

