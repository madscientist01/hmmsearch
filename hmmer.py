#!/usr/bin/python
import urllib, urllib2,os,re
import xml.etree.ElementTree as ET
import jsonpickle
import argparse
import re
import glob

def readFasta(filename) :
	#
	# Simple FASTA reader. return sequence name with sequence as a string
	#
	sequence = ''
	sequenceName = ''
	if os.path.exists(filename):		
		f = open(filename)
		while True:
			line = f.readline().strip()
			if not line: 
				break 
			match = re.match('^>(.*)',line)
			if match:
				sequenceName = match.group(1)
			else:
				sequence = sequence+line			
		return(sequenceName,sequence)	    	
	else:
		return()

class HmmerHit(object):
	#
	# Class for the Hmmer domain hits
	#
	def __init__(self, **kwargs):
		self.name=kwargs.get('name')
		self.acc=kwargs.get('acc')
		self.desc = kwargs.get('desc')
		self.evalue=kwargs.get('evalue')
		self.ievalue=kwargs.get('ievalue')
		self.cevalue=kwargs.get('cevalue')
		self.start = kwargs.get('start')
		self.end=kwargs.get('end')
		self.bitscore = kwargs.get('bitscore')

# install a custom handler to prevent following of redirects automatically.
class SmartRedirectHandler(urllib2.HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        return headers


class Hmmer(object):

	def __init__(self, **kwargs):

		self.file=kwargs.get('file')
		self.db=kwargs.get('db')
		self.cutoff=kwargs.get('evalue')
		(self.name, self.sequence) = readFasta(self.file)
		self.length=len(self.sequence)
		self.hits = []
		
	def run(self):

		opener = urllib2.build_opener(SmartRedirectHandler())
		urllib2.install_opener(opener)
		print self.name
		parameters = {
		              'hmmdb':self.db,
		              'seq':self.sequence
		             }
		enc_params = urllib.urlencode(parameters);
		#post the seqrch request to the server
		request = urllib2.Request('http://hmmer.janelia.org/search/hmmscan',enc_params)
		#get the url where the results can be fetched from
		results_url = urllib2.urlopen(request).getheader('location')
		# modify the range, format and presence of alignments in your results here
		res_params = {
		              'output':'xml'
		             }
		# add the parameters to your request for the results
		enc_res_params = urllib.urlencode(res_params)
		modified_res_url = results_url + '?' + enc_res_params
		# send a GET request to the server
		results_request = urllib2.Request(modified_res_url)
		data = urllib2.urlopen(results_request)
		# print out the results
		result = data.read()	
		if result:
			#
			# Parse using ElementTree Modules (http://docs.python.org/2/library/xml.etree.elementtree.html)
			#
			root = ET.fromstring(result)
			for child in root.iter('hits'):
				name = child.get('name')
				desc = child.get('desc')
				acc = child.get('acc')
				evalue=child.get('evalue')
				
				for element in child.iter('domains'):
					if (float(element.get('ievalue'))<float(self.cutoff)):
						cevalue = element.get('cevalue')
						ievalue = element.get('ievalue')
						start = element.get('iali')
						end = element.get('jali')
						bitscore = element.get('bitscore')
						hit = HmmerHit(	name=name, desc=desc, acc=acc,bitscore=bitscore,
										evalue=evalue, ievalue=ievalue, cevalue=cevalue, start=start, end=end)
						self.hits.append(hit)

			return(True)
		else:
		    print "Failed to retrieve results" 
		    return(False)

def main(argument):

	if not argument.files:
		files = glob.glob('*.fasta')
	else:
		files =argument.files

	hmmerResults=[]
	for file in files:
		hmmer = Hmmer(file=file,db=argument.db,evalue=argument.evalue)
		if hmmer.run():
			for hit in hmmer.hits:
				print hmmer.name, hmmer.length, hit.name, hit.desc, hit.acc, hit.start, hit.end, hit.cevalue, hit.ievalue, hit.bitscore
			hmmerResults.append(hmmer)

	pickled = jsonpickle.encode(hmmerResults)
	f=open("test.js","w")
	f.write("var data="+pickled)
	f.close()

if __name__ == "__main__":


	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--fasta', nargs='+', dest='files',default=[],
	                    help='Files to process')
	parser.add_argument('-d', '--database', dest='db',default='pfam',
						help='HMM database')
	parser.add_argument('-e', '--evalue_cutoff', dest='evalue',default=1e-5,
						help='E-value cutoff')
	results = parser.parse_args()
	main(results)

