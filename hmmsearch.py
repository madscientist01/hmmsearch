#!/usr/bin/python
#
# Hmmer3 based Protein Domain Visualizer 
#
# Written by Madscientist (http://madscientist.wordpress.com, https://github.com/madscientist01/)
#

# mysql query 
# mysql -N --user=genome --host=genome-mysql.cse.ucsc.edu -A -D hg19 -e 'select spDisplayID,geneSymbol,refseq from kgXref inner join kgProtAlias on kgXref.kgID=kgProtAlias.kgID where alias="NP_005883" or alias="NP_005210 g"' 
from __future__ import division
import urllib
import urllib2
import os
import xml.etree.ElementTree as ET
import argparse
import re
import glob
import sys
import subprocess

def extractMultiFasta(multifasta,searchId):
	
	capture = False
	fileNames=[]
	buffer=[]
	f=open(multifasta)
	header = re.compile('^>(.*)')
	gi = re.compile(:')
	# gid = ''
	refseq = ''
	while True:
		line = f.readline()
		if not line:
			break
		match = header.match(line)
		if match:
			if capture:
				print refseq+'.fasta'
				fw = open(refseq+'.fasta','w')
				fw.writelines(buffer)
				fw.close
				fileNames.append(refseq+'.fasta')
				capture=False
				refseq=''
				buffer=[]

			id = match.group(1).strip()
			if id in searchId:
				capture = True
				mat = gi.match(id)
				if mat:
					# gid = mat.group(1)
					refseq = mat.group(2)
		if capture:
			buffer.append(line)

	if capture:
		print refseq+'.fasta'
		fw = open(refseq+'.fasta','w')
		fw.writelines(buffer)
		fw.close
		capture=False
		refseq=''
		buffer=[]
		fileNames.append(refseq+'.fasta')

	return(fileNames)		


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
		return(None)

def fetchfromUniprot(uniprotId):
	#
	# fetch fasta file using Rest interface of uniprot 
	#  http://www.uniprot.org/faq/28
	#
	uniprotURL = "http://www.uniprot.org/uniprot/{0}.fasta".format(uniprotId)
	fileName = uniprotId+".fasta"
	try:
		req = urllib2.Request(uniprotURL)
		u = urllib2.urlopen(req)
		sequence = u.read()
		fileName = uniprotId+".fasta"
		f = open(fileName,'w')
		f.write(sequence)
		f.close()
		return(fileName)

	except urllib2.URLError:
		print "{0} is not valid uniprot id".format(uniprotId)
    	return(None)

	return(fileName)


class SmartRedirectHandler(urllib2.HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        return headers


class HmmFetch(object):
	#
	# Hmmscan wrapper class
	#
	def __init__(self, hmmdb, domainName):

		self.hmmdb=hmmdb
		self.domainName = domainName

	def run(self):
		hmmFileName = self.domainName+'.hmm'
		p = subprocess.Popen(['hmmfetch','-o',hmmFileName,self.hmmdb,self.domainName],stdout=subprocess.PIPE)
		p_stdout = p.stdout.read()
		regex = re.compile("Retrieved HMM")
		match = regex.search(p_stdout)
		if match:
			print "{0} is generated".format(hmmFileName)
			return(hmmFileName)
		else:
			print "Problems. {0} is not generated".format(hmmFileName)
			os.remove(hmmFileName)
			return(None)


class HmmerDomainHit(object):
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
		self.start = int(kwargs.get('start'))
		self.end=int(kwargs.get('end'))
		self.bitscore = kwargs.get('bitscore')
		self.exclude = False
		self.color = None
		self.label = True
		self.border = True
		self.startshow = True
		self.endshow = True
		self.gradient = False

class HmmerSearchHit(object):

	def __init__(self, **kwargs):

		self.target=kwargs.get('target')
		self.query=kwargs.get('query')
		self.desc = kwargs.get('desc')
		self.acc=kwargs.get('acc')
		self.evalue = kwargs.get('evalue')
		self.bitscore = kwargs.get('bitscore')
		self.domainNo = kwargs.get('domainNo')
		self.domains = []

class HmmerSearch(object):
	#
	# Hmmscan wrapper class
	#
	def __init__(self, **kwargs):

		self.file=kwargs.get('file')
		self.db=kwargs.get('db')
		self.cutoff=kwargs.get('evalue')
		self.localHmmDB = kwargs.get('localHmmDB')
		self.threshold = kwargs.get('threshold')
		self.hits = []

	def runLocal(self):
		#
		# Using Hmmscan in locally installed Hmmer3 packages, find locations of domains in the Fasta sequence and store 
		#
		if os.path.exists(self.db):
			proteinOutputFile = self.file+".tb"
			domainOutputFile = self.file+".dtb"
			hmmeroutFile = self.file+".hmmer"
			print "Processing {0}".format(self.file)
			#
			# hmmscan --domtbout outputFile --cut_ga, DBfile, queryFile
			#
			if self.threshold=='cut_ga':
				p = subprocess.Popen(['hmmsearch','--domtblout',domainOutputFile,'--tblout',proteinOutputFile, '--cut_ga',self.file, self.db],stdout=subprocess.PIPE)
			else:
				p = subprocess.Popen(['hmmsearch','--domtblout',domainOutputFile,'--tblout',proteinOutputFile, self.file, self.db],stdout=subprocess.PIPE)
			
			p_stdout = p.stdout.read()
			fw = open(hmmeroutFile,'w')
			fw.write(p_stdout)
			fw.close()

			f = open(proteinOutputFile,'r')
			data = f.readlines()
			f.close()
			#
			# Parse table formatted domain datafile generated by '--domtblout'.
			#
			for line in data:
				if line[0:1] != "#":
					splited = line.split()
					if (float(splited[4])<self.cutoff):
						target = splited[0]
						query = splited[2]
						acc = splited[3]
						evalue = splited[4]
						bitscore = splited[5]
						domainNo = splited[17]
						desc = line[156:]
						hit = HmmerSearchHit(target=target, query=query, desc=desc, acc=acc, bitscore=bitscore,
										evalue=evalue, domainNo=domainNo)
						self.hits.append(hit)		
			return(True)
		else:
			print "{0} file is not exist".format(self.db)
			sys.exit() 
			return(False)

def main(argument):

	hmmFiles = []
	if (len(argument.domains)>0):
		domainNames = argument.domains
		hmmdb = argument.hmmdb
		for domainName in domainNames:
			singleHmm = HmmFetch(hmmdb, domainName)
			hmmFile = singleHmm.run()
			if hmmFile:
				hmmFiles.append(hmmFile)

	if (len(argument.hmmFiles)>0):
		for hmmFile in argument.hmmFiles:
			if os.path.exists(hmmFile):
				hmmFiles.append(hmmFile)
	if argument.threshold:
		threshold="cut_ga"
	else:
		threshold="No"

	ids = []
	for hmmFile in hmmFiles:

		hmmSearch = HmmerSearch(file=hmmFile, db=argument.proteindb, evalue=argument.evalue, threshold=threshold)
		hmmSearch.runLocal()
		for hit in hmmSearch.hits:
			ids.append(hit.target+' '+hit.desc.strip())
	extractMultiFasta(argument.proteindb,ids)


if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--domains', nargs='+', dest='domains',default=[],
	                    help='Domains to search')
	parser.add_argument('-hd', '--hmm', dest='hmmFiles',default = [], help='HMM db')
	parser.add_argument('-p', '--proteindb', dest='proteindb', default = 'protein.fa',
						help='protein db')
	parser.add_argument('-b', '--hmmdb', dest='hmmdb', default = 'Pfam-A.hmm',
						help='protein db')
	parser.add_argument('-e', '--evalue_cutoff', dest='evalue',default=1e-5,
						help='E-value cutoff')
	parser.add_argument('-l', '--local',action='store_true', dest='local', default=False, 
						help='run local Hmmer')
	parser.add_argument('-t', '--no_threshold', dest='threshold',action='store_false', default=True,
						help='Turn of Pfam gathering threshold. Enable to look up more weak(unreliable) domains')
	results = parser.parse_args()
	main(results)

