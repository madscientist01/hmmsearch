#!/usr/bin/python
#
# Hmmer3 based Protein Domain Visualizer 
#
# Written by Madscientist (http://madscientist.wordpress.com, https://github.com/madscientist01/)
#

from __future__ import division
import urllib
import urllib2
from urllib2 import HTTPError
import os
import xml.etree.cElementTree as ET
import argparse
import glob
import sys
import subprocess
from hmmerhit import HmmerHit
from htmltable import HTMLTable
from inputfile import InputFile
from fetchutil import readFasta, readAccession
from svgdrawer import SVGDrawer
from annotation import miscAnnotation
from uniprotannotation import UniprotAnnotation
from phmmer import PhmmerSearch
from psipred import PsipredAnnotation

# install a custom handler to prevent following of redirects automatically.
class SmartRedirectHandler(urllib2.HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        return headers

class SVGList(object):

	def __init__(self):

		self.header = "<br>"
		self.footer = "<br>"
		self.svgEmbeddingTemplate = """
<div id="{0}">{1}</div>
"""
		self.svgEmbedContent=""
		
	def svgContentFill(self,svgContent):
		fill = self.svgEmbeddingTemplate.format(*svgContent)
		self.svgEmbedContent = self.svgEmbedContent+fill

class Hmmer(object):
	#
	# Hmmscan wrapper class
	#
	def __init__(self, **kwargs):

		self.file=kwargs.get('file')
		self.db=kwargs.get('db')
		self.cutoff=kwargs.get('evalue')
		self.localHmmDB = kwargs.get('localHmmDB')
		self.threshold = kwargs.get('threshold')
		fasta = readFasta(self.file)
		if fasta:
			(self.name, self.sequence) = fasta
		else:
			self.name = None
			self.sequence = None
		accession = readAccession(self.file)
		if accession:
			(self.source,self.accession)= accession
		else:
			self.source = None
			self.accession = None
		self.length=len(self.sequence)
		self.features = {}
		hits=[]
		self.features['domain']=hits
		self.tier = {}
		self.tier[0]='Pfam'
	
	def exclude(self):
		#
		# Exclude overlapped domain. If two domains are overlapped, the one have higher bitscore will be retained
		#
		if len(self.features['domain'])>1:
			for i in range (len(self.features['domain'])-1):			
				for j in range (i+1,len(self.features['domain'])):
					if (not self.features['domain'][i].exclude and not self.features['domain'][j].exclude):	
						if (self.features['domain'][i].end > self.features['domain'][j].start) and (self.features['domain'][i].end<self.features['domain'][j].end):
							if (self.features['domain'][i].bitscore > self.features['domain'][j].bitscore):
								self.features['domain'][j].exclude = True			
							else:
								self.features['domain'][i].exclude = True

						if (self.features['domain'][j].end > self.features['domain'][i].start) and (self.features['domain'][j].end<self.features['domain'][i].end):
							if (self.features['domain'][i].bitscore > self.features['domain'][j].bitscore):
								self.features['domain'][j].exclude = True			
							else:
								self.features['domain'][i].exclude = True
			
	def runRemote(self):
		#
		# Using Hmmscan in Hmmer3 web service, find locations of domains in the Fasta sequence and store into class
		#
		if os.path.exists(self.file+".xml") :
			print "{0} is already processed. Skipped.".format(self.file+".xml")
			f = open(self.file+".xml")
			read = f.readlines()
			self.parseHmmerScanXML(''.join(read))
			return(True)
		else:
			opener = urllib2.build_opener(SmartRedirectHandler())
			urllib2.install_opener(opener)
			print "Running Hmmscan of {0} using web service..".format(self.file)
			if not self.db in ['pfam','superfamily','tigrfam', 'gene3d']:
				print "{0} is not valid db. It should be pfam, superfamily, tigrfam or gene3d."
				print "search will be carried out with pfam"
				self.db = 'pfam'

			parameters = {
			              'hmmdb':self.db,
			              'threshold':self.threshold,
			              'seq':self.sequence
			             }
			enc_params = urllib.urlencode(parameters);
			
			try:
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
					f=open(self.file+".xml", "w")
					f.write(result)
					f.close()
					
					self.parseHmmerScanXML(result)
					return(True)
				else:
					print "Failed to retrieve results" 
					return(False)
			except HTTPError:
				print "Hmmscan error"
				return(False)

        def parseHmmerScanXML(self, result):
			"""
			Parse Hmmerscan XML output into 

			"""
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
						labellink = "http://pfam.sanger.ac.uk/family/{0}".format(acc)
						if self.source=='uniprot':
							numberlink = "http://www.uniprot.org/blast/?about={0}[{1}-{2}]".format(self.accession, start, end)

						hit = HmmerHit(	name=name, desc=desc, acc=acc,bitscore=bitscore,labellink=labellink, numberlink=numberlink,
										evalue=evalue, ievalue=ievalue, cevalue=cevalue, start=start, end=end)
						self.features['domain'].append(hit)

			self.exclude()
			return()

	def runLocal(self):
		#
		# Using Hmmscan in locally installed Hmmer3 packages, find locations of domains in the Fasta sequence and store 
		#
		if os.path.exists(self.db):
			outputFile = self.file+".tb"
			hmmeroutFile = self.file+".hmmer"
			print "Running Hmmscan of {0} using local hmmer3".format(self.file)
			#
			# hmmscan --domtbout outputFile --cut_ga, DBfile, queryFile
			#
			if self.threshold=='cut_ga':
				p = subprocess.Popen(['hmmscan','--domtblout',outputFile,'--cut_ga',self.db, self.file],stdout=subprocess.PIPE)
			else:
				p = subprocess.Popen(['hmmscan','--domtblout',outputFile,self.db, self.file],stdout=subprocess.PIPE)
			
			p_stdout = p.stdout.read()
			fw = open(hmmeroutFile,'w')
			fw.write(p_stdout)
			fw.close()

			f = open(outputFile,'r')
			data = f.readlines()
			f.close()
			#
			# Parse table formatted domain datafile generated by '--domtblout'.
			#
			for line in data:
				if line[0:1] != "#":
					splited = line.split()			
					if (float(splited[12])<float(self.cutoff)):
						cevalue = splited[11]
						ievalue = splited[12]
						start = splited[17]
						end = splited[18]
						bitscore = splited[7]
						acc = splited[1]
						desc = line[181:202]
						evalue = splited[6]
						name = splited[0]
						labellink = "http://pfam.sanger.ac.uk/family/{0}".format(acc)
						if self.source == "uniprot":
							positionlink = "http://www.uniprot.org/blast/?about={0}[{1}-{2}]".format(self.accession, start, end)
							hit = HmmerHit(	name=name, desc=desc, acc=acc,bitscore=bitscore,labellink=labellink, numberlink=positionlink,
											evalue=evalue, ievalue=ievalue, cevalue=cevalue, start=start, end=end)
						else:
							hit = HmmerHit(	name=name, desc=desc, acc=acc,bitscore=bitscore,labellink=labellink,
											evalue=evalue, ievalue=ievalue, cevalue=cevalue, start=start, end=end)
						self.features['domain'].append(hit)

			self.exclude()			
			return(True)
		else:
			print "{0} file is not exist".format(self.db)
			sys.exit() 
			return(False)


class HmmerScanRunner(object):
	#
	# Class for the Hmmer domain hits
	#
	def __init__(self, **kwargs):
		self.inputFileName = kwargs.get('inputFileName')
		self.inputFile = None
		self.files=kwargs.get('files')
		self.db = kwargs.get('db')
		self.evalue = kwargs.get('evalue')
		self.local = kwargs.get('local')
		self.outputHTML = kwargs.get('outputHTML')
		self.outputSVG = kwargs.get('outputSVG')
		self.threshold = kwargs.get('threshold')
		self.scaleFactor = kwargs.get('scaleFactor',1.5)
		self.hmmerResults = []
		self.titlemode = kwargs.get('titlemode',False)
	
	def processHmmerResults(self):

	#
	# inject information of input file into hmmerResults 
	#
		if self.inputFile:
			colors = ['aliceblue','antiquewhite', 'aqua', 'aquamarine', 'azure', 'beige', 'bisque', 'black', 'blanchedalmond', 
					'blue', 'blueviolet', 'brown', 'burlywood', 'cadetblue', 'chartreuse', 'chocolate', 'coral', 
					'cornflowerblue', 'cornsilk', 'crimson', 'cyan', 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray', 
					'darkgreen', 'darkgrey', 'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid', 
					'darkred', 'darksalmon', 'darkseagreen', 'darkslateblue', 'darkslategray', 'darkslategrey', 'darkturquoise', 
					'darkviolet', 'deeppink', 'deepskyblue', 'dimgray', 'dimgrey', 'dodgerblue', 'firebrick', 'floralwhite', 
					'forestgreen', 'fuchsia', 'gainsboro', 'ghostwhite', 'gold', 'goldenrod', 'gray', 'green', 'greenyellow',
			 		'grey', 'honeydew', 'hotpink', 'indianred', 'indigo', 'ivory', 'khaki', 'lavender', 'lavenderblush',
			  		'lawngreen', 'lemonchiffon', 'lightblue', 'lightcoral', 'lightcyan', 'lightgoldenrodyellow', 'lightgray',
			   		'lightgreen', 'lightgrey', 'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightslategray',
			    	'lightslategrey', 'lightsteelblue', 'lightyellow', 'lime', 'limegreen', 'linen', 'magenta', 'maroon',
			     	'mediumaquamarine', 'mediumblue', 'mediumorchid', 'mediumpurple', 'mediumseagreen', 'mediumslateblue', 
			     	'mediumspringgreen', 'mediumturquoise', 'mediumvioletred', 'midnightblue', 'mintcream', 'mistyrose', 
			     	'moccasin', 'navajowhite', 'navy', 'oldlace', 'olive', 'olivedrab', 'orange', 'orangered', 'orchid', 
			     	'palegoldenrod', 'palegreen', 'paleturquoise', 'palevioletred', 'papayawhip', 'peachpuff', 'peru', 
			     	'pink', 'plum', 'powderblue', 'purple', 'red', 'rosybrown', 'royalblue', 'saddlebrown', 
			     	'salmon', 'sandybrown', 'seagreen', 'seashell', 'sienna', 'silver', 'skyblue', 'slateblue', 
			     	'slategray', 'slategrey', 'snow', 'springgreen', 'steelblue', 'tan', 'teal', 'thistle', 
			     	'tomato', 'turquoise', 'violet', 'wheat', 'white', 'whitesmoke', 'yellow', 'yellowgreen']

			for i in range(len(self.inputFile.fileNames)):
				for hmmerResult in self.hmmerResults:
					if hmmerResult.file == self.inputFile.fileNames[i]:
						hmmerResult.name = self.inputFile.proteinNames[i]

			for domain in self.inputFile.domainDefinitions:
				for hmmerResult in self.hmmerResults:
					if (domain.proteinName == hmmerResult.name):
						if domain.action == 'remove':
							for (tier, hits) in hmmerResult.features.items():
								for hit in hits:
									if (hit.name == domain.domainName) and (hit.start == domain.start) and (hit.end == domain.end):
										hit.exclude = True	
						if domain.action == 'new':
							pseudoHit = HmmerHit(name=domain.domainName, start=domain.start, end=domain.end)
							hmmerResult.features['domain'].append(pseudoHit)
			
			for renameDef in self.inputFile.renameDefinitions:
				for hmmerResult in self.hmmerResults:
					for (tier, hits) in hmmerResult.features.items():
						for hit in hits:
							if (hit.name == renameDef.domainName):
								hit.name = renameDef.newName


			for colorDef in self.inputFile.colorDefinitions:
				gradient = False
				border = True
				label = True
				number = True
				startshow = True
				endshow = True
				if colorDef.option:
					options = colorDef.option.strip('+')
					if 'gradient' in options:
						gradient = True
					if 'noborder' in options:
						border = False
					if 'nolabel' in options:
						label = False
					if 'nonumber' in options:
						startshow = False
						endshow = False
					if 'nostart' in options:
						startshow = False
					if 'noend' in options:
						endshow = False

				for hmmerResult in self.hmmerResults:
					for (tier, hits) in hmmerResult.features.items():
						for hit in hmmerResult.hits:
							if (hit.name == colorDef.domainName and colorDef.color in colors):
								hit.color = colorDef.color
								hit.gradient = gradient
								hit.border = border
								hit.label = label
								hit.number = number
								hit.startshow = startshow
								hit.endshow = endshow


		return 
	

	def loadAnnotation(self):

		for hmmerResult in self.hmmerResults:
			disorder = miscAnnotation(method='disorder',tier=1,color='grey')
			disorder.readFile(hmmerResult.file)
			hmmerResult.features['disorder'] = disorder.hits
			hmmerResult.tier[1]='Disorder'
			coils = miscAnnotation(method='coils',tier=2,color='green')
			coils.readFile(hmmerResult.file)
			hmmerResult.features['coils'] = coils.hits
			hmmerResult.tier[2]='Coiled-coil'
			uniprot = UniprotAnnotation(tier=3)
			uniprot.readFile(hmmerResult.file)
			hmmerResult.features['uniprot'] = uniprot.hits
			hmmerResult.tier[3]='UniProt'
			phmmer = PhmmerSearch(tier=4)
			phmmer.readFile(hmmerResult.file)
			hmmerResult.features['PDB']=phmmer.hits
			hmmerResult.tier[4]='PDB'
			psipred = PsipredAnnotation(file=hmmerResult.file,tier=5)
			results = psipred.psipredReader()
			if results:
				hmmerResult.features['SS']=psipred.hits
				hmmerResult.tier[5]='SS'

		return

	def run(self):

		if self.inputFileName:
			self.inputFile=InputFile(inputFileName=self.inputFileName)
			self.inputFile.readInputFile()
			files = self.inputFile.fileNames
			self.db = self.inputFile.db
			self.threshold = self.inputFile.threshold
			self.local = self.inputFile.local
			self.evalue = self.inputFile.evalue
			self.outputHTML = self.inputFile.outputHTML
			self.outputSVG = self.inputFile.outputSVG
		else:
			if not self.files:
				self.files = glob.glob('*.fasta')
			files =self.files

		if self.threshold:
			threshold = "cut_ga"
		else:
			threshold = "No"

		for file in files:
			hmmer = Hmmer(file=file,db=self.db,evalue=self.evalue,threshold=threshold)
			if self.local:
				if hmmer.runLocal():
					self.hmmerResults.append(hmmer)	
			else:	
				if hmmer.runRemote():
					self.hmmerResults.append(hmmer)

		if len(self.hmmerResults)>0:
			self.processHmmerResults()
			self.loadAnnotation()
			self.hmmerResults.sort(key=lambda x:x.name)
			draw = SVGDrawer(outputSVG=self.outputSVG, scaleFactor = self.scaleFactor, hmmerResults = self.hmmerResults, outputHTML=self.outputHTML, titlemode = self.titlemode)
			draw.drawSVG()
			draw.sequencesDraw()
			(svgFileNames,svgContent) = draw.drawMultiSVG()
			header = ['Accession','Name','Domain','length']
			table = HTMLTable(header = header)
			svgList = SVGList()
		
			for hmmerResult in self.hmmerResults:
				print hmmerResult.name, hmmerResult.accession, hmmerResult.source
				if hmmerResult.source == 'refseq':
					linkAddress = "<a href='http://www.ncbi.nlm.nih.gov/protein/{0}'>{0}</a>"
				if hmmerResult.source == 'uniprot':
					linkAddress = "<a href='http://www.uniprot.org/uniprot/{0}'>{0}</a>"
				if hmmerResult.source in ['refseq', 'uniprot']:
					accessionLink = linkAddress.format(hmmerResult.accession)
					nameLink = "<a href='#{1}'>{0}</a>".format(hmmerResult.name,hmmerResult.accession)
				else:
					accessionLink = ""
					nameLink = ""

				domain = []
				for hit in hmmerResult.features['domain']:
					if hit.acc:
						domainName="<a href='http://pfam.sanger.ac.uk/family/{0}'>{1}</a>".format(hit.acc, hit.name)
					else:
						domainName = hit.name
					startEnd = "({0}-{1})".format(hit.start, hit.end)
					if hmmerResult.source =='uniprot':
						startEndLink = "<a href='http://www.uniprot.org/blast/?about={0}[{1}-{2}]'>{3}</a>".format(hmmerResult.accession, hit.start, hit.end, startEnd)
					else:
						startEndLink = startEnd
					domainName = domainName + startEndLink
					domain.append(domainName)
				domains=','.join(domain)
				length = hmmerResult.length			
				table.tableContentFill([accessionLink, nameLink, domains, length])
				svgList.svgContentFill([hmmerResult.accession,svgContent[hmmerResult.name]])
			svg = svgList.svgEmbedContent
			table.extra = "<br><div>"+svg+"</div>"
			table.tableGenerate(self.outputHTML)

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--fasta', nargs='+', dest='files',default=[],
	                    help='Files to process')
	parser.add_argument('-d', '--database', dest='db',default='pfam',
						help='HMM database')
	parser.add_argument('-e', '--evalue_cutoff', dest='evalue',default=1e-5,
						help='E-value cutoff')
	parser.add_argument('-l', '--local',action='store_true', dest='local', default=False, 
						help='run local Hmmer')
	parser.add_argument('-o', '--outputhtml', dest='outputHTML', default='output.html', 
						help='Output HTML filename')
	parser.add_argument('-s', '--outputsvg', dest='outputSVG', default='output.svg', 
						help='Output SVG filename')
	parser.add_argument('-t', '--no_threshold', dest='threshold',action='store_false', default=True,
						help='Turn of Pfam gathering threshold. Enable to look up more weak(unreliable) domains')
	parser.add_argument('-i', '--input_file', dest='inputFileName', default='hmmer.INP',
						help='Read configuration file')
	results = parser.parse_args()
	if not os.path.exists(results.inputFileName):
		results.inputFileName = None
	hmmerscan = HmmerScanRunner(files=results.files,inputFileName=results.inputFileName, db=results.db, evalue=results.evalue,
								local = results.local, outputHTML = results.outputHTML, outputSVG = results.outputSVG, threshold=results.threshold )
	hmmerscan.run()
	

