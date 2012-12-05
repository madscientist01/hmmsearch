#!/usr/bin/python
#
# Hmmer3 based Protein Domain Visualizer 
#
# Written by Madscientist (http://madscientist.wordpress.com, https://github.com/madscientist01/)
#

from __future__ import division
import urllib
import urllib2
import os
import xml.etree.cElementTree as ET
import argparse
import re
import glob
import sys
import subprocess
from htmltable import HTMLTable

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
		f.close()			
		return(sequenceName,sequence)	    	
	else:
		return(None)

def readAccession(filename) :
	#
	# From FASTA header, guess accession number and source
	#
	refseqRegex = re.compile('>gi\|(\S+)\|ref\|(\S+)\|')
	uniprotRegex = re.compile('>sp\|(\S+)\|(\S+)')
	db = ""
	accession = ""
	if os.path.exists(filename):		
		f = open(filename)
		line = f.readline().strip()
		if not line: 
			return(None) 
		refSeqMatch = refseqRegex.match(line)
		if refSeqMatch:
			db = "refseq"
			accession = refSeqMatch.group(2)
		else:
			uniprotMatch = uniprotRegex.match(line)
			if uniprotMatch:
				db = "uniprot"
				accession = uniprotMatch.group(1)	
		f.close()
	return(db, accession)	    	
	


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

class DomainDefinition(object):

	def __init__(self, **kwargs):
		self.proteinName = kwargs.get('proteinName')
		self.domainName = kwargs.get('domainName')
		self.accession = kwargs.get('acession')
		self.start = int(kwargs.get('start'))
		self.end = int(kwargs.get('end'))
		self.action = kwargs.get('action')

class ColorDefinition(object):
	def __init__(self, **kwargs):
		self.domainName = kwargs.get('domainName')
		self.accession = kwargs.get('acession')
		self.color = kwargs.get('color')
		self.option = kwargs.get('option')

class RenameDefinition(object):
	def __init__(self, **kwargs):
		self.domainName = kwargs.get('domainName')
		self.newName = kwargs.get('newName')

class InputFile(object):
	#
	# Input File (Hmmer.INP) processing class
	#
	def __init__(self, **kwargs):
		self.domainDefinitions = []
		self.colorDefinitions = []
		self.renameDefinitions = []
		self.inputFileName = kwargs.get('inputFileName')
		self.fileNames = []
		self.proteinNames = []
		self.searchDB = ""
		self.outputFileName = ""
		self.commands = {'FILE':'file','NEW':'new','RENAME':'rename','REMOVE':'remove', 'COLOR':'color', 'SIZE':'size', 'EVALUE':'evalue', 'DB':'db'}

	def file(self,content):
		data = content.split('\t')
		uniprot = re.compile("uniprot:(\S+)")
		if len(data)>1:
			name = data[0]
			fileName = data[1]
			match = uniprot.match(fileName)
			if match:
				fileName = fetchfromUniprot(match.group(1))
			if fileName:
				if os.path.exists(fileName):
					self.fileNames.append(fileName)
					self.proteinNames.append(name)
		return

	def new(self,content):
		data = content.split('\t')
		if len(data) > 1:
			proteinName = data[0]
			domainName = data[1]
			if data[2] < data[3]:
				start=data[2]
				end = data[3]
			else:
				start=data[3]
				end = data[2]

			if len(data) > 4:
				color = data[4]
				colorDef = ColorDefinition(domainName=domainName, color=color)
				self.colorDefinitions.append(colorDef)
			if len(data) > 5:
				colorDef.option = data[5]

			domain = DomainDefinition(proteinName=proteinName, domainName=domainName, start=start, end=end, action='new')
			self.domainDefinitions.append(domain)
		return

	def rename(self,content):
		data = content.split('\t')
		if len(data) > 1:
			domainName = data[0]
			newName = data[1]

			if len(data) > 2:
				color = data[2]
				colorDef = ColorDefinition(domainName=newName, color=color)
				self.colorDefinitions.append(colorDef)
			if len(data) > 3:
				colorDef.option = data[3]

			rename = RenameDefinition(domainName=domainName, newName=newName)
			self.renameDefinitions.append(rename)
		
		return

	def remove(self,content):
		data = content.split('\t')
		if len(data) > 1:
			proteinName = data[0]
			domainName = data[1]
			if data[2] < data[3]:
				start=data[2]
				end = data[3]
			else:
				start=data[3]
				end = data[2]
			domain = DomainDefinition(proteinName=proteinName, domainName=domainName, start=start, end=end, action='remove')
			self.domainDefinitions.append(domain)
		return

	def color(self,content):
		data = content.split('\t')
		if len(data) > 1:
			domainName = data[0]
			color = data[1]
			colorDef = ColorDefinition(domainName=domainName, color=color)
			self.colorDefinitions.append(colorDef)
			if len(data) > 2:
				colorDef.option = data[2]
		return

	def readInputFile(self):

		if os.path.exists(self.inputFileName):
			f = open (self.inputFileName, 'r')
			inputFileContent = f.readlines()
			regex = re.compile("^(\S+)\s+(.*)$")
			for line in inputFileContent:
				if line[0:1] != "#":
					match = regex.match(line)
					if match:
						keyword = match.group(1)
						content = match.group(2)
						if keyword in self.commands:
							process=getattr(self, self.commands[keyword])
							process(content)
			return(True)
		else:
			return(False)

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

# install a custom handler to prevent following of redirects automatically.
class SmartRedirectHandler(urllib2.HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        return headers


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
		(self.name, self.sequence) = readFasta(self.file)
		(self.source,self.accession)= readAccession(self.file)
		self.length=len(self.sequence)
		self.hits = []
	
	def exclude(self):
		#
		# Exclude overlapped domain. If two domains are overlapped, the one have higher bitscore will be retained
		#

		if len(self.hits)>1:

			# for hit in self.hits:
			# 	print hit.name, hit.start, hit.end, hit.bitscore, hit.exclude, "\n"

			for i in range (len(self.hits)-1):			
				for j in range (i+1,len(self.hits)):
					# if (self.hits[i].acc == self.hits[j].acc and not self.hits[i].exclude and not self.hits[j].exclude):
					if (not self.hits[i].exclude and not self.hits[j].exclude):	
						if (self.hits[i].end > self.hits[j].start) and (self.hits[i].end<self.hits[j].end):
							if (self.hits[i].bitscore > self.hits[j].bitscore):
								self.hits[j].exclude = True			
							else:
								self.hits[i].exclude = True

						if (self.hits[j].end > self.hits[i].start) and (self.hits[j].end<self.hits[i].end):
							if (self.hits[i].bitscore > self.hits[j].bitscore):
								self.hits[j].exclude = True			
							else:
								self.hits[i].exclude = True
			
	def runRemote(self):
		#
		# Using Hmmscan in Hmmer3 web service, find locations of domains in the Fasta sequence and store into class
		#
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

			self.exclude()
			return(True)
		else:
			print "Failed to retrieve results" 
			return(False)

	def runLocal(self):
		#
		# Using Hmmscan in locally installed Hmmer3 packages, find locations of domains in the Fasta sequence and store 
		#
		if os.path.exists(self.db):
			outputFile = self.file+".tb"
			hmmeroutFile = self.file+".hmmer"
			print "Running Hmmscan of {0}".format(self.file)
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
						hit = HmmerHit(	name=name, desc=desc, acc=acc,bitscore=bitscore,
										evalue=evalue, ievalue=ievalue, cevalue=cevalue, start=start, end=end)
						self.hits.append(hit)

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
		self.outputfile = kwargs.get('outputfile')
		self.threshold = kwargs.get('threshold')
		self.hmmerResults = []
		self.titlemode = kwargs.get('titlemode',False)
		
	def drawSVG(self):
		#
		# Draw SVG based on the hmmer domaim
		#
		x = 50
		y = 50
		leftMargin = 150
		rightMargin = 100
		fontSize = 16
		hmmcolors = {}
		maxLength = 0
		maxTitleLength = 0
		boxHeight = 20
		maxLength = 0
		for hmmerResult in self.hmmerResults:
			if hmmerResult.length>maxLength:
				maxlength=hmmerResult.length
			if len(hmmerResult.name)>maxTitleLength:
				maxTitleLength = len(hmmerResult.name)
			if maxLength< hmmerResult.length:
				maxLength= hmmerResult.length

		canvasWidth=int(maxlength*0.7)
		
		effectiveWidth = canvasWidth -leftMargin - rightMargin		
		#
		# In Python 2.x, int division by int is int. but it became float in Python 3.*. 
		#  from __future__ import division	was used 
		#
		conversion = effectiveWidth / maxLength
		#
		# If protein name width is bigger than leftMargin, wrote Label at top of proteins (self.titlemode=True)
		#
		if maxTitleLength*fontSize*conversion*0.7 > leftMargin:
			self.titlemode = True
		else:
			self.titlemode = False

		if self.titlemode:
			yDelta=120
		else:
			yDelta = 60	

		canvasHeight = len(self.hmmerResults)*yDelta+100
	
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
		colorIndex = 10

		gradientid = "test"
		doc = ET.Element('svg', width=str(canvasWidth), height=str(canvasHeight), version='1.1', xmlns='http://www.w3.org/2000/svg')
		doc.attrib["xmlns:xlink"]="http://www.w3.org/1999/xlink"
		defs = ET.Element('defs')
		gradientList = []
		for hmmerResult in self.hmmerResults:
			for hit in hmmerResult.hits:
				if hit.gradient:
					if not hit.color in gradientList:
						gradientid = "gradient_"+hit.color
						gradientList.append(gradientid)
						gradient = ET.Element('linearGradient', id=gradientid, x1="0%",y1="-20%",x2="0%",y2="120%")
						stop1 = ET.Element('stop', offset="0%", style="stop-color:white;stop-opacity:1")
						stop2 = ET.Element('stop', offset="40%", style="stop-color:"+hit.color+";stop-opacity:1")
						stop3 = ET.Element('stop', offset="60%", style="stop-color:"+hit.color+";stop-opacity:1")					
						stop4 = ET.Element('stop', offset="100%", style="stop-color:white;stop-opacity:1")
						gradient.append(stop1)
						gradient.append(stop2)
						gradient.append(stop3)
						gradient.append(stop4)	
						defs.append(gradient)

		doc.append(defs)
		
		for hmmer in self.hmmerResults:
			# Draw Protein Text
			if len(hmmer.source)>0 and len(hmmer.accession)>0:

				if hmmer.source == 'refseq':
					linkAddress = "http://www.ncbi.nlm.nih.gov/protein/{0}"
				if hmmer.source == 'uniprot':
					linkAddress = "http://www.uniprot.org/uniprot/{0}"


				anchor = ET.Element('a')
		 	 	anchor.attrib["xlink:href"]="#"+hmmer.accession
		 	 	doc.append(anchor)

				link = ET.Element('a')
		 		link.attrib["xlink:href"]=linkAddress.format(hmmer.accession)


			if self.titlemode:
				text = ET.Element('text', x=str(leftMargin), y=str(int(y-fontSize*2.5)), fill='black', 
									style='font-family:Sans-Serif;font-size:16px;text-anchor:left;dominant-baseline:middle')
	 		else:
	 			text = ET.Element('text', x=str(x), y=str(y), fill='black', 
									style='font-family:Sans-Serif;font-size:16px;text-anchor:left;dominant-baseline:middle')	 		
	 		text.text = hmmer.name
	 		text.attrib["id"]=hmmer.accession
		 		
	 		if len(hmmer.source)>0 and len(hmmer.accession)>0:
	 			link.append(text)
	 			doc.append(link)
	 		else:
	 			doc.append(text)

	 		# Draw Line
	 		line = ET.Element('line', x1=str(leftMargin), y1=str(y), x2=str(leftMargin+int(hmmer.length*conversion)),
	 							 y2=str(y), style='stroke:rgb(200,200,200);stroke-width:4')
	 		doc.append(line)
	 		# Start and End Amino Acid Number
	 		start = ET.Element('text', x=str(leftMargin-fontSize), y=str(y), fill='black', 
								style='font-family:Sans-Serif;font-size:13px;text-anchor:right;dominant-baseline:middle')
	 		start.text = '1'
	 		doc.append(start)
	 		end = ET.Element('text', x=str(leftMargin+int(hmmer.length*conversion)), y=str(y), fill='black', 
								style='font-family:Sans-Serif;font-size:13px;text-anchor:left;dominant-baseline:middle')
	 		end.text = str(hmmer.length)
	 		doc.append(end)
	 		#
	 		# Draw Domains
	 		#
	 		for hit in hmmer.hits:
	 			if not hit.exclude:
	 				#
	 				# if color of domain is not assigned, choose color in table and increase it.
	 				#
	 				if not hit.color:
			 			if hit.name in hmmcolors:
			 				color = hmmcolors[hit.name]
			 			else:
			 				hmmcolors[hit.name]=colors[colorIndex]
			 				color = colors[colorIndex]
			 				if colorIndex < len(colors)-1:
			 					colorIndex+=1
			 				else:
			 					colorIndex = 0
			 		else:
			 			color = hit.color

			 		if hit.border :
			 			border = ';stroke-width:1;stroke:black'
			 		else:
			 			border = ''
			 		if hit.gradient:
			 			style='fill:url(#'+'gradient_'+hit.color+')'+border
			 		else:
			 			style = 'fill:'+color+border
		 			#
		 			# Draw rectanglar domains
		 			#
		 			link = ET.Element('a')
		 			link.attrib["xlink:href"]="http://pfam.sanger.ac.uk/family/{0}".format(hit.acc)
		 			rect = ET.Element('rect', x=str(leftMargin+int(hit.start*conversion)), y=str(y-boxHeight/2),
		 								 width=str(int((hit.end - hit.start)*conversion)), 
		 								 height=str(boxHeight), style=style)
		 			
		 			if hit.acc:
		 				link.append(rect)
		 				doc.append(link)
		 			else:
		 				doc.append(rect)

		 			#
		 			# Draw Domain Label
		 			#
		 			if hit.label:
			 			font = 13
			 			if (len(str(hit.name))*font*0.6>(hit.end - hit.start)*conversion):
							delta = -1*boxHeight
						else:
							delta = 0
						link = ET.Element('a')
		 				link.attrib["xlink:href"]="http://pfam.sanger.ac.uk/family/{0}".format(hit.acc)
			 			textLabel = ET.Element('text',x=str(leftMargin+int((hit.start+(hit.end-hit.start)*0.5)*conversion)),y=str(y+delta),fill='black',
			 									style='font-family:Sans-Serif;font-size:'+str(font)+'px;text-anchor:middle;alignment-baseline:middle')
		 				textLabel.text = hit.name
		 				if hit.acc:
		 					link.append(textLabel)
		 					doc.append(link)
		 				else:
		 					doc.append(textLabel)
		 			#
		 			# Draw start and end aa numbers of the domain
		 			#
		 			# Adjust location of number based on the domain length
		 			#
		 			font=9
		 			if (leftMargin+hit.start*conversion+font*0.6*(len(str(hit.start))))>(leftMargin+hit.end*conversion-len(str(hit.end))*font*0.6):
		 				deltaStart = int(font*-0.5*len(str(hit.start)))
		 				deltaEnd = int(font*0.5*len(str(hit.end)))
		 			else:
		 				deltaStart = 0
		 				deltaEnd = 0

			 		if hit.startshow:	
			 			hitStart = ET.Element('text', x=str(leftMargin+int(hit.start*conversion)+deltaStart), y=str(y+boxHeight), fill='black', 
										style='font-family:Sans-Serif;font-size:'+str(font)+'px;text-anchor:left;dominant-baseline:top')
				 		hitStart.text = str(hit.start)
				 		doc.append(hitStart)

					if hit.endshow:		
				 		hitEnd = ET.Element('text', x=str(leftMargin+int(hit.end*conversion)-len(str(hit.end))*font*0.6+deltaEnd), y=str(y+boxHeight), fill='black', 
											style='font-family:Sans-Serif;font-size:'+str(font)+'px;text-anchor:right;dominant-baseline:top')
				 		hitEnd.text = str(hit.end)
				 		doc.append(hitEnd)

	 		y+=yDelta

	 	f = open(self.outputfile, 'w')
		f.write('<?xml version=\"1.0\" standalone=\"no\"?>\n')
		f.write('<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n')
		f.write('\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n')
		f.write('\n'.join(ET.tostringlist(doc)))
		f.close()
		return('\n'.join(ET.tostringlist(doc)))

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
							for hit in hmmerResult.hits:
							#	print hit.name, domain.domainName, hit.start, domain.start, hit.end, domain.end
								if (hit.name == domain.domainName) and (hit.start == domain.start) and (hit.end == domain.end):
							#		print "excluded"
									hit.exclude = True	
						if domain.action == 'new':
							pseudoHit = HmmerHit(name=domain.domainName, start=domain.start, end=domain.end)
							hmmerResult.hits.append(pseudoHit)
			
			for renameDef in self.inputFile.renameDefinitions:
				for hmmerResult in self.hmmerResults:
					for hit in hmmerResult.hits:
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

	def run(self):

		if self.inputFileName:
			self.inputFile=InputFile(inputFileName=self.inputFileName)
			self.inputFile.readInputFile()
			files = self.inputFile.fileNames
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
			svg=self.drawSVG()
			header = ['Accession','Name','Domain','length']
			table = HTMLTable(header = header)
			for hmmerResult in self.hmmerResults:
				if hmmerResult.source == 'refseq':
					linkAddress = "<a id='{0}' href='http://www.ncbi.nlm.nih.gov/protein/{0}'>{0}</a>"
				if hmmerResult.source == 'uniprot':
					linkAddress = "<a id='{0}' href='http://www.uniprot.org/uniprot/{0}'>{0}</a>"
				accession = linkAddress.format(hmmerResult.accession)
				name = "<a href='#{1}'>{0}</a>".format(hmmerResult.name,hmmerResult.accession)
				domain = []

				for hit in hmmerResult.hits:
					domainName="<a href='http://pfam.sanger.ac.uk/family/{0}'>{1}</a>".format(hit.acc, hit.name)
					domain.append(domainName)
				domains=','.join(domain)
				length = hmmerResult.length			
				table.tableContentFill([accession, name, domains, length])
			table.extra = "<div>"+svg+"</div>"
			table.tableGenerate('test.html')

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
	parser.add_argument('-o', '--output', dest='outputfile', default='output.svg', 
						help='Output svg filename')
	parser.add_argument('-t', '--no_threshold', dest='threshold',action='store_false', default=True,
						help='Turn of Pfam gathering threshold. Enable to look up more weak(unreliable) domains')
	parser.add_argument('-i', '--input_file', dest='inputFileName', default='hmmer.INP',
						help='Read configuration file')
	results = parser.parse_args()
	if not os.path.exists(results.inputFileName):
		results.inputFileName = None
	hmmerscan = HmmerScanRunner(files=results.files,inputFileName=results.inputFileName, db=results.db, evalue=results.evalue,
								local = results.local, outputfile = results.outputfile, threshold=results.threshold )
	hmmerscan.run()
	

