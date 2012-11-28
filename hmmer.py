#!/usr/bin/python
#
# Hmmer3 based Protein Domain Visualizer 
#

from __future__ import division
import urllib
import urllib2
import os
import xml.etree.ElementTree as ET
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

def fetchfromUniprot(uniprotId):
	print uniprotId
	fileName = ""
	return(fileName)

class DomainDefinition(object):

	def __init__(self, **kwargs):
		self.proteinName = kwargs.get('proteinName')
		self.domainName = kwargs.get('domainName')
		self.accession = kwargs.get('acession')
		self.start = kwargs.get('start')
		self.end = kwargs.get('end')
		self.domainName2 = kwargs.get('domainName2')
		#
		# action = ['new','rename', 'remove']
		#
		self.action = kwargs.get('action')

class ColorDefinition(object):
	def __init__(self, **kwargs):
		self.domainName = kwargs.get('domainName')
		self.accession = kwargs.get('acession')
		self.color = kwargs.get('color')
		self.gradient = kwargs.get('gradient')
	

class InputFile(object):

	def __init__(self, **kwargs):
		self.domainDefinitions = []
		self.colorDefinitions = []
		self.inputFileName = kwargs.get('inputFileName')
		self.fileNames = {}
		self.searchDB = ""
		self.outputFileName = ""
		self.commands = {'FILE':'file','NEW':'new','RENAME':'rename','REMOVE':'remove', 'DOMAIN':'domain', 'SIZE':'size', 'EVALUE':'evalue', 'DB':'db'}

	def file(self,content):
		data = content.split()
		uniprot = re.compile("uniprot:(\S+)")
		if len(data)>1:
			name = data[0]
			fileName = data[1]
			match = uniprot.match(fileName)
			if match:
				fileName = fetchfromUniprot(match.group(1))
			if os.path.exists(fileName):
				self.fileNames[name]=fileNames
		return

	def new(self,content):
		data = content.split()
		if len(data) > 3:
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
				colorDef.gradient = data[5]

			domain = DomainDefinition(proteinName=proteinName, domainName=domainName, start=start, end=end)
			self.domainDefinitions.append(domain)
		return

	def rename(self,content):
		# print "RENAME"
		# print content
		return

	def remove(self,content):
		# print "REMOVE"
		# print content
		return


	def domain(self,content):
		# print "DOMAIN"
		# print content
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
							print keyword
							process=getattr(self, self.commands[keyword])
							process(content)

			for domain in self.domainDefinitions:
				print domain.proteinName, domain.domainName, domain.start, domain.end

			for colorDef in self.colorDefinitions:
				print colorDef.domainName, colorDef.color, colorDef.gradient
			# 	print name,fileName

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
		self.exclude = False;

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
	
	def exclude(self):

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
			
	def run(self):
		#
		# Using Hmmscan in Hmmer3 web service, find locations of domains in the Fasta sequence and store 
		#

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

			self.exclude()

			return(True)
		else:
		    print "Failed to retrieve results" 
		    return(False)


def drawSVG(hmmerResults, filename):
	#
	# Draw SVG based on the hmmer domaim
	#
	hmmcolors = {}
	canvasHeight = 400
	canvasWidth=1200
	colors = ['aliceblue','antiquewhite','aqua','aquamarine','azure','beige','bisque','black','blanchedalmond','blue','blueviolet','brown','burlywood','cadetblue','chartreuse','chocolate','coral','cornflowerblue','cornsilk','crimson','cyan','darkblue','darkcyan','darkgoldenrod','darkgray','darkgreen','darkgrey','darkkhaki','darkmagenta','darkolivegreen','darkorange','darkorchid','darkred','darksalmon','darkseagreen','darkslateblue','darkslategray','darkslategrey','darkturquoise','darkviolet','deeppink','deepskyblue','dimgray','dimgrey','dodgerblue','firebrick','floralwhite','forestgreen','fuchsia','gainsboro','ghostwhite','gold','goldenrod','gray','grey','green','greenyellow','honeydew','hotpink','indianred','indigo','ivory','khaki','lavender','lavenderblush','lawngreen','lemonchiffon','lightblue','lightcoral','lightcyan','lightgoldenrodyellow','lightgray','lightgreen','lightgrey']
	colorIndex = 10
	doc = ET.Element('svg', width=str(canvasWidth), height=str(canvasHeight), version='1.1', xmlns='http://www.w3.org/2000/svg')
	x = 50
	y = 50
	leftMargin = 150
	rightMargin = 100
	fontSize = 16
	effectiveWidth = canvasWidth -leftMargin - rightMargin
	boxHeight = 20
	maxLength = 0
	for hmmer in hmmerResults:
		if maxLength< hmmer.length:
			maxLength= hmmer.length
	#
	# In Python 2.x, int division by int is int. but it became float in Python 3.*. 
	#  from __future__ import division	was used 
	#

	conversion = effectiveWidth / maxLength

	for hmmer in hmmerResults:
		text = ET.Element('text', x=str(x), y=str(y), fill='black', 
							style='font-family:Sans-Serif;font-size:16px;text-anchor:left;dominant-baseline:middle')
 		text.text = hmmer.name
 		doc.append(text)
 		line = ET.Element('line', x1=str(leftMargin), y1=str(y), x2=str(leftMargin+int(hmmer.length*conversion)),
 							 y2=str(y), style='stroke:rgb(200,200,200);stroke-width:4')
 		doc.append(line)
 		start = ET.Element('text', x=str(leftMargin-fontSize), y=str(y), fill='black', 
							style='font-family:Sans-Serif;font-size:13px;text-anchor:right;dominant-baseline:middle')
 		start.text = '1'
 		doc.append(start)
 		end = ET.Element('text', x=str(leftMargin+int(hmmer.length*conversion)), y=str(y), fill='black', 
							style='font-family:Sans-Serif;font-size:13px;text-anchor:left;dominant-baseline:middle')
 		end.text = str(hmmer.length)
 		doc.append(end)

 		for hit in hmmer.hits:
 			if not hit.exclude:
	 			if hit.acc in hmmcolors:
	 				color = hmmcolors[hit.acc]
	 			else:
	 				hmmcolors[hit.acc]=colors[colorIndex]
	 				color = colors[colorIndex]
	 				if colorIndex < len(colors)-1:
	 					colorIndex+=1
	 				else:
	 					colorIndex = 0

	 			rect = ET.Element('rect', x=str(leftMargin+int(hit.start*conversion)), y=str(y-boxHeight/2),
	 								 width=str(int((hit.end - hit.start)*conversion)), 
	 								 height=str(boxHeight), style='fill:'+color+';stroke-width:1;stroke:black')

	 			font = 13
	 			if (len(str(hit.name))*font*0.6>(hit.end - hit.start)*conversion):
					delta = -1*boxHeight
				else:
					delta = 0

	 			textLabel = ET.Element('text',x=str(leftMargin+int((hit.start+(hit.end-hit.start)*0.5)*conversion)),y=str(y+delta),fill='black',
	 									style='font-family:Sans-Serif;font-size:'+str(font)+'px;text-anchor:middle;alignment-baseline:middle')

	 			textLabel.text = hit.name
	 			doc.append(rect)
	 			doc.append(textLabel)
	 			font=9

	 			if (leftMargin+hit.start*conversion+font*0.6*(len(str(hit.start))))>(leftMargin+hit.end*conversion-len(str(hit.end))*font*0.6):
	 				deltaStart = int(font*-0.5*len(str(hit.start)))
	 				deltaEnd = int(font*0.5*len(str(hit.end)))
	 			else:
	 				deltaStart = 0
	 				deltaEnd = 0
	 			hitStart = ET.Element('text', x=str(leftMargin+int(hit.start*conversion)+deltaStart), y=str(y+boxHeight), fill='black', 
								style='font-family:Sans-Serif;font-size:'+str(font)+'px;text-anchor:left;dominant-baseline:top')
		 		hitStart.text = str(hit.start)
		 		doc.append(hitStart)
		 		hitEnd = ET.Element('text', x=str(leftMargin+int(hit.end*conversion)-len(str(hit.end))*font*0.6+deltaEnd), y=str(y+boxHeight), fill='black', 
									style='font-family:Sans-Serif;font-size:'+str(font)+'px;text-anchor:right;dominant-baseline:top')
		 		hitEnd.text = str(hit.end)
		 		doc.append(hitEnd)

 		y+=60

 	f = open(filename, 'w')
	f.write('<?xml version=\"1.0\" standalone=\"no\"?>\n')
	f.write('<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n')
	f.write('\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n')
	f.write(ET.tostring(doc))
	f.close()
	return

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
	drawSVG(hmmerResults,"data.svg")

if __name__ == "__main__":

	inputFile = InputFile(inputFileName='hmmer.inp')
	inputFile.readInputFile()

	# parser = argparse.ArgumentParser()
	# parser.add_argument('-f', '--fasta', nargs='+', dest='files',default=[],
	#                     help='Files to process')
	# parser.add_argument('-d', '--database', dest='db',default='pfam',
	# 					help='HMM database')
	# parser.add_argument('-e', '--evalue_cutoff', dest='evalue',default=1e-5,
	# 					help='E-value cutoff')
	# results = parser.parse_args()
	# main(results)

