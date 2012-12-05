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
import glob
import sys
import subprocess
from htmltable import HTMLTable
from inputfile import InputFile
from fetchutil import readFasta, readAccession

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
		self.outputHTML = kwargs.get('outputHTML')
		self.outputSVG = kwargs.get('outputSVG')
		self.threshold = kwargs.get('threshold')
		self.hmmerResults = []
		self.titlemode = kwargs.get('titlemode',False)
	
	def maxLength(self):
		maxTitleLength = 0
		maxLength = 0
		for hmmerResult in self.hmmerResults:
			if hmmerResult.length>maxLength:
				maxLength=hmmerResult.length
			if len(hmmerResult.name)>maxTitleLength:
				maxTitleLength = len(hmmerResult.name)
			if maxLength< hmmerResult.length:
				maxLength= hmmerResult.length
		return(maxLength,maxTitleLength)

	def saveSVG(self,filename,doc):
		#
		# Save SVG based on the ElementTreeDoc
		#
		f = open(filename, 'w')
		f.write('<?xml version=\"1.0\" standalone=\"no\"?>\n')
		f.write('<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n')
		f.write('\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n')
		f.write('\n'.join(ET.tostringlist(doc)))
		f.close()
		return()

	def drawSVG(self):
		
		# Draw SVG based on the hmmer domaim
		x = 50
		y = 70
		leftMargin = 150
		rightMargin = 100
		fontSize = 16
		boxHeight = 20
		maxLength,maxTitleLength = self.maxLength()
		canvasWidth=int(maxLength*0.7)	
		effectiveWidth = canvasWidth -leftMargin - rightMargin		
		
		# In Python 2.x, int division by int is int. but it became float in Python 3.*. 
		#  from __future__ import division	was used 
		
		conversion = effectiveWidth / maxLength
		# If protein name width is bigger than leftMargin, wrote Label at top of proteins (self.titlemode=True)
		if maxTitleLength*fontSize*conversion*0.7 > leftMargin:
			self.titlemode = True
		else:
			self.titlemode = False

		if self.titlemode:
			yDelta=120
			leftMargin = 20
			effectiveWidth = canvasWidth -leftMargin - rightMargin
			conversion = effectiveWidth / maxLength
		else:
			yDelta = 60	

		canvasHeight = len(self.hmmerResults)*yDelta+100
		# doc is elementTree container for svg
		doc = ET.Element('svg', width=str(canvasWidth), height=str(canvasHeight), version='1.2', xmlns='http://www.w3.org/2000/svg')
		doc.attrib["xmlns:xlink"]="http://www.w3.org/1999/xlink"	
		# Color and Gradient Assignment
		defs = self.colorDef()
		doc.append(defs)
		
		#Draw several domain architecture in single SVG
		for hmmer in self.hmmerResults:
			doc = self.singleSVG(hmmer,doc,x,y,leftMargin,fontSize,conversion,boxHeight)
	 		y+=yDelta

	 	self.saveSVG(self.outputSVG,doc)
		return('\n'.join(ET.tostringlist(doc)))
	

	def drawMultiSVG(self):
		
		# Draw SVG based on the hmmer domaim
		x = 50
		y = 60
		leftMargin = 150
		rightMargin = 100
		fontSize = 16
		boxHeight = 20
		maxLength,maxTitleLength = self.maxLength()
		canvasWidth=int(maxLength*0.7)	
		canvasHeight = 110
		effectiveWidth = canvasWidth -leftMargin - rightMargin		
		conversion = effectiveWidth / maxLength
		svgFileNames = {}
		svgContent = {}
		# If protein name width is bigger than leftMargin, wrote Label at top of proteins (self.titlemode=True)
	
		if maxTitleLength*fontSize*conversion*0.7 > leftMargin:
			self.titlemode = True
		else:
			self.titlemode = False

		if self.titlemode:
			leftMargin = 20
			effectiveWidth = canvasWidth -leftMargin - rightMargin
			conversion = effectiveWidth / maxLength
		#
		# Draw several SVGs per each protein
		#
		for hmmer in self.hmmerResults:
			doc = ET.Element('svg', width=str(canvasWidth), height=str(canvasHeight), version='1.2', xmlns='http://www.w3.org/2000/svg')
			doc.attrib["xmlns:xlink"]="http://www.w3.org/1999/xlink"	
			# Color and Gradient Assignment
			defs = self.colorDef()
			doc.append(defs)
			doc = self.singleSVG(hmmer,doc,x,y,leftMargin,fontSize,conversion,boxHeight)
			if hmmer.accession:
				svgFileName = hmmer.accession+".svg"
			else:
				svgFileName = hmmer.name+".svg"
	 		self.saveSVG(svgFileName,doc)
	 		svgFileNames[hmmer.name]=svgFileName
	 		svgContent[hmmer.name]=ET.tostring(doc)
		return(svgFileNames,svgContent)


	def colorDef(self):
		
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
		hmmcolors = {}
		colorIndex = 10
		defs = ET.Element('defs')
		gradientList = []
		for hmmerResult in self.hmmerResults:
			for hit in hmmerResult.hits:
				if not hit.color:
		 			if hit.name in hmmcolors:
		 				hit.color = hmmcolors[hit.name]
		 			else:
		 				hmmcolors[hit.name]=colors[colorIndex]
		 				hit.color = colors[colorIndex]
		 				if colorIndex < len(colors)-1:
		 					colorIndex+=1
		 				else:
		 					colorIndex = 0
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
		return(defs)	

	def singleSVG(self,hmmer,doc,x,y,leftMargin,fontSize,conversion,boxHeight):
			# Draw Protein Text
		if len(hmmer.source)>0 and len(hmmer.accession)>0:

			if hmmer.source == 'refseq':
				linkAddress = "http://www.ncbi.nlm.nih.gov/protein/{0}"
			if hmmer.source == 'uniprot':
				linkAddress = "http://www.uniprot.org/uniprot/{0}"
			
			link = ET.Element('a')
	 		link.attrib["xlink:href"]=linkAddress.format(hmmer.accession)


		if self.titlemode:
			text = ET.Element('text', x=str(leftMargin-fontSize), y=str(int(y-fontSize*2.5)), fill='black', 
								style='font-family:Sans-Serif;font-weight:bold;font-size:16px;text-anchor:left;dominant-baseline:middle')
 		else:
 			text = ET.Element('text', x=str(x), y=str(y), fill='black', 
								style='font-family:Sans-Serif;font-weight:bold;font-size:16px;text-anchor:left;dominant-baseline:middle')	 		
 		text.text = hmmer.name
 		text.attrib["id"]=hmmer.accession
 		text.attrib["name"]=hmmer.accession
 				 		
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
						delta = int (boxHeight/3.5)
					link = ET.Element('a')
	 				link.attrib["xlink:href"]="http://pfam.sanger.ac.uk/family/{0}".format(hit.acc)
		 		# 	textLabel = ET.Element('text',x=str(leftMargin+int((hit.start+(hit.end-hit.start)*0.5)*conversion)),y=str(y+delta),fill='black',
		 		# 							style='font-family:Sans-Serif;font-size:'+str(font)+'px;text-anchor:middle;alignment-baseline:middle')
					textLabel = ET.Element('text',x=str(leftMargin+int((hit.start+(hit.end-hit.start)*0.5)*conversion)),y=str(y+delta),fill='black',
		 									style='font-family:Sans-Serif;font-size:'+str(font)+'px;text-anchor:middle')
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


	 			startEndLink = ET.Element('a')
	 			startEndLink.attrib["xlink:href"]="http://www.uniprot.org/blast/?about={0}[{1}-{2}]".format(hmmer.accession, hit.start, hit.end)
	 			

		 		if hit.startshow:	
		 			hitStart = ET.Element('text', x=str(leftMargin+int(hit.start*conversion)+deltaStart), y=str(y+boxHeight), fill='black', 
									style='font-family:Sans-Serif;font-size:'+str(font)+'px;text-anchor:left;dominant-baseline:top')
			 		hitStart.text = str(hit.start)
			 		if hmmer.source == 'uniprot':
						startEndLink.append(hitStart)
						doc.append(startEndLink)
					else:
			 			doc.append(hitStart)

				if hit.endshow:		
			 		hitEnd = ET.Element('text', x=str(leftMargin+int(hit.end*conversion)-len(str(hit.end))*font*0.6+deltaEnd), y=str(y+boxHeight), fill='black', 
										style='font-family:Sans-Serif;font-size:'+str(font)+'px;text-anchor:right;dominant-baseline:top')
			 		hitEnd.text = str(hit.end)
			 		if hmmer.source == 'uniprot':
						startEndLink.append(hitEnd)
						doc.append(startEndLink)
					else:
			 			doc.append(hitEnd)
		return(doc)


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
			self.db = self.inputFile.db
			self.threshold = self.inputFile.threshold
			self.local = self.inputFile.local
			self.evalue = self.inputFile.evalue
			self.outputHTML = self.inputFile.outputHTML
			print self.outputHTML
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
			self.hmmerResults.sort(key=lambda x:x.name)		
			self.drawSVG()
			(svgFileNames,svgContent) = self.drawMultiSVG()
			header = ['Accession','Name','Domain','length']
			table = HTMLTable(header = header)
			svgList = SVGList()
		
			for hmmerResult in self.hmmerResults:
				print hmmerResult.name
				if hmmerResult.source == 'refseq':
					linkAddress = "<a href='http://www.ncbi.nlm.nih.gov/protein/{0}'>{0}</a>"
				if hmmerResult.source == 'uniprot':
					linkAddress = "<a href='http://www.uniprot.org/uniprot/{0}'>{0}</a>"
				accessionLink = linkAddress.format(hmmerResult.accession)
				nameLink = "<a href='#{1}'>{0}</a>".format(hmmerResult.name,hmmerResult.accession)
				domain = []

				for hit in hmmerResult.hits:
					domainName="<a href='http://pfam.sanger.ac.uk/family/{0}'>{1}</a>".format(hit.acc, hit.name)
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
	

