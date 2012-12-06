import xml.etree.cElementTree as ET

class SVGDrawer(object):
	#
	# Class for the Hmmer domain hits
	#
	def __init__(self, **kwargs):
		self.outputHTML = kwargs.get('outputHTML')
		self.outputSVG = kwargs.get('outputSVG')
		self.hmmerResults = kwargs.get('hmmerResults')
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
		
		
		conversion = float(effectiveWidth) / float(maxLength)
		# If protein name width is bigger than leftMargin, wrote Label at top of proteins (self.titlemode=True)
		if maxTitleLength*fontSize*conversion*0.7 > leftMargin:
			self.titlemode = True
		else:
			self.titlemode = False

		if self.titlemode:
			yDelta=120
			leftMargin = 20
			effectiveWidth = canvasWidth -leftMargin - rightMargin
			conversion = float(effectiveWidth) / float(maxLength)
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
		conversion = float(effectiveWidth) / float(maxLength)
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
			conversion = float(effectiveWidth) / float(maxLength)
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
					gradientid = "gradient_"+hit.color
					if not gradientid in gradientList:
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

	 			if (hit.desc and hit.desc=="feature"):
	 				rectYPos = int(y+boxHeight*1.3)
	 				rectHeight = int(boxHeight/3)
	 				labelYPos = int(y+boxHeight*2.2)
	 				
	 			else:
	 				rectYPos = y-boxHeight/2
	 				rectHeight = boxHeight
	 				labelfont = 13
		 			if (len(str(hit.name))*labelfont*0.5>(hit.end - hit.start)*conversion):
						labelYPos = int(y-0.8*boxHeight)
					else:
						labelYPos = y+int (boxHeight/3.5)
					numberYPos = y+boxHeight
	 			
	 			rect = ET.Element('rect', x=str(leftMargin+int(hit.start*conversion)), y=str(rectYPos),
	 								 width=str(int((hit.end - hit.start)*conversion)), 
	 								 height=str(rectHeight), style=style)
	 			
	 			if hit.acc:
	 				link = ET.Element('a')
	 				link.attrib["xlink:href"]="http://pfam.sanger.ac.uk/family/{0}".format(hit.acc)
	 				link.append(rect)
	 				doc.append(link)
	 			else:
	 				doc.append(rect)
	 			#
	 			# Draw Domain Label
	 			#
	 			if hit.label:
		 			textLabel = ET.Element('text',x=str(leftMargin+int((hit.start+(hit.end-hit.start)*0.5)*conversion)),y=str(labelYPos),fill='black',
		 									style='font-family:Sans-Serif;font-size:'+str(labelfont)+'px;text-anchor:middle')
	 				textLabel.text = hit.name
	 				if hit.acc:
	 					link = ET.Element('a')
	 					link.attrib["xlink:href"]="http://pfam.sanger.ac.uk/family/{0}".format(hit.acc)
	 					link.append(textLabel)
	 					doc.append(link)
	 				else:
	 					doc.append(textLabel)
	 			#
	 			# Draw start and end aa numbers of the domain
	 			#
	 			# Adjust location of number based on the domain length
	 			#
	 			numberFont=9
	 			if (leftMargin+hit.start*conversion+numberFont*0.6*(len(str(hit.start))))>(leftMargin+hit.end*conversion-len(str(hit.end))*numberFont*0.6):
	 				deltaStart = int(numberFont*-0.5*len(str(hit.start)))
	 				deltaEnd = int(numberFont*0.5*len(str(hit.end)))
	 			else:
	 				deltaStart = 0
	 				deltaEnd = 0
	 			startEndLink = ET.Element('a')
	 			startEndLink.attrib["xlink:href"]="http://www.uniprot.org/blast/?about={0}[{1}-{2}]".format(hmmer.accession, hit.start, hit.end)
	 			
		 		if hit.startshow:	
		 			hitStart = ET.Element('text', x=str(leftMargin+int(hit.start*conversion)+deltaStart), y=str(numberYPos), fill='black', 
									style='font-family:Sans-Serif;font-size:'+str(numberFont)+'px;text-anchor:left;dominant-baseline:top')
			 		hitStart.text = str(hit.start)
			 		if hmmer.source == 'uniprot':
						startEndLink.append(hitStart)
						doc.append(startEndLink)
					else:
			 			doc.append(hitStart)

				if hit.endshow:		
			 		hitEnd = ET.Element('text', x=str(leftMargin+int(hit.end*conversion)-len(str(hit.end))*numberFont*0.5+deltaEnd), y=str(numberYPos), fill='black', 
										style='font-family:Sans-Serif;font-size:'+str(numberFont)+'px;text-anchor:right;dominant-baseline:top')
			 		hitEnd.text = str(hit.end)
			 		if hmmer.source == 'uniprot':
						startEndLink.append(hitEnd)
						doc.append(startEndLink)
					else:
			 			doc.append(hitEnd)
		return(doc)


