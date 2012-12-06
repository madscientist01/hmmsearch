#!/usr/bin/python
import urllib2
import xml.etree.cElementTree as ET
from hmmer import HmmerHit
from fetchutil import readFasta, readAccession
from svgdrawer import SVGDrawer
import re
class SmartRedirectHandler(urllib2.HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        return headers

class UniprotAnnotation(object):

	def __init__(self, **kwargs):

		self.file=kwargs.get('file')
		self.name=kwargs.get('name')
		self.sequence = kwargs.get('sequence')
		self.source = kwargs.get('source')
		self.accession = kwargs.get('accession')
		self.length=kwargs.get('length')
		self.hits = []

	def readFile(self,fileName):
		(self.name, self.sequence) = readFasta(fileName)
		if self.sequence:
			self.length = len(self.sequence)
			self.file = fileName
		(self.source,self.accession)= readAccession(self.file)
		if self.source == "uniprot":
			self.uniprotXMLReader(self.accession)
		return

	def readAccession(self,accession):
		self.accession = accession
		self.uniprotXMLReader(self.accession)
		return

	def	uniprotXMLReader(self,uniprotId):
		opener = urllib2.build_opener(SmartRedirectHandler())
		urllib2.install_opener(opener)
		uniprotURL = "http://www.uniprot.org/uniprot/{0}.xml".format(uniprotId)
		# fileName = self.accession+".xml"
		stripnumber = re.compile('(.*)\s\d+$')
		try:
			req = urllib2.Request(uniprotURL)
			u = urllib2.urlopen(req)
			xml = u.read()
			# fileName = uniprotId+".xml"
			# f = open(fileName,'w')
			# f.write(xml)
			# f.close()
			root = ET.fromstring(xml)
			for child in root.iter('{http://uniprot.org/uniprot}entry'):
				gc =child.find('{http://uniprot.org/uniprot}sequence')
				self.length = int(gc.get('length'))
				self.sequence = gc.text

			for child in root.iter('{http://uniprot.org/uniprot}recommendedName'):
				gc =child.find('{http://uniprot.org/uniprot}fullName')
				self.name = gc.text	

			for child in root.iter('{http://uniprot.org/uniprot}feature'):
				# print child.tag, child.attrib
				featureType = child.get('type')
				description = child.get('description')

				if featureType == 'domain' or featureType == 'region of interest':
					
					match = stripnumber.match(description)
					if match:
						description = match.group(1)
					if featureType =='region of interest':
						featureType = 'feature'

					for begin in child.iter('{http://uniprot.org/uniprot}begin'):
						beginPos = begin.get('position')
					for end in child.iter('{http://uniprot.org/uniprot}end'):
						endPos = end.get('position')
					domainHit = HmmerHit(desc=featureType, name=description, start=beginPos, end=endPos)
					if int(endPos)-int(beginPos) < 35:
						domainHit.startshow=False
						domainHit.endshow=False
						domainHit.border=False
					if featureType == 'feature':
						domainHit.border = False
						domainHit.startshow=False
						domainHit.endshow=False
					self.hits.append(domainHit)	
			self.source = "uniprot"					
			return()

		except urllib2.URLError:
			print "{0} is not valid uniprot id".format(uniprotId)
	    	return()

if __name__ == "__main__":
	list = ["A1IGU5", "O00499", "P49418", "Q17R89", "Q5VWJ9", "Q68EM7", "Q6XZF7", "Q99961", "Q99962", "Q99963", "Q9NQY0", "Q9NR46", "Q9UBW5", "Q9UNH7", "Q9Y371", "Q9Y3L3"]
	uniprotResults=[]
	for item in list:
		uniprot = UniprotAnnotation()
		uniprot.readAccession(item)
		uniprotResults.append(uniprot)
		print uniprot.name, uniprot.accession, uniprot.length, uniprot.source
		for hit in uniprot.hits:
			print hit.desc, hit.name, hit.start, hit.end 
	drawer = SVGDrawer(outputSVG='uniprot.svg', hmmerResults = uniprotResults)
	drawer.drawSVG()



