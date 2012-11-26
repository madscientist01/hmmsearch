

#!/usr/bin/python
#
# Alignment Colorization using rate4site score 
# Require Biopython
#		  
#		  ReportLab
#
#
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.colors import *
from Bio import AlignIO
import re

alignment = AlignIO.read("CH2_Domains.aln", "clustal")
c=canvas.Canvas("Test.pdf", pagesize=landscape(letter))
c.setFont("Courier",13)
c.setStrokeColor(red)
titleregex = re.compile("^([a-zA-Z\-]+)")

#
#
# Read Rate4Site File
# parsed with awk '!/^#/ {print $1,"\t" $2 "\t" $3}' r4s.res >data.res
#

f = open("data.res")
filecontent = f.readlines()
amino = []
residues = []
sum = 0
max = 0
min = 0
#
# Normalize Rate4Site Values as  0..1 ranges
# 
#
for line in filecontent:
	if line.strip():
		res,am,nu = line.split()
		num = float(nu)
		amino.append(am)
		residues.append(num)
		if num>=max :
			max = num
		if num<=min :
			min = num
		sum+=num


for k in xrange(len(residues)):
	residues[k] = (residues[k]+(-1)*min)/(max-min)
	
width=70
cursor = 0
end =0
widthstep=1
heightstep=1
xstart=150
ystart=500
xtitlemargin = 80
xlmargin = 20
xrmargin = 5
xdel=8
ydel=15
cursor=0
fontsize = 14


length = alignment.get_alignment_length()

linenumber = []
for record in alignment:
	linenumber.append(1)

score = []
residue = 0
for singleAminoAcid in alignment[0].seq:
	if singleAminoAcid == "-":
		score.append(10)
	elif (singleAminoAcid == amino[residue]) :
		score.append (residues[residue])
		residue+=1


while (cursor < length):
	if (length-cursor) > width:
		end+=width
	else:
		end=length
	j=0
	for record in alignment:
		sequence = record.seq
		i=1
		
:
		match = titleregex.match(record.id)
		if match:
			title = match.group(0)
			c.setFillColorRGB(0,0,0)
			c.drawString(xstart-xlmargin, ystart+ydel,str(linenumber[j]))
			c.setFillColorRGB(0,0,0)
			c.drawString(xstart-xlmargin-xtitlemargin, ystart+ydel,title)
			
		
		for oneamino in sequence[cursor:end]:
			if oneamino!="-":
				print score[cursor+i-1], cursor+i-1	
				c.setFillColorRGB(0,0.3+score[cursor+i-1], 0.3+score[cursor+i-1])
				c.setStrokeColorRGB(0,0.3+score[cursor+i-1],0.3+score[cursor+i-1])
				c.rect(xstart+i*xdel,ystart+ydel,xdel,fontsize,fill=True)
				linenumber[j]+=1
			if score[cursor+i-1]>0.3:
				
				c.setFillColorRGB(0,0,0)
			else:
				c.setFillColorRGB(0.7,0.7,0.7)	
			c.drawString(xstart+i*xdel,ystart+ydel,oneamino)			
			i+=1

		c.setFillColorRGB(0,0,0)
		c.drawString(xstart+i*xdel+xrmargin, ystart+ydel,str(linenumber[j]-1))
		
		ydel-=15
		j+=1	
		
		if ystart+ydel <100:
			c.showPage()
			c.setFont("Courier",fontsize)
			ystart=720
			ydel=15

	cursor+=width
	ydel-=15
	if ystart+ydel <100:
			c.showPage()
			c.setFont("Courier",fontsize)
			ystart=720
			ydel=15

c.showPage()
c.save()

#	withoutgap = re.sub("-", "",str(record.seq))
#	print len(record.seq), len(withoutgap)
#	print record.seq
#	print withoutgap	