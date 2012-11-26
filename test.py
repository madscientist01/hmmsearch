#!/usr/bin/python
from reportlab.pdfgen import canvas
from reportlab.lib.colors import *
import re
import numpy as np


rate4site_filename="r4s.res"
c=canvas.Canvas("Test.pdf")
c.setFont("Courier",14)
c.setStrokeColor(red)
#
#
#
#
f=open(rate4site_filename)
filecontent=f.readlines()
searchstring="#The alpha parameter"
match = False
StartRegex = re.compile("^" + searchstring)
LineRegex = re.compile("^(\d+)\s+(\w)\s+([\w\.\-]+)")                
seq=[]
dat=[]
res=[]
for line in filecontent:

    if not match :
	match = StartRegex.match(line)
    else :
	linedata = LineRegex.match(line):
	if linedata:
	    res.append(int(linedata.group(1)))
	    seq.append(linedata.group(2))
	    dat.append(float(linedata.group(3)))

array_dat = np.array(dat)
diff=array_dat.max()-array_dat.min()
normalized = 1-(array_dat-array_dat.min())/diff
print normalized
        

width=50
cursor = 0
end =0
widthstep=1
heightstep=1
xstart=50
ystart=750
xdel=10
ydel=15

while (cursor < len(seq)):
	if (len(seq)-cursor) > width:
		end+=width
	else:
		end=len(seq)
	
	i=0
	c.setStrokeColorRGB(0,0,0)
	c.setFillColorRGB(0,0,0)
	c.drawString(xstart-30,ystart+ydel,str(res[cursor]))
	for oneamino in seq[cursor:end]:
                Strokecolor=Color(0,normalized[cursor+i]*0.8,(normalized[cursor+i]), alpha=1)
                TransparentColor=Color(0,normalized[cursor+i]*0.8,(normalized[cursor+i]),alpha=0.7)
                c.setStrokeColor(Strokecolor)
		c.setFillColor(Strokecolor)	
		c.drawString(xstart+i*xdel,ystart+ydel,oneamino)
                c.setFillColor(TransparentColor)
		c.rect(xstart+i*xdel,ystart+ydel+10,xdel,40*normalized[cursor+i],fill=True)
		i+=1		
	cursor+=width
	ydel-=60

c.showPage()
c.save()
