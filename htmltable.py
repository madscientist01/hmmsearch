#!/usr/bin/python
#
# Simple HTML Table Generator
#
# Written by Madscientist (http://madscientist.wordpress.com, https://github.com/madscientist01/)
#

class HTMLTable(object):

	def __init__(self,**kwargs):

		self.htmlHeader = """
<!DOCTYPE html>
<html>
	<head>
    	<title>List</title>
 		<meta charset='utf-8'>
	</head>
	<script type="text/javascript" charset="utf8" src="{0}"></script>
"""
		self.include = """
		<!-- DataTables CSS -->
		<link rel="stylesheet" type="text/css" href="http://ajax.aspnetcdn.com/ajax/jquery.dataTables/1.9.4/css/jquery.dataTables.css">
	 
		<!-- jQuery -->
		<script type="text/javascript" charset="utf8" src="http://ajax.aspnetcdn.com/ajax/jQuery/jquery-1.8.2.min.js"></script>
	 
		<!-- DataTables -->
		<script type="text/javascript" charset="utf8" src="http://ajax.aspnetcdn.com/ajax/jquery.dataTables/1.9.4/jquery.dataTables.min.js"></script>
"""	
		self.scriptheader = """
		<script>
		$(document).ready(function(){
			$('#listtable').dataTable();
"""
		self.scriptfooter = """ 		
		});
		</script>
"""
		self.scriptcontent = ""
		self.style = """
		<style>
		table{
		    font-family: "Arial",Sans-Serif;
		    font-size: 12px;
		    #margin: 45px;
		    width:1000px;
		    text-align: left;
		    border-collapse: collapse;  
			}
		div.head {
			width:800px;
			font-family: Sans-Serif;
			font-size: 14px;
			border:3px solid #EEEEEE;
			border-radius: 10px;
			padding: 10px;
		    align :center;
			background-color: #FFFFFF;
			}
		div {
			margin: 0px 0px 0px 0px
		}		
		</style>
		<body>
"""
		self.tablefooter = "</table>"
		self.htmlfooter = "</body></html>"
		self.tableContentTemplate = ""
		self.tableContent=""
		self.tableHeader=""
		self.extra=""
		self.header = kwargs.get('header')
		self.colNo = len(self.header)
		self.tableHeaderGenerate()

	def tableHeaderGenerate(self):
		self.tableHeader = "<br><br><table id='listtable'><thead><tr>" +\
						''.join(["<td>"+i+"</td>" for i in self.header])+\
						"</tr></thead><tbody>"	
		self.tableContentTemplate =	"<tr>"+\
					''.join(["<td>{"+str(i)+"}</td>" for i in range(self.colNo)])+\
					"</tr>"

	def tableContentFill(self,tableContent):
		fill = self.tableContentTemplate.format(*tableContent)
		self.tableContent = self.tableContent+fill


	def tableGenerate(self,filename):
		table = self.htmlHeader+self.include+self.scriptheader+self.scriptcontent+self.scriptfooter+self.style+self.tableHeader+self.tableContent+self.tablefooter+self.extra+self.htmlfooter
		f = open(filename,'w')
		f.write(table)
		f.close()
		return()
