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
		self.tier = 0
		self.exclude = False
		self.color = None
		self.label = True
		self.border = True
		self.startshow = True
		self.endshow = True
		self.gradient = False