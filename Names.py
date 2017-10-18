class Names:
	names = {}
	def __init__(self,names):
		self.names = names

	def getName(self,taxid):
		if self.names.get(taxid): 
			return self.names[taxid]
		else:
			return "No name"