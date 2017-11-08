class TTree:
	columns = ['RefID','TaxID','ParentTaxID','RankID']
	def __init__(self,mat):
		self.cols = {c:i for i,c in enumerate(TTree.columns)}
		self.ttree = mat
		self.leafs = [] #TTree with only leaf entries - set later
	def __iter__(self):
		for v in self.ttree:
			yield {c:v[i] for i,c in enumerate(TTree.columns)}
			
	def getSize(self):
		return self.ttree.shape[0]
	def getCol(self,col):
		return self.ttree[:,self.cols[col]]
	def getSubSet(self,ind):
		return TTree(self.ttree[ind])
	def setLeafs(self,ttree):
		self.leafs = ttree
	def getLeafs(self):
		return self.leafs