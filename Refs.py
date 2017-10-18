class Refs:
	columns = ['RefID','SeqLen']
	def __init__(self,mat):
		self.cols = {c:i for i,c in enumerate(Refs.columns)}
		self.refs = mat
	def getSize(self):
		return self.refs.shape[0]
	def getCol(self,col):
		return self.refs[:,self.cols[col]]
	def getSubSet(self,ind):
		return Refs(self.refs[ind])
	def __iter__(self):
		for v in self.refs:
			yield  {c:v[i] for i,c in enumerate(Refs.columns)}