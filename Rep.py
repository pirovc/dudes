import numpy as np
from collections import  defaultdict
from Util import *

class Rep:
	columns = ['TaxID','ParentTaxID','RankID','Abundance']
	def __init__(self,mat=None):
		self.cols = {c:i for i,c in enumerate(Rep.columns)}
		self.rep = []
		if not mat is None: self.rep = np.array(mat,ndmin=2)
	def __iter__(self):
		for v in self.rep:
			yield {c:v[i] for i,c in enumerate(Rep.columns)}
	def getSize(self):
		return len(self.rep)
	def getCol(self,col):
		return self.rep[:,self.cols[col]]
	def getSubSet(self,ind):
		return Rep(self.rep[ind])	

	def setRepresentatives(self,ident,last_rank):
		# Sum abundances iterations
		abundances_iter = defaultdict(float)
		for id in ident: abundances_iter[id['TaxID']] += id['CumulativeAbundance']
				
		# Set representatives for output
		for ids in ident.getSubSet(np.unique(ident.getCol('TaxID'),True)[1]):
			# Chose single abundance for strain level and summed cumulative abundance for lower levels (because of the intra-strains)
			if ids['RankID']>=getIndexRank("strain"):
				ab = ids['Abundance']
			else:
				ab = abundances_iter[ids['TaxID']]
			self.rep.append([ids['TaxID'],ids['ParentTaxID'],ids['RankID'],ab])
		self.rep = np.array(self.rep,dtype=int)
		
		if last_rank=="strain": 
			# Fix the multi-level strains (bind them to the species level)
			for r in self.rep[self.rep[:,2]>getIndexRank("strain")]:
				stack = [r[1]] # parent_taxid
				while stack:
					parent = stack.pop()
					#if reached the species level, adjust entry
					if self.rep[self.rep[:,0]==parent,2]==getIndexRank("species"):
						self.rep[self.rep[:,0]==r[0],1] = parent
						self.rep[self.rep[:,0]==r[0],2] = getIndexRank("strain")
					else:
						stack.append(self.rep[self.rep[:,0]==parent,1])
		
		# Sort by rankid and abundance
		self.rep = self.rep[np.lexsort((-self.rep[:,3],self.rep[:,2]))]
		
		# Filter entries with 0 abundance (intra-strains without direct matches)
		self.rep = self.rep[self.rep[:,3]>0]