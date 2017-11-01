from collections import defaultdict
import numpy as np

class Bins:
	def __init__(self,smap,bin_size):
		self.bins = defaultdict(lambda : defaultdict(lambda:(0,0)))
		#Sort smap by RefID
		order = np.argsort(smap.getCol('RefID'))
		sort_smap = smap.getSubSet(order)
		#Identify where the RefID changes
		change_idx = np.where(sort_smap.getCol('RefID')[:-1] != sort_smap.getCol('RefID')[1:])[0]
		#For each set of matches for each refid
		for ref in np.split(sort_smap.smap,change_idx+1):
			refid = ref[0,0]
			# Count how many matches each bin had
			count = np.bincount(ref[:,1]//bin_size)
			# Bin score - Sum of the match scores
			binscore = np.bincount(ref[:,1]//bin_size, weights=ref[:,2])
			self.bins[refid] = dict(enumerate(zip(binscore,count)))
	
	def getBins(self,refids):
		b = []
		m = 0
		for r in refids:
			bins_matches = list(zip(*list(self.bins[r].values())))
			if bins_matches:
				b.extend(bins_matches[0])
				m += sum(bins_matches[1])
		return np.array(b),m
