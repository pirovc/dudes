from dudes.Ranks import Ranks
import numpy as np
import sys

def printDebug(DEBUG, l):
	if DEBUG: sys.stderr.write(str(l) + "\n")
		
def group_max(groups, data, pre_order=None): 
	if pre_order is None:	
		order = np.lexsort((data, groups))
	else: 
		order = pre_order
	groups = groups[order] #this is only needed if groups is unsorted
	data = data[order]
	index = np.empty(len(groups), 'bool')
	index[-1] = True
	index[:-1] = groups[1:] != groups[:-1]
	if pre_order is None:
		return order, index
	else:
		return index # Return the data array in an orderer way (matching the output of np.unique(groups))

def getNameRank(rankid):
	# Returns the fixed ranks based on rankid
	if rankid<len(Ranks.ranks):
		return Ranks.ranks[rankid]
	else:
		return Ranks.ranks[-1] # more than one no_rank/strain

def getIndexRank(rank):
	# Returns the fixed ranks based on rankid
	return Ranks.ranks.index(rank)