import numpy as np
from collections import Counter, defaultdict
from dudes.Ranks import Ranks
from dudes.Util import *


class Ident:
    columns = [
        "Iter",
        "TaxID",
        "ParentTaxID",
        "RankID",
        "CumulativeMatches",
        "MatchScoreSum",
        "Abundance",
        "CumulativeAbundance",
    ]

    def __init__(self, mat=None, pv=None):
        """
        Hold identifications. Contains an array with the following columns:
        - Iter: iteration number
        - TaxID: taxon ID
        - ParentTaxID: taxon ID of the parent of TaxID
        - RankID: ID of the Rank
        - CumulativeMatches: no idea
        - MatchScoreSum: no idea
        - Abundance: no idea
        - CumulativeAbundance: no idea

        :param mat:
        :param pv:
        """
        self.cols = {c: i for i, c in enumerate(Ident.columns)}
        self.ident = []
        self.pvals = dict()
        if not mat is None:
            self.ident = np.array(mat, ndmin=2)
        if not pv is None:
            self.pvals = pv

    def __iter__(self):
        for v in self.ident:
            yield {c: v[i] for i, c in enumerate(Ident.columns)}

    def getIter(self):
        return self.ident[-1, 0]

    def getSize(self):
        return len(self.ident)

    def getCol(self, col):
        return self.ident[:, self.cols[col]]

    def getSubSet(self, ind):
        return Ident(self.ident[ind])

    def add(self, id):
        if id.getSize():
            if len(self.ident):
                self.ident = np.vstack((self.ident, id.ident))  # Add a Ident object
            else:
                self.ident = np.array(id.ident, ndmin=2)  # First entry
            self.pvals.update(id.pvals)  # Update p-val dict

    def remove(self, id):
        for r in id:
            # Reassign the parent node to the children
            self.ident[
                np.logical_and(
                    self.ident[:, 0] == r["Iter"], self.ident[:, 2] == r["TaxID"]
                ),
                2,
            ] = r["ParentTaxID"]
            # Remove entry
            self.ident = self.ident[
                ~np.logical_and(
                    self.ident[:, 0] == r["Iter"], self.ident[:, 1] == r["TaxID"]
                )
            ]
            # Remove pval
            if (r["Iter"], r["TaxID"]) in self.pvals:
                del self.pvals[(r["Iter"], r["TaxID"])]

    def setMatchScoreSum(self, iter, taxid, mssum):
        idx = np.logical_and(self.getCol("Iter") == iter, self.getCol("TaxID") == taxid)
        self.ident[idx, 5] = mssum

    def setAbundance(self, iter, taxid, ab):
        idx = np.logical_and(self.getCol("Iter") == iter, self.getCol("TaxID") == taxid)
        self.ident[idx, 6] = ab

    def setCumulativeAbundance(self, iter, taxid, ab):
        idx = np.logical_and(self.getCol("Iter") == iter, self.getCol("TaxID") == taxid)
        self.ident[idx, 7] = ab

    def getLeafs(self, taxid_start, iter=None):
        leaf_ident = Ident()
        if iter != None:
            sub_ident = self.ident[self.ident[:, 0] == iter]
        else:
            sub_ident = self.ident
        stack = list(sub_ident[sub_ident[:, 2] == taxid_start])
        while stack:
            node = stack.pop()
            children = list(sub_ident[sub_ident[:, 2] == node[1]])
            if children:
                stack.extend(children)
            else:
                leaf_ident.add(Ident(sub_ident[sub_ident[:, 1] == node[1]]))
        return leaf_ident

    def printTree(self, taxid_start, total_matches_filter, total_abundance_norm, names):
        tree = []
        if np.sum(self.ident) != 0:
            for it in range(self.getIter() + 1):
                tree.append(str(it))
                sub_ident = self.ident[self.ident[:, 0] == it]
                stack = list(sub_ident[sub_ident[:, 2] == taxid_start])
                while stack:
                    node = stack.pop()
                    tree.append(
                        "-%s %d %s %s- %d (%d) M:%d(%.6f) MSS:%d A:%.6f CA:%.6f %s %s"
                        % (
                            (node[3]) * "-",
                            node[3],
                            getNameRank(node[3]),
                            (len(Ranks.ranks) - node[3] - len(getNameRank(node[3])) + 9)
                            * "-",
                            node[1],
                            (node[2]),
                            node[4],
                            node[4] / total_matches_filter,
                            node[5],
                            node[6] / float(total_abundance_norm),
                            node[7] / float(total_abundance_norm),
                            Counter(self.pvals[(node[0], node[1])]),
                            names.getName(node[1]),
                        )
                    )
                    children = list(sub_ident[sub_ident[:, 2] == node[1]])
                    stack.extend(children)
        return tree
