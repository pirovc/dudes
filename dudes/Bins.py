from collections import defaultdict
from typing import Tuple

import numpy as np


class Bins:
    def __init__(self, smap, bin_size):
        """
        Create dictionary, where key: RefID, value: bin dict.
        Bin dict key: bin index no, value: tuple(binscore, matches in bin)
        The bin score is the sum the scores of the matches in this bin.

        :param smap: SMap object
        :param bin_size: int
        """
        self.bins = defaultdict(lambda: defaultdict(lambda: (0, 0)))
        # Sort smap by RefID
        order = np.argsort(smap.getCol("RefID"))
        sort_smap = smap.getSubSet(order)
        # Identify where the RefID changes (indices of RefIDs before change, e.g. change_idx([1,1,2]) = array([1])
        change_idx = np.where(
            sort_smap.getCol("RefID")[:-1] != sort_smap.getCol("RefID")[1:]
        )[0]
        # For each set of matches for each refid
        for ref in np.split(sort_smap.smap, change_idx + 1):
            refid = ref[0, 0]
            # Count how many matches each bin had:
            # array where index=bin index and value = how many matches in bin[i]
            count = np.bincount(ref[:, 1] // bin_size)
            # Bin score - Sum of the match scores
            binscore = np.bincount(ref[:, 1] // bin_size, weights=ref[:, 2])
            self.bins[refid] = dict(enumerate(zip(binscore, count)))

    def getBins(self, refids) -> Tuple[np.array, int]:
        """
        get array of bin scores and total number of matches for query refids

        :param refids: iterable of RefIDs
        :return: tuple(array(bin_scores), total number of matches in refids)
        """
        b = []
        m = 0
        for r in refids:
            bins_matches = list(zip(*list(self.bins[r].values())))
            if bins_matches:
                b.extend(bins_matches[0])
                m += sum(bins_matches[1])
        return np.array(b), m
