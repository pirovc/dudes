class SMap:
    columns = ["RefID", "MatchPosStart", "MatchScore", "ReadID"]

    def __init__(self, mat):
        self.cols = {c: i for i, c in enumerate(SMap.columns)}
        self.smap = mat
        self.index = 0

    def getSize(self):
        return self.smap.shape[0]

    def getCol(self, col):
        """
        get column by column name
        :param col: str, valid values: ['RefID','MatchPosStart','MatchScore','ReadID']
        :return: single column init data
        """
        return self.smap[:, self.cols[col]]

    def getSubSet(self, ind):
        """
        get new SMap object, initialized with subset rows 'ind' of self.smap

        :param ind: row numbers
        :return:
        """
        return SMap(self.smap[ind])

    def __iter__(self):
        """
        iter over rows

        :return: iter(dict of col names: values)
        """
        for v in self.smap:
            yield {c: v[i] for i, c in enumerate(SMap.columns)}

    def next(self):
        if self.index == self.getSize():
            raise StopIteration
        d = {c: self.smap[self.index, i] for i, c in enumerate(SMap.columns)}
        self.index = self.index + 1
        return d
