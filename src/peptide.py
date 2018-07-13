
class Peptide:

    def __init__(self, _id, _protein, _description, _charge, _seq, _dat):
        self.id = _id
        self.protein = _protein
        self.description = _description
        self.seq = _seq
        self.charge = _charge
        self.dat = [float(x) for x in _dat]
        self.datStr = _dat
        self.nPeptides = len([x for x in self.dat if x != 0])
        self.outliersList = list()

    def rmZeros(self, _na = "NA"):
        for i in range(0, len(self.dat)):
            if(self.dat[i] == 0):
                self.datStr[i] = _na

    def normRatios(self, _den, _mult = 1):
        for i in range(0, len(self.dat)):
            if self.dat[i] != 20:
                self.dat[i] = (self.dat[i] / _den[i]) * _mult
            if self.dat[i] != 0:
                self.datStr[i] = str(self.dat[i])

    def writeWide(self, outF, _includeDatType = True, _includePepSummary = True):
        if _includeDatType:
            outF.write("peptide\t")

        if _includePepSummary or not _includeDatType:
            outF.write(self.id +
                       '\t' + self.protein +
                       '\t' + self.description)
        else:
            outF.write('\t' * 2)

        outF.write('\t' + self.seq + '\t' + self.charge)

        for val in self.datStr:
            outF.write('\t' + val)

        ret = str()
        for i in range(0, len(self.datStr)):
            ret += '\t' + self.outliersList[i]

        return ret

    def writeLong(self, _index, outF, _sampleStr, _includeDatType = True):
        if _includeDatType:
            outF.write("peptide\t")

        outF.write(self.id +
                   '\t' + self.protein +
                   '\t' + self.description +
                   '\t' + self.seq +
                   '\t' + self.charge +
                   '\t' + _sampleStr + '\t' + self.datStr[_index])

        return self.outliersList[_index]
