
import peptide
import numpy
import utils
import constants as const

class Protein:
    colnames = list()
    splitSamples = list()
    headerIndexMap = dict()
    verbose = True
    nCols = 0

    def __init__(self, _id, _descrption, _protein):
        self.id = _id
        self.protein = _protein
        self.description = _descrption
        self.peptides = list()
        self.mrList = list()
        self.avgrList = list()
        self.sdList = list()
        self.nOutliersDetected = list()
        self.nPeptidesList = list()
        self.nPeptides = 0

    @staticmethod
    def setColnames(_colnames):
        Protein.colnames = _colnames
        Protein.nCols = len(Protein.colnames)

    def addPeptide(self, _id, _protein, _description, _charge, _seq, _dat):
        self.peptides.append(peptide.Peptide(_id, _protein, _description, _charge, _seq, _dat))

    def addPeptideLine(self, elems):
        self.addPeptide(self.id,
                        self.protein,
                        self.description,
                        elems[self.headerIndexMap[const.CHARGE_COLNAME]],
                        elems[self.headerIndexMap[const.SEQ_COLNAME]],
                        [elems[self.headerIndexMap[name]] for name in self.colnames])

    def getRatios(self, _ratios, skipNull = True):
        for i in range(0, len(_ratios)):
            for pep in self.peptides:
                _dat = pep.dat[i]
                if skipNull:
                    if _dat != 0:
                        _ratios[i].append(_dat)
                    else: continue
                else: _ratios[i].append(_dat)

    def normRatios(self, _den, _mult = 1):
        for pep in self.peptides:
            pep.normRatios(_den, _mult)

    def rmZeros(self, _na="NA"):
        for pep in self.peptides:
            pep.rmZeros(_na)

    def calcSummary(self, _outlierTest = 1, _na = "NA"):
        # init empty list of lists to hold ratios for each sample
        ratios = [[] for i in self.colnames]
        booList = list()

        self.getRatios(ratios, skipNull=True)
        #if _outlierTest != 0:
        for i in range(0, len(self.colnames)):
            preLen = len(ratios[i])
            if len(ratios[i]) > 1:
                ratios[i], booList = utils.rmOutliers(ratios[i], _outlierTest)
            self.nOutliersDetected.append(preLen - len(ratios[i]))

            booListCounter = 0
            for pep in self.peptides:
                if pep.dat[i] != 0 and len(booList) > 0:
                    pep.outliersList.append(str(int(not booList[booListCounter])))
                    booListCounter += 1
                else:
                    pep.outliersList.append(_na)
        #else:
        #    self.nOutliersDetected = [0 for x in self.colnames]

        for i in range(0, len(self.colnames)):
            if len(ratios[i]) == 0:
                self.mrList.append(_na)
                self.avgrList.append(_na)
                self.sdList.append(_na)
                self.nPeptidesList.append(_na)
            else:
                self.mrList.append(numpy.median(ratios[i], axis=0))
                self.avgrList.append(numpy.mean(ratios[i], axis=0))
                self.sdList.append(utils.meanAbsDev(ratios[i]))
                self.nPeptidesList.append((len([x for x in ratios[i] if x != 0])))

        for pep in self.peptides:
            if pep.nPeptides > 0:
                self.nPeptides += 1


    def writeWide(self, outF, _output, _peptideSummary, _protSummary = True, _skipNull=True):
        if _skipNull and self.nPeptides == 0:
            return

        assert(_output == "protein" or _output == "peptide" or _output == "both")
        includeProtein = (_output == "protein" or _output == "both")
        includePeptide = (_output == "peptide" or _output == "both")

        if includeProtein:
            if _output == "both":
                outF.write("protein\t")

            outF.write(self.id +
                   '\t' + self.protein +
                   '\t' + self.description + ('\t' * 2))

            # stuff which will already be in peptide dat
            for val in self.mrList:
                outF.write('\t' + str(val))

        # populate str with stuff to include after each peptide
        datTemp = str()
        outliersProt = str()
        for val in self.avgrList:
            datTemp += '\t' + str(val)

        for val in self.sdList:
            datTemp += '\t' + str(val)

        for val in self.nPeptidesList:
            datTemp += '\t' + str(val)

        for val in self.nOutliersDetected:
            outliersProt += '\t' + str(val)

        if includeProtein:
            outF.write(datTemp + outliersProt + '\n')

        if includePeptide:
            for pep in self.peptides:
                if pep.nPeptides == 0 and _skipNull:
                    continue
                else:
                    outliersPep = pep.writeWide(outF, (_output == "both"), _protSummary)
                    if _protSummary:
                        outF.write(datTemp)
                        if _peptideSummary == "byProtein":
                            outF.write(outliersProt)
                    elif includeProtein:
                        outF.write('\t' * (self.nCols * 3))
                    if _peptideSummary == "byPeptide":
                        outF.write(outliersPep)
                    outF.write('\n')

    def writeLong(self, outF, _output, _index, _parseReplicate, _skipNull = True):
        if _skipNull and self.nPeptides == 0:
            return

        # check if any ratios were found for the current index
        if (not any([x.dat[_index] > 0 for x in self.peptides])) and _skipNull:
            return

        assert (_output == "protein" or _output == "peptide" or _output == "both")
        includeProtein = (_output == "protein" or _output == "both")
        includePeptide = (_output == "peptide" or _output == "both")

        if includeProtein:
            if _output == "both":
                outF.write("protein\t")

            outF.write(self.id +
                       '\t' + self.protein +
                       '\t' + self.description)

            if includePeptide:
                outF.write('\t' * 2)

        if _parseReplicate:
            sTemp = (str(self.splitSamples[_index][0]) + "_" +
                     str(self.splitSamples[_index][1]) + '\t' +
                     str(self.splitSamples[_index][0]) + '\t' +
                     str(self.splitSamples[_index][1]))
        else:
            sTemp = str(self.splitSamples[_index])

        if includeProtein:
            outF.write('\t' + sTemp)
            if includePeptide:
                outF.write('\t')

            outF.write('\t' + str(self.mrList[_index]) +
                       '\t' + str(self.avgrList[_index]) +
                       '\t' + str(self.sdList[_index]) +
                       '\t' * 2 + str(self.nOutliersDetected[_index]) +
                       '\t' + str(self.nPeptidesList[_index]) + '\n')

        if includePeptide:
            for pep in self.peptides:
                if pep.dat[_index] == 0 and _skipNull:
                    continue
                outliersPep = pep.writeLong(_index, outF, sTemp, (_output == "both"))
                outF.write('\t' + str(self.mrList[_index]) +
                           '\t' + str(self.avgrList[_index]) +
                           '\t' + str(self.sdList[_index]) +
                           '\t' + str(outliersPep) +
                           '\t' + str(self.nOutliersDetected[_index]) +
                           '\t' + str(self.nPeptidesList[_index]) + '\n')




