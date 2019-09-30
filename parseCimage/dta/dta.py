
from .protein import Protein
import re
import numpy
from .constants import *
import sys
import copy

class DtaFile:
    headerIndexMap = {INDEX_COLNAME : None,
                      ID_COLNAME : None,
                      DESCRIPTION_COLNAME : None,
                      PROTEIN_COLNAME : None,
                      SEQ_COLNAME : None,
                      MASS_COLNAME : None,
                      CHARGE_COLNAME : None,
                      SEGMENT_COLNAME : None}
    colHeadersPre = ["ID", "protein", "description",
                     "sequence", "charge", "segment"]
    colHeadersLong = ["sample", "ratio", "median_ratio", "avg_ratio",
                      "median_avg_dev", "outlierBool", "nOutliers", "nPeptides"]
    verbose = True

    def __init__(self, _fname):
        self.fname = _fname
        self.proteins = dict()
        self.colnames = list()


    @staticmethod
    def setVerbosity(boo):
        DtaFile.verbose = boo
        Protein.verbose = boo


    def parseColHeaders(self, elems, _sampleRegex):
        # find indecies of static col headers
        for i, s in enumerate(elems):
            if s in self.headerIndexMap.keys():
                if self.verbose:
                    sys.stdout.write("Found " + s + " column index...\n")
                self.headerIndexMap[s] = i

        # check that required columns were found
        for name in self.headerIndexMap.keys():
            if self.headerIndexMap[name] is None:
                sys.stderr.write(name + " column index not found!\n")
                return False

        # find peptide ratio columns and put indecies in headerIndexMap
        assert isinstance(_sampleRegex, str)
        #pattern = re.compile(_sampleRegex)
        #self.colnames = filter(pattern.match, elems)
        self.colnames = [x for x in elems if re.search(_sampleRegex, x)]
        if len(self.colnames) == 0:
            sys.stdout.write("No data columns detected!\n")
            return False
        if self.verbose:
            sys.stdout.write("Found sample columns:\n")
        for name in self.colnames:
            if self.verbose:
                sys.stdout.write("\t" + name + '\n')
            self.headerIndexMap[name] = elems.index(name)

        return True


    def read(self, _sampleRegex, _inputFormat = "byProtein"):
        if _inputFormat == "byProtein":
            self.readByProtein(_sampleRegex)
        elif _inputFormat == "byPeptide":
            self.readByPeptide(_sampleRegex)
        else:
            raise RuntimeError("\nInvalid input format!\n")


    def readByPeptide(self, _sampleRegex):
        inF = open(self.fname, 'rU')
        lines = inF.readlines()

        headersFound = False
        indexColNum = 0

        for i in range(0, len(lines)):
            elems = lines[i].rstrip('\r\n').split('\t')
            if not headersFound:
                if elems[0] == INDEX_COLNAME:  # find header line by index column
                    headersFound = True

                    if not self.parseColHeaders(elems, _sampleRegex):
                        return False

                    Protein.setColnames(self.colnames)
                    Protein.headerIndexMap = self.headerIndexMap
                    indexColNum = self.headerIndexMap[INDEX_COLNAME]
                    continue

            # find peptide lines and ignore index lines by looking for blank ID column
            if elems[self.headerIndexMap[ID_COLNAME]].strip() != "":
                # check that file is not byProtein
                if elems[indexColNum].strip() != "":
                    sys.stderr.write("\nError!\nbyProtein input file detected!"
                                     " Use --input option to specify input format!\nExiting...\n\n")
                    return False

                curID = elems[self.headerIndexMap[ID_COLNAME]]
                if curID not in self.proteins.keys():
                    self.proteins[curID] = Protein(curID,
                                                   elems[self.headerIndexMap[DESCRIPTION_COLNAME]],
                                                   elems[self.headerIndexMap[PROTEIN_COLNAME]])
                    self.proteins[curID].addPeptideLine(elems)
                else:
                    self.proteins[curID].addPeptideLine(elems)

        if self.verbose:
            sys.stdout.write(str(len(self.proteins)) + " proteins read in...\n")


    def readByProtein(self, _sampleRegex):
        inF = open(self.fname, 'rU')
        lines = inF.readlines()
        headersFound = False

        for i in range(0, len(lines)):
            elems = lines[i].rstrip('\r\n').split('\t')
            if not headersFound:
                if elems[0] == INDEX_COLNAME: # find header line by index column
                    headersFound = True

                    if not self.parseColHeaders(elems, _sampleRegex):
                        return False

                    Protein.setColnames(self.colnames)
                    Protein.headerIndexMap = self.headerIndexMap
                    continue

            if elems[0].isdigit():
                curID = elems[self.headerIndexMap[ID_COLNAME]]
                if curID.strip() == "":
                    sys.stderr.write("\nError!\nbyPeptide input file detected!"
                                     " Use --input option to specify input format!\nExiting...\n\n")
                    return False
                self.proteins[curID] = Protein(curID,
                                               elems[self.headerIndexMap[DESCRIPTION_COLNAME]],
                                               elems[self.headerIndexMap[PROTEIN_COLNAME]])
            else:
                assert curID.strip() != ""
                self.proteins[curID].addPeptideLine(elems)

        if self.verbose:
            sys.stdout.write(str(len(self.proteins)) + " proteins read in...\n")

    def rmZeros(self, _na = "NA"):
        for key, prot in self.proteins.items():
            prot.rmZeros(_na)


    def normalizeBySample(self, _targetMedian = 1):
        # find median ratios for each sample
        # initialize dict of list to store ratios
        ratios = [[] for i in range(len(self.colnames))]

        for key, prot in self.proteins.items():
            prot.getRatios(ratios, skipNull=True)

        medians = list()
        for i in range(0, len(ratios)):
            medians.append(numpy.median(ratios[i]))

        # normalize ratios in each sample
        for key in self.proteins.keys():
            self.proteins[key].normRatios(medians, _targetMedian)


    def makeWideHeaders(self, _samplePrefix, _includeSummary = True):
        ret = copy.deepcopy(self.colHeadersPre)
        temp = [name.replace(_samplePrefix, "") for name in self.colnames]

        # add ratio names
        for name in temp:
            ret.append(name + "_mr")

        if _includeSummary:
            for name in temp:
                ret.append(name + "_avg")

            # add sd names
            for name in temp:
                ret.append(name + "_mad")

            # add pep count
            for name in temp:
                ret.append(name + "_nPeptides")

            for name in temp:
                ret.append(name + "_nOutliers")

        return ret


    def calcSummary(self, _outlierTest = 1, _na = "NA"):
        for prot in self.proteins.keys():
            self.proteins[prot].calcSummary(_outlierTest, _na)


    def writeWide(self, fname, _output, _peptideSummary, _samplePrefix,
                  _proteinSummary, _skipNull):
        outF = open(fname, 'w')

        headers = self.makeWideHeaders(_samplePrefix,
                                       not((not _proteinSummary) and (_output == "peptide")))
        if _output == 'both':
            headers = ["DatType"] + headers
        for i, it in enumerate(headers):
            if i == 0:
                outF.write(it)
            else: outF.write('\t' + it)
        outF.write('\n')

        for key in self.proteins.keys():
            self.proteins[key].writeWide(outF, _output, _peptideSummary, _proteinSummary, _skipNull)


    def writeLong(self, fname, _output, _samplePrefix, _parseReplicate, _skipNull):

        outF = open(fname, 'w')

        # make headers
        assert (_output == "protein" or _output == "peptide" or _output == "both")
        headers = self.colHeadersPre + self.colHeadersLong
        if _output == 'both':
            headers = ["DatType"] + headers
        elif _output == "protein":
            headers = [x for x in headers if x not in ["sequence", "charge", "ratio"]]
        tempNames = [name.replace(_samplePrefix, "") for name in self.colnames]

        if _parseReplicate:
            i = headers.index("sample")
            headers.insert(i, "long_sample_name")
            headers.insert(i + 2, "replicate")
            samples = [x.rsplit("_", 1) for x in tempNames]
            Protein.splitSamples = samples
        else:
            Protein.splitSamples = tempNames

        for i, it in enumerate(headers):
            if i == 0:
                outF.write(it)
            else: outF.write('\t' + it)
        outF.write('\n')

        for i, name in enumerate(tempNames):
            for key in self.proteins.keys():
                self.proteins[key].writeLong(outF, _output, i, _parseReplicate, _skipNull)







