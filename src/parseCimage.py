
import argparse
import sys
import os
import dta
import textwrap


PROG_VERSION = "1.0"
OFNAME_WIDE = "_fixed.tsv"
OFNAME_LONG = "_fixed_long.tsv"
DEFAULT_SAMPLE_RE = "^mr.set_[0-9]+$"

def main(argv):

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                            description = textwrap.dedent('''\
                            Format cimage combined_dta.tsv file.

                            By default the all columns with headers following \"mass\" with the prefix
                            \"mr.set_\" will be assumed to be columns containing peptide ratios.
                            All peptide ratios are normalized to a median ratio of 1 by sample.
                            Next median and mean ratios are calculated for each protein after removing
                            outliers using the interquartile range method.

                            '''),
                            prog = "parseCimage",
                            epilog="parseCimage was written by Aaron Maurais.\n"
                                   "Email questions or bugs to aaron.maurais@bc.edu")

    parser.add_argument("input_file",
                        help="file to read")

    parser.add_argument("--version", action='version', version='%(prog)s ' + PROG_VERSION)

    parser.add_argument('-b', '--na',
                         help="Specify how to represent ratios with no data in output\nfiles.",
                         default="")

    parser.add_argument("-i", "--input",
                        help = "Specify input format to expect. Default is byProtein",
                        choices = ["byPeptide", "byProtein"],
                        default = "byProtein")

    parser.add_argument("-o", "--output",
                        help = "Choose whether to include peptide, protein or both lines\nin output files.",
                        choices = ['peptide', 'protein', 'both'],
                        default = 'both')

    parser.add_argument("-n", "--normalize",
                        help = "Specify whether to normalize the median peptide ratio in\neach sample to 1.",
                        type=int,
                        choices = [0, 1],
                        default = 1)

    parser.add_argument("-z", "--nullP",
                        help = "Choose whether to include proteins with no ratios in\noutput file.",
                        type=int,
                        choices = [0, 1],
                        default = 0)

    parser.add_argument("-s", "--parentDat",
                        help = "Specify whether to include parent protein data on each\npeptide line. "
                               "Default is 1.",
                        type=int,
                        choices = [0, 1],
                        default = 1)

    parser.add_argument("-d", "--direction",
                        help="Output file direction.",
                        choices=['wide', 'long', 'both'],
                        default='wide')

    parser.add_argument("-c", "--peptideSummary",
                        help = "Choose how to calculate nOutlier column in output files.\nbyPeptide is the default. \n"
                               "  byPeptide : The value in the nOutlier will be 1 if the\n      peptide "
                               "ratio was used to calculate the median and\n      average ratio for that peptide, "
                               "and 0 if it was\n      found to be an outlier.\n"
                               "  byProtein : The value in nOutlier will be the total\n      "
                               "outliers found for that sample.\n",
                        choices = ["byPeptide", "byPeptide"],
                        default = "byPeptide")

    parser.add_argument("--outlierTest",
                        help= "Specify outlier test to use. 1 is the default.\n"
                        "  0 : do not remove outliers before calculating median\n"
                              "      and mean protein ratios.\n"
                        "  1 : use interquartile range outlier test.\n"
                        "  2 : use modified z score outlier test.\n",
                        type=int,
                        choices = [0, 1, 2],
                        default = 1)

    parser.add_argument("-r", "--sampleRegex",
                        help = "Specify regex to find peptide ratio column.\nDefault is " +
                               '\"' + DEFAULT_SAMPLE_RE + '\"',
                        default = DEFAULT_SAMPLE_RE)

    parser.add_argument("-p", "--samplePrefix",
                        help = "Specify prefix for sample names. The prefix string will \nbe removed from"
                        " all sample names where it is found\nin output files. Default is \"mr.\".",
                        default = "mr.")

    parser.add_argument("-f", "--parseReplicate",
                        help = "Choose whether to parse sample replicate from column\nheaders. "
                               "The sample replicate is taken as the text\nfollowing the last "
                               "underscore in the sample name.\nDefault is 1.",
                        type = int,
                        choices = [0, 1],
                        default = 1)

    parser.add_argument("-q", "--quiet",
                        help = "quiet output.",
                        action = "store_true",
                        default = False)

    args = parser.parse_args()

    if not args.quiet:
        sys.stdout.write("parseCimage v" + PROG_VERSION + '\n')

    fname = os.path.abspath(args.input_file)
    if args.sampleRegex != DEFAULT_SAMPLE_RE:
        args.samplePrefix = ""

    # initialize and read data file
    if not args.quiet:
        sys.stdout.write("\nReading " + args.input_file + "...\n")
    dtaFile = dta.DtaFile(fname)
    dtaFile.setVerbosity(not args.quiet)
    if not dtaFile.read(args.sampleRegex, args.input):
        sys.stderr.write("Error reading " + fname + '\n')
        exit()

    # process and perform statistics
    if not args.quiet:
        sys.stdout.write("\nCalculating summary data...\n")
    dtaFile.rmZeros(args.na)
    if args.normalize:
        dtaFile.normalizeBySample(1)
    dtaFile.calcSummary(args.outlierTest, args.na)

    # write output files
    if args.direction == "wide" or args.direction == "both":
        if not args.quiet:
            sys.stdout.write("\nWriting data...")
        ofname = os.path.splitext(fname)[0] + OFNAME_WIDE
        if not dtaFile.writeWide(ofname, args.output, args.peptideSummary,
                                 args.samplePrefix, args.parentDat, (not args.nullP)):
            sys.stderr.write("Failed to write" + ofname + "\nExiting...\n")
            exit()
        if not args.quiet:
            sys.stdout.write("done\n" + "Data written in wide format to: " + os.path.basename(ofname) + '\n')

    if args.direction == "long" or args.direction == "both":
        if not args.quiet:
            sys.stdout.write("\nWriting data...")
        ofname = os.path.splitext(fname)[0] + OFNAME_LONG
        if not dtaFile.writeLong(ofname, args.output, args.samplePrefix,
                                 args.parseReplicate, (not args.nullP)):
            sys.stderr.write("Failed to write" + ofname + "\nExiting...\n")
            exit()
        if not args.quiet:
            sys.stdout.write("done\n" + "Data written in long format to: " + os.path.basename(ofname) + '\n')

if __name__ == "__main__":
    main(sys.argv)

