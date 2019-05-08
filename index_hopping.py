#!/usr/bin/env python3

import sys
import os
import gzip
import pandas as pd
import argparse
import multiprocessing
import time
import numpy as np

# --------------------------------
# v05.2 - Accepts parameter "demuxed_reads" to compute stats directly
#         Single function to show stats to screen and file
#         Provides now complete stats on number of reads and index hopping rate
#
# v05.1 - Added Saving of files at the end
# v05.0 - Computes 0 and 1 mismatches - splits list of barcodes on 0m and 1m
#
# v04.0 - Removed parameter -m - Now it expects ALWAYS a undetermined.fastq.gz file from 0m demux
#                                It will always process up to 1 mismatch
#
# v03.4 - Added parameter -r (row) to indicate where the data starts in the Sample Sheet
# v03.2 - Added quick search for 0 mismatch using hash table with valid barcodes
#         Split has tables on two, one for i5 and one for i7
#         Added --max parameter
#         Reduced time to 11 sec / 1M reads (0 mismatches) and 13 sec / 1M reads (1 mismatch)
#
# v03.1 - Added -m parameter
#
# v03 - Added saving list of compatible (but invalid) barcode pairs
#
# Hashing  -  From https://interactivepython.org/runestone/static/pythonds/SortSearch/Hashing.html
# Reading  -  From https://darrenjw.wordpress.com/tag/biopython/
# --------------------------------


def line_to_record(line=None):
    record = dict()
    record["index1"] = str(line.split(':')[-1]).split("+")[0].rstrip()
    record["index2"] = str(line.split(':')[-1]).split("+")[1].rstrip()
    return record

def print_results_both(time_start,linesProcessed,MismatchesMatrixNoValid, MismatchesMatrixValid,myStatus,args):
    outputPrefix        = args.output
    logFileName    = outputPrefix + "_log.txt"
    f = open(logFileName, "w")
    print_results(time_start, linesProcessed, MismatchesMatrixNoValid, MismatchesMatrixValid, f, myStatus, args)
    f.close()
    f = sys.stdout
    print_results(time_start, linesProcessed, MismatchesMatrixNoValid, MismatchesMatrixValid, f, myStatus, args)

def print_results(time_start,linesProcessed, MismatchesMatrixNoValid, MismatchesMatrixValid,f,myStatus,args):
    time_end = time.time()
    time_million = 1000000 * (time_end - time_start) / linesProcessed

    demuxReads          = int(args.demuxReads)
    ignoreValid         = args.ignorevalid
    outputPrefix        = args.output
    fastqFileName       = args.input

    outputFileName = outputPrefix + "_reads.fastq"
    logFileName    = outputPrefix + "_log.txt"

    if ignoreValid:
        saveValidStr = "No"
    else:
        saveValidStr = "Yes"

    #Counts the amount of valid compatible barcodes (no swap, but removed by demux because of 1 mismatch)
    MMShape = MismatchesMatrixValid.shape
    mySumValid = 0
    for i in range(0, MMShape[0]):
        for j in range(0, MMShape[1]):
            mySumValid = mySumValid + MismatchesMatrixValid[i, j]

    #Counts the amount of non-valid compatible barcodes (swap)
    MMShape = MismatchesMatrixNoValid.shape
    mySumNoValid = 0
    for i in range(0,MMShape[0]):
        for j in range(0,MMShape[1]):
            mySumNoValid = mySumNoValid + MismatchesMatrixNoValid[i,j]

    # Reads that are compatible and no valid (wrong pair), with 0 mismatches
    mySumNoValid0m = MismatchesMatrixNoValid[0,0]

    # Compatible Reads = Compatible Valid + Compatible Non Valid
    linesCompatible = mySumNoValid + mySumValid

    f.write("\n")
    f.write("Status ................................ : " + myStatus + "\n")
    f.write("Input File ............................ : " + fastqFileName + "\n")
    f.write("Output File ........................... : " + outputFileName + "\n")
    f.write("Save Compatible Valid Barcodes ........ : " + saveValidStr + "\n")
    f.write("Seconds (Per Million barcodes) ........ : %.2f" % (time_million) + "\n")
    f.write("Seconds (Total Processing Time) ....... : %.2f" % (time_end - time_start) + "\n")
    f.write("\n")
    f.write("A = Compatible Valid Barcodes ......... : " + "{:>12.0f}".format(mySumValid) + "\n")
    f.write("B = Compatible Non Valid Barcodes ..... : " + "{:>12.0f}".format(mySumNoValid) + "\n")
    f.write("C = No Compatible Barcodes ............ : " + "{:>12.0f}".format(linesProcessed - linesCompatible) + "\n")
    f.write("Barcodes in Undetermined File = A+B+C . : " + "{:>12.0f}".format(linesProcessed) + "\n")
    f.write("Compatible Barcodes = A+B ............. : " + "{:>12.0f}".format(linesCompatible) + "\n")
    if demuxReads>0:
        f.write("D = Demuxed Valid Barcodes ............ : " + "{:>12.0f}".format(demuxReads)  + "\n")
        f.write("Total Barcodes = A+B+C+D .............. : " + "{:>12.0f}".format(linesProcessed + demuxReads) + "\n")
    f.write("\n")
    if demuxReads>0:
        f.write("Reads Demuxed with 0m ................. : " + "{:>12.0f}".format(demuxReads) + "\n")
        f.write("Reads Demuxed with 1m ................. : " + "{:>12.0f}".format(demuxReads+mySumValid) + "\n")
    f.write("Reads lost when using 0m at Demux ..... : " + "{:>12.0f}".format(mySumValid) + "\n")
    if demuxReads>0:
        f.write("Read lost rate when using 0m at Demux . : " + "{:>12.3f}".format( 100*mySumValid/(demuxReads+mySumValid)) + " %" + "\n")
    f.write("\n")
    f.write("0m - Reads with Index Hopping  ........ : " + "{:>12.0f}".format(mySumNoValid0m) + "\n")
    if demuxReads>0:
        f.write("0m - Reads Demuxed .................... : " + "{:>12.0f}".format(demuxReads) + "\n")
        f.write("0m - Reads Total ...................... : " + "{:>12.0f}".format(demuxReads + mySumNoValid0m) + "\n")
        f.write("0m - Index Hopping Rate ............... : " + "{:>12.3f}".format( 100*mySumNoValid0m/(demuxReads+mySumNoValid0m)) + " %" + "\n")
    f.write("\n")
    f.write("1m - Reads with Index Hopping ......... : " + "{:>12.0f}".format(mySumNoValid) + "\n")
    if demuxReads>0:
        f.write("1m - Reads Demuxed .................... : " + "{:>12.0f}".format(demuxReads+mySumValid) + "\n")
        f.write("1m - Reads Total ...................... : " + "{:>12.0f}".format(demuxReads + mySumValid + mySumNoValid) + "\n")
        f.write("1m - Index Hopping Rate ............... : " + "{:>12.3f}".format( 100*mySumNoValid/(demuxReads + mySumValid + mySumNoValid) )  + " %" + "\n")
    f.write("\n")
    MMShape = MismatchesMatrixValid.shape
    for i in range(0,MMShape[0]):
        for j in range(0,MMShape[1]):
            f.write("Compatible Valid (i7="+str(i)+"m,i5="+str(j)+"m) ........ : " + "{:>12.0f}".format(MismatchesMatrixValid[i,j]) + "\n")
    f.write("Compatible Valid Total ................ : " + "{:>12.0f}".format(mySumValid) + "\n")
    f.write("\n")
    MMShape = MismatchesMatrixNoValid.shape
    for i in range(0,MMShape[0]):
        for j in range(0,MMShape[1]):
            f.write("Compatible Non Valid (i7="+str(i)+"m,i5="+str(j)+"m) .... : " + "{:>12.0f}".format(MismatchesMatrixNoValid[i,j]) + "\n")
    f.write("Compatible Non Valid Total ............ : " + "{:>12.0f}".format(mySumNoValid) + "\n")


def save_results(outputPrefix, Prefix, ConfusionMatrix, index7list, index5list):

    cfFile = outputPrefix + "_"  + Prefix + "_confusion_list.csv"
    with open(cfFile, 'w') as f:
        f.write("i7Id,i7Seq,i5SeqPair,i5Id,i5Seq,i7SeqPair,Count\n")
        MMShape = ConfusionMatrix.shape
        for i in range(0,MMShape[0]):
            for j in range(0,MMShape[1]):
                f.write(str(i) + "," + str(index7list[i]) + "," + str(index5list[i]) + "," + str(j) + "," + str(index5list[j]) + "," + str(index7list[j]) + "," + str(ConfusionMatrix[i,j])  + "\n")
    f.close()

    cfFile = outputPrefix + "_"  + Prefix + "_confusion_matrix.csv"
    with open(cfFile, 'w') as f:
        MMShape = ConfusionMatrix.shape
        f.write("Index")
        for j in range(0, MMShape[1]):
            f.write(",")
            f.write(str(j))
        f.write("\n")
        for i in range(0,MMShape[0]):
            f.write(str(i))
            for j in range(0,MMShape[1]):
                f.write("," + str(ConfusionMatrix[i,j]))
            f.write("\n")

    f.close()

    summaryI7FileName = outputPrefix + "_"  + Prefix + "_summ_i7.csv"
    with open(summaryI7FileName, 'w') as f:
        f.write("Index,Sequence,Hits,Different\n")
        MMShape = ConfusionMatrix.shape
        for i in range(0,MMShape[0]):
            mySum = 0
            nonZero = 0
            for j in range(0,MMShape[1]):
                mySum = mySum + ConfusionMatrix[i,j]
                if ConfusionMatrix[i,j] > 0:
                    nonZero = nonZero + 1
            f.write(str(i) + "," + index7list[i] + "," + str(mySum) + "," + str(nonZero) + "\n")
    f.close()

    summaryI5FileName = outputPrefix + "_"  + Prefix + "_summ_i5.csv"
    with open(summaryI5FileName, 'w') as f:
        f.write("Index,Sequence,Hits,Different\n")
        MMShape = ConfusionMatrix.shape
        for j in range(0, MMShape[1]):
            mySum = 0
            nonZero = 0
            for i in range(0, MMShape[0]):
                mySum = mySum + ConfusionMatrix[i,j]
                if ConfusionMatrix[i, j] > 0:
                    nonZero = nonZero + 1
            f.write(str(j) + "," + index5list[j] + "," + str(mySum) + "," + str(nonZero) + "\n")
    f.close()

def readRead(s):
    return [s.readline(),s.readline(),s.readline(),s.readline()]

def writeRead(sread,s):
    for i in range(4):
        s.write(sread[i])


# Initializes Hash tables with compatible indices
def generateCompatibleTables0m(indexlist):
    myCollectionCompatible0m = {}
    position = 0
    for myIndex in indexlist:
        myCollectionCompatible0m[myIndex] = position
        position = position + 1
    return myCollectionCompatible0m


# Initializes Hash tables with compatible indices
# It generates all barcodes
# that are one mismatch from the barcodes in the list
# The value stored in each entry indicates the position
# of the barcode in the original list
def generateCompatibleTables1m(indexlist):
    myCollectionCompatible1m = {}
    position = 0
    for myIndex in indexlist:
        for i in range(0,len(myIndex)):
            myTemp = list(myIndex)
            myTemp[i] = 'A'
            myIndexB = ''.join(myTemp)
            if myIndexB!=myIndex:
                myCollectionCompatible1m[myIndexB] = position
            myTemp[i] = 'T'
            myIndexB = ''.join(myTemp)
            if myIndexB!=myIndex:
                myCollectionCompatible1m[myIndexB] = position
            myTemp[i] = 'C'
            myIndexB = ''.join(myTemp)
            if myIndexB!=myIndex:
                myCollectionCompatible1m[myIndexB] = position
            myTemp[i] = 'G'
            myIndexB = ''.join(myTemp)
            if myIndexB!=myIndex:
                myCollectionCompatible1m[myIndexB] = position
            myTemp[i] = 'N'
            myIndexB = ''.join(myTemp)
            if myIndexB!=myIndex:
                myCollectionCompatible1m[myIndexB] = position

        position = position + 1

    return myCollectionCompatible1m

def mainLoop():

    # command line arguments
    parser = argparse.ArgumentParser(description='Determines the number and percentage of compatible but invalid barcode pairs')
    parser.add_argument('-i',  '--input', required=True, help='Name of the fastq gzipped file with undemuxed reads - must be the result of demuxing with 0 mismatches')
    parser.add_argument('-s',  '--samples', required=True, help='Name of the sample sheet CSV file with barcodes as used for demux. Must contain ONLY valid barcodes')
    parser.add_argument('-r',  '--row', required=False, help='Row with columns header in sample sheet (Default = 17)', default=17)
    parser.add_argument('-dr', '--demuxReads', help="Number of reads demuxed with 0 mismatches. If provided will generate swap stats", default = 0)
    parser.add_argument('-o',  '--output', required=False, help='prefix for output files. It may include folder name. No need to end it with "_"', default="")
    parser.add_argument('-n',  '--printEach', required=False, help='Number of iters before printing', default=1000)
    parser.add_argument('-x',  '--max', required=False, help='Max Number of reads to process (Default = 0 - All)', default=0)
    parser.add_argument('-iv', '--ignore-valid', help="Do not save valid pairs (undemuxed because of 1 mismatch) in output fastq file", dest='ignorevalid', action='store_true')
    parser.add_argument('-v',  '--verbose', help="Makes verbose", dest='verbose', action='store_true')

    args = parser.parse_args()
    fastqFileName       = args.input
    sampleSheetFileName = args.samples
    isVerbose           = args.verbose
    printEach           = int(args.printEach)
    maxReads            = int(args.max)
    sampleSheetRow      = int(args.row)
    ignoreValid         = args.ignorevalid
    outputPrefix        = args.output
    demuxReads          = int(args.demuxReads)


    maxMismatchsAllowed = 1

    # Reads Barcodes List from Demux Sample Sheet (same used for demux)
    df = pd.read_csv(sampleSheetFileName, header = 0, skiprows=sampleSheetRow-1)
    index7list = df.loc[:,'index']
    index5list = df.loc[:,'index2']

    # Open undemuxed file
    myInput = gzip.open(fastqFileName, 'rt')

    # Open output file for saving
    outputFileName = outputPrefix + "_reads.fastq"
    if outputFileName:
        myOutput = open(outputFileName, 'w')
    else:
        myOutput = None

    numI7 = index7list.size
    numI5 = index5list.size
    ConfusionMatrixValid   = np.zeros( (numI7, numI5) )
    ConfusionMatrixNoValid = np.zeros( (numI7, numI5) )

    MismatchesMatrixNoValid = np.zeros( (maxMismatchsAllowed+1, maxMismatchsAllowed+1) )
    MismatchesMatrixValid   = np.zeros( (maxMismatchsAllowed+1, maxMismatchsAllowed+1) )

    # Initializes tables with valid 17 and 15 indices
    myCollectionValidi70m = generateCompatibleTables0m(index7list)
    myCollectionValidi50m = generateCompatibleTables0m(index5list)
    myCollectionValidi71m = generateCompatibleTables1m(index7list)
    myCollectionValidi51m = generateCompatibleTables1m(index5list)


    linesProcessed  = 0
    linesCompatible = 0
    time_start = time.time()
    with myInput:
        s = readRead(myInput)
        while s[0]:

            # Get record and barcodes
            record = line_to_record(s[0])

            index17 = record["index1"]
            index15 = record["index2"]

            isCompatible = False
            isValid      = False
            saveOutput   = False

            if index17 in myCollectionValidi70m:
                myPosiI7 = myCollectionValidi70m[index17]
                i7dist = 0
            elif index17 in myCollectionValidi71m:
                myPosiI7 = myCollectionValidi71m[index17]
                i7dist = 1
            else:
                i7dist = maxMismatchsAllowed+1

            if index15 in myCollectionValidi50m:
                myPosiI5 = myCollectionValidi50m[index15]
                i5dist = 0
            elif index15 in myCollectionValidi51m:
                myPosiI5 = myCollectionValidi51m[index15]
                i5dist = 1
            else:
                i5dist = maxMismatchsAllowed+1

            if i7dist<=maxMismatchsAllowed and i5dist<=maxMismatchsAllowed:
                isCompatible = True


            if isCompatible:
                linesCompatible += 1
                if  myPosiI7 == myPosiI5:
                    MismatchesMatrixValid[i7dist, i5dist] += 1
                    ConfusionMatrixValid[myPosiI7,myPosiI5] = ConfusionMatrixValid[myPosiI7,myPosiI5] + 1
                    isValid = True
                else:
                    MismatchesMatrixNoValid[i7dist, i5dist] += 1
                    ConfusionMatrixNoValid[myPosiI7,myPosiI5] = ConfusionMatrixNoValid[myPosiI7,myPosiI5] + 1

            if isCompatible and not (isValid and ignoreValid):
                saveOutput = True

            if saveOutput and myOutput:
                writeRead(s,myOutput)

            linesProcessed = linesProcessed + 1
            if linesProcessed % printEach == 0:
                print_results_both(time_start, linesProcessed, MismatchesMatrixNoValid, MismatchesMatrixValid, "Processing", args)

            if maxReads>0 and linesProcessed >= maxReads:
                break

            s = readRead(myInput)

    save_results(outputPrefix, "NoValid", ConfusionMatrixNoValid, index7list, index5list)
    save_results(outputPrefix, "Valid",   ConfusionMatrixValid,   index7list, index5list)
    print_results_both(time_start, linesProcessed, MismatchesMatrixNoValid, MismatchesMatrixValid,"Finished",args)

if __name__ == '__main__':
    mainLoop()

