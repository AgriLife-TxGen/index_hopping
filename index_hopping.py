#!/usr/bin/env python3

import sys
import os
import gzip
import pandas as pd
import argparse
import multiprocessing
import time

# --------------------------------
# v03.2 - Added quick search for 0 mismatch using hash table with valid barcodes
#         Split has tables on two, one for i5 and one for i7
#         Added --max parameter
#
# v03.1 - Added -m parameter
#
# v03 - Added saving list of compatible (but invalid) barcode pairs
#
# --------------------------------


def lines_to_record(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    record = {k: v for k, v in zip(ks, lines)}
    record["index1"] = str(record["name"].split(':')[-1]).split("+")[0]
    record["index2"] = str(record["name"].split(':')[-1]).split("+")[1]
    return record

def line_to_record(line=None):
    record = {}
    record["index1"] = str(line.split(':')[-1]).split("+")[0].rstrip()
    record["index2"] = str(line.split(':')[-1]).split("+")[1].rstrip()
    return record

def compute_distance(str1,str2,maxMismatchsAllowed):
    maxlen = len(str2) if len(str1) < len(str2) else len(str1)
    differences = 0
    for i in range(maxlen):
        letter1=str1[i:i+1]
        letter2=str2[i:i+1]
        if letter1 != letter2:
            differences +=1
            if differences>maxMismatchsAllowed:
                break
    return differences

def determine_valid(index17,index15,index7list,index5list,maxMismatchsAllowed):

    dist7 = len(index17)
    dist5 = len(index15)
    isValid = {'i7' : 0, 'i5' : 0}
    for index7 in index7list:
        dist7 = min(compute_distance(index17, index7, maxMismatchsAllowed), dist7)
        if dist7 <= maxMismatchsAllowed:
            isValid['i7'] = 1
            break
    for index5 in index5list:
        dist5 = min(compute_distance(index15, index5, maxMismatchsAllowed), dist5)
        if dist5 <= maxMismatchsAllowed:
            isValid['i5'] = 1
            break

    return isValid

def print_results(time_start,linesProcessed,fastqFileName,linesValid,logFileName,outputFileName,myStatus):
    time_end = time.time()
    time_million = 1000000 * (time_end - time_start) / linesProcessed
    print("Input File              : " + fastqFileName)
    print("Output File             : " + outputFileName)
    print("Lines Processed         : " + str(linesProcessed))
    print("Compatible Barcodes     : " + str(linesValid))
    print("No Compatible Barcodes  : " + str(linesProcessed - linesValid))
    print("Percentage (Compatible) : " + str(100 * (linesValid) / linesProcessed) + "%")
    print("Seconds (Total)         : %.2f" % (time_end - time_start))
    print("Seconds (Million)       : %.2f" % (time_million))
    print("Status                  : "+myStatus)
    with open(logFileName, 'w') as f:
        f.write("File\t" + fastqFileName + "\n")
        f.write("Lines Processed\t" + str(linesProcessed) + "\n")
        f.write("Compatible Barcodes\t" + str(linesValid) + "\n")
        f.write("No Compatible Barcodes\t" + str(linesProcessed - linesValid) + "\n")
        f.write("Percentage (Compatible)\t" + str(100 * (linesValid) / linesProcessed) + "\n")
        f.write("Seconds (Total)\t%.2f" % (time_end - time_start) + "\n")
        f.write("Seconds (Million)\t%.2f" % (time_million) + "\n")
        f.write("Status\t"+myStatus+"\n")
    f.close()

def readRead(s):
    return [s.readline(),s.readline(),s.readline(),s.readline()]

def writeRead(sread,s):
    for i in range(4):
        s.write(sread[i])

def mainLoop():

    # command line arguments
    parser = argparse.ArgumentParser(description='Determines the number and percentage of compatible but invalid barcode pairs')
    parser.add_argument('-i',  '--input', required=True, help='Name of the fastq gzipped file with undemuxed reads')
    parser.add_argument('-o',  '--output', required=False, help='Name of the fastq file where compatible reads will be stored (only stats if not defined)', default="")
    parser.add_argument('-lg', '--log', required=True, help='Name of the text file where stats will be saved')
    parser.add_argument('-s',  '--samples', required=True, help='Name of the CSV file with samples information')
    parser.add_argument('-v',  '--verbose', help="Makes verbose", dest='verbose', action='store_true')
    parser.add_argument('-n',  '--printEach', required=False, help='Number of iters before printing', default=1000)
    parser.add_argument('-m',  '--mismatch', required=False, help='Number of mismatches allowed (default = 1)', default=1)
    parser.add_argument('-x',  '--max', required=False, help='Max Number of reads to process (Default = 0 - All)', default=0)

    args = parser.parse_args()
    fastqFileName = args.input
    sampleSheetFileName = args.samples
    isVerbose = args.verbose
    outputFileName = args.output
    printEach = int(args.printEach)
    logFileName = args.log
    maxMismatchsAllowed = int(args.mismatch)
    maxReads = int(args.max)

    # Reads Barcodes List from Demux Sample Sheet (same used for demux)
    df = pd.read_csv(sampleSheetFileName, header = 0, skiprows=16)
    index7list = df.loc[:,'index']
    index5list = df.loc[:,'index2']

    if maxMismatchsAllowed==0:
        sizeCollection = 10**3
    else:
        sizeCollection = 10**5

    maxEntriesPerHash = 100

    myCollectionValidi7 = [None] * sizeCollection
    myCollectionValidi5 = [None] * sizeCollection

    myCollectionInvalidi7 = [None] * sizeCollection
    myCollectionInvalidi5 = [None] * sizeCollection

    myInput = gzip.open(fastqFileName, 'rt')
    if outputFileName:
        myOutput = open(outputFileName, 'w')
    else:
        myOutput = None

    # Initializes Hash tables with valid 17 and 15 indices
    for index17 in index7list:
        hashKeyi7 = hash(index17) % sizeCollection
        if myCollectionValidi7[hashKeyi7]:
            myCollectionValidi7[hashKeyi7].append(index17)
            if isVerbose==1:
                print("i7= " + index17 + ", Key= " + str(hashKeyi7) + ", Added to existing list")
        else:
            myCollectionValidi7[hashKeyi7] = [index17]
            if isVerbose==1:
                print("i7= " + index17 + ", Key= " + str(hashKeyi7) + ", Added new")

    for index15 in index5list:
        hashKeyi5 = hash(index15) % sizeCollection
        if myCollectionValidi5[hashKeyi5]:
            myCollectionValidi5[hashKeyi5].append(index15)
        else:
            myCollectionValidi5[hashKeyi5] = [index15]

    linesProcessed = 0
    linesValid = 0
    time_start = time.time()
    with myInput:
        s = readRead(myInput)
        while s[0]:

            # Get record and barcodes
            record = line_to_record(s[0])

            index17 = record["index1"]
            index15 = record["index2"]

            isValid = 0
            hashKeyi7 = hash(index17) % sizeCollection
            hashKeyi5 = hash(index15) % sizeCollection

            if ( (myCollectionValidi7[hashKeyi7] and index17 in myCollectionValidi7[hashKeyi7]) and
                (myCollectionValidi5[hashKeyi5] and index15 in myCollectionValidi5[hashKeyi5]) ):
                isValid = 1

            elif maxMismatchsAllowed>0:
                if ((myCollectionInvalidi7[hashKeyi7] and index17 in myCollectionInvalidi7[hashKeyi7]) or
                        (myCollectionInvalidi5[hashKeyi5] and index15 in myCollectionInvalidi5[hashKeyi5])):
                    isValid = 0

                else:
                    isValid2 = determine_valid(index17, index15, index7list, index5list, maxMismatchsAllowed)
                    if isValid2['i7']==1 and isValid2['i5']==1:
                        isValid = 1

                    if isValid == 0:
                        if isValid2['i7']==0:
                            if not (myCollectionInvalidi7[hashKeyi7] and index17 in myCollectionInvalidi7[hashKeyi7]):
                                if myCollectionInvalidi7[hashKeyi7]:
                                    if len(myCollectionInvalidi7[hashKeyi7]) < maxEntriesPerHash:
                                        myCollectionInvalidi7[hashKeyi7].append(index17)
                                else:
                                    myCollectionInvalidi7[hashKeyi7] = [index17]
                        if isValid2['i5']==0:
                            if not (myCollectionInvalidi5[hashKeyi5] and index15 in myCollectionInvalidi5[hashKeyi5]):
                                if myCollectionInvalidi5[hashKeyi5]:
                                    if len(myCollectionInvalidi5[hashKeyi5]) < maxEntriesPerHash:
                                        myCollectionInvalidi5[hashKeyi5].append(index15)
                                else:
                                    myCollectionInvalidi5[hashKeyi5] = [index15]

                    else:
                        if not (myCollectionValidi7[hashKeyi7] and index17 in myCollectionValidi7[hashKeyi7]) :
                            if myCollectionValidi7[hashKeyi7]:
                                if len(myCollectionValidi7[hashKeyi7]) < maxEntriesPerHash:
                                    myCollectionValidi7[hashKeyi7].append(index17)
                            else:
                                myCollectionValidi7[hashKeyi7] = [index17]
                        if not (myCollectionValidi5[hashKeyi5] and index15 in myCollectionValidi5[hashKeyi5]) :
                            if myCollectionValidi5[hashKeyi5]:
                                if len(myCollectionValidi5[hashKeyi5]) < maxEntriesPerHash:
                                    myCollectionValidi5[hashKeyi5].append(index15)
                            else:
                                myCollectionValidi5[hashKeyi5] = [index15]

            else:
                isValid = 0


            if isValid == 1:
                linesValid += 1
                if myOutput:
                    writeRead(s,myOutput)

            linesProcessed = linesProcessed + 1
            if linesProcessed % printEach == 0:
                print_results(time_start, linesProcessed, fastqFileName, linesValid, logFileName, outputFileName, "Processing")

            if maxReads>0 and linesProcessed >= maxReads:
                break

            s = readRead(myInput)

    print_results(time_start, linesProcessed, fastqFileName, linesValid,logFileName,outputFileName,"Finished")


if __name__ == '__main__':
    mainLoop()