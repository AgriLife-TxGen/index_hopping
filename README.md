# index_hopping
Tool to extract and count invalid index pairs

## Usage
```
index_hopping_v05.2.py [-h] -i INPUT -s SAMPLES [-r ROW]
                              [-dr DEMUXREADS] [-o OUTPUT] [-n PRINTEACH]
                              [-x MAX] [-iv] [-v]

Determines the number and percentage of compatible but invalid barcode pairs

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Name of the fastq gzipped file with undemuxed reads -
                        must be the result of demuxing with 0 mismatches
  -s SAMPLES, --samples SAMPLES
                        Name of the sample sheet CSV file with barcodes as
                        used for demux. Must contain ONLY valid barcodes
  -r ROW, --row ROW     Row with columns header in sample sheet (Default = 17)
  -dr DEMUXREADS, --demuxReads DEMUXREADS
                        Number of reads demuxed with 0 mismatches. If provided
                        will generate swap stats
  -o OUTPUT, --output OUTPUT
                        prefix for output files. It may include folder name.
                        No need to end it with "_"
  -n PRINTEACH, --printEach PRINTEACH
                        Number of iters before printing
  -x MAX, --max MAX     Max Number of reads to process (Default = 0 - All)
  -iv, --ignore-valid   Do not save valid pairs (undemuxed because of 1
                        mismatch) in output fastq file
  -v, --verbose         Makes verbose
```

## Example
In this case we assume that 1,000,000 reads were successfully demuxed with 0 mismatches, and the
reads that were not assigned to samples, are stored in the file undetermined.fastq.gz, while the
barcode information is stored in the sample sheet used for demuxing, named samplesheet.csv

```
index_hopping.py -i undetermined.fastq.gz -s samplesheet.csv -o output/results -dr 1000000
```

## Output
Several files are created with the specified prefix:
```
results_reads.fastq                  - fastq file containing the reads containing compatible barcodes
results_log.txt                      - Log File with statistics
results_NoValid_confusion_list.csv   - File containing count of observed invalid pairs, as a list
results_NoValid_confusion_matrix.csv - File containing count of observed invalid pairs, as a matrix 
results_NoValid_summ_i7.csv          - File containing count of observed invalid pairs for each i7 barcode
results_NoValid_summ_i5.csv          - File containing count of observed invalid pairs for each i5 barcode
results_Valid_confusion_list.csv     - File containing count of observed valid pairs, as a list
results_Valid_confusion_matrix.csv   - File containing count of observed valid pairs, as a matrix 
results_Valid_summ_i7.csv            - File containing count of observed valid pairs for each i7 barcode
results_Valid_summ_i5.csv            - File containing count of observed valid pairs for each i5 barcode
```

## Results Log 
The log file contains the following information, which is also displayed on screen

```
Status ................................ : Finished
Input File ............................ : undetermined.fastq.gz
Output File ........................... : output/results_reads.fastq
Save Compatible Valid Barcodes ........ : Yes
Seconds (Per Million barcodes) ........ : 
Seconds (Total Processing Time) ....... : 

A = Compatible Valid Barcodes ......... : 
B = Compatible Non Valid Barcodes ..... : 
C = No Compatible Barcodes ............ : 
Barcodes in Undetermined File = A+B+C . : 
Compatible Barcodes = A+B ............. : 
D = Demuxed Valid Barcodes ............ : 
Total Barcodes = A+B+C+D .............. : 

Reads Demuxed with 0m ................. : 
Reads Demuxed with 1m ................. : 
Reads lost when using 0m at Demux ..... : 
Read lost rate when using 0m at Demux . : 

0m - Reads with Index Hopping  ........ : 
0m - Reads Demuxed .................... : 
0m - Reads Total ...................... : 
0m - Index Hopping Rate ............... : 

1m - Reads with Index Hopping ......... : 
1m - Reads Demuxed .................... : 
1m - Reads Total ...................... : 
1m - Index Hopping Rate ............... : 

Compatible Valid (i7=0m,i5=0m) ........ : 
Compatible Valid (i7=0m,i5=1m) ........ : 
Compatible Valid (i7=1m,i5=0m) ........ : 
Compatible Valid (i7=1m,i5=1m) ........ : 
Compatible Valid Total ................ : 

Compatible Non Valid (i7=0m,i5=0m) .... : 
Compatible Non Valid (i7=0m,i5=1m) .... : 
Compatible Non Valid (i7=1m,i5=0m) .... : 
Compatible Non Valid (i7=1m,i5=1m) .... : 
Compatible Non Valid Total ............ : 
```

