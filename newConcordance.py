#!/user/bin/env python
#Created by Megan Neveau
#Last edited 6/9/17
import sys
import os.path
import argparse
import tempfile
from subprocess import Popen, PIPE
from collections import defaultdict
import FisherExact

#check to ensure files have corresponding index files
def checkForIndex(file):
    if os.path.isfile(file):
        print("Index file exists")
    else:
        print("Index file does not exist")
        sys.exit()

#execute bam-readcount
def bamReadcount(bamFile):
    bamReadcountCmd = ['/usr/bin/bam-readcount', '-f', refFasta, '-l', snpFile]
    print(refFasta, snpFile)
    outputFile = os.path.join(tempdir, bamFile + '_rc.tsv')
    print(outputFile)
    bamReadcountCmd.append(bamFile)
    print(bamReadcountCmd)
    execution = Popen(bamReadcountCmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = execution.communicate()
    if execution.returncode == 0:
        with open(outputFile, 'wb') as outputFh:
            outputFh.write(stdout)
            print("Output written")
        outputFh.close()
        #quits if output file is empty
        if os.stat(outputFile).st_size == 0:
            print("No output from bam-readcount")
            sys.exit()
    else:
        print("no output written")
        sys.exit(stderr)
    return outputFile 

def parseRc(rcOutputFile, i):
    parseFile = os.path.join(tempdir, "parse_file_" + str(i))
    print(parseFile)
    with open(parseFile, 'w') as parseFileFh:
        with open(rcOutputFile, 'r') as rcOutput:
            for line in rcOutput:
                #replace : with /t
                editedLine = line.replace(":","\t")
                #split at each tab
                field = editedLine.split("\t")
                #write the array elements you want to keep onto parse_file 
                #(chr pos refbase A # C # G # T #)
                parseFileFh.write(field[0] + "\t" + field[1] + "\t" + field [2] + "\t" + field[18] + "\t" + field[19] + "\t" + field[32] + "\t" + field[33] + "\t" + field[46] + "\t" + field[47] + "\t" + field[60] + "\t" + field[61] + "\t" + field[74] + "\t" + field[75] + " \n")
        rcOutput.close()
        parseFileFh.close()
        with open(parseFile, 'r') as f:
           print(f.read())
        f.close()
    return parseFile

#make genotype calls with fisher exact test
def geno_calc(parsed_file, i):
    deep_enough = 0
    output_file = os.path.join(tempdir, "geno_calc" + str(i)) 
    with open(parsed_file, 'r') as data_file:
        with open(output_file, 'w') as out_f:
            for line in data_file:
                counter, letter, alpha, num = 0, 3, [None, None], [0,0]
                print(line)
                field = line.rstrip('\n').split("\t")
                depth = int(field[4]) + int(field[6]) + int(field[8]) + int(field[10])
                if depth >= 10:
                    print("Deep enough! : " + str(depth))
                    deep_enough += 1

                    for val in field[4:11:2]:
                        print(val)
                        print(int(val))
                        if int(val) > 5:
                            print(int(val) > 5)
                            alpha[counter] = field[letter]
                            num[counter] = val
                            counter += 1
                            letter += 2
                            print(alpha + num)
                            print(counter)
                    if int(num[0]) >= 12 and int(num[1]) == 0:
                        field.append(alpha[0] + "\n")
                    elif counter == 2 and (int(num[0]) + int(num[1]) >= 12):
                                total = int(num[0]) + int(num[1])
                                half = total / 2
                                print(str(total) + " " + str(half))
                                p_hetero = FisherExact.fisher_exact([[num[0], num[1]], [half, half]])
                                p_homo_1 = FisherExact.fisher_exact([[num[0], num[1]], [total, 0]])
                                p_homo_2 = FisherExact.fisher_exact([[num[0], num[1]], [0, total]])
                                if p_hetero > p_homo_1 and p_hetero > p_homo_2:
                                    field.append(alpha[0] + "/" + alpha[1] + " \n")
                                elif p_homo_1 > p_homo_2 and p_homo_1 > p_hetero:
                                    field.append(alpha[0] + "\n")
                                else:
                                    field.append(alpha[1] + "\n")   
                    else:
                            field.append("NA\n")
                    print(field)
                    line = "\t".join(field)
                    print(line)
                    out_f.write(line)
                else:
                    field.append("NA\n")
                    out_f.write(line)
            #error out if no data with 10x coverage
            if deep_enough == 0:
                print("No site has greater than 10x coverage")
                sys.exit()
        out_f.close()
    data_file.close()
    return output_file
        
'''
#run R script to make genotype calls
def runR(parsedFile):
    rScriptCmd = ['Rscript']
    rScript = os.path.join(os.path.dirname(os.path.realpath(__file__)), "Concordance.R")
    rScriptCmd.append(rScript)
    rScriptCmd.append(parsedFile)
    output = os.path.join(tempdir, parsedFile +  ".r_file")
    rScriptCmd.append(output)
    print(rScriptCmd)
    execution = Popen(rScriptCmd)
    execution.communicate()
    if execution.returncode != 0:
        print("Error in R Script execution")
        sys.exit()
    return output
'''

#open the output file from R and get the genotype
def getGenotypes(rOutputFile,  n):
    total, snp = [], []
    validGenotypes = 0
    with open (rOutputFile, 'r') as rFile:
        for line in rFile:
            field = line.split('\t')
            print(field[11] + field[12])
            if field[13] != "NA":
                total.append(field[0] + '_' + field[1])
                snp.append(field[0] + '_' + field[1] + '_' + field[13])
                print(total)
                print(snp)
                validGenotypes = validGenotypes + 1
                print("Valid_genotypes: " + str(validGenotypes))
            #if full genotype output is wanted, creates dictionary to do so
            if outputGeno != None:
                print("Genotype file will be created")
                genotypesDict[field[0] + "\t" + field[1] + "\t" + field[2]]["samp" + str(n)] = field[13]
        #quit if no valid genotypes
        if validGenotypes <= 1:
            print("Not enough information to continue")
            sys.exit()
    print(total)
    print(snp)
    print(genotypesDict)
    rFile.close()
    return total, snp

#Program start:
#makes sure that all of the required arguments are present
parser = argparse.ArgumentParser()
parser.add_argument("bam_file_1")
parser.add_argument("bam_file_2")
parser.add_argument("ref_fasta")
parser.add_argument("snp_file")
parser.add_argument("output_file")
parser.add_argument("--output_geno")
args = parser.parse_args()
##unnecessary?
bam1 = args.bam_file_1 
bam2 = args.bam_file_2
refFasta = args.ref_fasta
snpFile = args.snp_file
output = args.output_file
outputGeno = args.output_geno

#checks for index for bam & reference FASTA files
checkForIndex(bam1 + ".bai") #will need + ".bai" & fasta will need + ".fa"
checkForIndex(bam2 + ".bai")
checkForIndex(refFasta + ".fai")

#create temp directory for munging
tempdir = tempfile.mkdtemp()
print(tempdir)

#execute bam-readcount on both bam files, returns rc .tsv files
rc1 = bamReadcount(bam1)
rc2 = bamReadcount(bam2)

#parse the rc return file for the data fields that we want
parsed1 = parseRc(rc1, 1)
parsed2 = parseRc(rc2, 2)

#run R on the fields to determine genotypes
geno_out_1 = geno_calc(parsed1, 1)
geno_out_2 = geno_calc(parsed2, 2)

#get genotype data from R results
genotypesDict = defaultdict(dict)
total1, snp1 = getGenotypes(geno_out_1, 1)
total2, snp2 = getGenotypes(geno_out_2, 2)

#write full genotype file in .bed format if specified in args
if outputGeno != None:
    with open(outputGeno, 'w') as genoFile:
        #write header
        genoFile.write("Chr" + "\t" + "Start" + "\t" + "Stop" + "\t" + "Genotype" + "\t" + "Sample1" + "\t" + "Sample2" + "\n")
        for key in genotypesDict:
            field = key.split('\t')
            print(field[0] + "\t" +  field[1] + "\t"  + field[1] + "\t" + genotypesDict[key]['samp1'] + "\t" + genotypesDict[key]['samp2'])
            genoFile.write(field[0] + "\t" +  field[1] + "\t"  + field[1] + "\t" + genotypesDict[key]['samp1'] + "\t" + genotypesDict[key]['samp2'] + "\n")
    genoFile.close()

#calculates matches
samplesWoComparison = 0
allSamples = total1 + total2
totalSamples = len(allSamples)
for samp in allSamples:
    if allSamples.count(samp) != 2:
        samplesWoComparison += 1
totalSamplesFixed = totalSamples - samplesWoComparison
totalMatchesPossible = totalSamplesFixed / 2
print(totalSamplesFixed)
print(totalMatchesPossible)

totalUniq = len(list(set(snp1 + snp2)))
print(totalUniq)

totalMatches = totalSamplesFixed - totalUniq

#write calculations to output file
with open(output, 'w') as outputFile:
    outputFile.write("%7s\t%22s\t%15s"%("Matches", "Total Matches Possible","Percent Matched"))
    outputFile.write("%7s\t%22s\t%15s"%(totalMatches, totalMatchesPossible, (totalMatches / totalMatchesPossible) * 100))

    print("%7s\t%22s\t%15s"%("Matches", "Total Matches Possible","Percent Matched"))
    print("%7s\t%22s\t%15s"%(totalMatches, totalMatchesPossible, (totalMatches / totalMatchesPossible) * 100))

outputFile.close()




