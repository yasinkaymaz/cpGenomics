#!/usr/bin/python

import os
import pandas  as pd
import sys
from Bio import AlignIO
import Bio.Align
dir = os.path.dirname(__file__)
print(os.getcwd())
print(dir)
alndata = []

RefName=snakemake.params.RefName
GenomeName = snakemake.params.species

InputAlign=snakemake.input.msaf
Repeat_inputFile = snakemake.params.repf
RepeatList = []
with open(Repeat_inputFile) as repeatFile:
    for line in repeatFile:
        repstart = int(line.strip().split("\t")[1])
        repend = int(line.strip().split("\t")[2])
        for i in range(repstart,repend):
            RepeatList.append(i)
RepeatList = set(RepeatList)
print(len(set(RepeatList)))

outfile1 = open(snakemake.output.varpos, "w")

tmpoutfile = open(snakemake.output.tmpaln, "w")
alignment = AlignIO.read(open(InputAlign), "fasta")
for record in alignment:
    tmpoutfile.write(str(record.id)+"\t"+str(record.seq)+"\n")
tmpoutfile.close()
#
#
UncoveredPos_1=[]
with open(snakemake.output.tmpaln, "r") as alnfile:
    df = pd.read_csv(alnfile, sep="\t",header=None)
    dfseq = []
    # take the names of sequences to index
    index = df[0]
    str_seq1 = df.iloc[1][1]
    print("len : " +  str(len(str_seq1)))
#     # create an empty list of length of sequence letters
    columns = list(range(len(str_seq1)))
#     #new data frame
    new_df = pd.DataFrame(index=index, columns= columns)
    #create new dataframe
    for i in range(len(df)):
#         # change a string of sequence to a list of sequence
        dfseq = list(df.iloc[i][1])
        new_df.loc[i:,:] = dfseq

    Ref_pos = 0
#    MatchCount = 0
    MissMatchCount = 0
    for i in range(len(str_seq1)):
#         # place the the sequence in a position to a set
        set_seq = set(new_df.loc[:,i])

        if new_df.loc[RefName, i] != '-':
            Ref_pos = Ref_pos +1
        else:
            pass
#         #put the bases at the location of the query sequence 1 and 2 into a set
        set_pair = set(new_df.loc[ [GenomeName, RefName ],i  ])

        if 'N' in set_pair or 'n' in set_pair:
            UncoveredPos_1.append(Ref_pos)
            print(set_pair)
        #if there is a mismatch error and this position is not in repeat regions, count as sequencing error.
        elif len(set_pair)>1 and Ref_pos not in RepeatList and '-' not in set_pair:
            print(set_pair, Ref_pos)
            MissMatchCount = MissMatchCount +1
            outfile1.write(str(RefName)+"\t"+str(Ref_pos)+"\t"+str(Ref_pos+1)+"\n")

outfile1.close()
print("Now smoothing")
outfile1 = open(snakemake.output.smvar, "w")

#smoothing:
#Define a window, typically 100 bases.
win=600
#Define a smoothing interval.
sm=200
genomeLen=Ref_pos

#Read the list of genomic positions that are mismatch between the aligned two genomes. Relative to Reference.
MismatchPositions =[]
with open(snakemake.output.varpos, "r") as snpfile:
    for line in snpfile:
        MismatchPositions.append(int(line.strip().split("\t")[1]))

outfile1.write("Chrom"+"\t"+"start"+"\t"+"end"+"\t"+"Sim2Ref"+"\n")
for i in range(0,genomeLen-win+sm,sm):
    windowPoss = range(i, i+win)
    dissim=len(set(windowPoss)&set(MismatchPositions))
    CoveredLen=win-len(set(windowPoss)&set(UncoveredPos_1))+1
    PercentDis = 100*float(dissim)/CoveredLen

    outfile1.write(RefName+"\t"+str(i)+"\t"+str(i+sm)+"\t"+str(PercentDis)+"\n")

outfile1.close()
