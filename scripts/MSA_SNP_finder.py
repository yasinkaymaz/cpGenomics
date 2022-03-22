#!/usr/bin/python
"""
Determining single nucleotide variant locations from multiple sequence alignment files.
"""
import os
import pandas  as pd
import sys
from Bio import AlignIO
import Bio.Align
dir = os.path.dirname(__file__)
print(os.getcwd())
print(dir)
from collections import defaultdict, Counter
#from collections import namedtuple
import collections

alndata = []

Ref = snakemake.params.Refname
Repeat_inputFile = snakemake.params.repf


#Find the genomic locations fall into miropeats repeat regions
RepeatList = []
with open(Repeat_inputFile) as repeatFile:
    for line in repeatFile:
        repstart = int(line.strip().split("\t")[1])
        repend = int(line.strip().split("\t")[2])
        for i in range(repstart,repend):
            RepeatList.append(i)
RepeatList = set(RepeatList)
print(len(set(RepeatList)))

outfile1 = open(snakemake.input.msaf+"_SNP_positions.bed", "w")

tmpoutfile = open(snakemake.input.msaf+".tmp.file.aln", "w")
alignment = AlignIO.read(open(snakemake.input.msaf), "fasta")
for record in alignment:
	#print(record.id, record.seq)
    tmpoutfile.write(str(record.id)+"\t"+str(record.seq)+"\n")
tmpoutfile.close()

#with open(sys.argv[3], "r") as alnfile:
with open(snakemake.input.msaf+".tmp.file.aln", "r") as alnfile:
    df = pd.read_csv(alnfile,sep="\t", header=None)
    dfseq = []
    # take the names of sequences to index
    index = df[0]
    str_seq1 = df.iloc[1][1]
    print("len : " +  str(len(str_seq1)))
    # create an empty list of length of sequence letters
    columns = list(range(len(str_seq1)))
    #new data frame
    new_df = pd.DataFrame(index=index, columns= columns)
    #print(new_df)
    #create new dataframe
    for i in range(len(df)):
        # change a string of sequence to a list of sequence
        dfseq = list(df.iloc[i][1])
        # You can edit a subset of a dataframe by using loc:
        # df.loc[<row selection>, <column selection>]
        # place the list of sequence to the row that it belongs to in the new data frame
        new_df.loc[i:,:] = dfseq

    # for each row:
    Ref_pos = 0
    MatchCount = 0
    MissMatchCount = 0
    SeqCount=len(df)
    for i in range(len(str_seq1)):
        # place the the sequence in a position to a set
        #print(collections.Counter(new_df.loc[:,i]), collections.Counter(new_df.loc[:,i]).most_common(4))
        nuccompos = collections.Counter(new_df.loc[:,i])
    #    print("n_count:", nuccompos['n'], "gap_count:", nuccompos['-'])
        del nuccompos['n']
        del nuccompos['-']
        #print(nuccompos.most_common(4))

        set_seq = set(new_df.loc[:,i])

        if new_df.loc[Ref,i] != '-':
            Ref_pos = Ref_pos +1
        else:
            pass
    	#if there is a mismatch error and this position is not in repeat regions, count as sequencing error.
        if len(nuccompos.most_common(4))>1 and Ref_pos not in RepeatList and float(100*nuccompos.most_common(4)[1][1]/SeqCount) > 10.00:
#    		print(nuccompos.most_common(4), Ref_pos, SeqCount, 100*nuccompos.most_common(4)[0][1]/SeqCount,100*nuccompos.most_common(4)[1][1]/SeqCount)
            MissMatchCount = MissMatchCount +1
            outfile1.write(str(Ref)+"\t"+str(Ref_pos)+"\t"+str(Ref_pos+1)+"\t"+str(1)+"\n")

    	#If there is a perfect match and the position is not in the repeat regions, count as match.
        elif len(nuccompos.most_common(4)) == 1 and Ref_pos not in RepeatList:
            MatchCount = MatchCount +1
#outfile2.write(sys.argv[1]+"\t"+sys.argv[2]+"\t"+str(MissMatchCount)+"\t"+str(MatchCount)+"\t"+str(float(MissMatchCount)/(MissMatchCount+MatchCount)))
print("Match Count:", MatchCount, Ref_pos)
outfile1.close()
outfile2 = open(snakemake.input.msaf+"_SNP_positions.hist","w")

snp_list = []
with open(snakemake.input.msaf+"_SNP_positions.bed") as snpBedfile:
    for line in snpBedfile:
	#	chrom = str(line.strip().split("\t")[0])
        snpstart = int(line.strip().split("\t")[1])
		#snpend = int(line.strip().split("\t")[2])
		#snp = int(line.strip().split("\t")[3])
        snp_list.append(snpstart)
		# if snpstart%10 != 0:
		# 	snp_histcount += 1
		# else:
		# 	outfile2.write(str(chrom)+"\t"+str(snpstart-10)+"\t"+str(snpstart)+"\t"+str(snp_histcount)+"\n")
		# #	print(chrom, snpstart-10, snpstart, snp_histcount)
		# 	snp_histcount = 0

snp_histcount = 0
for i in range(Ref_pos):
    if i in set(snp_list):
        snp_histcount += 1
    else:
        snp_histcount = snp_histcount
    if i%100 == 0 and i != 0:
        print(Ref, i-0, i, snp_histcount)
        outfile2.write(str(Ref)+"\t"+str(i-100)+"\t"+str(i)+"\t"+str(snp_histcount)+"\n")
        snp_histcount = 0

outfile2.close()


#import numpy as np
#import matplotlib.pyplot as plt
