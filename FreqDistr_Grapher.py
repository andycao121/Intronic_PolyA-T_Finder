# Andy Cao

#import packages
import sys
import matplotlib.pyplot as plt
import numpy as np


# creates the normalized polyA/T frequency distribution graphs for specified chromosomes
def makeGraph(chr, polylist, type, bin):
    n, bins, patches = plt.hist(polylist, bins=bin, density=True, alpha=0.75)
    plt.xlabel("Sequence Length (bp)")
    plt.ylabel("Normalized Frequency")
    plt.title("Normalized Frequency Distribution of Poly" + type + " Sequences for " + chr)
    plt.xticks(np.arange(10, max(polylist)+1, 10.0))
    plt.savefig(chr + "_" + type + ".png")
    plt.close()


# converts string of numbers (separated by spaces) into list of integers representing sequence lengths
def convertList(lista, listt):
    lista, listt = [int(n) for n in lista], [int(n) for n in listt]
    tempa = []
    tempt = []
    for i, a in enumerate(lista):
        if a > 0:
            for j in range(a):
                tempa.append(i+1)
    for i, t in enumerate(listt):
        if t > 0:
            for j in range(t):
                tempt.append(i+1)
    return tempa, tempt


# determines number of bins in histogram
def findBins(polylist):
    temp = polylist[::-1]
    for i, t in enumerate(temp):
        if t != "0":
            return 186 - i


# "whole_genome" arg will output a frequency distribution of the whole genome
# chromosome labels (e.g. 1, 2, X, Y) will output multiple frequency distributions for each chromosome
file = open(sys.argv[1], "r")
if sys.argv[2] == "whole_genome":
    lines = [x.strip("\n") for x in file if x != "\n"]
    freqa = []
    freqt = []
    binsa = []
    binst = []
    for i, l in enumerate(lines):
        if "chr" in l:
            lista = lines[i+1].split()
            listt = lines[i+2].split()
            binsa.append(findBins(lista))
            binst.append(findBins(listt))
            tempa, tempt = convertList(lista, listt)
            freqa += tempa
            freqt += tempt
    makeGraph("Whole Genome", freqa, "A", max(binsa))
    makeGraph("Whole Genome", freqt, "T", max(binst))
else:
    chromosomes = sys.argv[2:]
    lines = [x.strip("\n") for x in file if x != "\n"]
    for i, l in enumerate(lines):
        if "chr" in l:
            chr = l
            lista = lines[i+1].split()
            listt = lines[i+2].split()
            if chr[3:] in chromosomes:
                bina = findBins(lista)
                bint = findBins(listt)
                lista, listt = convertList(lista, listt)
                makeGraph(chr, lista, "A", bina)
                makeGraph(chr, listt, "T", bint)
