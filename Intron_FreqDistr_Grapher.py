# Andy Cao

#import packages
import sys
import matplotlib.pyplot as plt
import numpy as np


# creates the normalized intronic polyA/T frequency distribution graphs for specified chromosomes
def makeGraph(polylist, type, bin):
    legend = ["Inside", "Partial"]
    n, bins, patches = plt.hist(polylist, bins=bin, density=True, alpha=0.75, stacked=True, label=legend)
    plt.xlabel("Sequence Length (bp)")
    plt.ylabel("Normalized Frequency")
    plt.title("Normalized Frequency Distribution of Intronic Poly" + type + " in Whole Genome")
    plt.xticks(np.arange(10, max(polylist[0])+1, 10.0))
    plt.legend(loc="upper right")
    plt.savefig("Intronic_Poly_" + type + ".png")
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


file = open(sys.argv[1], "r")
lines = [x.strip("\n") for x in file if x != "\n"]
freqi = []
freqp = []
bins = []
for i, l in enumerate(lines):
    if "Poly" in l:
        listi = lines[i+1].split()
        listp = lines[i+2].split()
        bins.append(findBins(listi))
        bins.append(findBins(listp))
        freqi, freqp = convertList(listi, listp)
        makeGraph([freqi, freqp], l[4], max(bins))
        freqi = []
        freqp = []
        bins = []
