# Andy Cao

# import packages
import sys

# create output files
within = open(sys.argv[3] + "_within.txt", "w")
partial = open(sys.argv[3] + "_partial_within.txt", "w")
summary = open(sys.argv[3] + "_intron_summary.txt", "w")
freq = open(sys.argv[3] + "_intron_freq.txt", "w")


# creates two output files: one for polyA/T sequences fully in introns and one for those partly in introns
def findOverlap(atdict, introndict, atlines):
    summdict = {"A": [0, 0], "T": [0, 0]}
    freqdict = {"A": [[0]*200, [0]*200], "T": [[0]*200, [0]*200]}
    for chr in introndict.keys():
        if chr in atdict.keys():
            maxcount = len(atdict[chr])
            count = 0
            for seq in introndict[chr]:
                while count < maxcount:
                    seqstart = atdict[chr][count][0]
                    seqend = atdict[chr][count][1]
                    seqtype = atdict[chr][count][2]
                    if seq[0] <= seqstart and seqend <= seq[1]:
                        within.write(atlines[chr][count])
                        summdict[seqtype][0] += 1
                        freqdict[seqtype][0][seqend-seqstart] += 1
                        count += 1
                    elif seq[0] > seqend:
                        count += 1
                    elif seq[1] < seqstart:
                        break
                    else:
                        partial.write(atlines[chr][count])
                        summdict[seqtype][1] += 1
                        freqdict[seqtype][1][seqend - seqstart] += 1
                        count += 1
    summary.write("PolyA\t" + str(summdict["A"][0]) + "\t" + str(summdict["A"][1]) + "\t" + str(sum(summdict["A"])) + "\n")
    summary.write("PolyT\t" + str(summdict["T"][0]) + "\t" + str(summdict["T"][1]) + "\t" + str(sum(summdict["T"])) + "\n")
    freq.write("PolyA\n")
    for i in freqdict["A"][0]:
        freq.write(str(i) + " ")
    freq.write("\n")
    for p in freqdict["A"][1]:
        freq.write(str(p) + " ")
    freq.write("\n\n")
    freq.write("PolyT\n")
    for i in freqdict["T"][0]:
        freq.write(str(i) + " ")
    freq.write("\n")
    for p in freqdict["T"][1]:
        freq.write(str(p) + " ")


# returns a dictionary of lists of start and end locations of the sequences in BED file
def organize(file):
    dict = {}
    seq = open(file, "r")
    lines = [x.strip("\n") for x in seq if x != "\n"]
    for l in lines:
        temp = l.split()
        if temp[0] not in dict.keys():
            dict[temp[0]] = []
        seqtype = "A"*(temp[3] == "polyA") or "T"
        dict[temp[0]].append((int(temp[1]), int(temp[2]), seqtype))
    for chr in dict.keys():
        dict[chr] = sorted(dict[chr])
    return dict


# returns the BED file lines in a dict
def lines(file):
    dict = {}
    seq = open(file, "r")
    lines = [x for x in seq if x != "\n"]
    for l in lines:
        temp = l.split()
        if temp[0] not in dict.keys():
            dict[temp[0]] = []
        dict[temp[0]].append(l)
    return dict


polyat = sys.argv[1]
introns = sys.argv[2]
atdict = organize(polyat)
introndict = organize(introns)
atlines = lines(polyat)
findOverlap(atdict, introndict, atlines)
