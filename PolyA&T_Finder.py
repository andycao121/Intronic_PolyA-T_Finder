# Andy Cao

# import packages
from Bio import SeqIO
import sys
import glob


# create output BED file and summary file
output = open(sys.argv[4] + "_output.txt", "w")
summary = open(sys.argv[4] + "_summary.txt", "w")
freq = open(sys.argv[4] + "_freq.txt", "w")


# takes a fasta/fastq file, sequence length, # mismatches
# prints out priming ranges in the form of x-y inclusive
# first nucleotide is labelled 1
def findInternalPrimingSites(file, length, mismatches):
    filetype = file[file.index(".") + 1:]
    for seq in SeqIO.parse(file, filetype):
        polyAmm, freqA = locatePolyA(seq.seq.upper(), length, mismatches, seq.id)
        polyTmm, freqT = locatePolyT(seq.seq.upper(), length, mismatches, seq.id)
        summary.write(seq.id + "\n# Poly-A\t")
        for i, n in enumerate(polyAmm):
            summary.write(str(n) + " (" + str(i) + ")\t")
        summary.write(str(sum(polyAmm)) + "\n# Poly-T\t")
        for i, n in enumerate(polyTmm):
            summary.write(str(n) + " (" + str(i) + ")\t")
        summary.write(str(sum(polyTmm)) + "\n\n")
        freq.write(seq.id + "\n")
        for a in freqA:
            freq.write(str(a) + " ")
        freq.write("\n")
        for t in freqT:
            freq.write(str(t) + " ")
        freq.write("\n\n")


# finds all ranges of polyA within a specific sequence
def locatePolyA(sequence, length, mismatches, id):
    mmlist = [0] * (mismatches + 2)
    freqA = [0] * 200
    numA = length - mismatches
    start = 0
    while start <= len(sequence) - length:
        tempseq = sequence[start:start + length]
        if "N" in tempseq:
            start += tempseq.rindex("N") + 1
        else:
            countA = tempseq.count("A")
            if countA < numA:
                start = start + newStart(tempseq, mismatches, "A")
            else:
                countMM = length - countA
                end = start + length
                while end < len(sequence) - 1 and sequence[end] == "A":
                    end += 1
                if sequence[start] != "A":
                    if end - start > length:
                        start += 1
                        countMM -= 1
                    if end - start == length and end < len(sequence) - 1 and "N" != sequence[end] != "A":
                        end += 1
                        countMM += 1
                output.write(id + "\t" + str(start + 1) + "\t" + str(end) + "\tpolyA\t" + str(countMM) + "\n")
                mmlist[countMM] += 1
                freqA[end-start-1] += 1
                start = end
    return mmlist, freqA


# finds all ranges of polyT within a specific sequence
def locatePolyT(sequence, length, mismatches, id):
    mmlist = [0] * (mismatches + 2)
    freqT = [0] * 200
    numT = length - mismatches
    start = 0
    while start <= len(sequence) - length:
        tempseq = sequence[start:start + length].upper()
        if "N" in tempseq:
            start += tempseq.rindex("N") + 1
        else:
            countT = tempseq.count("T")
            if countT < numT:
                start = start + newStart(tempseq, mismatches, "T")
            else:
                countMM = length - countT
                end = start + length
                while end < len(sequence) - 1 and sequence[end] == "T":
                    end += 1
                if sequence[start] != "T":
                    if end - start > length:
                        start += 1
                        countMM -= 1
                    if end - start == length and end < len(sequence) - 1 and "N" != sequence[end] != "T":
                        end += 1
                        countMM += 1
                output.write(id + "\t" + str(start + 1) + "\t" + str(end) + "\tpolyT\t" + str(countMM) + "\n")
                mmlist[countMM] += 1
                freqT[end-start-1] += 1
                start = end
    return mmlist, freqT


# determines location of first A/T after (n+1)th to last non-A/T in shortened sequence (n = # mismatches)
def newStart(sequence, mismatches, base):
    count = 1
    for index, bp in enumerate(sequence):
        if bp != base:
            count += 1
            if count >= mismatches:
                return index + 1


directory = sys.argv[1]
length = int(sys.argv[2])
mm = int(sys.argv[3])
for file in glob.glob(directory + "\*.fasta"):
    print(file)
    findInternalPrimingSites(file, length, mm)
