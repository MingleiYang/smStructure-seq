# from cbam2m5 import bam2m5, mapScore
# from BioUtil import samFile, xzopen, cachedFasta
import argparse
from Bio import SeqIO
import re

def parserM5(ref,m5File,output):
    seqs = SeqIO.read(open(ref),format = "fasta")
    outfile = open(output,"w")
    seqlen = len(seqs.seq)
    m5Reads = open(m5File,"r")
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    for read in m5Reads: ### the issue here is that the m5 file sometimes has more than one spaces, therefore the string.split() doesn't work well
        qName,qLength,qStart,qEnd,qStrand,tName,tLength,tStart,tEnd,tStrand,score,numMatch,numMismatch,numIns,numDel,mapQV,qAlignedSeq,matchPattern,tAlignedSeq = re.split(r'\s{1,}',read.strip())
        # print(qAlignedSeq + "\n" + matchPattern + "\n" +tAlignedSeq)
        if tStrand == "+":
            patL = "-" * int(tStart)
            patR = "-" * (int(seqlen) - int(tEnd))
            tseq = "".join(check_mutation(qAlignedSeq, matchPattern, tAlignedSeq))
            print(qName + "\t" + tName + "\t" + tStart + "\t" + tEnd + "\t" + patL + tseq + patR, file = outfile)
        else: # tStrand == "-"
            patR = "-" * int(tStart)
            patL = "-" * (int(seqlen) - int(tEnd))
            reverse_tseq = "".join(check_mutation(qAlignedSeq, matchPattern,tAlignedSeq))
            tseq = "".join(complement.get(base, base) for base in reversed(reverse_tseq))
            # print(tseq)
            print(qName + "\t" + tName + "\t" + tStart + "\t" + tEnd + "\t" + patL + tseq + patR, file=outfile)
    outfile.close()




def check_mutation(query,matchPat,reference):
    final = []
    for item in zip(query,matchPat,reference):
        if item[1] == "|":
            final.append(item[0])
        elif item[1] == "*" and item[2] == "-":
            continue
        elif item[1] == "*" and item[0] == "-":
            final.append("D")
        else:
            final.append("E")
    return(final)



def main():
    parser = argparse.ArgumentParser(description = "Convert pacbio m5 format file to simple m5 for ploting")
    parser.add_argument("fasta", metavar="ref.fa", help="reference file")
    parser.add_argument("inM5", metavar = "in.m5", help = "input M5 file")
    parser.add_argument("outM5", metavar="out.m5.fixed", help = "output simple m5 file ")
    args = parser.parse_args()
    parserM5(args.fasta, args.inM5, args.outM5)

if __name__ == '__main__':
    main()
