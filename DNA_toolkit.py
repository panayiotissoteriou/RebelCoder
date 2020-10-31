import random
randDNAstr = ''.join([random.choice(Nucleotides)
                        for nuc in range(20)])          # random DNA seq with 20 letters

# Validating a DNA seq
Nucleotides = ['A', 'C', 'G', 'T']

def validateSeq(dna_seq):
    tmpseq = dna_seq
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

# countFrequency:
def countNucFrequency(seq):
    tmpFreqDict = {'A':0, "C":0, 'G':0, 'T':0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict

# Transcritpion: replacing T with U
def transcription(seq):
    return seq.replace('T', 'U')

# Reverse Complement - using dictionary
def RevComplement(seq):
    empty_string = ''
    DNA_reverseComplement = {"A":"T", "C":"G", "G":"C", "T":"A"}
    for i in seq[::-1]:
        empty_string += DNA_reverseComplement[i]
    return empty_string


# GC content:
def gc_content(seq):
    GCcount = 0
    for i in seq:
        if i == "C" or i == "G":
            GCcount += 1
    gcContent = float(GCcount)/float(len(seq)) * 100
    return gcContent

# GC content of subsections of length k
def GCcontent_subsect(seq, k):
    res = []
    for i in range(0, len(seq)-k + 1, k):
        subsec = seq[i : i+k]
        res.append(gc_content(subsec))
    return res

# DNA codons
DNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
}

# RNA codons
RNA_Codons = {
    # 'M' - START, '_' - STOP
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UGU": "C", "UGC": "C",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "UUU": "F", "UUC": "F",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAU": "H", "CAC": "H",
    "AUA": "I", "AUU": "I", "AUC": "I",
    "AAA": "K", "AAG": "K",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUG": "M",
    "AAU": "N", "AAC": "N",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UGG": "W",
    "UAU": "Y", "UAC": "Y",
    "UAA": "_", "UAG": "_", "UGA": "_"
}


# Translation
def translate_seq(seq, init_pos=0):
    prot_seq = []
    for pos in range(init_pos, len(seq) -2, 3):         # range(start at pos 0, end at pos 2 from end, 3-letter steps)
        prot_seq.append(DNA_Codons[seq[pos:pos+3]])
    return prot_seq

# Codon usage
import collections

def codon_usage(seq, aminoacid):
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i : i+3]] == aminoacid:
            tmpList.append(seq[i:i+3])

    freqDict = dict(collections.Counter(tmpList))
    TotalWight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / TotalWight, 2)
    return freqDict

# Open reading frames:
def gen_reading_frame(seq):
    frames = []
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq, 1))
    frames.append(translate_seq(seq, 2))
    frames.append(translate_seq(RevComplement(seq), 0))
    frames.append(translate_seq(RevComplement(seq), 1))
    frames.append(translate_seq(RevComplement(seq), 2))
    return frames

# Part 6: looking for a start codon
# computes all possible proteins in an AA seq and returns a list of possible proteins
def prot_from_RF(aa_seq):
    current_prot = []
    proteins = []                           # accumulates all proteins in the sequence
    for aa in aa_seq:
        if aa == '_':                       # stop accumulating AAs if '_' is found
            if current_prot:                # checks if current_prot contains anything - needed if multiple 'M' are present
                for p in current_prot:
                    proteins.append(p)
                current_prot = []           # empties current_prot, so that it can start accumulate for 2nd possible prot seq
        else:
            if aa == "M":                           # start accumulating AAs if "M" is found
                 current_prot.append("")            # starts a string with "M" as first letter
            for i in range(len(current_prot)):
                current_prot[i] += aa               #adds the rest of AAs into that string, after M, until "_" is encountered
    return proteins

#Part 7: Function that finds proteins from DNA seqs
#must do a few things:
# Generate all ORFs
# Extract all proteins
# Return a sorted/unsorted list
# compute all possible proteins for all ORFs
def all_prots_from_ORFs(seq, startReadPos=0, endReadPos=0, ordered=False):
    if endReadPos > startReadPos:                           # this is to specify if we want RFs from the whole seq, or from a range
        rfs = gen_reading_frame(seq[startRead:endRead])
    else:
        rfs = gen_reading_frame(seq)

    res = []
    for rf in rfs:
        prots = prot_from_RF(rf)        #from each RF, a prot is generated
        for p in prots:
            res.append(p)               #protein generated is added to list res

    if ordered:
        return sorted(res, key=len, reverse=True)   #sorts the list res, by length of elements, from largest to smallest
    return res

print("ORFs:", all_prots_from_ORFs("ACGTGTTAGGAGTGCAGGTTACTCGACTGTCATAAGTGATCGCCCTCGATAACGAGGTAGAGGAGG", 0, 0, True))


# Utilities - EXTRA useful code

# Adding colour
def coloured(seq):
    bcolours = {
            'A' : '\033[92m',
            'C' : '\033[94m',
            'G' : '\033[93m',
            'T' : '\033[91m',
            'U' : '\033[91m',
            'reset' : '\033[0;0m'
    }

    tempString = ""
    for i in seq:
        if i in bcolours:
            tempString += bcolours[i] + i
        else:
            tempString += bcolours[reset] + i

    return tempString + '\033[0;0m]'


def readTextFile(filePath):
    with open(filePath, 'r') as f:
        return "".join([l.strip() for l in f.readlines()])


def writeTextFile(filePath, seq, mode='w'):
    with open(filePath, mode) as f:
        f.write(seq + '\n')


def read_FASTA(filePath):
    with open(filePath, 'r') as f:
        FASTAFile = [l.strip() for l in f.readlines()]

    FASTADict = {}
    FASTALabel = ""

    for line in FASTAFile:
        if '>' in line:
            FASTALabel = line
            FASTADict[FASTALabel] = ""
        else:
            FASTADict[FASTALabel] += line

    return FASTADict
