from Bio import SeqIO

def nucleotides_to_codons(s,offset=0):
    nucleotides_to_codons = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L"
        , "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M", "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V"
        , "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P"
        , "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A"
        , "UAU": "Y", "UAC": "Y", "UAA": "_", "UAG": "_", "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q"
        , "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E"
        , "UGU": "C", "UGC": "C", "UGA": "_", "UGG": "W", "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R"
        , "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
    codons = []
    for i in range(0+offset, len(s)+offset, 3):
        nucle = s[i:i + 3]
        if len(nucle) < 3:
            break
        codon = nucleotides_to_codons[nucle]
        # if codon == "_":
        #     break
        codons.append(codon)

    return "".join(codons)

def reverse_complement(s):
    complement = str.maketrans("ATGC", "TACG")
    return s.translate(complement)[::-1]


def findORFs(s):
    start_indexes =[i for i,x in enumerate(s) if x =="M"]
    end_indexes = [i for i,x in enumerate(s) if x =="_"]
    pairs_indexes = set()
    for i in start_indexes:
        for j in range(0,len(end_indexes)):
            if(i>=end_indexes[j]):
                continue
            if(j == 0):
                pairs_indexes.add((i,end_indexes[j]))
            elif(i > end_indexes[j-1]):
                pairs_indexes.add((i, end_indexes[j]))
    return {s[start_index:end_index] for start_index, end_index in pairs_indexes}


sequences = list(SeqIO.parse("06 rosalind_orf.txt", 'fasta'))
s = str(sequences[0].seq)
all_options = []

# 1
RNA = s.replace("T", "U")
for i in range(3):
    proteins = nucleotides_to_codons(RNA, i)
    all_options.append(proteins)

# 2
rev_RNA = reverse_complement(s).replace("T", "U")
for i in range(3):
    proteins = nucleotides_to_codons(rev_RNA, i)
    all_options.append(proteins)

result_set = set()
for option in all_options:
    result_set |= findORFs(option)

for protein in result_set:
    print(protein)



