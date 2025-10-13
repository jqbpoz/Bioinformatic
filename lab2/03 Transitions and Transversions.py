from Bio import SeqIO

sequences = SeqIO.parse("03 rosalind_tran.txt", 'fasta')
all_substrings = []
for sequence_data in sequences:
    all_substrings.append(str(sequence_data.seq))

transitions = ["AG","GA","CT","TC"]
tranvertions = ["AC","CA","AT","TA","CG","GC","GT","TG"]

first,second = all_substrings
transitions_sum =0
tranvertions_sum =0
for f,s in zip(first,second):
    mutation = f+s
    if mutation in tranvertions:
        tranvertions_sum +=1
    if mutation in transitions:
        transitions_sum +=1

print(transitions_sum/tranvertions_sum)
