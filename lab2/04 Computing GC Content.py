from Bio import SeqIO

sequences = SeqIO.parse("04 rosalind_gc.txt", 'fasta')
sequences_data = []
for sequence_data in sequences:
    sample = {"name":sequence_data.name,"len":len(sequence_data.seq),
            "CG_counts":sequence_data.seq.count("C") + sequence_data.seq.count("G")}
    sample["CG_percentage"]=(int(sample["CG_counts"])/int(sample["len"])*100)

    sequences_data.append(sample)

max_sample = max(sequences_data,key=lambda x :x['CG_percentage'])
print(max_sample['name'])
print(max_sample['CG_percentage'])