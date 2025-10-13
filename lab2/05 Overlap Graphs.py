from Bio import SeqIO

sequences = SeqIO.parse("05 rosalind_grph.txt", 'fasta')

adjacency_set = set()
samples = []
k = 3
for idx, sequence_data in enumerate(sequences):
    sample = {"idx":idx,"name":sequence_data.name,"seq":sequence_data.seq}
    sequence_data.idx = idx
    sample['sufix'] = (sample['seq'][-k:])
    sample['prefix'] = (sample['seq'][:k])
    samples.append(sample)

for sample1 in samples:
    for sample2  in samples:
        id1 = sample1['idx']
        id2 = sample2['idx']
        if(id1 == id2):
            continue
        if(sample1['sufix'] == sample2['prefix']):
            adjacency_set.add((id1,id2))

# print(samples)
# print(adjacency_set)


for x ,y in adjacency_set:
    print(f"{samples[x]['name']} {samples[y]['name']}")
    