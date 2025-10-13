from Bio import SeqIO

sequences = SeqIO.parse("02 rosalind_lcsm.txt", 'fasta')
all_substrings = []
common_substrings = set()
for sequence_data in sequences:
    all_substrings.append(str(sequence_data.seq))

shortest_len = len(min(all_substrings, key=len))

def substring(s,maxlen):
    l = []
    for i in range(0,len(s)+1):
        for j in range(i+1,len(s)+1):
            if(abs(i-j)>maxlen):
                continue
            l.append(s[i:j])
    return set(l)

for string in all_substrings:
    if len(common_substrings) == 0:
        common_substrings = substring(string,shortest_len)
        continue
    common_substrings = common_substrings.intersection(substring(string,shortest_len))

print(max(common_substrings,key=len))

