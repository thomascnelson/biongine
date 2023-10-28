import biongine as bn

# Create Sequence objects
seq1 = bn.Sequence("ATCGTAGCTAGCATGCTAGCTAGC", name="Sample DNA", identifier="ID1")
seq2 = bn.Sequence("ATCGTAGCTAGCATGCTAGCTAGC", name="Sample DNA2", identifier="ID2")
seq3 = bn.Sequence("AUCGUAGCUAGCAUGCUAGCUAGC", name="Sample RNA", identifier="ID3")

# print prints a sane representation
print(seq1)
print(seq2)
print(seq3)

# comparison == operator overridden to compare sequence only
seq1 == seq2
seq1 is seq2
seq1.name == seq2.name
seq1 == seq3 #autoconverts RNA sequene to DNA for direct comparison

# attributes and methods
print("Original Sequence:", seq1.sequence)
print("Reverse Complement:", seq1.reverse_complement())
print("Reverse Complement:", seq3.reverse_complement())
print("Translated Protein Sequence:", seq1.translate())
print("GC Content:", seq1.gc_content())
print("Nucleotide Counts:", seq1.nucleotide_counts())
print("Nucleotide Counts:", seq3.nucleotide_counts())
print("Transcribed RNA Sequence:", seq1.transcribe())
print("Transcribed RNA Sequence:", seq3.transcribe())

# Adding sequences to the start and end
seq1.prepend("AAAAAAAA")
seq1.append("GGGGGGGG")
print(seq1)

# concatenation
concatenated_sequence = seq1 + seq2
print(concatenated_sequence.sequence)
print(concatenated_sequence.name)
concatenated_sequence2 = seq1 + seq3 # can't concatenate DNA and RNA

# Computing similarity with another sequence
similarity = seq1.similarity(seq2.sequence)
print("Similarity to Other Sequence:", similarity)

# Example usage for the fasta reader
# Reads the file and returns a pandas dataframe
file_name = 'transcripts.fasta'
df = bn.fasta_to_dataframe(file_name)
print(df)

df.loc[df['Identifier'] == 'NM_017409'].Sequence.values[0]
df.loc[df['Identifier'] == 'NM_017409'].Identifier.values[0]

new_seq = bn.Sequence(df.loc[df['Identifier'] == 'NM_017409'].Sequence.values[0], 
                      name="Some gene name", 
                      identifier=df.loc[df['Identifier'] == 'NM_017409'].Identifier.values[0])
print(new_seq)
new_seq.gc_content()
new_seq.reverse_complement()

new_seq2 = bn.Sequence(df.loc[df['Identifier'] == 'NM_153693'].Sequence.values[0], 
                      name="Some other gene name", 
                      identifier=df.loc[df['Identifier'] == 'NM_153693'].Identifier.values[0])

new_seq.similarity(new_seq2)

new_seq.find('tccAA')

new_seq.sequence[257:277]


