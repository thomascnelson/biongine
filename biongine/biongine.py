# welcome to the biongine module
''' This is a useful set of biology related python functionality.
It is an ongoing project, please contact me if you are interested
in contributing.'''

import re
import pandas as pd

# class to represent nucleotide sequences
class Sequence:
    def __init__(self, sequence, name="", identifier=""):
        self.sequence = str(sequence).upper()
        self.name = name
        self.identifier = identifier
        if 'U' in sequence:
            self.sequence_type = "RNA"
        else:
            self.sequence_type = "DNA"
        
    def reverse_complement(self):
        if self.sequence_type == "DNA":
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            reverse_comp_seq = [complement[base] for base in self.sequence[::-1]]
            return "".join(reverse_comp_seq)
        else:
            complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
            reverse_comp_seq = [complement[base] for base in self.sequence[::-1]]
            return "".join(reverse_comp_seq)

    def translate(self):
        if self.sequence_type == "RNA":
            sequence = self.sequence.replace('U', 'T')
        else:
            sequence = self.sequence
        codon_table = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        protein_seq = [codon_table[sequence[i:i+3]] for i in range(0, len(sequence), 3)]
        return "".join(protein_seq)

    def gc_content(self):
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        return (gc_count / len(self.sequence)) * 100 if len(self.sequence) > 0 else 0

    def nucleotide_counts(self):
        if self.sequence_type == "RNA":
            counts = {
                'A': self.sequence.count('A'),
                'U': self.sequence.count('U'),
                'C': self.sequence.count('C'),
                'G': self.sequence.count('G')
            }
            return counts
        else:
            counts = {
                'A': self.sequence.count('A'),
                'T': self.sequence.count('T'),
                'C': self.sequence.count('C'),
                'G': self.sequence.count('G')
            }
            return counts

    def prepend(self, to_add):
        self.sequence = to_add + self.sequence

    def append(self, to_add):
        self.sequence += to_add

    def transcribe(self):
        if self.sequence_type == "DNA":
            return self.sequence.replace('T', 'U')
        else:
            return self.sequence  # If it's already RNA, no change

    def similarity(self, other_sequence):
        # compute similarity with simple dynamic programming method
        # sequences can be of different lengths
        '''This method uses dynamic programming to compute the edit distance 
        between the two sequences and then calculates the similarity based on 
        the edit distance and the length of the longer sequence.'''
        if self.sequence_type == "RNA":
            sequence =  self.sequence.replace('U', 'T')
        else:
            sequence = self.sequence
        
        if isinstance(other_sequence, Sequence):
            other_sequence = other_sequence.sequence
        
        m = len(sequence)
        n = len(other_sequence)

        # Create a 2D matrix to store alignment scores
        dp = [[0] * (n + 1) for _ in range(m + 1)]

        # Initialize the first row and column
        for i in range(m + 1):
            dp[i][0] = i
        for j in range(n + 1):
            dp[0][j] = j

        # Fill in the matrix using dynamic programming
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if sequence[i - 1] == other_sequence[j - 1]:
                    cost = 0
                else:
                    cost = 1
                dp[i][j] = min(dp[i - 1][j - 1] + cost, dp[i - 1][j] + 1, dp[i][j - 1] + 1)

        # Compute similarity based on the alignment score
        max_length = max(m, n)
        similarity = 1.0 - (dp[m][n] / max_length)
        return similarity
    
    def find(self, pattern):
        matches = [match.start() for match in re.finditer(pattern.upper(), self.sequence)]
        if len(matches) != 0:
            return matches
        else:
            return None        
    
    # override the == operator to compare sequences themselves
    def __eq__(self, other):
        if isinstance(other, Sequence):
            return self.sequence.replace('U', 'T') == other.sequence.replace('U', 'T')
        return False
        
    def __add__(self, other_sequence):
        if isinstance(other_sequence, Sequence):
            if self.sequence_type != other_sequence.sequence_type:
                raise TypeError("Sequences are not the same type, cant concatenate DNA and RNA")
            else:
                new_sequence = self.sequence + other_sequence.sequence
                new_name = f"{self.name} + {other_sequence.name}"
            return Sequence(new_sequence, new_name)
        else:
            raise TypeError("Unsupported operand type. You can only concatenate two DNASequence objects.")

    def __str__(self):
        return f"Name: {self.name}\nIdentifier: {self.identifier}\nSequence Type: {self.sequence_type}\nSequence: {self.sequence}"



# function to read a fasta file and return it as a dataframe
def fasta_to_dataframe(file_name):
    # Create empty lists to store identifiers and sequences
    identifiers = []
    sequences = []

    # Read the FASTA file and parse it
    with open(file_name, 'r') as file:
        identifier = None
        sequence = ""

        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if identifier is not None:
                    identifiers.append(identifier)
                    sequences.append(sequence)
                identifier = line[1:]
                sequence = ""
            else:
                sequence += line

        # Append the last sequence
        if identifier is not None:
            identifiers.append(identifier)
            sequences.append(sequence)

    # Create a Pandas DataFrame
    df = pd.DataFrame({'Identifier': identifiers, 'Sequence': sequences})

    return df










