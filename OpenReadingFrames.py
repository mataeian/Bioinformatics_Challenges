# Problem
# Either strand of a DNA double helix can serve as the coding strand for RNA transcription. Hence, a given DNA string implies six total reading frames, or ways in which the same region of DNA can be translated into amino acids: three reading frames result from reading the string itself, whereas three more result from reading its reverse complement.

# An open reading frame (ORF) is one which starts from the start codon and ends by stop codon, without any other stop codons in between. Thus, a candidate protein string is derived by translating an open reading frame into amino acids until a stop codon is reached.

# Given: A DNA string s of length at most 1 kbp in FASTA format.

# Return: Every distinct candidate protein string that can be translated from ORFs of s. Strings can be returned in any order.


from Bio.Seq import Seq
from Bio import SeqIO

# Define the codon table
codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}

def translate_dna(dna):
    protein = ""
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        if codon in codon_table:
            protein += codon_table[codon]
    return protein


def find_orfs(dna):
    orfs = set() #set to remove duplicates
    for i in range(3):
        for j in range(i, len(dna) - 2, 3):
            codon = dna[j:j+3]
            if codon == "ATG":
                protein = ""
                for k in range(j, len(dna) - 2, 3):
                    codon = dna[k:k+3]
                    if codon in ["TAA", "TAG", "TGA"]:
                        orfs.add(protein)
                        break
                    protein += codon_table.get(codon, '')
    return orfs


with open("rosalind_orf.txt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        dna = str(record.seq)


reverse_complement = str(Seq(dna).reverse_complement())
orfs = find_orfs(dna).union(find_orfs(reverse_complement))


for protein in orfs:
    print(protein)
    

    
    
    
    
    
    
    
