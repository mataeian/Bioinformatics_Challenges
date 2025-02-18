# Problem
# After identifying the exons and introns of an RNA string, we only need to delete the introns and concatenate the exons to form a new string ready for translation.

# Given: A DNA string s  (of length at most 1 kbp) and a collection of substrings of s acting as introns. All strings are given in FASTA format.

# Return: A protein string resulting from transcribing and translating the exons of s. (Note: Only one solution will exist for the dataset provided.)

#Solution:
# Extract the DNA sequence and the intron sequences from the FASTA input.
# Remove the introns from the DNA sequence to obtain the exons.
# Transcribe the exons into RNA.
# Translate the RNA into a protein sequence.
from Bio import SeqIO
from Bio.Seq import Seq

def remove_introns(dna_sequence, introns):
    for intron in introns:
        dna_sequence = dna_sequence.replace(intron, "")
    return dna_sequence

def transcribe_and_translate(dna_sequence):
    rna_sequence = dna_sequence.transcribe()
    protein_sequence = rna_sequence.translate(to_stop=True)
    return protein_sequence

def main():

    sequences = list(SeqIO.parse("rosalind_splc.txt", "fasta"))
    
    # Extract the DNA sequence and introns
    dna_sequence = str(sequences[0].seq)
    introns = [str(record.seq) for record in sequences[1:]]
    
    # Remove introns to get exons
    exons = remove_introns(dna_sequence, introns)
    
    # Transcribe and translate the exons
    protein_sequence = transcribe_and_translate(Seq(exons))
    
    # Output the protein sequence
    print(protein_sequence)

if __name__ == "__main__":
    main()