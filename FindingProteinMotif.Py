# Problem
# To allow for the presence of its varying forms, a protein motif is represented by a shorthand as follows: [XY] means "either X or Y" and {X} means "any amino acid except X." For example, the N-glycosylation motif is written as N{P}[ST]{P}.

# You can see the complete description and features of a particular protein by its access ID "uniprot_id" in the UniProt database, by inserting the ID number into

# http://www.uniprot.org/uniprot/uniprot_id
# Alternatively, you can obtain a protein sequence in FASTA format by following

# http://www.uniprot.org/uniprot/uniprot_id.fasta
# For example, the data for protein B5ZC00 can be found at http://www.uniprot.org/uniprot/B5ZC00.

# Given: At most 15 UniProt Protein Database access IDs.
# Return: For each protein possessing the N-glycosylation motif, output its given access ID followed by a list of locations in the protein string where the motif can be found.
# N{P}[ST]{P}


import re
import requests


# motif_pattern = r'N[^P][ST][^P]'

# def fetch_protein_sequence(uniprot_id):
#     """
#     Fetch the protein sequence from UniProt in FASTA format and return the sequence as a string.
#     """
#     url = f"http://www.uniprot.org/uniprot/{uniprot_id}.fasta"
#     response = requests.get(url)
#     fasta_data = response.text
#     # Remove the header
#     sequence = ''.join(fasta_data.split('\n')[1:]).strip()
#     return sequence

# def find_motif_locations(sequence, pattern):
#     """
#     Find all positions (1-based index) of the N-glycosylation motif in the protein sequence.
#     """
#     return [m.start() + 1 for m in re.finditer(pattern, sequence)]

# def find_proteins_with_motif(uniprot_ids):
#     """
#     For each UniProt ID, find the locations of the N-glycosylation motif and return the results.
#     """
#     results = {}
#     for uniprot_id in uniprot_ids:
#         sequence = fetch_protein_sequence(uniprot_id)
#         motif_positions = find_motif_locations(sequence, motif_pattern)
#         if motif_positions:  # Only include proteins with at least one motif match
#             results[uniprot_id] = motif_positions
#     return results
motif_pattern = re.compile(r"N[^P][ST][^P]")

def fetch_protein_sequence(uniprot_id):
    """
    Fetch the protein sequence from UniProt using the given UniProt ID.
    """
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        # Parse the FASTA file to extract the sequence
        fasta_lines = response.text.split("\n")
        sequence = "".join(fasta_lines[1:])  # Skip the header line
        return sequence
    else:
        print(f"Failed to fetch sequence for {uniprot_id}")
        return None

def find_motif_locations(sequence, motif_pattern):
    """
    Find all locations of the motif in the protein sequence.
    """
    locations = []
    for match in motif_pattern.finditer(sequence):
        # Positions are 1-based in the problem description
        locations.append(match.start() + 1)
    return locations

def main(uniprot_ids):
    """
    Process each UniProt ID, fetch its sequence, and find motif locations.
    """
    for uniprot_id in uniprot_ids:
        sequence = fetch_protein_sequence(uniprot_id)
        if sequence:
            locations = find_motif_locations(sequence, motif_pattern)
            if locations:
                print(uniprot_id)
                print(" ".join(map(str, locations)))

with open("rosalind_mprt.txt", "r") as file:
    uniprot_ids = file.read().splitlines()  # Read all lines into a list, each line as one UniProt ID


results = find_proteins_with_motif(uniprot_ids)


with open("output.txt", "w") as output_file:
    for uniprot_id, positions in results.items():
        output_file.write(f"{uniprot_id}\n")
        output_file.write(' '.join(map(str, positions)) + '\n')