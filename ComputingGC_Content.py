#Problem
# The GC-content of a DNA string is given by the percentage of symbols in the string that are 'C' or 'G'. For example, the GC-content of "AGCTATAG" is 37.5%. Note that the reverse complement of any DNA string has the same GC-content.
# DNA strings must be labeled when they are consolidated into a database. A commonly used method of string labeling is called FASTA format. In this format, the string is introduced by a line that begins with '>', followed by some labeling information. Subsequent lines contain the string itself; the first line to begin with '>' indicates the label of the next string.
# In Rosalind's implementation, a string in FASTA format will be labeled by the ID "Rosalind_xxxx", where "xxxx" denotes a four-digit code between 0000 and 9999.
# Given: At most 10 DNA strings in FASTA format (of length at most 1 kbp each).
# Return: The ID of the string having the highest GC-content, followed by the GC-content of that string. Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated; please see the note on absolute error below.

with open("rosalind_gc.txt", "r") as file:
    fasta_input = file.read()

def calculate_gc_content(dna_string):
    g_count = dna_string.count('G')
    c_count = dna_string.count('C')
    gc_content = (g_count + c_count)/ len(dna_string) * 100
    return gc_content

def find_max_gc_content(fasta_input):
    lines = fasta_input.strip().split('\n') #This splits the cleaned input string into a list of lines, where each line is a separate string.
    max_gc_content = -1 
    max_gc_id = ""

    current_id = ""
    current_sequence = ""

    for line in lines:
        if line.startswith('>'):
            if current_sequence:
                gc_content = calculate_gc_content(current_sequence)
                if gc_content > max_gc_content:
                    max_gc_content = gc_content
                    max_gc_id = current_id
        
            current_id = line[1:]  
            current_sequence = ""
        else:
            current_sequence += line.strip() 


    if current_sequence:
        gc_content = calculate_gc_content(current_sequence)
        if gc_content > max_gc_content:
            max_gc_content = gc_content
            max_gc_id = current_id

    return max_gc_id, max_gc_content

max_id, max_gc_content = find_max_gc_content(fasta_input)
print(f"{max_id}\n{max_gc_content:.6f}")