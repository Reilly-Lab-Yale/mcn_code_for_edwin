from collections import Counter
from Bio import SeqIO

def generate_consensus(fasta_file, threshold=0.1):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    seq_length = len(sequences[0].seq)
    
    # Initialize a list of Counters for each position
    position_counts = [Counter() for _ in range(seq_length)]
    
    # Count bases at each position
    for record in sequences:
        for i, base in enumerate(record.seq):
            position_counts[i][base] += 1
    
    # Generate consensus sequence
    consensus = []
    num_sequences = len(sequences)
    for count in position_counts:
        # Find the most common base and its frequency
        most_common_base, most_common_count = count.most_common(1)[0]
        if most_common_count / num_sequences > threshold:
            consensus.append(most_common_base)
        else:
            consensus.append('N')
    
    return ''.join(consensus)

# Specify the input file
fasta_file = "input.fasta"

# Generate and print the consensus sequence
consensus_sequence = generate_consensus("/home/er793/palmer_scratch/Edwin/Edwin_oligos.fasta")
print("Consensus Sequence:")
print(consensus_sequence)
