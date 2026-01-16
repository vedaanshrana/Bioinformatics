import sys
import gzip
import os
from collections import defaultdict


def parse_fasta_gz(file_path):
    """
    Parses a gzipped FASTA file and returns the genome sequence as a single string.
    """
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
        sys.exit(1)

    print(f"Loading sequence from {os.path.basename(file_path)}...")
    sequence = []
    try:
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                if not line.startswith(">"):
                    sequence.append(line.strip().upper())
    except Exception as e:
        print(f"Error reading compressed file: {e}")
        sys.exit(1)

    return "".join(sequence)


def find_clumps(genome, k, L, t):
    """
    Finds k-mers that form (L, k, t)-clumps.
    A k-mer forms a clump if it appears at least 't' times within a window of length 'L'.
    """
    print(f"Scanning for clumps (k={k}, L={L}, t={t})...")

    clump_kmers = set()  # Stores unique k-mers that satisfy the condition
    window_counts = defaultdict(int)

    n = len(genome)
    if n < L:
        print("Error: Genome is shorter than the window length L.")
        return []

    # 1. Initialize the first window (from index 0 to L)
    # The k-mers in this window start at indices 0 to L-k
    for i in range(L - k + 1):
        kmer = genome[i: i + k]
        window_counts[kmer] += 1

    # Check the first window for any clumps
    for kmer, count in window_counts.items():
        if count >= t:
            clump_kmers.add(kmer)

    # 2. Slide the window across the genome by 1 position at a time
    # Range: from 1 up to the point where the window fits
    for i in range(1, n - L + 1):
        # REMOVE the k-mer that is sliding OUT (the one at i-1)
        first_kmer = genome[i - 1: i - 1 + k]
        window_counts[first_kmer] -= 1

        # ADD the k-mer that is sliding IN (the one at i + L - k)
        # Why this index? The window ends at i+L, so the last k-mer starts at (i+L) - k
        new_kmer = genome[i + L - k: i + L]
        window_counts[new_kmer] += 1

        # Check if the newly added k-mer is now a clump
        if window_counts[new_kmer] >= t:
            clump_kmers.add(new_kmer)

    return list(clump_kmers)


def main():
    # --- Configuration ---
    K = 8  # Length of k-mer
    T = 3  # Minimum occurrences
    L = 1000  # Window length

    # --- File Handling ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    fasta_file = "genomic.fa.gz"
    fasta_path = os.path.join(script_dir, fasta_file)

    # --- Execution ---

    # 1. Load Genome
    genome_seq = parse_fasta_gz(fasta_path)

    # 2. Find Clumps
    # Note: Clump finding requires checking *every* single position (step=1),
    # unlike the previous visualization which skipped by 500.
    clumps = find_clumps(genome_seq, K, L, T)

    # 3. Output Results
    print("-" * 40)
    print(f"Total distinct (L, k, t)-clumps found: {len(clumps)}")
    print("-" * 40)

    if len(clumps) > 0:
        print(f"First 10 detected clump k-mers:")
        for i, kmer in enumerate(clumps[:10]):
            print(f"{i + 1}. {kmer}")

        if len(clumps) > 10:
            print(f"... and {len(clumps) - 10} more.")
    else:
        print("No k-mers found meeting the criteria.")


if __name__ == "__main__":
    main()