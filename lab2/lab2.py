import sys
import gzip
import os
import collections
import matplotlib.pyplot as plt


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
        with gzip.open(file_path, 'rt') as f:  # 'rt' mode for reading text
            for line in f:
                if not line.startswith(">"):
                    sequence.append(line.strip().upper())
    except Exception as e:
        print(f"Error reading compressed file: {e}")
        sys.exit(1)

    return "".join(sequence)


def find_top_kmer(sequence, k):
    """
    Counts all k-mers in the entire sequence to find the most frequent one (Global Max).
    This serves as our 'overrepresented' target to check for ORI signals.
    """
    print(f"Calculating global frequencies for k={k}...")
    counts = collections.Counter()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        counts[kmer] += 1

    # Get the most common k-mer
    top_kmer, count = counts.most_common(1)[0]
    print(f"Most overrepresented k-mer globally: {top_kmer} (Count: {count})")
    return top_kmer


def sliding_window_analysis(sequence, k, window_size, step, target_kmer):
    """
    Slides across the genome and counts the target k-mer in each window.
    """
    print(f"Running sliding window analysis (Window={window_size}, Step={step})...")
    x_positions = []
    y_scores = []

    seq_len = len(sequence)

    for i in range(0, seq_len - window_size + 1, step):
        window_seq = sequence[i: i + window_size]

        # Count occurrences of the target k-mer in this window
        # Note: We use string count which is non-overlapping for simplicity in standard labs,
        # but for strict k-mer overlapping, manual iteration is preferred.
        # Here we use standard count for speed.
        count = window_seq.count(target_kmer)

        # X position is the center of the window
        x_positions.append(i + window_size // 2)
        y_scores.append(count)

    return x_positions, y_scores


def plot_enrichment(x, y, kmer, filename):
    """
    Visualizes the enrichment pattern.
    """
    plt.figure(figsize=(12, 6))
    plt.plot(x, y, color='blue', linewidth=1, label=f'{kmer} count')
    plt.title(f'ORI Signal Checker: Enrichment of {kmer} (Window=5000bp)', fontsize=14)
    plt.xlabel('Genomic Position (bp)', fontsize=12)
    plt.ylabel('K-mer Count per Window', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()

    # Save the plot
    output_img = filename.replace('.fa.gz', '_plot.png')
    plt.savefig(output_img)
    print(f"Plot saved as: {output_img}")
    plt.show()


def main():
    # --- Configuration ---
    K = 8
    WINDOW_SIZE = 5000
    STEP = 500

    # --- File Handling using sys and os ---
    # Determine the script's directory to find files relatively
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Define file names as per requirements
    fasta_file = "genomic.fa.gz"
    gff_file = "genomic.gff.gz"  # Not used for plotting, but requested in context

    # Construct absolute paths
    fasta_path = os.path.join(script_dir, fasta_file)

    # --- Execution Flow ---

    # 1. Load Genome
    genome_seq = parse_fasta_gz(fasta_path)

    # 2. Identify Overrepresented K-mer (The "Signal")
    # In a real ORI scenario, you might look for DnaA boxes (e.g., TTATCCACA).
    # Here, we auto-detect the most frequent k-mer to plot.
    target_kmer = find_top_kmer(genome_seq, K)

    # 3. Sliding Window Counting
    x, y = sliding_window_analysis(genome_seq, K, WINDOW_SIZE, STEP, target_kmer)

    # 4. Visualization
    plot_enrichment(x, y, target_kmer, fasta_file)


if __name__ == "__main__":
    main()