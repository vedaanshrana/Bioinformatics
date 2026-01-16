import sys
import gzip
import os
import itertools
import collections
import matplotlib.pyplot as plt


def parse_fasta_gz(file_path):
    """
    Parses a gzipped FASTA file and returns the genome sequence.
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
        print(f"Error reading file: {e}")
        sys.exit(1)

    return "".join(sequence)


def calculate_cumulative_skew(sequence):
    """
    Calculates the Cumulative GC Skew across the entire genome using the
    formula: Cumulative sum of (G-C)/(G+C).
    Returns the cumulative skew array and the index of the minimum value.
    """
    print("Step 1: Calculating Cumulative GC Skew...")

    # We calculate skew step-by-step for the cumulative plot
    # To save memory/time on huge genomes, we can just do G-C (simplified skew)
    # But here we adhere to the standard (G-C)/(G+C)

    skew_values = []
    current_skew = 0

    # processing in chunks for speed, or char by char.
    # For accuracy of "Cumulative Skew Diagrams", we usually sum (G-C) across the genome.
    # Let's use the sliding window approach as requested to generate data points.

    # Note: Standard cumulative skew diagrams are often just walking the genome.
    # Given the prompt asks for "window size = 5000", we will calculate skew per window
    # and accumulate those window scores.

    return []  # Logic moved to main loop for combined processing


def sliding_window_analysis(sequence, k, window_size, step):
    """
    Performs simultaneous GC Skew calculation and K-mer counting.
    """
    print(f"Step 2: Analyzing genome in windows (Size={window_size}, Step={step})...")

    x_positions = []
    raw_skew_values = []

    seq_len = len(sequence)

    for i in range(0, seq_len - window_size + 1, step):
        window = sequence[i: i + window_size]

        # --- GC Skew Calculation ---
        g = window.count('G')
        c = window.count('C')
        if g + c == 0:
            skew = 0
        else:
            skew = (g - c) / (g + c)

        raw_skew_values.append(skew)
        x_positions.append(i + window_size // 2)

    # Calculate Cumulative Skew
    cumulative_skew = list(itertools.accumulate(raw_skew_values))

    # Find the position of the Minimum Cumulative Skew (Predicted ORI)
    min_skew_val = min(cumulative_skew)
    min_skew_index = cumulative_skew.index(min_skew_val)
    predicted_ori_pos = x_positions[min_skew_index]

    print(f"  -> GC Skew Minimum found at approx {predicted_ori_pos} bp")

    return x_positions, cumulative_skew, predicted_ori_pos


def find_ori_specific_kmer(sequence, ori_pos, k, search_radius=2000):
    """
    Instead of finding the global top k-mer, this looks at the sequence
    specifically around the Predicted ORI to find the local top k-mer (DnaA box candidate).
    """
    print(f"Step 3: Identifying consensus signal at ORI region (+/- {search_radius}bp)...")

    start = max(0, ori_pos - search_radius)
    end = min(len(sequence), ori_pos + search_radius)
    ori_region = sequence[start:end]

    # Count k-mers in this specific region
    counts = collections.Counter()
    for i in range(len(ori_region) - k + 1):
        kmer = ori_region[i: i + k]
        counts[kmer] += 1

    top_kmer, count = counts.most_common(1)[0]
    print(f"  -> Most frequent k-mer at ORI: {top_kmer} (Found {count} times in region)")
    return top_kmer


def get_kmer_distribution(sequence, target_kmer, window_size, step):
    """
    Scans the whole genome for the specific ORI k-mer to see if it peaks at the ORI.
    """
    print(f"Step 4: Mapping distribution of {target_kmer} across genome...")
    y_counts = []

    for i in range(0, len(sequence) - window_size + 1, step):
        window = sequence[i: i + window_size]
        y_counts.append(window.count(target_kmer))

    return y_counts


def plot_combined_results(x, cum_skew, kmer_counts, kmer_seq, ori_pos, filename):
    """
    Plots GC Skew and K-mer enrichment on two aligned subplots.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    # --- Plot 1: Cumulative GC Skew ---
    ax1.plot(x, cum_skew, color='green', label='Cumulative GC Skew')
    ax1.axvline(x=ori_pos, color='red', linestyle='--', alpha=0.8, label='Predicted ORI')
    ax1.set_title('Metric 1: Cumulative GC Skew (Min = ORI)', fontsize=14)
    ax1.set_ylabel('Cumulative Skew')
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend()

    # --- Plot 2: K-mer Enrichment ---
    ax2.plot(x, kmer_counts, color='blue', label=f'Freq of {kmer_seq}')
    ax2.axvline(x=ori_pos, color='red', linestyle='--', alpha=0.8)
    ax2.set_title(f'Metric 2: Enrichment of ORI Signal ({kmer_seq})', fontsize=14)
    ax2.set_xlabel('Genomic Position (bp)', fontsize=12)
    ax2.set_ylabel('Count per Window')
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend()

    plt.tight_layout()

    output_img = filename.replace('.fa.gz', '_ORI_finder.png')
    plt.savefig(output_img)
    print(f"Consolidated ORI Plot saved as: {output_img}")
    plt.show()


def main():
    # --- Configuration ---
    K = 8
    WINDOW_SIZE = 5000
    STEP = 500

    # --- File Handling ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    fasta_file = "genomic.fa.gz"
    fasta_path = os.path.join(script_dir, fasta_file)

    # --- Execution Pipeline ---

    # 1. Load Genome
    genome_seq = parse_fasta_gz(fasta_path)

    # 2. GC Skew Analysis
    # We find the "valley" in the skew plot which suggests the ORI location
    x_data, skew_data, predicted_ori = sliding_window_analysis(genome_seq, K, WINDOW_SIZE, STEP)

    # 3. Identify Candidate K-mer
    # We look closely at that predicted location to find the repeating signal (DnaA box)
    target_kmer = find_ori_specific_kmer(genome_seq, predicted_ori, K)

    # 4. K-mer Distribution
    # We count this signal across the whole genome for verification
    kmer_data = get_kmer_distribution(genome_seq, target_kmer, WINDOW_SIZE, STEP)

    # 5. Combined Visualization
    plot_combined_results(x_data, skew_data, kmer_data, target_kmer, predicted_ori, fasta_file)


if __name__ == "__main__":
    main()