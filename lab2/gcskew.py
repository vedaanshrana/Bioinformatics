import sys
import gzip
import os
import matplotlib.pyplot as plt

# Check if itertools is available (standard library), if not we do manual sum
import itertools


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


def calculate_gc_skew(sequence, window_size, step):
    """
    Calculates GC skew (G-C)/(G+C) for sliding windows.
    Returns the raw skew values and the genomic positions.
    """
    print(f"Calculating GC Skew (Window={window_size}, Step={step})...")

    skew_values = []
    x_positions = []
    seq_len = len(sequence)

    for i in range(0, seq_len - window_size + 1, step):
        window = sequence[i: i + window_size]

        g_count = window.count('G')
        c_count = window.count('C')

        # Avoid division by zero if a window has no G or C (e.g., all Ns or As/Ts)
        if g_count + c_count == 0:
            skew = 0
        else:
            skew = (g_count - c_count) / (g_count + c_count)

        skew_values.append(skew)
        # Store the center position of the window
        x_positions.append(i + window_size // 2)

    return x_positions, skew_values


def plot_cumulative_skew(x, raw_skew_values, filename):
    """
    Calculates the cumulative sum of the skew and plots it.
    The minimum of this plot often indicates the Origin of Replication (Ori).
    """
    # Calculate Cumulative Skew: Running total of the skew values
    cumulative_skew = list(itertools.accumulate(raw_skew_values))

    # Identify the Minimum point (Predicted Ori)
    min_skew = min(cumulative_skew)
    min_index = cumulative_skew.index(min_skew)
    ori_position = x[min_index]

    # Plotting
    plt.figure(figsize=(12, 6))

    # Plot the cumulative line
    plt.plot(x, cumulative_skew, color='green', label='Cumulative GC Skew')

    # Mark the minimum point
    plt.axvline(x=ori_position, color='red', linestyle='--', alpha=0.7, label=f'Predicted Ori (~{ori_position} bp)')

    plt.title('Cumulative GC Skew Analysis', fontsize=14)
    plt.xlabel('Genomic Position (bp)', fontsize=12)
    plt.ylabel('Cumulative Skew', fontsize=12)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)

    # Save the plot
    output_img = filename.replace('.fa.gz', '_gc_skew_plot.png')
    plt.savefig(output_img)
    print(f"Plot saved as: {output_img}")
    print(f"Predicted Origin of Replication (Minimum Skew) at approx: {ori_position} bp")
    plt.show()


def main():
    # --- Configuration ---
    WINDOW_SIZE = 5000
    STEP = 500

    # --- File Handling ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    fasta_file = "genomic.fa.gz"
    fasta_path = os.path.join(script_dir, fasta_file)

    # --- Execution ---

    # 1. Load Genome
    genome_seq = parse_fasta_gz(fasta_path)

    # 2. Calculate Windowed Skew
    x_pos, raw_skews = calculate_gc_skew(genome_seq, WINDOW_SIZE, STEP)

    # 3. Calculate Cumulative Skew and Plot
    plot_cumulative_skew(x_pos, raw_skews, fasta_file)


if __name__ == "__main__":
    main()