import sys
import os

def count_fasta_records(fasta_file):
    """
    Counts the number of FASTA records in a multi-FASTA file.
    """
    record_count = 0

    try:
        with open(fasta_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    record_count += 1

    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.")
        sys.exit(1)

    except IOError:
        print(f"Error: Unable to read file '{fasta_file}'.")
        sys.exit(1)

    if record_count == 0:
        print("Error: No FASTA records found.")
        sys.exit(1)

    return record_count


if __name__ == "__main__":

    fasta_file = "seq.mfa"   # multi-FASTA file name

    if not os.path.isfile(fasta_file):
        print(f"Error: '{fasta_file}' does not exist.")
        sys.exit(1)

    num_records = count_fasta_records(fasta_file)
    print(f"Number of FASTA records: {num_records}")
