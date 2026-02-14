import os

class SequenceAligner:
    def __init__(self, match, mismatch, gap_open, gap_ext):
        # Scoring parameters
        self.match = match
        self.mismatch = mismatch
        self.gap_open = gap_open
        self.gap_ext = gap_ext

    def read_fasta(self, filepath):
        """Reads a FASTA file and returns the sequence."""
        if not os.path.exists(filepath):
            print(f"Error: {filepath} not found.")
            return ""

        seq = []
        with open(filepath, 'r') as f:
            for line in f:
                if not line.startswith(">"):
                    seq.append(line.strip().upper())
        return "".join(seq)

    def align(self, seq1, seq2):
        """
        Performs Local Alignment using Smith-Waterman with Affine Gap Penalties (Gotoh).
        """
        n, m = len(seq1), len(seq2)

        # Initialize matrices with zeros
        M = [[0] * (m + 1) for _ in range(n + 1)]  # Match/Mismatch state
        X = [[0] * (m + 1) for _ in range(n + 1)]  # Gap in Seq1 state
        Y = [[0] * (m + 1) for _ in range(n + 1)]  # Gap in Seq2 state

        max_score = 0
        max_pos = (0, 0)

        # Fill the matrices
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                score = self.match if seq1[i - 1] == seq2[j - 1] else self.mismatch

                X[i][j] = max(M[i][j - 1] + self.gap_open, X[i][j - 1] + self.gap_ext)
                Y[i][j] = max(M[i - 1][j] + self.gap_open, Y[i - 1][j] + self.gap_ext)
                M[i][j] = max(0, M[i - 1][j - 1] + score, X[i][j], Y[i][j])

                if M[i][j] > max_score:
                    max_score = M[i][j]
                    max_pos = (i, j)

        # Pass all three matrices directly to the traceback
        return self._traceback(seq1, seq2, M, X, Y, max_pos, max_score)

    def _traceback(self, seq1, seq2, M, X, Y, max_pos, max_score):
        """Traces back dynamically using state transitions (Gotoh)."""
        align1, align2, match_str = [], [], []
        i, j = max_pos

        # We always start the traceback in the Match/Mismatch matrix (M)
        state = 'M'

        while i > 0 and j > 0:
            if state == 'M':
                if M[i][j] == 0:
                    break  # Local alignment ends when score hits 0

                score = self.match if seq1[i - 1] == seq2[j - 1] else self.mismatch

                # Where did this M score come from?
                if M[i][j] == M[i - 1][j - 1] + score:
                    align1.append(seq1[i - 1])
                    align2.append(seq2[j - 1])
                    match_str.append('|' if seq1[i - 1] == seq2[j - 1] else '.')
                    i -= 1
                    j -= 1
                    # State remains 'M'
                elif M[i][j] == X[i][j]:
                    state = 'X'  # Switch to Gap in Seq1 state
                elif M[i][j] == Y[i][j]:
                    state = 'Y'  # Switch to Gap in Seq2 state

            elif state == 'X':
                # Moving Left (Gap in Seq1)
                align1.append('-')
                align2.append(seq2[j - 1])
                match_str.append(' ')

                # Did this gap extend an existing one, or open a new one?
                if X[i][j] == X[i][j - 1] + self.gap_ext:
                    state = 'X'  # Continue extending
                else:
                    state = 'M'  # Gap was opened here, return to M
                j -= 1

            elif state == 'Y':
                # Moving Up (Gap in Seq2)
                align1.append(seq1[i - 1])
                align2.append('-')
                match_str.append(' ')

                if Y[i][j] == Y[i - 1][j] + self.gap_ext:
                    state = 'Y'
                else:
                    state = 'M'
                i -= 1

        # The sequences were built backwards, so reverse them
        align1 = "".join(align1[::-1])
        align2 = "".join(align2[::-1])
        match_str = "".join(match_str[::-1])

        return align1, match_str, align2, max_score


def main():
    # 1. Take File Inputs
    file1 = "seq1.fa"
    file2 = "seq2.fa"

    # 2. User Input Parameters
    # (Using standard LALIGN default magnitudes: Match +5, Mismatch -4, Open -10, Ext -1)
    print("--- Nucleotide Sequence Alignment Tool ---")
    try:
        match = float(input("Enter Match Score (e.g., 5): "))
        mismatch = float(input("Enter Mismatch Penalty (negative, e.g., -4): "))
        gap_open = float(input("Enter Gap Opening Penalty (negative, e.g., -10): "))
        gap_ext = float(input("Enter Gap Extension Penalty (negative, e.g., -1): "))
    except ValueError:
        print("Invalid input. Using default LALIGN parameters (5, -4, -10, -1).")
        match, mismatch, gap_open, gap_ext = 5, -4, -10, -1

    # Initialize aligner
    aligner = SequenceAligner(match, mismatch, gap_open, gap_ext)

    # Read sequences
    seq1 = aligner.read_fasta(file1)
    seq2 = aligner.read_fasta(file2)

    if not seq1 or not seq2:
        print("Please ensure both 'seq1.fa' and 'seq2.fa' exist in the same directory as this script.")
        return

    print("\nAligning sequences...")
    align1, match_str, align2, score = aligner.align(seq1, seq2)

    # 3. Output
    print("\n--- Best Local Alignment ---")
    print(f"Alignment Score: {score}")
    print(f"Seq1: {align1}")
    print(f"      {match_str}")
    print(f"Seq2: {align2}")


if __name__ == "__main__":
    main()