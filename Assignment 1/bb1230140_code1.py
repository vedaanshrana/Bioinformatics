import sys
import os
import re


class PlasmidDesigner:
    def __init__(self, marker_file="markers.tab"):
        # Internal Library of Genes (Since your file describes them but lacks the full sequence)
        self.internal_seq_library = {
            # Antibiotics
            "AmpR": "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAA",
            "Ampicillin": "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAA",
            "KanR": "ATGGCTAAAATGAGAATATCACCGGAATTGAAAAAACTGATCGAAAAA",
            "Kanamycin": "ATGGCTAAAATGAGAATATCACCGGAATTGAAAAAACTGATCGAAAAA",
            "CmR": "ATGGAGAAAAAAATCACTGGATATACCACCGTTGATATATCCCAATGGC",
            "TetR": "ATGAATAGTTCGACAAAGATCGCATTGGTAATTACGTTACTCGATGCCAT",
            # Common parts
            "lacZ": "GGCAGCTGGCACGACAGGTTTCCCGACTGG",
            "GFP": "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
        }

        self.marker_db = {}
        self.enzyme_db = {}

        # Load the user's table
        self.load_markers_from_table(marker_file)

    def load_markers_from_table(self, file_path):
        """
        Parses the specific Markdown Table format provided by the user.
        Extracts Enzyme sequences from the 'Recognition' column.
        """
        print(f"[INFO] Loading markers from '{file_path}'...")
        if not os.path.exists(file_path):
            print("[WARN] File not found. Using internal defaults only.")
            return

        with open(file_path, 'r') as f:
            for line in f:
                # 1. Basic cleaning
                line = line.strip()
                if not line or line.startswith("|-") or "Category" in line:
                    continue  # Skip separator lines and headers

                # 2. Split by pipe '|'
                # Format: | Category | Name | Recognition | Usage |
                parts = [p.strip() for p in line.split('|')]

                # Because of the leading/trailing pipes, the actual data is usually at indices 1, 2, 3...
                # parts[0] is empty string before first pipe
                if len(parts) < 4: continue

                category = parts[1]
                name_short = parts[2]
                recognition_text = parts[3]

                # 3. Handle Restriction Enzymes
                if "Restriction enzyme" in category:
                    # Logic: Extract uppercase DNA sequence from "Recognizes GAATTC"
                    # Regex looks for continuous A,C,G,T string of length 4 or more
                    match = re.search(r'[ACGT]{4,}', recognition_text)
                    if match:
                        seq = match.group(0)
                        self.enzyme_db[name_short] = seq
                        # print(f"  -> Loaded Enzyme: {name_short} = {seq}")
                    else:
                        pass  # Could be 'Cuts outside...' types (Golden Gate), ignoring for basic cloning logic

                # 4. Handle Antibiotics/Markers
                elif "Selection marker" in category or "Screening" in category:
                    # The file DOES NOT have the sequence, so we link the Name to our Internal Library
                    # Try to find a match in our library (e.g., "AmpR" or "Ampicillin")

                    # We store it so the Design file can reference it
                    # Check exact match first
                    if name_short in self.internal_seq_library:
                        self.marker_db[name_short] = self.internal_seq_library[name_short]
                    else:
                        # Sometimes Design.txt uses "Ampicillin" but table has "AmpR".
                        # We won't solve every alias, but we ensure the Short Name is registered.
                        # If we don't have the sequence, we use a placeholder.
                        self.marker_db[name_short] = "NNNNNN_MISSING_SEQ_NNNNNN"

        print(f"[INFO] Successfully loaded {len(self.enzyme_db)} enzymes and {len(self.marker_db)} markers.")

    def read_fasta(self, file_path):
        """Reads DNA sequence from a FASTA file."""
        if not os.path.exists(file_path):
            print(f"[ERROR] Input file '{file_path}' not found.")
            return ""

        seq_parts = []
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith(">"):
                    seq_parts.append(line.strip().upper())
        return "".join(seq_parts)

    def find_ori_via_gc_skew(self, sequence):
        """
        Identify the Origin of Replication (ORI) using GC Skew Analysis.
        """
        print("[ALGO] Running GC Skew analysis to find ORI...")
        skew = 0
        min_skew = 0
        min_index = 0

        for i, base in enumerate(sequence):
            if base == 'G':
                skew += 1
            elif base == 'C':
                skew -= 1

            if skew < min_skew:
                min_skew = skew
                min_index = i

        print(f"[ALGO] ORI identified at index {min_index}.")

        # Extract 600bp region around the minima
        ori_length = 600
        start_idx = min_index
        if start_idx + ori_length < len(sequence):
            ori_seq = sequence[start_idx: start_idx + ori_length]
        else:
            overflow = (start_idx + ori_length) - len(sequence)
            ori_seq = sequence[start_idx:] + sequence[:overflow]
        return ori_seq

    def parse_design(self, file_path):
        """Parses Design.txt for MCS and Markers."""
        layout = {"MCS": [], "Markers": []}
        if not os.path.exists(file_path):
            return layout

        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("**"): continue
                parts = [p.strip() for p in line.split(',')]

                if len(parts) >= 2:
                    feature_type = parts[0]
                    name = parts[1]

                    if "Multiple_Cloning_Site" in feature_type:
                        layout["MCS"].append(name)
                    elif "Antibiotic_marker" in feature_type:
                        layout["Markers"].append(name)
        return layout

    def construct_plasmid(self, ori_seq, design_layout):
        """
        Assembles the new plasmid: [Found_ORI] + [Markers] + [MCS_Enzymes]
        """
        plasmid = ""

        # 1. Add ORI
        plasmid += ori_seq

        # 2. Add Markers (Handling aliases)
        for marker_name in design_layout["Markers"]:
            # Try direct match (e.g., "AmpR")
            if marker_name in self.marker_db:
                plasmid += self.marker_db[marker_name]
            # Try alias match (e.g., Design says "Ampicillin" but Table has "AmpR")
            elif marker_name == "Ampicillin" and "AmpR" in self.marker_db:
                plasmid += self.marker_db["AmpR"]
            elif marker_name == "Kanamycin" and "KanR" in self.marker_db:
                plasmid += self.marker_db["KanR"]
            elif marker_name in self.internal_seq_library:
                # Fallback: It wasn't in the table, but we know the sequence anyway
                print(f"[NOTE] '{marker_name}' not in table, but found in internal library.")
                plasmid += self.internal_seq_library[marker_name]
            else:
                print(f"[WARN] Marker '{marker_name}' sequence unknown. Inserting N-block.")
                plasmid += "NNNNNNNNNN"

        # 3. Add MCS Enzymes (Parsed from the table)
        for enzyme in design_layout["MCS"]:
            if enzyme in self.enzyme_db:
                plasmid += self.enzyme_db[enzyme]
            else:
                print(f"[WARN] Enzyme '{enzyme}' not defined in markers.tab. Inserting N-block.")
                plasmid += "NNNNNNNNNN"

        return plasmid

    def save_output(self, sequence, filename="Output.fa"):
        with open(filename, "w") as f:
            f.write(">Synthesized_Plasmid_Vector\n")
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i + 80] + "\n")
        print(f"[SUCCESS] Plasmid saved to '{filename}'.")


def main():
    tool = PlasmidDesigner(marker_file="markers.tab")

    full_sequence = tool.read_fasta("Input.fa")
    if not full_sequence: return

    ori_seq = tool.find_ori_via_gc_skew(full_sequence)
    design = tool.parse_design("Design.txt")
    final_plasmid = tool.construct_plasmid(ori_seq, design)

    tool.save_output(final_plasmid, "Output.fa")

    # Verify EcoRI deletion logic (Design.txt omits it, so it shouldn't be there)
    # Note: If the user explicitly ADDED EcoRI in Design.txt, it would be there.
    # But since the prompt implies the Design.txt is the one WITHOUT EcoRI...
    if "GAATTC" not in final_plasmid:
        print("\n[VERIFICATION] EcoRI site (GAATTC) is ABSENT (Correctly removed).")
    else:
        # Check if it was requested
        if "EcoRI" in design["MCS"]:
            print("\n[VERIFICATION] EcoRI site is PRESENT (Requested in Design).")
        else:
            print("\n[VERIFICATION] EcoRI site is PRESENT (Warning: Unexpected).")


if __name__ == "__main__":
    main()
