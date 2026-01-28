import sys
import os
import re


class PlasmidDesigner:
    def __init__(self, marker_file="markers.tab"):
        # 1. Comprehensive Internal Enzyme Database (Matches your markers.tab)
        self.enzyme_db = {
            "EcoRI": "GAATTC",
            "BamHI": "GGATCC",
            "HindIII": "AAGCTT",
            "PstI": "CTGCAG",
            "KpnI": "GGTACC",
            "SacI": "GAGCTC",
            "SalI": "GTCGAC",
            "XbaI": "TCTAGA",
            "NotI": "GCGGCCGC",
            "SmaI": "CCCGGG",
            "SphI": "GCATGC",
            "BsaI": "GGTCTC",
            "BbsI": "GAAGAC",
            "BsmBI": "CGTCTC"
        }

        # 2. Internal Library for Complex Sequences (Antibiotics/Reporters)
        self.internal_seq_library = {
            # Antibiotics & Aliases
            "AmpR": "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAA",
            "Ampicillin": "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAA",
            "KanR": "ATGGCTAAAATGAGAATATCACCGGAATTGAAAAAACTGATCGAAAAA",
            "Kanamycin": "ATGGCTAAAATGAGAATATCACCGGAATTGAAAAAACTGATCGAAAAA",

            # Functional Modules
            "lacZ_alpha": "GGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCACACAGGAAACAGCTATGACCATGATTAC",
            "Blue_White_Selection": "GGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCACACAGGAAACAGCTATGACCATGATTAC",
        }

        self.marker_db = {}

        # Load external definitions (overwrites defaults if found)
        self.load_markers_from_table(marker_file)

    def load_markers_from_table(self, file_path):
        """Parses the user's Markdown table to extract sequences."""
        print(f"[INFO] Loading markers from '{file_path}'...")
        if not os.path.exists(file_path):
            print("[WARN] markers.tab not found. Using internal defaults.")
            return

        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("|-") or "Category" in line: continue

                parts = [p.strip() for p in line.split('|')]
                if len(parts) < 4: continue

                # Table Columns: | Category | Name | Recognition | Usage |
                category = parts[1]
                name_short = parts[2]
                recognition_text = parts[3]

                # Extract Restriction Enzymes
                if "Restriction enzyme" in category:
                    match = re.search(r'[ACGT]{4,}', recognition_text)
                    if match:
                        self.enzyme_db[name_short] = match.group(0)

                # Extract Markers (Linking names to internal library)
                elif "Selection marker" in category or "Screening" in category:
                    if name_short in self.internal_seq_library:
                        self.marker_db[name_short] = self.internal_seq_library[name_short]

    def read_fasta(self, file_path):
        if not os.path.exists(file_path): return ""
        seq = []
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith(">"): seq.append(line.strip().upper())
        return "".join(seq)

    def find_ori_via_gc_skew(self, sequence):
        """Finds ORI using Cumulative GC Skew Minima."""
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

        ori_length = 600
        if min_index + ori_length < len(sequence):
            return sequence[min_index: min_index + ori_length]
        else:
            return sequence[min_index:] + sequence[:(min_index + ori_length) - len(sequence)]

    def parse_design(self, file_path):
        """
        Parses Design.txt. Format: Custom_Name, Lookup_Key
        """
        features = []
        if not os.path.exists(file_path): return features

        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("**"): continue

                parts = [p.strip() for p in line.split(',')]
                if len(parts) >= 2:
                    features.append((parts[0], parts[1]))
        return features

    def construct_plasmid(self, ori_seq, design_features):
        """Assembles plasmid based on the ORDER in Design file."""
        plasmid = ""
        ori_added = False

        print("\n--- Assembling Plasmid Features ---")

        for name, key in design_features:
            # 1. Handle ORI Request
            if key == "High_Copy_Replication" or "ori" in key.lower():
                print(f"  + Adding Algorithmically Found ORI ({key})")
                plasmid += ori_seq
                ori_added = True

            # 2. Handle Enzymes
            elif key in self.enzyme_db:
                print(f"  + Adding Enzyme: {key}")
                plasmid += self.enzyme_db[key]

            # 3. Handle Markers (External DB)
            elif key in self.marker_db:
                print(f"  + Adding Marker: {key}")
                plasmid += self.marker_db[key]

            # 4. Handle Markers (Internal Fallback)
            elif key in self.internal_seq_library:
                print(f"  + Adding Feature (Internal Lib): {key}")
                plasmid += self.internal_seq_library[key]

            else:
                print(f"  [WARN] Feature '{key}' not found. Inserting Placeholder.")
                plasmid += "NNNNNN"

        if not ori_added:
            print("  [NOTE] ORI not requested in Design file. Prepending default.")
            plasmid = ori_seq + plasmid

        return plasmid

    def save_output(self, sequence, filename="Output.fa"):
        with open(filename, "w") as f:
            f.write(">Synthesized_Plasmid_Vector\n")
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i + 80] + "\n")
        print(f"\n[SUCCESS] Plasmid saved to '{filename}'.")


def main():
    tool = PlasmidDesigner(marker_file="markers.tab")

    full_genome = tool.read_fasta("Input.fa")
    if not full_genome:
        print("Input.fa missing. Run generator first.")
        return

    ori_seq = tool.find_ori_via_gc_skew(full_genome)
    design_features = tool.parse_design("Design.txt")
    final_plasmid = tool.construct_plasmid(ori_seq, design_features)

    tool.save_output(final_plasmid, "Output.fa")

    # Verification: EcoRI (GAATTC) should be ABSENT
    if "GAATTC" not in final_plasmid:
        print("[CHECK] EcoRI (GAATTC) is ABSENT (Correctly Deleted).")
    else:
        print("[CHECK] EcoRI (GAATTC) is PRESENT.")


if __name__ == "__main__":
    main()

