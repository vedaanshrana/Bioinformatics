import sys
import os


class PlasmidDesigner:
    def __init__(self):
        # Database of standard sequences (Mock sequences for demonstration)
        # In a real scenario, these would be fetched from NCBI/GenBank.
        self.sequence_db = {
            # --- DEFAULT REPLICATION GENES (Based on Jain & Srivastava, 2013) ---
            # Essential for BHR plasmid replication
            "oriV": "TTATCCGCTCACAATTCCACACA",  # Vegetative Origin
            "repA": "ATGAGCCCGAAAGCCAGCCCG",  # Helicase
            "repB": "ATGGCGACCGAGTTGCTCTTG",  # Primase
            "repC": "ATGTCCGAGAGCCGGATCGTG",  # Initiator

            # --- ANTIBIOTIC MARKERS ---
            "Ampicillin": "TTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCAT",
            "Kanamycin": "ATGGCTAAAATGAGAATATCACCGGAATTGAAAAAACTGATCGAAAAA",
            "Tetracycline": "ATGAATAGTTCGACAAAGATCGCATTGGTAATTACGTTACTCGATGCCAT",
            "Chloramphenicol": "ATGGAGAAAAAAATCACTGGATATACCACCGTTGATATATCCCAATGGC",

            # --- RESTRICTION ENZYME SITES (Palindromic sequences) ---
            "EcoRI": "GAATTC",
            "BamHI": "GGATCC",
            "HindIII": "AAGCTT",
            "NotI": "GCGGCCGC",
            "XhoI": "CTCGAG",
            "PstI": "CTGCAG"
        }

    def read_fasta(self, file_path):
        """Reads the Input DNA (Insert) from a FASTA file."""
        if not os.path.exists(file_path):
            print(f"Error: {file_path} not found.")
            return ""

        sequence = ""
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith(">"):
                    sequence += line.strip().upper()
        return sequence

    def parse_design(self, file_path):
        """Parses the Design.txt file to extract MCS and Markers."""
        layout = {
            "MCS": [],
            "Markers": []
        }

        if not os.path.exists(file_path):
            print(f"Error: {file_path} not found.")
            return layout

        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("**"): continue

                parts = [p.strip() for p in line.split(',')]

                if len(parts) == 2:
                    feature_type, name = parts

                    if "Multiple_Cloning_Site" in feature_type:
                        # Pair: (MCS_ID, Enzyme_Name)
                        layout["MCS"].append(name)
                    elif "Antibiotic_marker" in feature_type:
                        # Pair: (Marker_ID, Antibiotic_Name)
                        layout["Markers"].append(name)

        return layout

    def get_gene_sequence(self, gene_name):
        """Retrieves sequence from DB or returns a placeholder N-string."""
        # Simple fuzzy matching for the demo
        for key, seq in self.sequence_db.items():
            if key.lower() in gene_name.lower():
                return seq
        return "NNNNNNNNNN"  # Placeholder if unknown

    def construct_plasmid(self, insert_dna, design_layout):
        """
        Assembles the plasmid:
        [Default Rep Module] + [Antibiotic Markers] + [MCS/Restriction 5'] + [INSERT] + [MCS/Restriction 3']
        """

        print("--- Constructing Plasmid ---")

        # 1. Add Default Replication Genes (Based on Paper)
        # "Each plasmid will have certain genes necessary for its replication"
        # Using RepA, RepB, RepC logic from IncQ plasmids
        plasmid_seq = ""
        print("Adding Default Replication Module (oriV, repA, repB, repC)...")
        plasmid_seq += self.sequence_db["oriV"]
        plasmid_seq += self.sequence_db["repA"]
        plasmid_seq += self.sequence_db["repB"]
        plasmid_seq += self.sequence_db["repC"]

        # 2. Add Antibiotic Markers
        for marker in design_layout["Markers"]:
            print(f"Adding Marker: {marker}...")
            plasmid_seq += self.get_gene_sequence(marker)

        # 3. Add Cloning Sites and Insert
        # Logic: We treat the first half of MCS list as 5' flank and rest as 3' flank
        # or simple flanking. Let's place the insert BETWEEN the enzymes.

        mcs_enzymes = design_layout["MCS"]

        # Add 5' Restriction Site (First Enzyme)
        if len(mcs_enzymes) > 0:
            enzyme_5 = mcs_enzymes[0]
            print(f"Adding 5' Restriction Site: {enzyme_5}...")
            plasmid_seq += self.get_gene_sequence(enzyme_5)

        # Add The Input DNA (The Insert)
        print("Ligating Input DNA Insert...")
        plasmid_seq += insert_dna

        # Add 3' Restriction Site (Second Enzyme if exists, else same as first for non-directional)
        if len(mcs_enzymes) > 1:
            enzyme_3 = mcs_enzymes[1]
            print(f"Adding 3' Restriction Site: {enzyme_3}...")
            plasmid_seq += self.get_gene_sequence(enzyme_3)
        elif len(mcs_enzymes) == 1:
            # If only one enzyme, it flanks both sides
            enzyme_3 = mcs_enzymes[0]
            print(f"Adding 3' Restriction Site: {enzyme_3}...")
            plasmid_seq += self.get_gene_sequence(enzyme_3)

        return plasmid_seq

    def save_to_file(self, sequence, filename="Output.Fa"):
        """Saves the final sequence in FASTA format."""
        with open(filename, "w") as f:
            f.write(">Constructed_Plasmid_Vector_BHR\n")
            # Write in lines of 80 chars
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i + 80] + "\n")
        print(f"--- Plasmid saved to {filename} ---")


# --- Execution ---
def main():
    # File names
    input_file = "Input.Fa"
    design_file = "Design.txt"
    output_file = "Output.Fa"

    # Check if files exist (For the sake of this code running immediately,
    # I will generate dummy files if they don't exist)
    if not os.path.exists(input_file):
        with open(input_file, "w") as f:
            f.write(">Target_Gene\nATGCGTACGTACGTACGTACGTACGTACGTTAG")

    if not os.path.exists(design_file):
        with open(design_file, "w") as f:
            f.write("Multiple_Cloning_Site1, EcoRI\n")
            f.write("Multiple_Cloning_Site2, BamHI\n")
            f.write("Antibiotic_marker1, Ampicillin\n")

    # Run Process
    designer = PlasmidDesigner()

    # 1. Read Insert
    insert_seq = designer.read_fasta(input_file)

    # 2. Read Design
    design = designer.parse_design(design_file)

    # 3. Build
    full_plasmid = designer.construct_plasmid(insert_seq, design)

    # 4. Output
    designer.save_to_file(full_plasmid, output_file)


if __name__ == "__main__":
    main()