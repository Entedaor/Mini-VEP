from Bio import SeqIO
from Bio.Seq import Seq

def fasta_reader(file_path):
    """Reads a FASTA file and returns the DNA sequence in uppercase."""
    try:
        record = next(SeqIO.parse(file_path, "fasta"))
        print(f"File Read Successfully")
        print(f"ID: {record.id}")
        return record.seq.upper()
    except FileNotFoundError:
        print("Error: The FASTA file could not be found.")
    except ValueError:
        print("Error: The file may contain multiple sequences.")

def highlight_protein_changes(original, mutant):
    print("Protein Comparison")

    min_len = min(len(original), len(mutant))

    result = []

    for i in range(min_len):
        if original[i] == mutant[i]:
            result.append(original[i])
        else:
            result.append(f"[{i+1}:{original[i]}→{mutant[i]}]")

    if len(original) != len(mutant):
        result.append("...")

    print("Result :", " ".join(result))

def point_mutation(coding_dna, original_protein):
    # 1. Get Mutation Details from User
    print("Point Mutation Input")
    pos = int(input(f"Enter mutation position (1-{len(coding_dna)}): ")) - 1
    # Validate that the mutation position is within the valid DNA sequence range
    if pos < 0 or pos >= len(coding_dna):
        print("Invalid position!")
        return

    new_nucleotide = input("Enter new nucleotide (A, T, G, C): ").upper()
    # Ensure that the input nucleotide is valid (only A, T, G, C are allowed)
    if new_nucleotide not in ["A", "T", "G", "C"]:
        print("Invalid nucleotide!")
        return
    
    # 2. Apply Point Mutation
    dna_list = list(coding_dna)
    dna_list[pos] = new_nucleotide
    mutant_dna = Seq("".join(dna_list))

    # 3. Mutant Translation
    mutant_protein = mutant_dna.translate(to_stop=True)

    # 4. Classification Logic
    print("Variant Effect Prediction")
        
    if original_protein == mutant_protein:
        print("Classification: SILENT MUTATION")
        print("Effect: No change in the amino acid sequence.")
            
    elif "*" in mutant_dna[pos//3 * 3 : pos//3 * 3 + 3].translate():
        # Note: This is a simplified check for a Stop Codon appearing
        print("Classification: NONSENSE MUTATION")
        print("Effect: Premature stop codon created.")
            
    else:
        print("Classification: MISSENSE MUTATION")
        print("Effect: Single amino acid change detected.")

    # 5. Visual Comparison (Extra feature)
    print(f"\nOriginal: {original_protein}")
    print(f"Mutated : {mutant_protein}")

    highlight_protein_changes(original_protein, mutant_protein)

def insertion_mutation(coding_dna, original_protein):
    """Handles nucleotide insertions and predicts frameshift effects."""
    print("Insertion Input")
    pos = int(input(f"Enter insertion position (1-{len(coding_dna)}): ")) - 1
    # Validate that the mutation position is within the valid DNA sequence range
    if pos < 0 or pos > len(coding_dna):
        print("Invalid position!")
        return

    inserted_bases = input("Enter nucleotides to insert (e.g., A, CC, GAT): ").upper()
    for base in inserted_bases:
        # Check that all inserted bases are valid nucleotides
        if base not in ["A", "T", "G", "C"]:
            print("Invalid nucleotide in insertion!")
            return
        
    # Apply Insertion
    # We take the DNA up to 'pos', add the new bases, and then the rest of the DNA
    mutant_dna_str = str(coding_dna[:pos]) + inserted_bases + str(coding_dna[pos:])
    mutant_dna = Seq(mutant_dna_str)

    # Mutant Translation
    # Frameshift mutations often change the protein length or result in early stop codons
    mutant_protein = mutant_dna.translate(to_stop=True)
    
    print("Variant Effect Prediction")
    
    # Logic: If the number of inserted bases is NOT a multiple of 3, it's a Frameshift
    if len(inserted_bases) % 3 != 0:
        print("Classification: FRAMESHIFT MUTATION (Insertion)")
        print(f"Effect: Reading frame shifted by {len(inserted_bases)} base(s).")
    else:
        print("Classification: IN-FRAME INSERTION")
        print(f"Effect: {len(inserted_bases)//3} amino acid(s) added without shifting the frame.")

    # Visual Comparison
    print(f"Original Protein: {original_protein}")
    print(f"Mutant Protein  : {mutant_protein}")
    
    # Extra: Show how the DNA length changed
    print(f"\nDNA Length Change: {len(coding_dna)} bp -> {len(mutant_dna)} bp")

    highlight_protein_changes(original_protein, mutant_protein)

def deletion_mutation(coding_dna, original_protein):
    """Handles nucleotide deletions and predicts frameshift or in-frame effects."""
    print("Deletion Input")
    pos = int(input(f"Enter starting position for deletion (1-{len(coding_dna)}): ")) - 1
    # Validate that the mutation position is within the valid DNA sequence range
    if pos < 0 or pos >= len(coding_dna):
        print("Invalid position!")
        return
    
    length = int(input("How many nucleotides do you want to delete?: "))
    # Ensure deletion length is valid and does not exceed sequence bounds
    if length <= 0 or pos + length > len(coding_dna):   
        print("Invalid deletion length!")
        return
    
    # Apply Deletion
    # We take the DNA up to 'pos', and skip the 'length' amount of bases
    mutant_dna_str = str(coding_dna[:pos]) + str(coding_dna[pos + length:])
    mutant_dna = Seq(mutant_dna_str)

    # Mutant Translation
    mutant_protein = mutant_dna.translate(to_stop=True)
    
    print("Variant Effect Prediction")
    
    # Logic: If the number of deleted bases is NOT a multiple of 3, it's a Frameshift
    if length % 3 != 0:
        print("Classification: FRAMESHIFT MUTATION (Deletion)")
        print(f"Effect: Reading frame shifted by deleting {length} base(s).")
    else:
        print("Classification: IN-FRAME DELETION")
        print(f"Effect: {length//3} amino acid(s) removed without shifting the frame.")

    # Visual Comparison
    print(f"\nOriginal Protein: {original_protein}")
    print(f"Mutant Protein  : {mutant_protein}")
    
    # Extra: Show how the DNA length changed
    print(f"DNA Length Change: {len(coding_dna)} bp -> {len(mutant_dna)} bp")

    highlight_protein_changes(original_protein, mutant_protein)

def main_menu():
    """Main program loop for Variant Effect Predictor."""
    filename = input("Enter the FASTA filename (e.g., gen.fasta):")
    dna_input = fasta_reader(filename)

    if not dna_input:
        print("Could not proceed without a valid sequence.")
        return

    # Find Start Codon and Set Frame
    start_index = dna_input.find("ATG")
    if start_index != -1:
        coding_dna = dna_input[start_index:]
        print(f"Sequence trimmed. Start index: {start_index}")
    else:
        print("Warning: No start codon (ATG) found. Using raw sequence.")
        coding_dna = dna_input

    original_protein = coding_dna.translate(to_stop=True)

    while True:
            print("MINI VEP - MAIN MENU")
            print("1. Substitution (Point Mutation)")
            print("2. Insertion (To be implemented)")
            print("3. Deletion (To be implemented)")
            print("4. Exit")

            choice = input("\nSelect an operation (1-4): ")

            if choice == '1':
                point_mutation(coding_dna, original_protein)
            elif choice == '2':
                insertion_mutation(coding_dna, original_protein)
            elif choice == '3':
                deletion_mutation(coding_dna, original_protein)
            elif choice == '4':
                print("\nExiting program. Goodbye!")
                break
            else:
                print("\n[ERROR] Invalid choice. Please select 1, 2, 3, or 4.")

main_menu()