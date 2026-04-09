**Mini VEP** is a Python tool that predicts the effects of genetic variants (point mutations, insertions, and deletions) on coding DNA sequences. It translates DNA to protein and visually highlights amino acid changes.

---

## 🔹 Features
- **Point Mutation Analysis:** Silent, Missense, or Nonsense classification.
- **Insertion & Deletion Handling:** Detects Frameshift and In-frame effects.
- **Protein Comparison:** Highlights amino acid differences between original and mutated proteins.
- **User-Friendly:** Interactive command-line menu.

---

## 🔹 Installation

1. Make sure you have **Python 3** installed.
2. Install **Biopython** library:

```bash
pip install biopython
```
3. Clone the repository:
```bash
git clone https://github.com/Entedaor/Mini-VEP.git
cd Mini_VEP
```

## 🔹 Usage
1. Prepare a DNA sequence in FASTA format (e.g., sample.fasta).
2. Run the program:
```bash
python mini_vep.py
```
3. Follow the interactive menu:
```bash
MINI VEP - MAIN MENU
1. Substitution (Point Mutation)
2. Insertion
3. Deletion
4. Exit
```
4. Enter mutation details and see predicted effects along with protein changes.

## 🔹 Example Output

Original Protein:
M A R K T H E P R O T E I N

Mutated Protein:
M A [5:T→I] K T H E P R O T E I N

[5:T→I] shows the amino acid change at position 5 (Threonine → Isoleucine).

## 🔹 Example Data

The repository includes two example CDS sequences to test Mini VEP:

- `example_data/E_coli_lacZ_example.fna` — LacZ gene from *Escherichia coli*
- `example_data/H_sapiens_HBB_example.fna` — Beta-globin gene from *Homo sapiens*

You can use these files with Mini VEP to try point mutations, insertions, and deletions.
Frameshift or nonsense mutations are also highlighted automatically.
