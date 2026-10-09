# Lysine-residue-contact-analysis

## Overview

Lysine-residue-contact-analysis is a bioinformatic tool designed to automate the identification of lysine (K) residues located within a user-defined distance cutoff of a target chain, using protein–nucleic acid complex models predicted by AlphaFold3. Developed alongside our recent research paper, this tool streamlines the analysis process by leveraging Biopython's KD-tree neighbour search to compute minimum atom-to-atom distances directly from your structure files.

## Key Features

- **Automated Contact Detection:** Reports every lysine residue whose atoms come within the cutoff of the target chain, sorted by distance.
- **AlphaFold3-Friendly Input:** Accepts standard `.cif` structure files, and `.pdb` files as well.
- **KD-Tree Distance Search:** Uses `Bio.PDB.NeighborSearch` for fast minimum-distance computation instead of scanning every atom pair.
- **Configurable Query:** Target chain (`--chain`), residue identity as one-letter or three-letter code (`--resname K` or `--resname LYS`), and cutoff (`--cutoff`) are all set from the command line.
- **Detailed Output:** Each hit reports the minimum distance plus the closest atom pair (for example `LYS784 NZ vs DC13 OP1`), which identifies the chemical contact that forms it.
- **Multi-Model Aware:** Scans only the first model by default, or every model with `--models all`.
- **CSV Export:** Writes the same contacts to a spreadsheet-ready CSV with `--csv`.
- **Scripted Workflow:** Runs entirely through Python, requiring minimal user intervention.

## Technologies Used

- Python
- Biopython
- NumPy

## Installation

Install Python 3.8 or newer, then:

Clone this repository:

```
git clone https://github.com/ihuskyyyy/Lysine-residue-contact-analysis.git
cd Lysine-residue-contact-analysis
```

(Optional) Create a virtual environment:

```
python -m venv venv
source venv/bin/activate  # On Windows use: venv\Scripts\activate
```

Install Python dependencies:

```
pip install -r requirements.txt
```

## Usage

Prepare your structure file:

Ensure your predicted complex model from AlphaFold3 is saved as a `.cif` file (`.pdb` input works too).

Run the tool:

```
python main.py --input your_model.cif --chain C --resname K --cutoff 4.0
```

Replace `your_model.cif` with your actual structure file. `--chain`, `--resname` and `--cutoff` are optional and default to chain `C`, residue `LYS` and `4.0` Å.

Results:

The script prints the chain composition of the model, then one line per residue within the cutoff, sorted by minimum distance:

```
LYS residues within <= 4 A of chain C: 4
A:LYS784  min=2.16 A  (NZ vs C:DC13:OP1)
A:LYS768  min=2.97 A  (C vs C:DT10:OP1)
A:LYS361  min=3.29 A  (NZ vs C:DG6:OP1)
A:LYS673  min=3.51 A  (NZ vs C:DA1:OP3)
```

Distances are minimum atom-to-atom distances between the residue and the target chain; hydrogen atoms are ignored.

| Option | Description |
| --- | --- |
| `--input`, `-i` | Input structure file (`.cif`/`.mmcif` or `.pdb`); a positional path works too |
| `--chain` | Target chain id (default: `C`) |
| `--resname` | Residue to search, one-letter or three-letter (default: `LYS`) |
| `--cutoff` | Distance cutoff in Å (default: `4.0`) |
| `--top` | Print only the N closest residues (default: `0`, meaning all) |
| `--models` | `first` (default) or `all` |
| `--csv` | Also write the contacts to this CSV file |

## Example

```
python main.py --input fold_2026_08_05_18_23_wt_brcas12b_pam_6_14nt_model_0.pdb --chain C --resname K --cutoff 4.0 --csv contacts.csv
```

```
=== overview ===
chain A: 8955 atoms, 1090 residues, LYS=7,9,29,44,56,74,...
chain B: 456 atoms, 22 residues, LYS=none
chain C: 287 atoms, 14 residues, LYS=none
chain D: 3190 atoms, 150 residues, LYS=none

LYS residues within <= 4 A of chain C: 4
A:LYS784  min=2.16 A  (NZ vs C:DC13:OP1)
A:LYS768  min=2.97 A  (C vs C:DT10:OP1)
A:LYS361  min=3.29 A  (NZ vs C:DG6:OP1)
A:LYS673  min=3.51 A  (NZ vs C:DA1:OP3)

CSV written to contacts.csv (4 contact(s))
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for bug fixes or feature suggestions.

## License

This project is licensed under the MIT License.

## Contact

For questions or collaboration opportunities, please contact XxinGuan.
