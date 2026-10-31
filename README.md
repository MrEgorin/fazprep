# FazPrep: 3D Molecular Conformer Generator and Structure Preparation

FazPrep is an open-source chemoinformatics platform designed for automated three-dimensional (3D) molecular structure preparation, conformer ensemble searching, and force field energy minimization.

The application operates either via an interactive Streamlit graphical web interface or integrated directly as a Python library into automated high-throughput computational chemistry pipelines.

---

## Free Academic and Industrial License

FazPrep is distributed under the permissive **MIT License**.

The software is completely free and open-source for:
- **Academic Research**: unrestricted use by universities, educational institutions, research laboratories, and for scientific publications.
- **Industrial and Commercial Applications**: full permission for commercial development, biotechnology workflows, drug discovery pipelines, and proprietary integration without royalties or licensing fees.

---

## Key Features

- **Input Flexibility**:
  - Direct SMILES string input.
  - Structure file upload: SDF, MOL, PDB, and MOL2 formats.

- **Conformer Generation**:
  - Distance-geometry stochastic sampling using the **ETKDGv3** algorithm (Experimental-Torsion Knowledge Distance Geometry version 3).
  - Configurable ensemble size (5 to 200+ conformers) and random seeds for strict reproducibility.
  - Automated hydrogen addition and valence verification.

- **Force Field Energy Minimization**:
  - **MMFF94s** (Merck Molecular Force Field 94s) default parameterization for organic drug-like molecules.
  - Automated fallback to the **UFF** (Universal Force Field) if unparameterized elements or functional groups are encountered.
  - Automatic identification and extraction of the lowest-energy conformer.

- **Standardized Multi-Format Export**:
  - **PDB**: optimized for molecular docking (AutoDock Vina, Glide, LeDock) and visualization (PyMOL, ChimeraX).
  - **SDF V3000**: full chemoinformatics specification retaining precise 3D atom coordinates and stereochemistry.
  - **MOL2**: Tripos MOL2 format generated via OpenBabel (with Gasteiger partial charges) or an internal standalone serializer.

- **Organized Workspace Archiving**:
  - Automated run versioning saved into sequential directories (`results/run_1`, `results/run_2`, etc.).
  - Direct file download buttons within the web interface.

---

## Computational Methodology (for Scientific Manuscripts)

You can adapt the following concise paragraph for the **Computational Methods** or **Materials and Methods** section of your research paper:

> "Three-dimensional molecular conformers were prepared and optimized using FazPrep, implemented on top of the RDKit chemoinformatics framework (v2023.3+). Explicit hydrogen atoms were added to the molecular structures, and a stochastic conformer ensemble was generated using the Experimental-Torsion Knowledge Distance Geometry algorithm (ETKDGv3; Riniker and Landrum, 2015). For each compound, an ensemble of conformers was subjected to geometric energy minimization using the Merck Molecular Force Field 94s (MMFF94s; Halgren, 1999), with an automated fallback to the Universal Force Field (UFF; Rappé et al., 1992) for unparameterized elements. The conformer corresponding to the global energy minimum was extracted and saved in PDB, SDF (V3000 specification), and Tripos MOL2 formats for subsequent molecular docking and molecular dynamics investigations."

---

## System Requirements

- **Python**: 3.10 or newer (Conda environment recommended).
- **RDKit**: 2023.3.1 or newer.
- **Streamlit**: 1.22.0 or newer.
- **OpenBabel** *(optional, recommended)*: for enhanced SYBYL atom typing and Gasteiger charges in MOL2 files.

---

## Installation

### Option 1: Conda Environment (Recommended)

```bash
conda create -n fazprep python=3.11 -y
conda activate fazprep
conda install -c conda-forge rdkit openbabel streamlit -y
```

### Option 2: Pip

```bash
pip install -r requirements.txt
```

---

## Usage

### 1. Launch the Streamlit Web Interface

From the repository directory:

```bash
streamlit run app.py
```

Or using the Python module launcher:

```bash
python -m streamlit run app.py
```

On Windows, you can also simply double-click or run:

```cmd
run.bat
```

Open your web browser at `http://localhost:8501`.

### 2. Python API Integration

FazPrep functions can be imported into your custom Python scripts and batch scripts:

```python
from rdkit import Chem
from optimizer import optimize_flexible_molecule, get_next_run_dir

# Parse input ligand SMILES
smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
mol = Chem.MolFromSmiles(smiles)

# Prepare sequential output folder
run_dir, run_num = get_next_run_dir("results")

# Execute conformer search and energy minimization
results = optimize_flexible_molecule(
    mol=mol,
    output_dir=run_dir,
    base_filename="aspirin",
    num_confs=50,
    random_seed=42
)

print(f"Results directory: {run_dir} (Run #{run_num})")
print(f"Minimum energy: {results['min_energy']:.2f} kcal/mol ({results['force_field']})")
print(f"PDB file: {results['pdb_path']}")
print(f"SDF file: {results['sdf_path']}")
print(f"MOL2 file: {results['mol2_path']}")
```

---

## Directory Structure

```
fazprep/
├── __init__.py        # Public package interface
├── app.py             # Streamlit graphical application
├── optimizer.py       # Core conformer generation, minimization, and export routines
├── requirements.txt   # Python dependency specifications
├── LICENSE            # MIT License text
├── README.md          # Complete documentation and publication methodology
└── .gitignore         # Version control ignore definitions
```

---

## License

Copyright (c) 2026 Yehor Franchuk / FAZENA Project.  
Released under the MIT License. Permitted for free academic and commercial use.
