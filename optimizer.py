import os
import glob
import tempfile
import shutil
import subprocess
from typing import Optional, Dict, Any, Tuple
from rdkit import Chem
from rdkit.Chem import AllChem


def get_next_run_dir(base_results_dir: str = "results") -> Tuple[str, int]:
    """
    Creates and returns the path to the next sequential run directory (e.g., results/run_1).
    """
    if not os.path.exists(base_results_dir):
        os.makedirs(base_results_dir, exist_ok=True)

    existing_runs = glob.glob(os.path.join(base_results_dir, "run_*"))
    run_numbers = []
    for r in existing_runs:
        folder_name = os.path.basename(r)
        try:
            num = int(folder_name.split("_")[1])
            run_numbers.append(num)
        except (IndexError, ValueError):
            pass

    next_run_num = max(run_numbers) + 1 if run_numbers else 1
    run_dir = os.path.join(base_results_dir, f"run_{next_run_num}")
    os.makedirs(run_dir, exist_ok=True)
    return run_dir, next_run_num


def load_molecule_from_file(uploaded_file, filename: Optional[str] = None) -> Optional[Chem.Mol]:
    """
    Reads an RDKit Mol object from an uploaded file (SDF, MOL2, PDB, MOL).
    Supports Streamlit UploadedFile objects as well as standard file objects.
    """
    fname = filename or getattr(uploaded_file, "name", "molecule.mol")
    file_ext = os.path.splitext(fname)[1].lower()

    if hasattr(uploaded_file, "read"):
        file_bytes = uploaded_file.read()
        if hasattr(uploaded_file, "seek"):
            uploaded_file.seek(0)
    else:
        with open(uploaded_file, "rb") as f:
            file_bytes = f.read()

    mol = None

    if file_ext in [".sdf", ".mol"]:
        mol = Chem.MolFromMolBlock(file_bytes.decode("utf-8", errors="ignore"))
    elif file_ext == ".pdb":
        mol = Chem.MolFromPDBBlock(file_bytes.decode("utf-8", errors="ignore"))
    elif file_ext == ".mol2":
        mol = Chem.MolFromMol2Block(file_bytes.decode("utf-8", errors="ignore"))

    if mol is not None:
        return mol

    # Fallback via temporary file on disk
    with tempfile.NamedTemporaryFile(delete=False, suffix=file_ext) as tmp:
        tmp.write(file_bytes)
        tmp_path = tmp.name

    try:
        if file_ext == ".mol2":
            mol = Chem.MolFromMol2File(tmp_path)
        elif file_ext in [".sdf", ".mol"]:
            mol = Chem.MolFromMolFile(tmp_path)
        elif file_ext == ".pdb":
            mol = Chem.MolFromPDBFile(tmp_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    return mol


def save_mol_to_mol2(mol: Chem.Mol, mol2_path: str, pdb_source_path: Optional[str] = None) -> bool:
    """
    Saves a molecule in Tripos MOL2 format.
    Uses OpenBabel if available on system, or falls back to an internal Tripos MOL2 serializer.
    """
    # 1. Attempt conversion using OpenBabel CLI
    obabel_candidates = [
        shutil.which("obabel"),
        r"C:\Program Files\OpenBabel-3.1.1\obabel.exe",
        r"C:\Program Files (x86)\OpenBabel-3.1.1\obabel.exe",
    ]
    obabel_bin = next((path for path in obabel_candidates if path and os.path.exists(path)), None)

    if obabel_bin and pdb_source_path and os.path.exists(pdb_source_path):
        try:
            res = subprocess.run(
                [obabel_bin, "-ipdb", pdb_source_path, "-omol2", "-O", mol2_path],
                capture_output=True,
                text=True,
                check=False,
            )
            if res.returncode == 0 and os.path.exists(mol2_path) and os.path.getsize(mol2_path) > 0:
                return True
        except Exception:
            pass

    # 2. Standalone fallback Tripos MOL2 writer without external dependencies
    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()
    mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "FazPrep_Molecule"

    lines = [
        "@<TRIPOS>MOLECULE",
        mol_name,
        f" {num_atoms} {num_bonds} 0 0 0",
        "SMALL",
        "NO_CHARGES",
        "",
        "@<TRIPOS>ATOM",
    ]

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        sym = atom.GetSymbol()
        lines.append(
            f"{idx + 1:>7} {sym + str(idx + 1):<8} {pos.x:>10.4f} {pos.y:>10.4f} {pos.z:>10.4f} {sym:<5} 1 <0> 0.0000"
        )

    lines.append("@<TRIPOS>BOND")
    for bond in mol.GetBonds():
        b_idx = bond.GetIdx()
        b_type = "1"
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            b_type = "2"
        elif bond.GetBondType() == Chem.BondType.TRIPLE:
            b_type = "3"
        elif bond.GetBondType() == Chem.BondType.AROMATIC:
            b_type = "ar"
        lines.append(
            f"{b_idx + 1:>6} {bond.GetBeginAtomIdx() + 1:>5} {bond.GetEndAtomIdx() + 1:>5} {b_type:>4}"
        )

    lines.append("")
    with open(mol2_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    return True


def optimize_flexible_molecule(
    mol: Chem.Mol,
    output_dir: str,
    base_filename: str = "molecule_optimized",
    num_confs: int = 50,
    random_seed: int = 42,
) -> Dict[str, Any]:
    """
    Generates 3D conformers (ETKDGv3), performs energy minimization (MMFF94s / UFF),
    and saves the lowest-energy conformer in PDB, SDF V3000, and MOL2 formats.
    """
    if mol is None:
        raise ValueError("Failed to process molecule. Please verify input data.")

    os.makedirs(output_dir, exist_ok=True)

    # 1. Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # 2. Conformer ensemble generation (ETKDGv3)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(random_seed)
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params))

    if len(conf_ids) == 0:
        conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, useRandomCoords=True))
        if len(conf_ids) == 0:
            raise RuntimeError("Failed to generate 3D conformers for this molecule.")

    # 3. Energy minimization with force field
    best_conf_id = -1
    min_energy = float("inf")
    use_mmff = AllChem.MMFFHasAllMoleculeParams(mol)
    used_ff_name = "MMFF94s" if use_mmff else "UFF"

    for cid in conf_ids:
        if use_mmff:
            ff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
            ff = AllChem.MMFFGetMoleculeForceField(mol, ff_props, confId=cid)
        else:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)

        if ff is not None:
            ff.Minimize()
            energy = ff.CalcEnergy()
            if energy < min_energy:
                min_energy = energy
                best_conf_id = cid

    if best_conf_id == -1:
        raise RuntimeError("Force field energy minimization failed.")

    # Retain only the lowest-energy conformer
    best_mol = Chem.Mol(mol)
    best_mol.RemoveAllConformers()
    best_mol.AddConformer(mol.GetConformer(best_conf_id))

    # 4. Save results to output files
    pdb_path = os.path.join(output_dir, f"{base_filename}.pdb")
    sdf_path = os.path.join(output_dir, f"{base_filename}_v3000.sdf")
    mol2_path = os.path.join(output_dir, f"{base_filename}.mol2")

    # Save PDB
    Chem.MolToPDBFile(best_mol, pdb_path)

    # Save SDF V3000
    writer = Chem.SDWriter(sdf_path)
    writer.SetForceV3000(True)
    writer.write(best_mol)
    writer.close()

    # Save MOL2
    save_mol_to_mol2(best_mol, mol2_path, pdb_source_path=pdb_path)

    return {
        "pdb_path": pdb_path,
        "sdf_path": sdf_path,
        "mol2_path": mol2_path,
        "min_energy": min_energy,
        "force_field": used_ff_name,
        "best_conf_id": best_conf_id,
    }
