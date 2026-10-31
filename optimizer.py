import os
import glob
import tempfile
import shutil
import subprocess
import zipfile
import csv
from typing import Optional, Dict, Any, Tuple, List, Callable
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


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


def parse_smiles_lines(content: str) -> List[Tuple[str, str]]:
    """
    Parses multi-line SMILES text (or .smi file content) into list of (smiles, molecule_name) tuples.
    Supports whitespace-separated or tab-separated formats: SMILES [Molecule_Name].
    """
    molecules = []
    lines = content.strip().splitlines()
    for i, line in enumerate(lines, start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        smiles = parts[0]
        name = parts[1].strip() if len(parts) > 1 else f"molecule_{i}"
        molecules.append((smiles, name))
    return molecules


def load_molecules_from_file(uploaded_file, filename: Optional[str] = None) -> List[Tuple[Chem.Mol, str]]:
    """
    Reads molecules from uploaded files.
    Supports multi-molecule files (.smi, .txt, .sdf) and single-molecule files (.mol, .pdb, .mol2).
    Returns a list of (RDKit Mol, molecule_name) tuples.
    """
    fname = filename or getattr(uploaded_file, "name", "molecule.mol")
    file_ext = os.path.splitext(fname)[1].lower()
    base_name = os.path.splitext(fname)[0]

    if hasattr(uploaded_file, "read"):
        file_bytes = uploaded_file.read()
        if hasattr(uploaded_file, "seek"):
            uploaded_file.seek(0)
    else:
        with open(uploaded_file, "rb") as f:
            file_bytes = f.read()

    results: List[Tuple[Chem.Mol, str]] = []

    # 1. SMILES format (.smi, .txt)
    if file_ext in [".smi", ".txt"]:
        text = file_bytes.decode("utf-8", errors="ignore")
        parsed_entries = parse_smiles_lines(text)
        for smiles, name in parsed_entries:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                mol.SetProp("_Name", name)
                results.append((mol, name))
        return results

    # 2. SDF / MOL multi-supplier
    if file_ext in [".sdf", ".mol"]:
        with tempfile.NamedTemporaryFile(delete=False, suffix=file_ext) as tmp:
            tmp.write(file_bytes)
            tmp_path = tmp.name

        try:
            suppl = Chem.SDMolSupplier(tmp_path, removeHs=False)
            for i, mol in enumerate(suppl, start=1):
                if mol is not None:
                    name = mol.GetProp("_Name") if mol.HasProp("_Name") and mol.GetProp("_Name").strip() else f"{base_name}_{i}"
                    mol.SetProp("_Name", name)
                    results.append((mol, name))
        finally:
            suppl = None
            if os.path.exists(tmp_path):
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass

        if results:
            return results

    # 3. PDB format
    if file_ext == ".pdb":
        mol = Chem.MolFromPDBBlock(file_bytes.decode("utf-8", errors="ignore"), removeHs=False)
        if mol is not None:
            mol.SetProp("_Name", base_name)
            results.append((mol, base_name))
            return results

    # 4. MOL2 format
    if file_ext == ".mol2":
        mol = Chem.MolFromMol2Block(file_bytes.decode("utf-8", errors="ignore"), removeHs=False)
        if mol is not None:
            mol.SetProp("_Name", base_name)
            results.append((mol, base_name))
            return results

    # Fallback via disk read
    with tempfile.NamedTemporaryFile(delete=False, suffix=file_ext) as tmp:
        tmp.write(file_bytes)
        tmp_path = tmp.name

    try:
        if file_ext == ".mol2":
            mol = Chem.MolFromMol2File(tmp_path, removeHs=False)
        elif file_ext in [".sdf", ".mol"]:
            mol = Chem.MolFromMolFile(tmp_path, removeHs=False)
        elif file_ext == ".pdb":
            mol = Chem.MolFromPDBFile(tmp_path, removeHs=False)
        else:
            mol = None

        if mol is not None:
            mol.SetProp("_Name", base_name)
            results.append((mol, base_name))
    finally:
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass

    return results


def save_mol_to_mol2(mol: Chem.Mol, mol2_path: str, pdb_source_path: Optional[str] = None) -> bool:
    """
    Saves a molecule in Tripos MOL2 format.
    Uses OpenBabel if available, or falls back to an internal standalone Tripos MOL2 serializer.
    """
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

    # Standalone fallback Tripos MOL2 writer
    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()
    mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else "Molecule"

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
    Adds explicit hydrogens, generates 3D conformers (ETKDGv3), minimizes energy (MMFF94s / UFF),
    and saves the single best conformer in PDB, SDF V3000, and MOL2 formats.
    """
    if mol is None:
        raise ValueError("Failed to process molecule. Please verify input data.")

    os.makedirs(output_dir, exist_ok=True)

    # 1. Add explicit hydrogen atoms (All-atom 3D model)
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
    best_mol.SetProp("_Name", base_filename)

    # 4. Save results to individual files
    pdb_path = os.path.join(output_dir, f"{base_filename}.pdb")
    sdf_path = os.path.join(output_dir, f"{base_filename}_v3000.sdf")
    mol2_path = os.path.join(output_dir, f"{base_filename}.mol2")

    Chem.MolToPDBFile(best_mol, pdb_path)

    writer = Chem.SDWriter(sdf_path)
    writer.SetForceV3000(True)
    writer.write(best_mol)
    writer.close()

    save_mol_to_mol2(best_mol, mol2_path, pdb_source_path=pdb_path)

    # Molecular formula and atom counts
    formula = rdMolDescriptors.CalcMolFormula(best_mol)
    total_atoms = best_mol.GetNumAtoms()
    heavy_atoms = best_mol.GetNumHeavyAtoms()

    return {
        "pdb_path": pdb_path,
        "sdf_path": sdf_path,
        "mol2_path": mol2_path,
        "min_energy": min_energy,
        "force_field": used_ff_name,
        "best_conf_id": best_conf_id,
        "best_mol": best_mol,
        "formula": formula,
        "total_atoms": total_atoms,
        "heavy_atoms": heavy_atoms,
    }


def optimize_batch(
    molecules_data: List[Tuple[Chem.Mol, str]],
    output_dir: str,
    num_confs: int = 50,
    random_seed: int = 42,
    progress_callback: Optional[Callable[[int, int, str], None]] = None,
) -> Dict[str, Any]:
    """
    Optimizes a batch of molecules, selecting the lowest-energy conformer for each.
    Exports individual structures (PDB, SDF, MOL2), a unified multi-molecule SDF file,
    a comprehensive CSV summary report, and a consolidated ZIP archive.
    """
    os.makedirs(output_dir, exist_ok=True)
    total_count = len(molecules_data)

    successful_mols: List[Chem.Mol] = []
    summary_records: List[Dict[str, Any]] = []

    for idx, (mol, name) in enumerate(molecules_data, start=1):
        clean_name = "".join(c if (c.isalnum() or c in "-_") else "_" for c in name).strip("_")
        if not clean_name:
            clean_name = f"molecule_{idx}"

        if progress_callback:
            progress_callback(idx, total_count, clean_name)

        try:
            smiles = Chem.MolToSmiles(mol) if mol is not None else "N/A"
            res = optimize_flexible_molecule(
                mol=mol,
                output_dir=output_dir,
                base_filename=clean_name,
                num_confs=num_confs,
                random_seed=random_seed,
            )
            successful_mols.append(res["best_mol"])
            summary_records.append({
                "Molecule Name": clean_name,
                "SMILES": smiles,
                "Formula": res["formula"],
                "Heavy Atoms": res["heavy_atoms"],
                "Total Atoms (with Hydrogens)": res["total_atoms"],
                "Minimum Energy (kcal/mol)": round(res["min_energy"], 3),
                "Force Field": res["force_field"],
                "Best Conformer Index": res["best_conf_id"],
                "Status": "Success",
                "PDB File": os.path.basename(res["pdb_path"]),
                "SDF File": os.path.basename(res["sdf_path"]),
                "MOL2 File": os.path.basename(res["mol2_path"]),
            })
        except Exception as err:
            smiles = Chem.MolToSmiles(mol) if mol is not None else "N/A"
            summary_records.append({
                "Molecule Name": clean_name,
                "SMILES": smiles,
                "Formula": "N/A",
                "Heavy Atoms": mol.GetNumHeavyAtoms() if mol is not None else 0,
                "Total Atoms (with Hydrogens)": 0,
                "Minimum Energy (kcal/mol)": None,
                "Force Field": "N/A",
                "Best Conformer Index": -1,
                "Status": f"Error: {err}",
                "PDB File": "",
                "SDF File": "",
                "MOL2 File": "",
            })

    # 1. Create combined multi-molecule SDF file
    combined_sdf_path = os.path.join(output_dir, "all_best_conformers_v3000.sdf")
    writer = Chem.SDWriter(combined_sdf_path)
    writer.SetForceV3000(True)
    for bmol in successful_mols:
        writer.write(bmol)
    writer.close()

    # 2. Create CSV summary report
    csv_report_path = os.path.join(output_dir, "summary_report.csv")
    fieldnames = [
        "Molecule Name",
        "SMILES",
        "Formula",
        "Heavy Atoms",
        "Total Atoms (with Hydrogens)",
        "Minimum Energy (kcal/mol)",
        "Force Field",
        "Best Conformer Index",
        "Status",
        "PDB File",
        "SDF File",
        "MOL2 File",
    ]
    with open(csv_report_path, "w", newline="", encoding="utf-8") as f:
        csv_writer = csv.DictWriter(f, fieldnames=fieldnames)
        csv_writer.writeheader()
        csv_writer.writerows(summary_records)

    # 3. Create unified ZIP archive with all structures and reports
    zip_path = os.path.join(output_dir, "all_optimized_structures.zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
        for fname in os.listdir(output_dir):
            fpath = os.path.join(output_dir, fname)
            if os.path.isfile(fpath) and fname != "all_optimized_structures.zip":
                zipf.write(fpath, arcname=fname)

    return {
        "total_count": total_count,
        "success_count": len(successful_mols),
        "failed_count": total_count - len(successful_mols),
        "combined_sdf_path": combined_sdf_path,
        "csv_report_path": csv_report_path,
        "zip_path": zip_path,
        "records": summary_records,
    }
