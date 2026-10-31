import os
import streamlit as st
import pandas as pd
from rdkit import Chem

# Support running directly or as an imported package
try:
    from .optimizer import (
        get_next_run_dir,
        parse_smiles_lines,
        load_molecules_from_file,
        optimize_flexible_molecule,
        optimize_batch,
    )
except ImportError:
    from optimizer import (
        get_next_run_dir,
        parse_smiles_lines,
        load_molecules_from_file,
        optimize_flexible_molecule,
        optimize_batch,
    )


def main():
    st.set_page_config(
        page_title="FazPrep — Molecular Conformer Generator",
        layout="wide",
    )

    st.title("FazPrep — 3D Conformer & Structure Preparation")
    st.markdown(
        "Automated platform for 3D molecular structure generation, conformer ensemble searching, "
        "and force field energy minimization. Explicit hydrogen atoms are added to all molecules. "
        "Supports single or batch processing (.smi, .sdf, multiline SMILES) with one-click export."
    )

    st.sidebar.header("Parameter Settings")
    num_confs = st.sidebar.slider(
        "Number of 3D Conformers (ETKDGv3)",
        min_value=5,
        max_value=200,
        value=50,
        step=5,
    )
    random_seed = st.sidebar.number_input("Random Seed", value=42, step=1)
    base_filename = st.sidebar.text_input("Default Output Filename / Prefix", value="optimized_molecule")

    st.sidebar.markdown("---")
    st.sidebar.subheader("Chemical Modeling Options")
    explicit_h_info = st.sidebar.checkbox(
        "Add Explicit Hydrogen Atoms",
        value=True,
        disabled=True,
        help="All structures are automatically protonated with full 3D coordinates generated for all hydrogen atoms.",
    )

    tab1, tab2, tab3 = st.tabs([
        "Single SMILES",
        "Multiple SMILES (Batch)",
        "Structure File (.smi / .sdf / .mol2 / .pdb / .mol)",
    ])

    molecules_to_process = []
    is_batch_mode = False

    with tab1:
        single_default_smiles = (
            "C1C=CC=C2C3C=CC=CC=3N(C3C=CC(N4N=NC(C(O)=O)=C4C(F)(F)F)=CC=3)C=12"
        )
        single_smiles = st.text_input("Enter single SMILES string:", value=single_default_smiles)
        single_name = st.text_input("Molecule Name:", value="faz_ligand_1")
        if single_smiles.strip():
            m = Chem.MolFromSmiles(single_smiles.strip())
            if m is not None:
                m.SetProp("_Name", single_name.strip() or "molecule")
                molecules_to_process = [(m, single_name.strip() or "molecule")]
                is_batch_mode = False
            else:
                st.error("Invalid SMILES string.")

    with tab2:
        default_batch_smiles = (
            "CC(=O)Oc1ccccc1C(=O)O Aspirin\n"
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O Ibuprofen\n"
            "CC(=O)Nc1ccc(O)cc1 Paracetamol\n"
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C Caffeine"
        )
        batch_text = st.text_area(
            "Paste multiple SMILES strings (one per line, format: SMILES [Molecule_Name]):",
            value=default_batch_smiles,
            height=160,
        )
        if batch_text.strip():
            parsed_entries = parse_smiles_lines(batch_text)
            batch_mols = []
            for s, n in parsed_entries:
                mol = Chem.MolFromSmiles(s)
                if mol is not None:
                    mol.SetProp("_Name", n)
                    batch_mols.append((mol, n))
            if batch_mols:
                # If tab2 text is selected and contains entries, set batch
                if st.session_state.get("active_input") == "batch_text" or len(batch_mols) > 1:
                    molecules_to_process = batch_mols
                    is_batch_mode = True
                st.info(f"Detected {len(batch_mols)} valid molecule(s) in SMILES list.")

    with tab3:
        uploaded_file = st.file_uploader(
            "Upload molecule structure file (.smi, .txt, .sdf, .mol2, .pdb, .mol)",
            type=["smi", "txt", "sdf", "mol2", "pdb", "mol"],
        )
        if uploaded_file is not None:
            try:
                file_mols = load_molecules_from_file(uploaded_file)
                if file_mols:
                    molecules_to_process = file_mols
                    is_batch_mode = len(file_mols) > 1
                    st.success(f"File '{uploaded_file.name}' loaded successfully: {len(file_mols)} molecule(s) ready.")
                else:
                    st.error("Failed to parse molecule(s) from uploaded file.")
            except Exception as e:
                st.error(f"Error reading file: {e}")

    st.divider()

    st.write(f"**Ready for processing:** {len(molecules_to_process)} molecule(s)")

    if st.button("Run Optimization", type="primary"):
        if not molecules_to_process:
            st.warning("Please provide at least one valid molecule via SMILES or file upload.")
            return

        run_dir, run_num = get_next_run_dir("results")

        if len(molecules_to_process) == 1 and not is_batch_mode:
            # Single molecule processing flow
            mol, mol_name = molecules_to_process[0]
            with st.spinner(f"Generating conformers and minimizing energy for {mol_name}..."):
                try:
                    result = optimize_flexible_molecule(
                        mol=mol,
                        output_dir=run_dir,
                        base_filename=mol_name or base_filename,
                        num_confs=int(num_confs),
                        random_seed=int(random_seed),
                    )

                    st.success(
                        f"Optimization completed. Results saved to: {run_dir} (Run #{run_num})"
                    )

                    col1, col2, col3, col4 = st.columns(4)
                    col1.metric("Minimum Energy", f"{result['min_energy']:.2f} kcal/mol")
                    col2.metric("Force Field", result["force_field"])
                    col3.metric("Formula (with H)", result["formula"])
                    col4.metric("Total Atoms", f"{result['total_atoms']} ({result['heavy_atoms']} heavy)")

                    st.subheader("Generated Files:")
                    c1, c2, c3 = st.columns(3)
                    with c1:
                        with open(result["pdb_path"], "rb") as f:
                            st.download_button(
                                "Download PDB",
                                f,
                                file_name=os.path.basename(result["pdb_path"]),
                                use_container_width=True,
                            )
                    with c2:
                        with open(result["sdf_path"], "rb") as f:
                            st.download_button(
                                "Download SDF (V3000)",
                                f,
                                file_name=os.path.basename(result["sdf_path"]),
                                use_container_width=True,
                            )
                    with c3:
                        with open(result["mol2_path"], "rb") as f:
                            st.download_button(
                                "Download MOL2",
                                f,
                                file_name=os.path.basename(result["mol2_path"]),
                                use_container_width=True,
                            )

                except Exception as e:
                    st.error(f"Execution error: {e}")

        else:
            # Batch optimization flow
            progress_bar = st.progress(0)
            status_text = st.empty()

            def update_progress(current, total, current_name):
                fraction = current / total
                progress_bar.progress(fraction)
                status_text.text(f"Optimizing molecule {current} of {total}: {current_name}...")

            with st.spinner("Processing batch conformer search and energy minimization..."):
                try:
                    batch_res = optimize_batch(
                        molecules_data=molecules_to_process,
                        output_dir=run_dir,
                        num_confs=int(num_confs),
                        random_seed=int(random_seed),
                        progress_callback=update_progress,
                    )

                    progress_bar.progress(1.0)
                    status_text.text("Batch processing finished.")

                    st.success(
                        f"Batch completed successfully. {batch_res['success_count']} of {batch_res['total_count']} "
                        f"molecules processed. Saved to directory: {run_dir} (Run #{run_num})"
                    )

                    # Display Summary Metrics
                    records = batch_res["records"]
                    df = pd.DataFrame(records)

                    valid_energies = [r["Minimum Energy (kcal/mol)"] for r in records if r["Minimum Energy (kcal/mol)"] is not None]
                    avg_energy = sum(valid_energies) / len(valid_energies) if valid_energies else 0.0

                    m1, m2, m3, m4 = st.columns(4)
                    m1.metric("Total Molecules", batch_res["total_count"])
                    m2.metric("Successful", batch_res["success_count"])
                    m3.metric("Failed", batch_res["failed_count"])
                    m4.metric("Mean Energy", f"{avg_energy:.2f} kcal/mol")

                    st.subheader("Single-Click Batch Exports:")
                    st.markdown("Export the lowest-energy conformer for every molecule with a single click:")

                    col_exp1, col_exp2, col_exp3 = st.columns(3)

                    with col_exp1:
                        with open(batch_res["combined_sdf_path"], "rb") as f:
                            st.download_button(
                                "Download All Best Conformers (Combined SDF)",
                                f,
                                file_name="all_best_conformers_v3000.sdf",
                                mime="chemical/x-mdl-sdfile",
                                use_container_width=True,
                            )

                    with col_exp2:
                        with open(batch_res["zip_path"], "rb") as f:
                            st.download_button(
                                "Download Complete Archive (ZIP)",
                                f,
                                file_name="all_optimized_structures.zip",
                                mime="application/zip",
                                use_container_width=True,
                            )

                    with col_exp3:
                        with open(batch_res["csv_report_path"], "rb") as f:
                            st.download_button(
                                "Download Summary Report (CSV)",
                                f,
                                file_name="summary_report.csv",
                                mime="text/csv",
                                use_container_width=True,
                            )

                    st.subheader("Optimization Summary Table:")
                    st.dataframe(df, use_container_width=True)

                except Exception as e:
                    st.error(f"Batch execution error: {e}")


if __name__ == "__main__":
    main()
