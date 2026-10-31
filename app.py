import os
import streamlit as st
from rdkit import Chem

# Support running directly or as an imported package
try:
    from .optimizer import (
        get_next_run_dir,
        load_molecule_from_file,
        optimize_flexible_molecule,
    )
except ImportError:
    from optimizer import (
        get_next_run_dir,
        load_molecule_from_file,
        optimize_flexible_molecule,
    )


def main():
    st.set_page_config(
        page_title="FazPrep — Molecular Conformer Generator",
        layout="wide",
    )

    st.title("FazPrep — 3D Conformer & Structure Preparation")
    st.markdown(
        "Automated platform for 3D molecular structure generation, conformer ensemble searching, "
        "and force field energy minimization with automatic run archiving (`results/run_X`)."
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
    base_filename = st.sidebar.text_input("Base Output Filename", value="optimized_molecule")

    tab1, tab2 = st.tabs(["SMILES Input", "File Upload (SDF / MOL2 / PDB / MOL)"])

    mol_input = None

    with tab1:
        default_smiles = (
            "C1C=CC=C2C3C=CC=CC=3N(C3C=CC(N4N=NC(C(O)=O)=C4C(F)(F)F)=CC=3)C=12"
        )
        smiles_str = st.text_input("Enter SMILES string:", value=default_smiles)
        if smiles_str.strip():
            mol_input = Chem.MolFromSmiles(smiles_str.strip())
            if mol_input is None:
                st.error("Invalid SMILES string.")

    with tab2:
        uploaded_file = st.file_uploader(
            "Upload molecule structure file",
            type=["sdf", "mol2", "pdb", "mol"],
        )
        if uploaded_file is not None:
            try:
                mol_input = load_molecule_from_file(uploaded_file)
                if mol_input is not None:
                    st.success(f"File '{uploaded_file.name}' loaded successfully.")
                else:
                    st.error("Failed to parse molecule file.")
            except Exception as e:
                st.error(f"Error reading molecule file: {e}")

    st.divider()

    if st.button("Run Optimization", type="primary"):
        if mol_input is None:
            st.warning("Please provide a valid SMILES string or upload a valid molecule file.")
            return

        with st.spinner("Generating conformers and minimizing energy..."):
            try:
                # Save runs into results/run_X
                run_dir, run_num = get_next_run_dir("results")

                results = optimize_flexible_molecule(
                    mol=mol_input,
                    output_dir=run_dir,
                    base_filename=base_filename,
                    num_confs=int(num_confs),
                    random_seed=int(random_seed),
                )

                st.success(
                    f"Process completed successfully. Results saved to: {run_dir} (Run #{run_num})"
                )

                col1, col2, col3 = st.columns(3)
                col1.metric("Minimum Energy", f"{results['min_energy']:.2f} kcal/mol")
                col2.metric("Force Field", results["force_field"])
                col3.metric("Best Conformer ID", results["best_conf_id"])

                st.subheader("Generated Files:")

                # Download buttons
                with open(results["pdb_path"], "rb") as f:
                    st.download_button(
                        "Download PDB",
                        f,
                        file_name=os.path.basename(results["pdb_path"]),
                    )

                with open(results["sdf_path"], "rb") as f:
                    st.download_button(
                        "Download SDF (V3000)",
                        f,
                        file_name=os.path.basename(results["sdf_path"]),
                    )

                with open(results["mol2_path"], "rb") as f:
                    st.download_button(
                        "Download MOL2",
                        f,
                        file_name=os.path.basename(results["mol2_path"]),
                    )

            except Exception as e:
                st.error(f"Execution error: {e}")


if __name__ == "__main__":
    main()
