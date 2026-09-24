from pathlib import Path

project_root = Path(r"C:\Users\sevco\Documents\phd\proj_2026_Zimcik_Xenium\Xenium_Analysis_SATB2_UMG_2026")


# RAW DATA
cell_boundaries_path = project_root / "data" / "cell_boundaries.parquet"
cell_feature_matrix_path = project_root / "data" / "cell_feature_matrix.parquet"
cells = project_root / "data" / "cells.parquet"
gene_panel_path = project_root / "data" / "gene_panel.json"
nucleus_boundaries_path = project_root / "data" / "nucleus_boundaries.parquet"
transcripts_path = project_root / "data" / "transcripts.parquet"

# PROSEG DATA
proseg_dir = project_root / "data_raw" / "proseg_out"

# 00_preprocess
preprocessed_dir = project_root / "data" / "00_preprocess"
preprocessed_figdir = preprocessed_dir / "figures" / "00_preprocess"

# 07_instant_ST_analysis
instant_dir = project_root / "results" / "07_instant"
instant_figdir = project_root / "figures" / "07_instant"
