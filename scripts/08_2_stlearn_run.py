#08_2_stlearn_run
"""
Run stLearn cell-cell interaction (CCI) analysis on one or more samples.

Example:
    python 08_2_stlearn_run.py \
        --samples wt_1 wt_2 wt_4 \
        --database connectomeDB2020_lit \
        --input-dir /mnt/scratch/sevcoviz/proj_2026_Zimcik_Xenium/data/08_stlearn \
        --output-dir /mnt/scratch/sevcoviz/proj_2026_Zimcik_Xenium/data/08_stlearn/cci \
        --fig-dir /mnt/scratch/sevcoviz/proj_2026_Zimcik_Xenium/figures/08_stlearn \
        --min-spots 20 --distance 30 --n-pairs 10000 --n-cpus 16 \
        --annotation annotation --cci-min-spots 5 --n-perms 1000
"""

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import stlearn as st


DEFAULT_DIR = "/mnt/scratch/sevcoviz/proj_2026_Zimcik_Xenium/data/08_stlearn"


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "1", "y"):
        return True
    if v.lower() in ("no", "false", "f", "0", "n"):
        return False
    raise argparse.ArgumentTypeError("Boolean value expected.")


def parse_samples(values):
    # Accept both "wt_1 wt_2" and "wt_1,wt_2"
    samples = []
    for v in values:
        samples.extend(s.strip() for s in v.split(",") if s.strip())
    return samples


def parse_args():
    parser = argparse.ArgumentParser(
        description="stLearn CCI analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # I/O
    io = parser.add_argument_group("Input / output")
    io.add_argument("--samples", nargs="+", required=True,
                    help="Samples to process, space- or comma-separated (e.g. wt_1 wt_2 wt_4)")
    io.add_argument("--input-dir", default=DEFAULT_DIR,
                    help="Directory containing <sample>.h5ad")
    io.add_argument("--output-dir", default=DEFAULT_DIR,
                    help="Directory to save results (<sample>_<database>_cci.h5ad)")
    io.add_argument("--fig-dir", default=DEFAULT_DIR,
                    help="Directory to save figures")

    # LR database
    db = parser.add_argument_group("LR database")
    db.add_argument("--database", nargs="+", default=["connectomeDB2020_lit"],
                    choices=["connectomeDB2020_lit", "connectomeDB2020_put"],
                    help="Ligand-receptor database(s)")
    db.add_argument("--species", default="mouse", choices=["mouse", "human"])

    # st.tl.cci.run
    run = parser.add_argument_group("st.tl.cci.run parameters")
    run.add_argument("--min-spots", type=int, default=20,
                     help="Filter out LR pairs with scores in fewer than min_spots")
    run.add_argument("--distance", type=float, default=30,
                     help="Single-cell mode: neighbourhood radius in coordinate units (um); "
                          "0 = within-cell. Ignored with --grid-n (see --grid-neighbours).")
    run.add_argument("--n-pairs", type=int, default=10000,
                     help="Number of random pairs for background")
    run.add_argument("--n-cpus", type=int, default=None,
                     help="Number of CPUs (default: all available)")

    # st.tl.cci.run_cci
    cci = parser.add_argument_group("st.tl.cci.run_cci parameters")
    cci.add_argument("--annotation", default="annotation",
                     help="obs column with cell type labels")
    cci.add_argument("--cci-min-spots", type=int, default=5,
                     help="min_spots for run_cci")
    cci.add_argument("--cell-prop-cutoff", type=float, default=0.1,
                     help="Min cell type proportion per grid spot; only used with --grid-n")
    cci.add_argument("--sig-spots", type=str2bool, default=True)
    cci.add_argument("--n-perms", type=int, default=1000)

    # st.tl.cci.grid
    grid = parser.add_argument_group("Gridding (st.tl.cci.grid)")
    grid.add_argument("--grid-n", type=int, default=None,
                      help="Aggregate cells into an N x N grid before CCI (sets spot_mixtures=True). "
                           "Omit to run at single-cell resolution (spot_mixtures=False).")
    grid.add_argument("--grid-neighbours", type=float, default=1.0,
                      help="Grid mode: neighbourhood radius in grid-spot spacings; distance = "
                           "value x spot spacing. 1 = 4 direct neighbours, 1.5 = 8 incl. diagonals, "
                           "0 = within-spot.")

    # Plotting
    pl = parser.add_argument_group("Plotting")
    pl.add_argument("--n-top", type=int, default=30,
                    help="Number of top LR pairs in lr_summary plot")

    args = parser.parse_args()
    args.samples = parse_samples(args.samples)
    # Gridded spots hold cell type proportions (in .uns[annotation]); single cells hold one label
    args.spot_mixtures = args.grid_n is not None
    return args


def prepare_adata(adata):
    adata = adata.copy()

    # Copy spatial coordinates from the same AnnData object
    adata.obs["imagecol"] = adata.obsm["spatial"][:, 0]
    adata.obs["imagerow"] = adata.obsm["spatial"][:, 1]

    # Create minimal spatial metadata required by stLearn
    adata.uns["spatial"] = {
        "library": {
            "scalefactors": {
                "spot_diameter_fullres": 1.0,
                "tissue_hires_scalef": 1.0,
                "tissue_lowres_scalef": 1.0
            },
            # No tissue image; empty dict (not {"hires": None}) so it survives
            # write_h5ad and stLearn plots find the key but no image
            "images": {},
            "metadata": {}
        }
    }

    return adata


def add_safe_labels(adata, annotation):
    """Copy obs[annotation] to obs[f"{annotation}_safe"] with h5py-safe category names.

    stLearn stores DataFrames with cell types as columns in .uns (grid proportions,
    per-LR CCI matrices); anndata writes each column as an HDF5 key, and h5py forbids
    "/" in keys. Returns the new column name and a {safe: original} mapping.
    """
    labels = adata.obs[annotation].astype("category")
    safe_name = {c: re.sub(r"\s*/\s*", "_", str(c)).strip() for c in labels.cat.categories}

    if len(set(safe_name.values())) != len(safe_name):
        raise ValueError(f"Sanitised cell type names are not unique: {safe_name}")

    safe_col = f"{annotation}_safe"
    # Keeps category order (and therefore matching colours)
    adata.obs[safe_col] = labels.cat.rename_categories(safe_name)
    if f"{annotation}_colors" in adata.uns:
        adata.uns[f"{safe_col}_colors"] = adata.uns[f"{annotation}_colors"]

    for orig, safe in safe_name.items():
        if orig != safe:
            print(f"  renamed cell type: '{orig}' -> '{safe}'")
    return safe_col, {safe: orig for orig, safe in safe_name.items()}


def save_adata(adata, path):
    adata = adata.copy()
    # None values cannot be written to h5ad; stLearn plotting needs the "images" key
    # to exist after reloading, so store it empty rather than deleting it
    for lib in adata.uns.get("spatial", {}).values():
        lib["images"] = {}
    adata.write_h5ad(path)


def process_sample(sample, lrs, run_tag, args):
    in_path = Path(args.input_dir) / f"{sample}.h5ad"
    print(f"\n=== {sample} ===\nReading {in_path}")
    adata = sc.read_h5ad(in_path)
    adata = prepare_adata(adata)
    # Cell type names end up as HDF5 keys, so "/" must go; use the safe column throughout.
    # It is categorical, which grid() requires (it calls .cat on the label column).
    label, safe_to_orig = add_safe_labels(adata, args.annotation)

    distance = args.distance
    if args.grid_n is not None:
        n_ = args.grid_n
        print(f"Gridding {adata.n_obs} cells into {n_} x {n_} = {n_ * n_} spots")
        # Returns a new AnnData: expression summed per bin, cell type proportions in
        # .uns[annotation], dominant type in .obs[annotation], coords = bin centres
        adata = st.tl.cci.grid(adata, n_row=n_, n_col=n_,
                               use_label=label, n_cpus=args.n_cpus)
        print(f"Grid has {adata.n_obs} non-empty spots")

        # Bins are rectangular (x and y ranges differ); use the larger side so that
        # grid_neighbours=1.5 always reaches the adjacent spot in both directions.
        # Small tolerance because cKDTree compares float distances to the radius.
        #      1 = direct neighbours       1.5 = + diagonals        2 = two spots away
                    #   . N .                      N N N                    . . N . .
                    #   N X N                      N X N                    . N N N .
                    #   . N .                      N N N                    N N X N N
                    #                                                       . N N N .
                    #                                                       . . N . .

        spacing = max(np.diff(adata.uns["grid_xedges"]).mean(),
                      np.diff(adata.uns["grid_yedges"]).mean())
        distance = args.grid_neighbours * spacing * 1.001
        print(f"Grid spot spacing {spacing:.2f}; distance = {args.grid_neighbours} x spacing = {distance:.2f}")
    else:
        print(f"Single-cell mode; distance = {distance}")

    st.tl.cci.run(adata, lrs,
                  min_spots=args.min_spots,  # Filter out any LR pairs with no scores for less than min_spots
                  distance=distance,         # Must be numeric: distance=None needs Visium metadata we don't have; 0 = within-spot
                  n_pairs=args.n_pairs,      # Number of random pairs to generate; recommend ~10,000
                  n_cpus=args.n_cpus,        # Number of CPUs for parallel. If None, detects & use all available.
                  )

    st.pl.lr_summary(adata, n_top=args.n_top, figsize=(6, 4))
    fig_path = Path(args.fig_dir) / f"{sample}_{run_tag}_lr_summary.png"
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"Saved {fig_path}")

    st.tl.cci.run_cci(
        adata,
        label,
        min_spots=args.cci_min_spots,
        spot_mixtures=args.spot_mixtures,
        cell_prop_cutoff=args.cell_prop_cutoff,
        sig_spots=args.sig_spots,
        n_perms=args.n_perms,
        n_cpus=args.n_cpus,
    )

    # Keep the mapping back to the original names (keys are safe, values may contain "/")
    adata.uns[f"{label}_to_original"] = safe_to_orig

    out_path = Path(args.output_dir) / f"{sample}_{run_tag}_cci.h5ad"
    save_adata(adata, out_path)
    print(f"Saved {out_path}")

def main():
    args = parse_args()
    print("Arguments:")
    for k, v in vars(args).items():
        print(f"  {k}: {v}")

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    Path(args.fig_dir).mkdir(parents=True, exist_ok=True)

    lrs = st.tl.cci.load_lrs(args.database, species=args.species)
    print(f"Loaded {len(lrs)} LR pairs from {args.database}")
    run_tag = "_".join(d.replace("connectomeDB2020_", "") for d in args.database)
    run_tag += f"_grid{args.grid_n}" if args.grid_n is not None else "_sc"

    for sample in args.samples:
        process_sample(sample, lrs, run_tag, args)


if __name__ == "__main__":
    main()
