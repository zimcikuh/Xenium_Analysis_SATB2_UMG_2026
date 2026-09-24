import logging as log
import numpy as np
import anndata as ad
from typing import Optional
import scipy.sparse as sp
import scanpy as sc


def adata_preprocessing(
    adata: ad.AnnData,
    min_cells: int = 3,
    min_genes: int = 200,
    target_sum: float = 1e4,
    log1p: bool = True,
    scale_max_value: float = 10,
    npcs: int = 50,
    num_highly_variable_genes: int = 2000,
    leiden_res: Optional[list] = None,
    n_neighbors: int = 15,
    min_dist: float = 0.3,
    spread: float = 0.8,
) -> ad.AnnData:

    if leiden_res is None:
        leiden_res = [0.5, 1.0, 1.5, 2.0]

    # ── 1. Preserve raw counts + spatial before any filtering ─────────────────
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
        log.info("  stored raw counts in layers['counts']")

    # ── 2. Filter ─────────────────────────────────────────────────────────────
    n_cells_before, n_genes_before = adata.n_obs, adata.n_vars
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    log.info("  filter: %d → %d cells,  %d → %d genes",
                n_cells_before, adata.n_obs, n_genes_before, adata.n_vars)

    # ── 3. Normalize ──────────────────────────────────────────────────────────
    sc.pp.normalize_total(adata, target_sum=target_sum)

    # ── 4. Log1p ──────────────────────────────────────────────────────────────
    if log1p:
        sc.pp.log1p(adata)

    # ── 5. Store lognorm layer (needed for sq.pl.spatial_scatter) ─────────────
    adata.layers["lognorm"] = adata.X.copy()

    # ── 6. Highly variable genes ──────────────────────────────────────────────
    if "counts" in adata.layers:
        print("Using seurat_v3 flavor for HVG selection (counts layer found)")
        sc.pp.highly_variable_genes(
            adata, n_top_genes=num_highly_variable_genes, flavor="seurat_v3", layer="counts",
        )
        log.info("  HVG: seurat_v3, n_top=%d", num_highly_variable_genes)
    else:
        print("No counts layer found, using seurat flavor for HVG selection")
        sc.pp.highly_variable_genes(adata, n_top_genes=num_highly_variable_genes, flavor="seurat")


    log.info("  HVG: %d / %d genes selected",
                adata.var["highly_variable"].sum(), adata.n_vars)

    # ── 7. Save .raw (full lognorm matrix, all genes) ─────────────────────────
    adata.raw = adata.copy()

    # ── 8. Subset to HVGs ─────────────────────────────────────────────────────
    adata = adata[:, adata.var["highly_variable"]].copy()

    # ── 9. Scale ──────────────────────────────────────────────────────────────
    sc.pp.scale(adata, zero_center=True, max_value=scale_max_value)

    # Sanitise any NaN / Inf introduced by scaling zero-variance genes
    if sp.issparse(adata.X):
        adata.X.data = np.nan_to_num(adata.X.data, nan=0.0, posinf=0.0, neginf=0.0)
    else:
        adata.X = np.nan_to_num(np.array(adata.X, dtype=np.float32),
                                nan=0.0, posinf=0.0, neginf=0.0)

    # ── 10. PCA → neighbors → UMAP → Leiden ──────────────────────────────────
    n_comps = min(npcs, adata.n_obs - 1, adata.n_vars - 1)
    sc.pp.pca(adata, n_comps=n_comps)
    log.info("  PCA: %d components", n_comps)
    # write explained variance for first 5 components
    for i in range(1, min(6, n_comps + 1)):
        log.info("  PCA: component %d explains %.2f%% of variance", i, adata.uns["pca"]["variance_ratio"][i - 1] * 100)

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps)
    sc.tl.umap(adata, min_dist=min_dist, spread=spread, random_state=42)

    for res in leiden_res:
        sc.tl.leiden(adata, flavor="igraph", n_iterations=-1,
                        resolution=res, key_added=f"leiden_{res}")
    log.info("  Leiden: resolutions %s", leiden_res)
    log.info("  preprocessing complete  (%d cells × %d HVGs)", adata.n_obs, adata.n_vars)

    return adata



### ----------------------------------------------
# Solve the problem of missing spatialdata_attrs in the table model of C5 (Satb2-/-_1) and C3 (Satb2+/+_3)
### ------------------------------------------------
# import gc
# import spatialdata as sd
# from spatialdata.models import TableModel
# from spatialdata.models import PointsModel, TableModel
# from spatialdata.transformations import get_transformation

# # Resume with C3; C1 and C2 are already repaired.
# for sample in ["C3", "C4", "C5", "C6", "C7", "C8"]:
#     src = RAW_DIR / sample / "proseg-output.zarr"
#     dst = FIXED_DIR / sample / "proseg-output.zarr"

#     print(f"Repairing {sample}...")
#     sdata = sd.read_zarr(src)

#     points = sdata.points["transcripts"]

#     # Preserve the spatial transformation before Dask operations.
#     transformations = get_transformation(points, get_all=True)

#     points["gene"] = points["gene"].astype("string")
#     points["assignment"] = points["assignment"].astype("UInt32")
#     points = points.repartition(partition_size="128MB")

#     # Restore SpatialData metadata removed by repartition().
#     points = PointsModel.parse(
#         points,
#         transformations=transformations,
#     )

#     sdata.points["transcripts"] = points

#     # Standardize table metadata and repair C5.
#     table = sdata.tables["table"].copy()

#     if "spatial" not in table.obsm:
#         table.obsm["spatial"] = table.obs[
#             ["centroid_x", "centroid_y"]
#         ].to_numpy()

#     table.uns.pop("spatialdata_attrs", None)

#     sdata.tables["table"] = TableModel.parse(
#         table,
#         region="cell_boundaries",
#         region_key="region",
#         instance_key="cell",
#     )

#     dst.parent.mkdir(parents=True, exist_ok=True)

#     # Replaces the partial C3 output from the failed validation.
#     sdata.write(dst, overwrite=True)

#     repaired = sd.read_zarr(dst)

#     assert repaired.is_self_contained()
#     assert "spatialdata_attrs" in repaired.tables["table"].uns
#     assert "spatial" in repaired.tables["table"].obsm

#     print(
#         f"{sample}: OK — "
#         f"{repaired.tables['table'].n_obs:,} cells, "
#         f"{repaired.points['transcripts'].npartitions} transcript partitions"
#     )

#     del sdata, repaired, table, points
#     gc.collect()


from utils.utils import PALETTE, LEVEL1_COLORS, L2_COLORS, L2_TO_L1
import pandas as pd
import matplotlib.pyplot as plt

def apply_palette(adata, key='final_annotation_code'):
    present = set(adata.obs[key].astype(str))
    if (miss := present - set(L2_COLORS)): print('Not in palette (-> grey):', sorted(miss))
    order = [t for t in L2_COLORS if t in present] + sorted(miss)
    adata.obs[key] = pd.Categorical(adata.obs[key].astype(str), categories=order)
    adata.uns[f'{key}_colors'] = [L2_COLORS.get(t, '#BDBDBD') for t in order]
    l1 = adata.obs[key].astype(str).map(lambda t: L2_TO_L1.get(t, 'Other'))
    g_order = [g for g in PALETTE if g in set(l1)] + (['Other'] if miss else [])
    adata.obs['level1'] = pd.Categorical(l1, categories=g_order)
    adata.uns['level1_colors'] = [LEVEL1_COLORS.get(g, '#BDBDBD') for g in g_order]

def show_palette():
    fig, axes = plt.subplots(1, len(PALETTE), figsize=(2.6 * len(PALETTE), 3))
    for ax, (g, d) in zip(axes, PALETTE.items()):
        ax.set_title(g, color=LEVEL1_COLORS[g], weight='bold', fontsize=9)
        for i, (t, c) in enumerate(d.items()):
            ax.add_patch(plt.Rectangle((0, -i), .12, .8, color=c))
            ax.text(.16, -i + .4, t.replace('_', ' '), va='center', fontsize=6.5)
        ax.set_xlim(0, 1); ax.set_ylim(-6.5, 1); ax.axis('off')
    plt.tight_layout(); plt.show()

######### COLOCALIZATION helper functions ##############
# ---- anchor and data -------------------------------------------------------
WAVE_GENE   = "Col1a1"
LAYER       = "counts"        # layer used for smoothing
SPATIAL_KEY = "spatial"        # adata.obsm key with micrometre coordinates
CELLTYPE_KEY = None            # auto-detected below; set explicitly to override
N_ROUNDS = 2

# ---- inference -------------------------------------------------------------
BLOCK_UM          = 300        # side length of spatial blocks (bootstrap / split-half)
N_BOOT            = 500        # block-bootstrap replicates
N_CONTROL_ANCHORS = 200        # matched random anchor genes forming the null
FDR_ALPHA         = 0.05
N_TOP             = 25         # genes carried into heatmap / enrichment

# ---- sensitivity sweep -----------------------------------------------------
SWEEP_NEIGHS = (6, 10, 20)
SWEEP_ROUNDS = (1, 2, 4, 8)



def corr_cols(X, y):
    """Pearson r between every column of X and the vector y."""
    Xc = X - X.mean(axis=0, keepdims=True)
    yc = y - y.mean()
    xs = np.sqrt(np.einsum("ij,ij->j", Xc, Xc))
    ys = np.sqrt(float(yc @ yc))
    denom = xs * ys
    out = np.full(X.shape[1], np.nan, dtype=np.float64)
    ok = denom > 0
    out[ok] = (Xc[:, ok].T @ yc) / denom[ok]
    return out


def corr_matrix(X, Y):
    """Pearson r between every column of X and every column of Y -> (nY, nX)."""
    Xc = X - X.mean(axis=0, keepdims=True)
    Yc = Y - Y.mean(axis=0, keepdims=True)
    xs = np.sqrt(np.einsum("ij,ij->j", Xc, Xc))
    ys = np.sqrt(np.einsum("ij,ij->j", Yc, Yc))
    num = Yc.T @ Xc
    denom = np.outer(ys, xs)
    return np.divide(num, denom, out=np.full_like(num, np.nan), where=denom > 0)

def get_coords(adata, key=SPATIAL_KEY):
    return np.asarray(adata.obsm[key], dtype=np.float64)[:, :2]



def row_stochastic(W):
    n = W.shape[0]
    W = (W + sp.eye(n, format="csr", dtype=np.float32)).tocsr()
    rs = np.asarray(W.sum(axis=1)).ravel()
    rs[rs == 0] = 1.0
    return (sp.diags(1.0 / rs).astype(np.float32) @ W).tocsr()


def smooth(P, X, n_rounds=N_ROUNDS):
    for _ in range(n_rounds):
        X = P @ X
    return X


def densify(M):
    return M.toarray() if sp.issparse(M) else np.asarray(M)


def effective_kernel_radius(P, coords, n_rounds=N_ROUNDS, n_sample=300, rng=None):
    """Weight-averaged Euclidean radius reached after `n_rounds` propagations."""
    rng = np.random.default_rng(42) if rng is None else rng
    n = P.shape[0]
    pick = rng.choice(n, size=min(n_sample, n), replace=False)
    E = sp.csr_matrix((np.ones(pick.size, np.float32),
                       (np.arange(pick.size), pick)), shape=(pick.size, n))
    for _ in range(n_rounds):
        E = (E @ P).tocsr()

    radii = np.empty(pick.size)
    for i, c in enumerate(pick):
        row = E.getrow(i).tocoo()
        d = np.linalg.norm(coords[row.col] - coords[c], axis=1)
        radii[i] = float((row.data * d).sum() / row.data.sum())
    return radii


def morans_I(W, X):
    """Moran's I for every column of X on the binary symmetrised graph W."""
    Wsym = ((W + W.T) > 0).astype(np.float32)
    n = X.shape[0]
    S0 = float(Wsym.sum())
    Z = X - X.mean(axis=0, keepdims=True)
    num = np.einsum("ij,ij->j", Z, Wsym @ Z)
    den = np.einsum("ij,ij->j", Z, Z)
    return (n / S0) * num / np.where(den == 0, np.nan, den)


def corr_cols(X, y):
    """Pearson r between every column of X and the vector y."""
    Xc = X - X.mean(axis=0, keepdims=True)
    yc = y - y.mean()
    denom = np.sqrt(np.einsum("ij,ij->j", Xc, Xc)) * np.sqrt(float(yc @ yc))
    out = np.full(X.shape[1], np.nan, dtype=np.float64)
    ok = denom > 0
    out[ok] = (Xc[:, ok].T @ yc) / denom[ok]
    return out

import scipy.stats as stats
def rank_cols(X):
    """Column-wise average ranks (ties -> mean rank); 1-D input ranked as a vector."""
    return stats.rankdata(X, axis=0).astype(np.float32)


def rank_covariates(Cov):
    """Rank every covariate column except the intercept (column 0)."""
    return np.hstack([Cov[:, :1], rank_cols(Cov[:, 1:])])


def partial_corr_cols(X, y, Cov):
    """Partial Pearson r of every column of X with y, both residualized on Cov."""
    X_res = residualize(X, Cov)
    y_res = residualize(y[:, None], Cov).ravel()
    return corr_cols(X_res, y_res)

def total_counts(adata, exclude=(WAVE_GENE,), layer=LAYER):
    """Total counts per cell, optionally excluding a gene."""
    src = adata.layers[layer] if layer in adata.layers else adata.X
    keep = ~np.isin(np.asarray(adata.var_names, dtype=str), list(exclude))
    return np.asarray(densify(src[:, keep].sum(axis=1)), dtype=np.float32).ravel()


def local_density(adata, k=10, key=SPATIAL_KEY):
    """-2 log10 of the mean distance to the k nearest cells (no graph is built or stored)."""
    from scipy.spatial import cKDTree
    coords = get_coords(adata, key)
    dist, _ = cKDTree(coords).query(coords, k=k + 1)
    return -2.0 * np.log10(dist[:, 1:].mean(axis=1) + 1e-6)


def build_covariates(adata, g=None, celltype_key=CELLTYPE_KEY, n_rounds=N_ROUNDS):
    """Design matrix (intercept + nuisance covariates); smoothed only when a graph `g` is given."""
    n = adata.n_obs

    tc = np.log10(total_counts(adata) + 1.0)
    dens = (-2.0 * np.log10(g["knn_dist"].mean(axis=1) + 1e-6)) if g is not None else local_density(adata)

    blocks = [tc[:, None], dens[:, None]]
    names = ["log_total_counts", "log_density"]

    if celltype_key is not None:
        dummies = pd.get_dummies(adata.obs[celltype_key].astype(str))
        dummies = dummies.loc[:, dummies.sum(axis=0) >= 10]
        # drop one level: the fractions sum to 1 and would be collinear with the intercept
        frac = dummies.to_numpy(dtype=np.float32)[:, 1:]
        blocks.append(frac)
        names += [f"frac_{c}" for c in dummies.columns[1:]]

    Cov = np.asarray(np.hstack(blocks), dtype=np.float32)
    if g is not None:
        Cov = np.asarray(smooth(g["P"], Cov, n_rounds), dtype=np.float32)
    Cov = (Cov - Cov.mean(0)) / np.where(Cov.std(0) == 0, 1.0, Cov.std(0))
    Cov = np.hstack([np.ones((n, 1), np.float32), Cov])
    return Cov, ["intercept"] + names

def residualize(X, Cov):
    """Residualize every column of X against the covariates in Cov."""
    beta, *_ = np.linalg.lstsq(Cov, X, rcond=None)
    return X - Cov @ beta

def run_colocalization(name, adata, g=None, wave_gene=WAVE_GENE, layer=LAYER,
                       n_rounds=N_ROUNDS, celltype_key=CELLTYPE_KEY, rng=None):
    """Per-cell correlation of every gene with `wave_gene`.
    g=None: raw (unsmoothed) expression; g=graph dict: expression smoothed over the kNN graph."""
    rng = np.random.default_rng(42) if rng is None else rng

    genes = np.asarray(adata.var_names, dtype=str)
    gidx = {gene: i for i, gene in enumerate(genes)}

    if wave_gene not in gidx:
        raise ValueError(f"Anchor gene {wave_gene} not found in adata.var_names")
    gi = gidx[wave_gene]

    raw = np.asarray(densify(adata.layers[layer]), dtype=np.float32)
    S = raw if g is None else np.asarray(smooth(g["P"], raw, n_rounds=n_rounds), dtype=np.float32)
    y = S[:, gi].copy()

    # partial pearson - density / depth / composition as covariates
    Cov, cov_names = build_covariates(adata, g, celltype_key=celltype_key, n_rounds=n_rounds)
    S_rk, y_rk, Cov_rk = rank_cols(S), rank_cols(y), rank_covariates(Cov)


    out = pd.DataFrame({
        "sample": name,
        "anchor gene": wave_gene,
        "gene": genes,
        "pearson_r": corr_cols(S, y),
        "pearson_r_partial": partial_corr_cols(S, y, Cov),
        "spearman_r": corr_cols(S_rk, y_rk),
        "spearman_r_partial": partial_corr_cols(S_rk, y_rk, Cov_rk),
    })


    return out


def orient_samples(adata, transforms, basis="spatial", sample_key="sample"):
    """Mirror / rotate individual samples around their own centre.

    transforms: {sample: {"flip": "x" | "y" | "xy" | None, "rotate": degrees}}
      flip "x" = mirror left-right, "y" = mirror up-down; rotate is counter-clockwise, applied after the flip.
    The original coordinates are kept in obsm[f"{basis}_orig"], and every call starts from them,
    so re-running the cell does not stack transformations.
    """
    orig = f"{basis}_orig"
    if orig not in adata.obsm:
        adata.obsm[orig] = adata.obsm[basis].copy()
    xy = np.asarray(adata.obsm[orig], dtype=float).copy()

    for smp, t in transforms.items():
        m = (adata.obs[sample_key] == smp).to_numpy()
        assert m.any(), f"sample {smp!r} not found in obs[{sample_key!r}]"
        c = xy[m, :2].mean(axis=0)
        p = xy[m, :2] - c
        flip = t.get("flip") or ""
        if "x" in flip:
            p[:, 0] *= -1
        if "y" in flip:
            p[:, 1] *= -1
        a = np.deg2rad(t.get("rotate", 0))
        R = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
        xy[m, :2] = p @ R.T + c

    adata.obsm[basis] = xy
    return adata