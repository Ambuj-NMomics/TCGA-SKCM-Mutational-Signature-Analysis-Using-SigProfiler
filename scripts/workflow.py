"""
TCGA-SKCM Mutational Signature Analysis Pipeline
=================================================
Author  : Ambuj Narayan Maurya
Project : TCGA Melanoma (SKCM) — COSMIC Signature Fitting
GitHub  : https://github.com/Ambuj-NMomics

Pipeline Steps
--------------
1. Install reference genome (GRCh38)
2. Generate SBS96 mutational matrix from MAF files
3. COSMIC signature fitting using SigProfilerAssignment
4. Exposure heatmap visualization
5. Cosine similarity analysis between samples

Requirements
------------
    pip install SigProfilerMatrixGenerator SigProfilerAssignment SigProfilerExtractor scikit-learn matplotlib pandas

Usage
-----
    python scripts.py

Note: Update PATH_CONFIG at the top with your actual file paths before running.
"""

# ──────────────────────────────────────────────
# IMPORTS
# ──────────────────────────────────────────────

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

# ──────────────────────────────────────────────
# PATH CONFIGURATION  ← edit these before running
# ──────────────────────────────────────────────

PATH_CONFIG = {
    "maf_dir"        : "/home/ambuj/Desktop/cancer/maf_files",
    "output_dir"     : "/home/ambuj/Desktop/cancer/maf_files/output",
    "assignment_dir" : "/home/ambuj/Desktop/cancer/SKCMAssignment",
    "figures_dir"    : "/home/ambuj/Desktop/cancer/figures",
}

# Auto-derived paths (don't touch these)
SBS96_MATRIX   = os.path.join(PATH_CONFIG["output_dir"], "SBS", "SKCM.SBS96.all")
ACTIVITIES_TXT = os.path.join(
    PATH_CONFIG["assignment_dir"],
    "Assignment_Solution", "Activities", "Assignment_Solution_Activities.txt"
)

# ──────────────────────────────────────────────
# HELPER UTILITIES
# ──────────────────────────────────────────────

def make_dir(path):
    """Create directory if it doesn't already exist."""
    os.makedirs(path, exist_ok=True)
    print(f"  [OK] Directory ready: {path}")


def check_file(path, label="File"):
    """Exit early with a clear message if a required file is missing."""
    if not os.path.exists(path):
        print(f"\n  [ERROR] {label} not found: {path}")
        print("  Make sure previous steps completed successfully.\n")
        sys.exit(1)
    print(f"  [OK] Found: {path}")


def section(title):
    """Print a section header to keep terminal output readable."""
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")


# ──────────────────────────────────────────────
# STEP 1 — INSTALL REFERENCE GENOME
# ──────────────────────────────────────────────

def install_genome(genome="GRCh38"):
    """
    Download and install the reference genome for SigProfilerMatrixGenerator.
    Only needs to run once — safe to call repeatedly (checks if already installed).
    """
    section("STEP 1 — Installing Reference Genome")

    try:
        from SigProfilerMatrixGenerator import install as genInstall
        print(f"  Installing {genome} ... (skip if already done, this takes a few minutes)")
        genInstall.install(genome)
        print(f"  [OK] {genome} ready.")
    except Exception as e:
        print(f"  [ERROR] Genome install failed: {e}")
        sys.exit(1)


# ──────────────────────────────────────────────
# STEP 2 — GENERATE SBS96 MUTATIONAL MATRIX
# ──────────────────────────────────────────────

def generate_matrix(project="SKCM", genome="GRCh38"):
    """
    Parse all MAF files in maf_dir and generate the SBS96 mutational matrix.

    Output lands in:
        maf_dir/output/SBS/SKCM.SBS96.all

    Also generates DBS and indel matrices as a bonus — ignore them for now
    unless you want to extend the analysis later.
    """
    section("STEP 2 — Generating SBS96 Matrix")

    check_file(PATH_CONFIG["maf_dir"], label="MAF directory")

    try:
        from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

        print(f"  Project   : {project}")
        print(f"  Genome    : {genome}")
        print(f"  MAF dir   : {PATH_CONFIG['maf_dir']}")
        print("  Running matrix generator ...\n")

        matGen.SigProfilerMatrixGeneratorFunc(
            project,
            genome,
            PATH_CONFIG["maf_dir"],
            plot=True        # generates plots alongside the matrix
        )

        print("\n  [OK] Matrix generation complete.")
        print(f"  Output: {SBS96_MATRIX}")

    except Exception as e:
        print(f"  [ERROR] Matrix generation failed: {e}")
        sys.exit(1)


# ──────────────────────────────────────────────
# STEP 3 — COSMIC SIGNATURE FITTING
# ──────────────────────────────────────────────

def run_cosmic_fitting(genome="GRCh38", cosmic_version=3.4):
    """
    Decompose the SBS96 matrix into COSMIC v3.4 signatures using NNLS.

    Why NNLS? Mutation counts can't go negative — so the solver is
    constrained accordingly. It reconstructs your observed profile
    as a weighted combination of known COSMIC signatures.

    Output lands in:
        assignment_dir/Assignment_Solution/Activities/
    """
    section("STEP 3 — COSMIC Signature Fitting")

    check_file(SBS96_MATRIX, label="SBS96 matrix")
    make_dir(PATH_CONFIG["assignment_dir"])

    try:
        from SigProfilerAssignment import Analyzer as Analyze

        print(f"  Loading matrix from: {SBS96_MATRIX}")
        matrix = pd.read_csv(SBS96_MATRIX, sep="\t", index_col=0)

        print(f"  Matrix shape: {matrix.shape[0]} contexts × {matrix.shape[1]} samples")
        print(f"  COSMIC version: {cosmic_version}")
        print("  Running assignment ...\n")

        Analyze.cosmic_fit(
            samples=matrix,
            output=PATH_CONFIG["assignment_dir"],
            input_type="matrix",
            context_type="96",
            genome_build=genome,
            cosmic_version=cosmic_version
        )

        print("\n  [OK] Signature fitting complete.")
        print(f"  Results: {PATH_CONFIG['assignment_dir']}")

    except Exception as e:
        print(f"  [ERROR] COSMIC fitting failed: {e}")
        sys.exit(1)


# ──────────────────────────────────────────────
# STEP 4 — EXPOSURE HEATMAP
# ──────────────────────────────────────────────

def plot_exposure_heatmap():
    """
    Visualize per-sample signature exposures as a heatmap.

    Rows  = COSMIC signatures
    Cols  = tumor samples
    Color = mutation count attributed to each signature

    Saved to: figures_dir/exposure_heatmap.png
    """
    section("STEP 4 — Exposure Heatmap")

    check_file(ACTIVITIES_TXT, label="Assignment activities file")
    make_dir(PATH_CONFIG["figures_dir"])

    df = pd.read_csv(ACTIVITIES_TXT, sep="\t", index_col=0)

    print(f"  Signatures detected : {df.shape[0]}")
    print(f"  Samples             : {df.shape[1]}")

    # Filter out signatures with zero contribution across all samples
    df = df.loc[df.sum(axis=1) > 0]
    print(f"  Active signatures   : {df.shape[0]}")

    fig, ax = plt.subplots(figsize=(18, 6))

    im = ax.imshow(
        df.values,
        aspect="auto",
        cmap="YlOrRd",
        norm=mcolors.PowerNorm(gamma=0.4)   # gamma < 1 pulls out low-exposure signals
    )

    plt.colorbar(im, ax=ax, label="Mutation Count (Exposure)", pad=0.01)

    ax.set_xticks(range(df.shape[1]))
    ax.set_xticklabels(df.columns, rotation=90, fontsize=5)
    ax.set_yticks(range(df.shape[0]))
    ax.set_yticklabels(df.index, fontsize=8)

    ax.set_xlabel("Tumor Samples", fontsize=11)
    ax.set_ylabel("COSMIC Signatures", fontsize=11)
    ax.set_title("Mutational Signature Exposures — TCGA-SKCM", fontsize=13, fontweight="bold")

    plt.tight_layout()

    out_path = os.path.join(PATH_CONFIG["figures_dir"], "exposure_heatmap.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"\n  [OK] Heatmap saved: {out_path}")


# ──────────────────────────────────────────────
# STEP 5 — COSINE SIMILARITY ANALYSIS
# ──────────────────────────────────────────────

def run_cosine_similarity():
    """
    Compute pairwise cosine similarity between tumor mutation profiles.

    Values:
        > 0.85  — nearly identical mutation spectra (strong shared process)
        0.5–0.8 — moderate overlap
        < 0.3   — distinct profiles, likely different secondary drivers

    Saves:
        figures_dir/cosine_similarity_heatmap.png
        figures_dir/cosine_similarity_matrix.csv
    """
    section("STEP 5 — Cosine Similarity Analysis")

    check_file(SBS96_MATRIX, label="SBS96 matrix")
    make_dir(PATH_CONFIG["figures_dir"])

    df = pd.read_csv(SBS96_MATRIX, sep="\t", index_col=0)
    samples = df.columns.tolist()

    print(f"  Samples: {len(samples)}")
    print("  Computing pairwise cosine similarity ...\n")

    sim_matrix = cosine_similarity(df.T)
    sim_df = pd.DataFrame(sim_matrix, index=samples, columns=samples)

    # Save CSV
    csv_path = os.path.join(PATH_CONFIG["figures_dir"], "cosine_similarity_matrix.csv")
    sim_df.to_csv(csv_path)
    print(f"  [OK] Similarity matrix saved: {csv_path}")

    # Quick stats
    upper = sim_matrix[np.triu_indices_from(sim_matrix, k=1)]
    print(f"\n  Similarity stats (upper triangle):")
    print(f"    Mean   : {upper.mean():.3f}")
    print(f"    Median : {np.median(upper):.3f}")
    print(f"    Min    : {upper.min():.3f}")
    print(f"    Max    : {upper.max():.3f}")

    # Plot heatmap
    n = len(samples)
    fig_size = max(10, n * 0.15)   # scale with sample count
    fig, ax = plt.subplots(figsize=(fig_size, fig_size * 0.8))

    im = ax.imshow(sim_matrix, cmap="coolwarm", vmin=0, vmax=1)
    plt.colorbar(im, ax=ax, label="Cosine Similarity", shrink=0.6)

    ax.set_xticks(range(n))
    ax.set_xticklabels(samples, rotation=90, fontsize=4)
    ax.set_yticks(range(n))
    ax.set_yticklabels(samples, fontsize=4)
    ax.set_title("Pairwise Cosine Similarity — TCGA-SKCM Mutation Profiles",
                 fontsize=12, fontweight="bold")

    plt.tight_layout()

    fig_path = os.path.join(PATH_CONFIG["figures_dir"], "cosine_similarity_heatmap.png")
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()

    print(f"  [OK] Similarity heatmap saved: {fig_path}")


# ──────────────────────────────────────────────
# MAIN — RUN FULL PIPELINE
# ──────────────────────────────────────────────

def main():
    print("\n" + "="*60)
    print("  TCGA-SKCM Mutational Signature Analysis Pipeline")
    print("  Author: Ambuj Narayan Maurya")
    print("="*60)

    # Comment out any step you've already run
    install_genome()          # Step 1 — only needed once
    generate_matrix()         # Step 2 — parses MAF files → SBS96 matrix
    run_cosmic_fitting()      # Step 3 — COSMIC signature decomposition
    plot_exposure_heatmap()   # Step 4 — heatmap visualization
    run_cosine_similarity()   # Step 5 — sample similarity analysis

    section("PIPELINE COMPLETE")
    print("  Figures  →", PATH_CONFIG["figures_dir"])
    print("  Results  →", PATH_CONFIG["assignment_dir"])
    print()


if __name__ == "__main__":
    main()
