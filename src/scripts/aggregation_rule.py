import argparse
import anndata as ad
import scanpy as sc
from codexanalyzer.data import load_zarr_to_adatas
from codexanalyzer.process import normalize_single_codex_slide, equal_sampling
import pandas as pd
import spatialproteomics as sp


def main(args):
    """Run the processing and unsupervised analysis pipeline."""
    # Load marker list
    # Assuming it's a file path if it exists
    try:
        marker_list = pd.read_csv(args.marker_list, header=None)[0].tolist()
    except Exception:
        # Fallback if it's a direct list string (unlikely for file path arg)
        marker_list = args.marker_list.split(',')

    adatas, _ = load_zarr_to_adatas(args.input_files, marker_list, verbose=True)

    # -- Parameters
    sc.settings.verbosity = 1  # verbosity level
    sc.settings.set_figure_params(dpi=80, facecolor='white')

    # 1. Load and normalize all AnnData objects
    normalized_adatas = []
    for i, adata in enumerate(adatas):
        print(f"Processing slide {i+1}/{len(adatas)}")
        adata_norm = normalize_single_codex_slide(adata.copy(), i, True)
        normalized_adatas.append(adata_norm)

    # 2. Merge normalized datasets
    print("Merging normalized datasets...")
    adata_merged = ad.concat(normalized_adatas, join="outer", label="batch", index_unique="-")

    # Calculate number of cells to sample per dataset
    n_samples = len(adata_merged.obs['sample_id'].unique())
    # Ensure target_cells is an integer
    target_cells = int(args.target_cells)
    cells_per_sample = target_cells // n_samples

    print(f"Sampling {cells_per_sample} cells from each of the {n_samples} samples...")
    print(f"Total target cells: {cells_per_sample * n_samples}")

    # 3. Downsampling (uniform random sampling)
    adata_merged = equal_sampling(adata_merged,
                                  cells_per_sample=cells_per_sample,
                                  sample_key='sample_id',
                                  random_state=args.random_seed)

    print("\nCells per sample after downsampling:")
    print(adata_merged.obs['sample_id'].value_counts())

    print(f"Final merged dataset: {adata_merged.n_obs} cells, {adata_merged.n_vars} markers")
    print(f"Batches: {adata_merged.obs['batch'].value_counts().to_dict()}")

    # 4. Principal Component Analysis (equivalent to RunPCA)
    print("Running PCA...")
    sc.tl.pca(adata_merged, n_comps=20, svd_solver='arpack')

    #Batch correction BBKNN
    print("Running Batch Correction ...")
    sc.external.pp.bbknn(
        adata_merged,
        batch_key='sample_id',
        neighbors_within_batch=3,
        computation="annoy",
    )

    # 5. UMAP embedding (equivalent to RunUMAP)
    print("Running UMAP...")
    sc.tl.umap(adata_merged)

    # 6. Find clusters (equivalent to FindNeighbors + FindClusters)
    print("Finding clusters...")
    sc.tl.leiden(adata_merged, resolution=0.25, random_state=1, flavor='igraph', n_iterations=2, directed=False, key_added='leiden_0.25')
    sc.tl.leiden(adata_merged, resolution=0.65, random_state=1, flavor='igraph', n_iterations=2, directed=False,key_added='leiden_0.65')
    sc.tl.leiden(adata_merged, resolution=0.85, random_state=1, flavor='igraph', n_iterations=2, directed=False,key_added='leiden_0.85')
    sc.tl.leiden(adata_merged, resolution=0.1, random_state=1, flavor='igraph', n_iterations=2, directed=False,key_added='leiden_0.1')

    # Save the processed AnnData object
    print(f"Saving processed AnnData object to {args.output_file}...")
    adata_merged.write_h5ad(args.output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate and process CODEX data.")
    
    parser.add_argument("--input_files", nargs='+', required=True, help="List of input quantified zarr files")
    parser.add_argument("--output_file", required=True, help="Output h5ad file path")
    parser.add_argument("--marker_list", required=True, help="Path to marker list file")
    parser.add_argument("--target_cells", type=int, default=600000, help="Target total cells for downsampling")
    parser.add_argument("--random_seed", type=int, default=42, help="Random seed")

    args = parser.parse_args()
    main(args)