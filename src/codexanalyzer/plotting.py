import matplotlib.pyplot as plt
import scanpy as sc
import scimap as sm
import pandas as pd
import numpy as np

def plot_umap(adata, color, title, **kwargs):
    """Generate a UMAP plot."""
    sc.pl.umap(adata, color=color, title=title, **kwargs)

def plot_marker_dotplot(adata, n_genes, **kwargs):
    """Generate a dotplot of marker genes."""
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=n_genes, **kwargs)

def plot_marker_heatmap(adata, marker_list, **kwargs):
    """Generate a heatmap of marker genes."""
    sc.pl.matrixplot(adata, var_names=marker_list, groupby='leiden', **kwargs)

def plot_umaps_by_sample(adata, sample_key='sample_name', **kwargs):
    """Generate UMAP plots for each sample."""
    sample_ids = adata.obs[sample_key].unique()
    n_samples = len(sample_ids)
    fig, axes = plt.subplots(nrows=(n_samples + 2) // 3, ncols=3,
                             figsize=(15, 5 * ((n_samples + 2) // 3)))
    axes = axes.flatten()

    for idx, sample in enumerate(sample_ids):
        mask = adata.obs[sample_key] == sample
        sc.pl.umap(adata[mask], ax=axes[idx], show=False, title=f'Sample: {sample}', **kwargs)

    for idx in range(n_samples, len(axes)):
        fig.delaxes(axes[idx])

    plt.tight_layout()
    plt.show()

def plot_spatial_by_sample(adata, color, sample_key='sample_name', **kwargs):
    """Generate spatial embedding plots for each sample."""
    sample_ids = adata.obs[sample_key].unique()
    n_samples = len(sample_ids)
    fig, axes = plt.subplots(nrows=(n_samples + 2) // 3, ncols=3,
                             figsize=(15, 5 * ((n_samples + 2) // 3)))
    axes = axes.flatten()

    for idx, sample in enumerate(sample_ids):
        mask = adata.obs[sample_key] == sample
        sc.pl.embedding(adata[mask], basis='spatial', color=color, ax=axes[idx],
                        show=False, title=f'Sample: {sample}', **kwargs)

    for idx in range(n_samples, len(axes)):
        fig.delaxes(axes[idx])

    plt.tight_layout()
    plt.show()

def plot_spatial_distance(adata, **kwargs):
    """Plot spatial distance using scimap."""
    sm.pl.spatial_distance(adata, **kwargs)

def plot_spatial_interaction(adata, **kwargs):
    """Plot spatial interaction using scimap."""
    sm.pl.spatial_interaction(adata, **kwargs)

def plot_stacked_barplot(adata, x_axis, y_axis, **kwargs):
    """Plot stacked barplot using scimap."""
    sm.pl.stacked_barplot(adata, x_axis=x_axis, y_axis=y_axis, **kwargs)