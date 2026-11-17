import numpy as np
import scanpy as sc

def clr_normalization(x):
    """
    Perform Centered Log-Ratio normalization on a vector or array.

    This is a Python implementation of the CLR function from Seurat:
    clr_function <- function(x) {
      return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
    }

    Parameters:
    -----------
    x : array-like
        Input data to normalize

    Returns:
    --------
    array-like
        CLR normalized data
    """
    x = np.asarray(x)

    # Create mask for positive values
    positive_mask = x > 0

    # Calculate the geometric mean of positive values
    # In log space: exp(sum(log(1+x)) / n)
    if np.any(positive_mask):
        log1p_sum = np.sum(np.log1p(x[positive_mask]))
        n = len(x)  # Use full length including zeros as in the R function
        geometric_mean = np.exp(log1p_sum / n)
    else:
        # Handle the case where there are no positive values
        geometric_mean = 1

    # Divide by geometric mean and apply log1p
    return np.log1p(x / geometric_mean)

def normalize_single_codex_slide(adata, sample_id, verbose=True):
    """
    Normalize a single CODEX slide before merging
    """
    if verbose:
        print(f"Normalizing slide {sample_id}...")

    # Add slide identifier
    adata.obs["batch"] = sample_id

    # 2. Filter Minimum count
    sc.pp.filter_cells(adata, min_counts=1)

    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)

    # Apply CLR normalization to the expression matrix
    X_clr = np.zeros_like(adata.X)     # Create a copy of the raw data

    # Apply CLR normalization per cell (row-wise)
    for i in range(adata.n_obs):
        X_clr[i, :] = clr_normalization(adata.X[i, :])

    # Replace the data with CLR-normalized values
    adata.X = X_clr

    # 3. Scale the data (important for CODEX due to different antibody intensities)
    sc.pp.scale(adata, max_value=10)

    # 4. Set all features as highly variable (since CODEX panels are targeted)
    adata.var['highly_variable'] = True
    adata.var['highly_variable_rank'] = range(len(adata.var))

    return adata

def stratified_downsample(adata, target_cells, batch_key='batch', random_state=42):
    """
    Perform stratified downsampling to maintain batch proportions
    """
    np.random.seed(random_state)

    # Calculate proportions for each batch
    batch_props = adata.obs[batch_key].value_counts(normalize=True)

    # Calculate number of cells to sample from each batch
    batch_numbers = (batch_props * target_cells).round().astype(int)

    # Adjust if total doesn't match target (due to rounding)
    diff = target_cells - batch_numbers.sum()
    if diff != 0:
        batch_numbers[0] += diff

    # Sample from each batch
    keep_idx = []
    for batch in batch_props.index:
        batch_idx = np.where(adata.obs[batch_key] == batch)[0]
        n_sample = min(len(batch_idx), batch_numbers[batch])
        sampled_idx = np.random.choice(batch_idx, size=n_sample, replace=False)
        keep_idx.extend(sampled_idx)

    return adata[keep_idx].copy()

def equal_sampling(adata, cells_per_sample, sample_key='sample_id', random_state=42):
    """
    Sample equal number of cells from each sample

    Parameters:
    -----------
    adata : AnnData object
        The annotated data matrix
    cells_per_sample : int
        Number of cells to sample from each sample
    sample_key : str
        Key in adata.obs containing sample identifiers
    random_state : int
        Random seed for reproducibility

    Returns:
    --------
    AnnData object with equal numbers of cells from each sample
    """
    np.random.seed(random_state)

    # Get unique samples
    samples = adata.obs[sample_key].unique()

    # Sample equal numbers from each sample
    keep_idx = []
    for sample in samples:
        sample_mask = adata.obs[sample_key] == sample
        sample_idx = np.where(sample_mask)[0]

        # Check if we have enough cells in this sample
        n_sample = min(len(sample_idx), cells_per_sample)
        if n_sample < cells_per_sample:
            print(f"Warning: Sample {sample} has only {len(sample_idx)} cells, "
                  f"sampling {n_sample} instead of {cells_per_sample}")

        # Sample cells
        sampled_idx = np.random.choice(sample_idx, size=n_sample, replace=False)
        keep_idx.extend(sampled_idx)

    return adata[keep_idx].copy()
