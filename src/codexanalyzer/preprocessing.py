import xarray as xr
from scipy.signal import medfilt2d
import pandas as pd
import numpy as np
import spatialproteomics as sp

def preprocess_codex_sample(ds, sample_quantiles, patch_size=500):
    """
    Performs the core preprocessing steps for a single CODEX sample.
    
    Steps:
    1. Thresholding and median filtering.
    2. Cell size filtering (50 < area < 570).
    3. Mask expansion.
    4. Quantification (mean intensity) and Arcsinh transformation.
    
    Parameters:
    -----------
    ds : xr.Dataset
        The input segmented Zarr dataset.
    sample_quantiles : list
        List of quantile thresholds for each marker.
    patch_size : int
        Size of the patch for plotting/QC (not used for processing itself, but useful for returning a slice).
        
    Returns:
    --------
    ds_processed : xr.Dataset
        The fully processed and quantified dataset.
    df_area : pd.DataFrame
        DataFrame of cell areas before filtering (for QC histogram).
    coord_slice : str
        String representation of the center slice for plotting.
    """
    
    # --- 1. Calculate Plotting Slice (Center Crop) ---
    try:
        x_dim = ds.sizes['x']
        y_dim = ds.sizes['y']
    except KeyError:
        # Fallback
        x_dim = ds.sizes.get('x', 5000)
        y_dim = ds.sizes.get('y', 5000)

    center_x = x_dim // 2
    center_y = y_dim // 2

    x_start = max(0, center_x - patch_size // 2)
    x_end = min(x_dim, center_x + patch_size // 2)
    y_start = max(0, center_y - patch_size // 2)
    y_end = min(y_dim, center_y + patch_size // 2)
    
    coord_slice = f"{x_start}:{x_end}, {y_start}:{y_end}"

    # --- 2. Thresholding & Filtering ---
    print("Applying thresholding and median filter...")
    # Using the pp accessor from spatialproteomics
    ds_processed = ds.pp.threshold(sample_quantiles).pp.apply(medfilt2d, kernel_size=3)

    # --- 3. Cell Size Filtering ---
    ds_processed = ds_processed.pp.add_observations("area")
    df_area = ds_processed.pp.get_layer_as_df("_obs") # Keep for histogram

    print(f"Filtering cells by area (50 < area < 570)... Initial: {len(df_area)}")
    ds_processed = ds_processed.pp.filter_by_obs("area", func=lambda x: (x > 50) & (x < 570))
    # df_filtered = ds_processed.pp.get_layer_as_df("_obs") 

    # --- 4. Mask Expansion ---
    print("Expanding cell masks...")
    ds_processed = ds_processed.pp.grow_cells(iterations=2)

    # --- 5. Quantification & Transformation ---
    print("Adding quantification (intensity_mean) and arcsinh transformation...")
    ds_processed = ds_processed.pp.add_quantification(func="intensity_mean")
    ds_processed = ds_processed.pp.transform_expression_matrix(method="arcsinh")

    return ds_processed, df_area, x_start, x_end, y_start, y_end
