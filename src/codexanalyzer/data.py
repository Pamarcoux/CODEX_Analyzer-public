import os
import glob
import xarray as xr
import anndata as ad
import spatialproteomics as sp


def load_zarr_to_adatas(input_source, marker_list, verbose=True):
    import spatialproteomics as sp
    """
    Load all .zarr files from a folder or a list of paths and convert them to AnnData objects

    Parameters:
    -----------
    input_source : str or list
        Path to folder containing .zarr files, or a list of file paths
    marker_list : list
        List of markers to select from each dataset
    verbose : bool
        Whether to print progress information

    Returns:
    --------
    adatas : list
        List of AnnData objects, one per .zarr file
    zarr_files : list
        List of .zarr file paths (for reference)
    """

    # Determine zarr files based on input type
    if isinstance(input_source, str):
        # Find all .zarr files in the folder
        zarr_pattern = os.path.join(input_source, "*.zarr")
        zarr_files = glob.glob(zarr_pattern)
        if not zarr_files:
            raise ValueError(f"No .zarr files found in {input_source}")
    elif isinstance(input_source, list):
        zarr_files = input_source
    else:
        raise ValueError("input_source must be a folder path (str) or a list of file paths")

    if verbose:
        print(f"Found {len(zarr_files)} .zarr files:")
        for f in zarr_files:
            print(f"  - {os.path.basename(f)}")

    adatas = []
    failed_files = []

    for i, zarr_path in enumerate(zarr_files):
        try:
            if verbose:
                print(f"\nProcessing {i + 1}/{len(zarr_files)}: {os.path.basename(zarr_path)}")

            # Load the zarr dataset
            ds = xr.open_zarr(zarr_path)

            # ds est ton xarray.Dataset chargé
            available_channels = list(ds['channels'].values)

            missing_markers = [m for m in marker_list if m not in available_channels]
            if missing_markers:
                print("Markers missing from dataset channels:", missing_markers)
            else:
                print("All markers are present in dataset.")

            # Select the markers using .pp
            ds_filtered = ds.pp[marker_list]

            # Préparer les dicts conditionnellement
            additional_layers = {}
            if "_percentage_positive" in ds_filtered.data_vars:
                additional_layers["percentage_positive"] = "_percentage_positive"

            additional_uns = {}
            if "_la_properties" in ds_filtered.data_vars:
                additional_uns["label_colors"] = "_la_properties"

            # Convert to AnnData
            adata = ds_filtered.tl.convert_to_anndata(
                expression_matrix_key="_intensity",
                additional_layers=additional_layers,
                additional_uns=additional_uns,
            )

            # Add metadata about the source file
            file_name = os.path.basename(zarr_path)
            adata.obs["source_file"] = file_name

            # Extract sample_id as the first part of filename (before first underscore)
            sample_id = file_name.split('_')[0] if '_' in file_name else file_name.replace('.zarr', '')
            adata.obs["sample_id"] = sample_id

            # Store the original file path in uns
            adata.uns["source_path"] = zarr_path

            adatas.append(adata)

            if verbose:
                print(f"  Successfully loaded: {adata.n_obs} cells, {adata.n_vars} markers")

        except Exception as e:
            print(f"  Error processing {os.path.basename(zarr_path)}: {str(e)}")
            failed_files.append(zarr_path)
            continue

    if verbose:
        print(f"\nSummary:")
        print(f"  Successfully loaded: {len(adatas)} files")
        print(f"  Failed: {len(failed_files)} files")
        if failed_files:
            print("  Failed files:")
            for f in failed_files:
                print(f"    - {os.path.basename(f)}")

    return adatas, zarr_files
