import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib
# Import from the refactored library module
from codexanalyzer.preprocessing import preprocess_codex_sample

# Use non-interactive backend for plots
matplotlib.use('Agg')

def main(args):
    # 1. ARGUMENT READING
    full_sample_name = args.sample
    # Clean sample name if needed
    suffix_to_remove = "_segmented"
    sample_name = full_sample_name[:-len(suffix_to_remove)] if full_sample_name.endswith(suffix_to_remove) else full_sample_name

    print(f"--- Starting Preprocessing for Sample: {sample_name} ---")

    # Ensure output directory for plots exists
    plot_dir = os.path.dirname(args.threshold_plot)
    os.makedirs(plot_dir, exist_ok=True)

    # Load config files
    markers = pd.read_csv(args.markers, header=None)[0].tolist()
    threshold_df = pd.read_csv(args.thresholds, index_col=0)

    # Retrieve quantiles
    try:
        sample_quantiles = threshold_df.loc[sample_name].tolist()
        if len(sample_quantiles) != len(markers):
            raise ValueError("Quantile list length does not match marker list length.")
    except KeyError:
        print(f"ERROR: Could not find quantiles for sample '{sample_name}' in the threshold file.")
        raise

    # 2. LOAD DATA
    print(f"Loading segmented data from: {args.zarr}")
    ds = xr.open_zarr(args.zarr)
    
    # 3. RUN CORE PREPROCESSING
    ds_processed, df_area, x_start, x_end, y_start, y_end = preprocess_codex_sample(ds, sample_quantiles)

    # 4. GENERATE PLOTS
    print(f"Generating and saving plots to {plot_dir}...")
    print(f"Plotting coordinates: x={x_start}:{x_end}, y={y_start}:{y_end}")

    # Create a spatially sliced version of ds_processed for plotting
    ds_plot_slice = ds_processed.isel(x=slice(x_start, x_end), y=slice(y_start, y_end))

    # Plot 1: Full Marker Check
    fig, ax = plt.subplots(9, 6, figsize=(24, 36))
    for i, marker in enumerate(markers):
        if i >= len(ax.flat): break
        try:
            # Use the sliced dataset for plotting
            (ds_plot_slice.pp[marker]
             .pl.colorize("gold")
             .pl.show(ax=ax.flat[i]))
            ax.flat[i].set_title(f"{marker} (q={sample_quantiles[i]})")
        except KeyError:
            print(f"SKIPPING: Marker '{marker}' not found in channels.")
            continue
    plt.tight_layout()
    fig.savefig(args.threshold_plot)
    plt.close(fig)

    # Plot 2: T-cell Visualization (if markers exist)
    try:
        fig, ax = plt.subplots(1, 3, figsize=(30, 15))
        _ = ds_plot_slice.pp["CD8"].pl.colorize().pl.show(render_segmentation=True, ax=ax[0])
        _ = ds_plot_slice.pp["CD4"].pl.colorize().pl.show(render_segmentation=True, ax=ax[1])
        _ = ds_plot_slice.pp["CD3"].pl.colorize().pl.show(render_segmentation=True, ax=ax[2])
        fig.savefig(args.tcell_plot)
        plt.close(fig)
    except KeyError as e:
         print(f"Skipping T-cell plot: missing marker {e}")

    # Plot 3: Area Histogram
    fig, ax = plt.subplots(1, 1)
    _ = ax.hist(df_area["area"], bins=100)
    ax.set_title(f"Cell Area Distribution for {sample_name}")
    fig.savefig(args.area_hist)
    plt.close(fig)

    # Plot 4: Dapi Segmentation Plot
    try:
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        _ = ds_plot_slice.pp["DAPI"].pl.colorize("gold").pl.show(render_segmentation=True, ax=ax[0])
        _ = ds_plot_slice.pp[["CD4", "CD8", "CD3", "BCL6"]].pl.show(render_segmentation=True,ax=ax[1])
        fig.savefig(args.segmentation_plot)
        plt.close(fig)
    except KeyError as e:
         print(f"Skipping Segmentation plot: missing marker DAPI")


    # 5. SAVE OUTPUT
    print(f"Saving final quantified Zarr to: {args.quant_output}")
    ds_processed.to_zarr(args.quant_output, mode='w')
    print("--- Pipeline Complete ---")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Full CODEX processing pipeline.")
    parser.add_argument("--zarr", required=True)
    parser.add_argument("--quant_output", required=True)
    parser.add_argument("--thresholds", required=True)
    parser.add_argument("--markers", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--threshold_plot", required=True)
    parser.add_argument("--tcell_plot", required=True)
    parser.add_argument("--area_hist", required=True)
    parser.add_argument("--segmentation_plot", required=True)


    args = parser.parse_args()
    main(args)