import argparse
import pandas as pd
from skimage.io import imread
import spatialproteomics as sp
import torch
import gc


def run_segmentation(input_path, output_path, marker_file, channel, gpu):
    """Processes a single image file for segmentation."""

    if gpu:
        print("Pre-clearing CUDA cache...")
        torch.cuda.empty_cache()
        gc.collect()

    # 1. Load Markers
    markers = pd.read_csv(marker_file, header=None)[0].tolist()

    # 2. Process Single Image
    print(f"Reading image: {input_path}")
    image = imread(input_path)

    print("Loading image into spatial object with markers...")
    sp_object = sp.load_image_data(image, channel_coords=markers)

    print(f"Running Cellpose on '{channel}' channel (GPU={gpu})...")
    sp_object = sp_object.tl.cellpose(channel=channel, gpu=gpu)

    print(f"Saving output to: {output_path}")
    sp_object.to_zarr(output_path)

    print("Segmentation complete for this file.")

    # 5. Clean up after yourself
    if gpu:
        torch.cuda.empty_cache()
        gc.collect()
        print("Segmentation complete. CUDA cache cleared.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Cellpose segmentation on a single CODEX image.")
    # Arguments passed directly from Snakemake
    parser.add_argument("--input", required=True, help="Path to the single input .tif image.")
    parser.add_argument("--output", required=True, help="Path to the output segmented .zarr file.")
    parser.add_argument("--markers", required=True, help="Path to the marker list CSV file.")
    parser.add_argument("--channel", required=True, help="Channel to use for segmentation.")
    parser.add_argument("--gpu", type=str, required=True, help="Use GPU (True/False).")

    args = parser.parse_args()

    # Convert the string 'True'/'False' from Snakemake to boolean
    gpu_bool = args.gpu.lower() == 'true'

    run_segmentation(
        input_path=args.input,
        output_path=args.output,
        marker_file=args.markers,
        channel=args.channel,
        gpu=gpu_bool
    )