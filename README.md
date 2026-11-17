# CODEX Analyzer

A Snakemake pipeline for batch processing of CODEX imaging data: raw multi-channel
TIFF → Cellpose segmentation → per-cell marker quantification → a single
batch-corrected `AnnData` object ready for clustering and spatial analysis.

Built and used on 54-marker CODEX panels of human lymph node tissue (B-cell
lymphoma cohort), but the panel is configurable.

## Pipeline

| Step | Script | Input | Output |
|------|--------|-------|--------|
| `segmentation` | `src/scripts/segmentation_rule.py` | `.tif` (one per ROI) | `*_segmented.zarr` |
| `preprocessing` | `src/scripts/preprocessing_rule.py` | `*_segmented.zarr` | `*_quantified.zarr` + 4 QC plots |
| `aggregation` | `src/scripts/aggregation_rule.py` | all `*_quantified.zarr` | `Processed_adata_<PROJECT>.h5ad` |

What each step does:

1. **Segmentation** — loads the image with the marker list as channel coordinates
   and runs `spatialproteomics.tl.cellpose(channel=...)` (default `DAPI`).
   Optionally GPU-accelerated.
2. **Preprocessing** — per-marker quantile thresholding + 3×3 median filter, cell
   size filtering (`50 < area < 570` px), 2-iteration mask expansion, mean
   intensity quantification, `arcsinh` transform. Also writes QC figures:
   all markers, CD3/CD4/CD8 overlay, cell area histogram, DAPI + mask overlay.
3. **Aggregation** — converts every zarr to `AnnData`, applies a Seurat-style CLR
   normalization per slide, scales, concatenates, samples an equal number of cells
   per sample, then PCA (20 PCs) → BBKNN batch correction (`sample_id`) → UMAP →
   Leiden at resolutions 0.1 / 0.25 / 0.65 / 0.85.

## Requirements

Conda (or mamba). The environment also serves as the `conda:` directive for every
rule, so Snakemake must be run with `--use-conda`.

```bash
conda env create -f environment.yml
conda activate codex_analysis
snakemake --version
```

`environment.yml` pins `pulp<2.8` (Snakemake 7 dependency resolver) and pulls
Cellpose, `spatialproteomics`, `scimap`, `scanpy`, `anndata`, `xarray`.

On first run Cellpose downloads its model weights — this needs network access once.

## Input data contract

- One multi-channel `.tif` per ROI, named `<sample>.tif`, under
  `<input_dir>/<project_name>/`.
- Channel order **must** match `config/marker_list.txt` exactly, one marker per
  line. Both are read positionally: the i-th line names the i-th channel.
- Marker names are used verbatim throughout (`DAPI`, `CD3`, `BCL6`, …), so renaming
  a marker renames it in the output `AnnData` too.

## Configuration

`config/snakemake_config.yaml`:

| Key | Meaning |
|-----|---------|
| `project_name` | Sub-directory of `input_dir` holding the TIFFs. Also names every output directory. |
| `input_dir` | Root of raw images, default `data/01_Raw_ROI_clean`. |
| `marker_file` | Channel/name mapping, default `config/marker_list.txt`. |
| `thresholds_file` | Per-sample quantile table, default `config/thresholds.csv`. |
| `channel` | Channel used for segmentation, default `DAPI`. |
| `gpu` | `True`/`False`; also drives the `gpu=1` Snakemake resource. |
| `target_cells` | Total cells kept by the stratified/equal sampling step. |
| `random_seed` | Seed for sampling. |

### `config/thresholds.csv` is not tracked

Thresholds are dataset-specific (tuned per sample and per marker) and are therefore
not included in this repository. You must provide the file yourself at the path set
by `thresholds_file`, with this schema:

```csv
,CD3,CD4,CD8,...
<SAMPLE_NAME>,0.90,0.85,0.80,...
```

- First column: sample name = TIFF basename without extension
  (`SAMPLE_01_ROI_01`, i.e. whatever you named the `.tif`).
- Remaining columns: marker names, in the same order as `marker_list.txt`.
- Values: quantiles in `[0, 1)` used as `ds.pp.threshold(...)` cut-offs. `0` means
  "no threshold" (typically used for `DAPI`).

The pipeline fails loudly if a sample is missing from the file or if the number of
quantiles does not match the number of markers.

## Running

```bash
snakemake -n                                  # dry run
snakemake --cores 4 --use-conda                # full pipeline
snakemake --cores 4 --use-conda --resources gpu=1   # if gpu: True
snakemake --cores 4 --use-conda --ts           # timestamps in the console
```

The default target (`rule all`) is
`data/04_Analysis/<PROJECT>/Processed_adata_<PROJECT>.h5ad`. Individual rules can be
run by naming their output, e.g.

```bash
snakemake --cores 4 --use-conda --until preprocessing
```

**CUDA out-of-memory during segmentation**: a single ROI is too large for the GPU.
Set `gpu: False` in the config and rerun; segmentation is then CPU-only.

## Output

Layout (all under the project root, `project_name` substituted for `<PROJECT>`):

```
data/01_Raw_ROI_clean/<PROJECT>/<sample>.tif      # input
data/02_Segmented/<PROJECT>/<sample>_segmented.zarr
data/03_Processed/<PROJECT>/<sample>_quantified.zarr
output/plots/<PROJECT>/Preprocessing/<sample>_{threshold,tcell,area_hist,segmentation}.png
data/04_Analysis/<PROJECT>/Processed_adata_<PROJECT>.h5ad
logs/<PROJECT>/segmentation|preprocessing/<sample>.log
logs/<PROJECT>/aggregation.log
```

The final `AnnData` has one row per cell, one column per marker, and:

| Slot | Content |
|------|---------|
| `X` | CLR-normalized, scaled mean intensities |
| `layers["percentage_positive"]` | fraction of positive pixels per cell/marker, when present in the source zarr |
| `obs` | `source_file`, `sample_id` (filename prefix before the first `_`), `batch`, per-cell metadata carried from the zarr (centroids, area), `leiden_0.1/0.25/0.65/0.85` |
| `obsm` | `X_pca`, `X_umap`, `spatial` |
| `uns` | `source_path`, `leiden_*` cluster colors |

Cell-type annotation is **not** automated here: the Leiden clusters
(`leiden_0.65` in practice) have to be labelled from a marker heatmap before any
phenotype-level analysis.

## Repository layout

```
Snakefile                     # 3 rules, config loaded from config/
config/                       # pipeline config, marker/channel list
src/codexanalyzer/            # library: data loading, preprocessing, normalization, plotting helpers
src/scripts/                  # argparse CLIs called by the Snakefile rules
environment.yml               # conda env, also used by the conda: directives
```

`src/` is not an installed package: the rules export `PYTHONPATH=src` and import
`codexanalyzer.*`. Run Snakemake from the repository root.

## What is deliberately not in this repository

- **`notebooks/`** — the exploratory analysis (cluster annotation, `scimap` /
  `squidpy` spatial statistics, figure generation) lives in per-dataset notebooks
  and is not published, as it is tied to identifiable cohort samples.
- **Raw and intermediate data** — `data/`, `output/`, `logs/`, `*.h5ad`,
  `*.zarr/`, `*.tif` are gitignored, as are large `AnnData` exports.
- **`config/thresholds.csv`** — see above.

## License

MIT — see `LICENSE`. Applies to the code only; no data is distributed.

## Contact

Paul Marcoux — Inserm / CRCT, Toulouse — `paul.marcoux@inserm.fr`
