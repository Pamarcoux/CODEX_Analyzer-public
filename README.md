# CODEX Analyzer

Snakemake workflow for batch processing of CODEX imaging data: raw multi-channel
TIFF → Cellpose segmentation → per-cell marker quantification → one batch-corrected
`AnnData` object ready for clustering and spatial analysis.

**This is a pipeline, not a package.** There is nothing to `pip install` and no
public Python API: the repo glues together existing community tools
(spatialproteomics, Cellpose, anndata/scanpy, scimap) behind three Snakemake rules.
The `codexanalyzer` module is internal plumbing, imported via `PYTHONPATH=src`.

Developed on 54-marker CODEX panels of human lymph node tissue (B-cell lymphoma
cohort); the panel is configurable.

## Pipeline

| Rule | Script | Input → output |
|------|--------|----------------|
| `segmentation` | `src/scripts/segmentation_rule.py` | `.tif` → `*_segmented.zarr` |
| `preprocessing` | `src/scripts/preprocessing_rule.py` | `*_segmented.zarr` → `*_quantified.zarr` + QC plots |
| `aggregation` | `src/scripts/aggregation_rule.py` | all `*_quantified.zarr` → `Processed_adata_<PROJECT>.h5ad` |

- **Segmentation** — `spatialproteomics.tl.cellpose()` on one channel (`DAPI` by
  default), optionally on GPU.
- **Preprocessing** — quantile thresholding + 3×3 median filter, cell size filter
  (`50 < area < 570` px), 2-iteration mask dilation, mean-intensity
  quantification, `arcsinh`.
- **Aggregation** — CLR normalization per slide, scaling, concatenation, equal
  sampling per sample, PCA → BBKNN → UMAP → Leiden (0.1/0.25/0.65/0.85).

## Install

```bash
conda env create -f environment.yml
conda activate codex_analysis
```

The env doubles as the `conda:` directive of every rule, so always run Snakemake
with `--use-conda`. Cellpose downloads its model weights on first use.

## Input

- One multi-channel `.tif` per ROI, named `<sample>.tif`, under
  `<input_dir>/<project_name>/`.
- Channel order **must** match `config/marker_list.txt` (one marker per line):
  both are read positionally.
- Marker names are used verbatim downstream.

## Config

`config/snakemake_config.yaml`: `project_name` (input sub-directory, also names
every output directory), `input_dir`, `marker_file`, `thresholds_file`, `channel`
(segmentation channel), `gpu`, `target_cells` (cells kept after sampling),
`random_seed`.

### `config/thresholds.csv` is not tracked

Thresholds are tuned per dataset and are not distributed. Provide the file at
`thresholds_file`:

```csv
,CD3,CD4,CD8,...
SAMPLE_01_ROI_01,0.90,0.85,0.80,...
```

Sample name = TIFF basename without extension; then one column per marker, in
`marker_list.txt` order; values are quantiles in `[0, 1)` (`0` = no threshold,
typically for `DAPI`). The pipeline stops with an explicit error on a missing
sample or a marker/quantile count mismatch.

## Run

```bash
snakemake -n                                    # dry run
snakemake --cores 4 --use-conda                 # default target: the aggregated h5ad
snakemake --cores 4 --use-conda --resources gpu=1
snakemake --cores 4 --use-conda --until preprocessing
```

CUDA OOM during segmentation means one ROI is too large for the GPU: set
`gpu: False` and rerun.

## Output

```
data/01_Raw_ROI_clean/<PROJECT>/<sample>.tif
data/02_Segmented/<PROJECT>/<sample>_segmented.zarr
data/03_Processed/<PROJECT>/<sample>_quantified.zarr
output/plots/<PROJECT>/Preprocessing/<sample>_{threshold,tcell,area_hist,segmentation}.png
data/04_Analysis/<PROJECT>/Processed_adata_<PROJECT>.h5ad
logs/<PROJECT>/{segmentation,preprocessing}/<sample>.log
```

The final `AnnData`: one row per cell, one column per marker.

| Slot | Content |
|------|---------|
| `X` | CLR-normalized, scaled mean intensities |
| `layers["percentage_positive"]` | fraction of positive pixels, when present in the source zarr |
| `obs` | `source_file`, `sample_id`, `batch`, centroids, `leiden_0.1/0.25/0.65/0.85` |
| `obsm` | `X_pca`, `X_umap`, `spatial` |

Cell types are **not** assigned automatically: clusters have to be labelled from a
marker heatmap before phenotype-level analysis.

## Repository layout

```
Snakefile            # 3 rules, config loaded from config/
config/              # pipeline config, marker/channel list
src/codexanalyzer/   # internal library: loading, preprocessing, normalization, plots
src/scripts/         # argparse CLIs called by the rules
environment.yml      # conda env, also used by the conda: directives
```

## Not included

`notebooks/` (per-dataset analysis: cluster annotation, scimap/squidpy spatial
statistics, figures) is not published, as it is tied to identifiable cohort
samples. Data, `output/`, `logs/`, `*.h5ad`, `*.zarr/` and `config/thresholds.csv`
are gitignored.

## Citations

The pipeline rests on spatialproteomics:

> Meyer-Bender M., Vöhringer H., Schniederjohann C., et al. Spatialproteomics: an
> interoperable toolbox for analyzing highly multiplexed fluorescence image data.
> *Nature Methods* 23, 1564–1573 (2026). doi:10.1038/s41592-026-03155-1

and was used in:

> Marcoux P., et al. *Blood* (2026).

## License

MIT — see `LICENSE`. Applies to the code only; no data is distributed.

## Contact

Paul Marcoux — Inserm / CRCT, Toulouse — `paul.marcoux@inserm.fr`
