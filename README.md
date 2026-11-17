# CODEX Analyzer

Snakemake workflow for batch processing of CODEX imaging data. It takes raw
multi-channel TIFF files, runs Cellpose segmentation, quantifies markers per cell
and writes one batch-corrected `AnnData` object for clustering and spatial
analysis.

The repository is a workflow, not an installable package. Clone it and run
Snakemake; there is no pip package and no public API. The scientific work is done
by spatialproteomics, Cellpose, anndata/scanpy and scimap, wrapped by three rules.
The `codexanalyzer` module is internal code loaded through `PYTHONPATH=src`.

It was developed on 54-marker CODEX panels of human lymph node tissue (B-cell
lymphoma cohort). The panel is configurable.

## Pipeline

| Rule | Script | Input, output |
|------|--------|---------------|
| `segmentation` | `src/scripts/segmentation_rule.py` | `.tif` to `*_segmented.zarr` |
| `preprocessing` | `src/scripts/preprocessing_rule.py` | `*_segmented.zarr` to `*_quantified.zarr` and QC plots |
| `aggregation` | `src/scripts/aggregation_rule.py` | all `*_quantified.zarr` to `Processed_adata_<PROJECT>.h5ad` |

Segmentation calls `spatialproteomics.tl.cellpose()` on one channel, `DAPI` by
default, with optional GPU use.

Preprocessing applies quantile thresholding and a 3x3 median filter, removes cells
outside 50 to 570 px of area, dilates masks by 2 iterations, quantifies mean
intensity and applies `arcsinh`.

Aggregation normalizes each slide with CLR, scales, concatenates, samples the same
number of cells per sample, then runs PCA, BBKNN, UMAP and Leiden clustering at
resolutions 0.1, 0.25, 0.65 and 0.85.

## Install

```bash
conda env create -f environment.yml
conda activate codex_analysis
```

The environment is also the `conda:` directive of each rule, so run Snakemake with
`--use-conda`. Cellpose downloads its model weights on first use.

## Input

- One multi-channel `.tif` per ROI, named `<sample>.tif`, under
  `<input_dir>/<project_name>/`.
- Channel order must match `config/marker_list.txt`, one marker per line. Both are
  read positionally.
- Marker names are used unchanged in every output.

## Config

`config/snakemake_config.yaml` holds `project_name` (input sub-directory, also used
to name output directories), `input_dir`, `marker_file`, `thresholds_file`,
`channel` (segmentation channel), `gpu`, `target_cells` (cells kept after sampling)
and `random_seed`.

### `config/thresholds.csv` is not tracked

Thresholds are tuned per dataset and are not distributed. Provide the file at the
path set in `thresholds_file`:

```csv
,CD3,CD4,CD8,...
SAMPLE_01_ROI_01,0.90,0.85,0.80,...
```

The sample name is the TIFF basename without extension. Columns are markers in
`marker_list.txt` order. Values are quantiles in `[0, 1)`, where `0` means no
threshold, as used for `DAPI`. The pipeline stops with an error if a sample is
missing or if the number of quantiles does not match the number of markers.

## Run

```bash
snakemake -n                                    # dry run
snakemake --cores 4 --use-conda                 # default target: the aggregated h5ad
snakemake --cores 4 --use-conda --resources gpu=1
snakemake --cores 4 --use-conda --until preprocessing
```

If segmentation runs out of CUDA memory, one ROI is too large for the GPU. Set
`gpu: False` and run again.

## Output

```
data/01_Raw_ROI_clean/<PROJECT>/<sample>.tif
data/02_Segmented/<PROJECT>/<sample>_segmented.zarr
data/03_Processed/<PROJECT>/<sample>_quantified.zarr
output/plots/<PROJECT>/Preprocessing/<sample>_{threshold,tcell,area_hist,segmentation}.png
data/04_Analysis/<PROJECT>/Processed_adata_<PROJECT>.h5ad
logs/<PROJECT>/{segmentation,preprocessing}/<sample>.log
```

The final `AnnData` has one row per cell and one column per marker.

| Slot | Content |
|------|---------|
| `X` | CLR-normalized, scaled mean intensities |
| `layers["percentage_positive"]` | fraction of positive pixels, when present in the source zarr |
| `obs` | `source_file`, `sample_id`, `batch`, centroids, `leiden_0.1/0.25/0.65/0.85` |
| `obsm` | `X_pca`, `X_umap`, `spatial` |

Cell types are not assigned automatically. Clusters have to be labelled from a
marker heatmap before phenotype-level analysis.

## Repository layout

```
Snakefile            # 3 rules, config loaded from config/
config/              # pipeline config, marker/channel list
src/codexanalyzer/   # internal code: loading, preprocessing, normalization, plots
src/scripts/         # argparse CLIs called by the rules
environment.yml      # conda env, also used by the conda: directives
```

## Not included

`notebooks/` contains the per-dataset analysis (cluster annotation, scimap and
squidpy spatial statistics, figures). It is not published because it is tied to
identifiable cohort samples. Data, `output/`, `logs/`, `*.h5ad`, `*.zarr/` and
`config/thresholds.csv` are gitignored.

## Citations

The pipeline uses spatialproteomics:

> Meyer-Bender M., Vöhringer H., Schniederjohann C., et al. Spatialproteomics: an
> interoperable toolbox for analyzing highly multiplexed fluorescence image data.
> *Nature Methods* 23, 1564-1573 (2026). doi:10.1038/s41592-026-03155-1

It was used for the CODEX analysis in:

> Marcoux P., Gava F., Tosolini M., et al. Patient-derived lymphoma spheroids
> reveal predictive markers of glofitamab resistance in relapsed/refractory
> B-NHL. *Blood* 147(23), 2770 (2026). doi:10.1182/blood.2025031309

## License

MIT, see `LICENSE`. It covers the code only. No data is distributed.

## Contact

Paul Marcoux, Inserm / CRCT, Toulouse, `paul.marcoux@inserm.fr`
