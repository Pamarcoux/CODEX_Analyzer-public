import os

configfile: "config/snakemake_config.yaml"

PROJECT = config["project_name"]
INPUT_DIR = os.path.join(config["input_dir"], PROJECT)

# Dynamic discovery of samples in the project directory
SAMPLES, = glob_wildcards(os.path.join(INPUT_DIR, "{sample}.tif"))

if not SAMPLES:
    raise ValueError(
        f"No .tif found in {INPUT_DIR}. "
        "Check 'input_dir' and 'project_name' in config/snakemake_config.yaml."
    )

rule all:
    input:
        f"data/04_Analysis/{PROJECT}/Processed_adata_{PROJECT}.h5ad"

rule segmentation:
    input:
        tif = os.path.join(config["input_dir"], "{Project}/{sample}.tif"),
        markers = config["marker_file"] 
    output:
        zarr = directory("data/02_Segmented/{Project}/{sample}_segmented.zarr")
    log:
        "logs/{Project}/segmentation/{sample}.log"
    conda:
        "environment.yml"
    params:
        channel = config["channel"],
        gpu = str(config["gpu"])
    resources:
        gpu = 1 if config["gpu"] else 0
    shell:
        """
        export PYTHONPATH=src:${{PYTHONPATH:-}}
        python src/scripts/segmentation_rule.py \
            --input \"{input.tif}\" \
            --output \"{output.zarr}\" \
            --markers \"{input.markers}\" \
            --channel \"{params.channel}\" \
            --gpu \"{params.gpu}\" > {log} 2>&1
        """

rule preprocessing:
    input:
        zarr = "data/02_Segmented/{Project}/{sample}_segmented.zarr",
        thresholds = config["thresholds_file"],
        markers = config["marker_file"] 
    output:
        quant_zarr = directory("data/03_Processed/{Project}/{sample}_quantified.zarr"),
        plot_threshold = "output/plots/{Project}/Preprocessing/{sample}_threshold.png",
        plot_tcell = "output/plots/{Project}/Preprocessing/{sample}_tcell.png",
        plot_hist = "output/plots/{Project}/Preprocessing/{sample}_area_hist.png",
        segmentation_plot = "output/plots/{Project}/Preprocessing/{sample}_segmentation.png"
    log:
        "logs/{Project}/preprocessing/{sample}.log"
    conda:
        "environment.yml"
    shell:
        """
        export PYTHONPATH=src:${{PYTHONPATH:-}}
        python src/scripts/preprocessing_rule.py \
            --zarr \"{input.zarr}\" \
            --quant_output \"{output.quant_zarr}\" \
            --thresholds \"{input.thresholds}\" \
            --markers \"{input.markers}\" \
            --sample \"{wildcards.sample}\" \
            --threshold_plot \"{output.plot_threshold}\" \
            --tcell_plot \"{output.plot_tcell}\" \
            --area_hist \"{output.plot_hist}\" \
            --segmentation_plot \"{output.segmentation_plot}\" > {log} 2>&1
        """

rule aggregation:
    input:
        zarrs = expand("data/03_Processed/{Project}/{sample}_quantified.zarr", Project=PROJECT, sample=SAMPLES),
        marker_list = config["marker_file"] 
    output:
        h5ad = "data/04_Analysis/{Project}/Processed_adata_{Project}.h5ad"
    log:
        "logs/{Project}/aggregation.log"
    conda:
        "environment.yml"
    params:
        target_cells = config["target_cells"],
        seed = config["random_seed"]
    shell:
        """
        export PYTHONPATH=src:${{PYTHONPATH:-}}
        python src/scripts/aggregation_rule.py \
            --input_files {input.zarrs} \
            --output_file \"{output.h5ad}\" \
            --marker_list \"{input.marker_list}\" \
            --target_cells {params.target_cells} \
            --random_seed {params.seed} > {log} 2>&1
        """
