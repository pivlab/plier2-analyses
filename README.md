# PLIER2 <img src="man/figures/plier2.png" width="121px" height="140px" align="right" style="padding-left:10px;background-color:white;" />

**PLIER2 (Pathway-Level Information ExtractoR 2)** is a scalable framework for extracting biologically meaningful latent variables from gene expression data using prior knowledge of gene sets. It builds upon and extends the original [PLIER](https://github.com/wgmao/PLIER) method, offering improvements in performance and reproducibility for large-scale datasets like GTEx and ARCHS4.

## 🔧 Installation

Set up the Conda environment using the provided environment file:

```bash
conda env create -f envs/plier2.yaml
conda activate plier2
```

## 📘 Vignettes

This repository includes several RMarkdown-based vignettes organized by analysis stage and dataset. They serve as both documentation and executable workflows for building and evaluating PLIER2 models.

### `vignettes/00_setup/`
This folder is reserved for initial setup scripts or supporting files. *(Currently empty.)*

### `vignettes/01_model_building/`
This section contains scripts to build PLIER2 models on different datasets.

#### **ARCHS4**
Located in `01_model_building/archs4/`, this subdirectory includes a complete pipeline:
- `archs4-01-data_download.Rmd`: Downloads ARCHS4 data.
- `archs4-02-data_preprocess.Rmd`: Prepares and normalizes the expression data.
- `archs4-03-data_filtering.Rmd`: Applies filters to remove low-variance genes and samples.
- `archs4-04-svd.Rmd`: Performs singular value decomposition (SVD) to estimate dimensionality.
- `archs4-05-simple_decomp.Rmd`: Applies a simple matrix factorization.
- `archs4-06-PLIERv2.Rmd`: Fits the final PLIER2 model.

Each `.Rmd` file has a corresponding `.html` output for visualization and review.

#### **GTEx**
Located in `01_model_building/gtex/`, this vignette:
- `gtex.Rmd`: Downloads, filters, and models the GTEx v8 RNA-Seq dataset using both PLIER and PLIER2.
Includes generated plots (e.g., `figure-html/unnamed-chunk-10-1.png`).

#### **recount2**
Located in `01_model_building/recount2/`, this vignette:
- `recount2.Rmd`: Processes recount2 data and builds a PLIER2 model.

### `vignettes/02_model_comparisons/`
This section compares model performance between PLIER and PLIER2.

#### **GTEx comparisons**
- `gtex_comparisons.Rmd`: Compares latent variable structure, pathway enrichments, and variance explained between PLIER and PLIER2 models.
- `gtex_PLIER2_params.Rmd`: Summarizes and evaluates parameter choices and outputs of the PLIER2 model.

### `vignettes/03_model_performance/`
Reserved for benchmarking and performance evaluation of PLIER2 on various datasets. *(Currently empty.)*

### `vignettes/04_model_biology/`
Planned for biological interpretation and downstream functional analysis using latent variables. *(Currently empty.)*

## ▶️ Rendering Vignettes

To reproduce the full PLIER2 modeling pipelines, render the corresponding RMarkdown files using the following commands:

### GTEx and Model Comparisons

```bash
Rscript -e "rmarkdown::render('vignettes/01_model_building/gtex/gtex.Rmd')"
Rscript -e "rmarkdown::render('vignettes/02_model_comparisons/gtex/gtex_comparisons.Rmd')"
```

### recount2

```bash
Rscript -e "rmarkdown::render('vignettes/01_model_building/recount2/recount2.Rmd')"
```

### ARCHS4 Pipeline

To enable efficient matrix operations when rendering ARCHS4 vignettes, you can control parallelization by setting the following environment variables:

```bash
n_jobs=30
export NUMBA_NUM_THREADS=$n_jobs
export MKL_NUM_THREADS=$n_jobs
export OPENBLAS_NUM_THREADS=$n_jobs
export NUMEXPR_NUM_THREADS=$n_jobs
export OMP_NUM_THREADS=$n_jobs
```

Then run each step of the ARCHS4 pipeline:

```bash
Rscript -e "rmarkdown::render('vignettes/01_model_building/archs4/archs4-01-data_download.Rmd')"
Rscript -e "rmarkdown::render('vignettes/01_model_building/archs4/archs4-02-data_preprocess.Rmd')"
Rscript -e "rmarkdown::render('vignettes/01_model_building/archs4/archs4-03-data_filtering.Rmd')"
Rscript -e "rmarkdown::render('vignettes/01_model_building/archs4/archs4-04-svd.Rmd')"
Rscript -e "rmarkdown::render('vignettes/01_model_building/archs4/archs4-05-simple_decomp.Rmd')"
Rscript -e "rmarkdown::render('vignettes/01_model_building/archs4/archs4-06-PLIERv2.Rmd')"
```

## 📂 Input Data

Each RMarkdown handles downloading and preprocessing of input data as needed. Data types include:

- **Gene expression matrix**: Gene × Sample matrix (e.g., TPM or counts normalized and log-transformed)
- **Pathway matrix**: Gene × Pathway binary matrix

## 📤 Output

Model objects created by PLIER2 contain:

- `Z`: Gene loadings matrix (`Genes × Latent Variables`)
- `B`: Latent variable activity matrix (`Latent Variables × Samples`)
- `U`: Prior information coefficients (`Pathways × Latent Variables`)
- `summary`: Metadata describing each latent variable (e.g., pathway enrichment, variance explained)

## 💡 Citation

If you use this software in your research, please cite:

> Subirana-Granés, M., Chikina, M., Pividori, M. *A scalable pathway-level information extractor for massive gene expression datasets*. (In preparation)

## 📫 Contact

For questions or issues, please open a GitHub [issue](https://github.com/your-username/PLIER2/issues).
