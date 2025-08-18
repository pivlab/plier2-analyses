# PLIER2: Bigger, Better, and Faster

**PLIER2** is a high-performance, scalable reimplementation of the [Pathway-Level Information Extractor (PLIER)](https://doi.org/10.1038/s41592-019-0456-1) framework for latent variable analysis of large-scale transcriptomic datasets.  
It integrates prior biological knowledge into matrix factorization to generate interpretable latent variables (LVs) that capture meaningful biological pathways while mitigating technical noise.  

This package implements major algorithmic and computational improvements over PLIER, enabling the analysis of modern transcriptomic compendia such as **GTEx**, **recount2**, and **ARCHS4**.

## Key Features

- **Two-phase modeling**:
  - **PLIERbase**: Unsupervised factorization without priors for fast initialization.
  - **PLIERfull**: Integration of prior knowledge via regularized regression (glmnet).
- **Scalable to massive datasets** using memory-mapped matrices from `bigstatsr`.
- **Automated per-LV regularization tuning** via internal cross-validation.
- **Improved biological specificity** in LV–pathway associations.
- **7×–41× speed improvements** over PLIER.
- **Full support for very large compendia** (e.g., successful ARCHS4 modeling).

## Dependencies

Before running, ensure you have:

1. **Conda environment**  
   Create the environment using:
   ```bash
   conda env create -f envs/plier2-analyses.yaml
   conda activate plier2-analyses
   ```

2. **Clone the PLIER2 repository**  
   ```bash
   git clone https://github.com/chikinalab/PLIER2.git
   ```

3. **Install PLIER2 in R**  
   ```r
   install.packages("devtools")
   devtools::install("PLIER2")
   ```

All required input data will be downloaded automatically by the provided scripts.

## How to Run

Once dependencies are installed and the environment is set up:

```bash
conda activate plier2-analyses

# Example: running a notebook or analysis script
jupyter nbconvert --to notebook --execute nbs/01_model_building/archs4/00_archs4.ipynb --inplace
```

The scripts automatically handle downloading all necessary datasets and generating intermediate results.

## Citation

If you use **PLIER2**, please cite:

> Subirana-Granés M, Nandi S, Zhang H, Chikina M, Pividori M.  
> *PLIER2: bigger, better and faster*. bioRxiv, 2025.  
> doi: [10.1101/2025.06.05.658122](https://doi.org/10.1101/2025.06.05.658122)

## License

This project is licensed under the [CC-BY 4.0 License](http://creativecommons.org/licenses/by/4.0/).

## Acknowledgments

Supported by the National Human Genome Research Institute,  
The Eunice Kennedy Shriver National Institute of Child Health and Human Development,  
the National Science Foundation, and the National Eye Institute.
