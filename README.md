# PLIER2: Bigger, Better, and Faster

**PLIER2** is a high-performance, scalable reimplementation of the [Pathway-Level Information Extractor (PLIER)](https://doi.org/10.1038/s41592-019-0456-1) framework for latent variable analysis of large-scale transcriptomic datasets.

It integrates prior biological knowledge into matrix factorization to generate interpretable latent variables (LVs) that capture meaningful biological pathways while mitigating technical noise.  

This package implements major algorithmic and computational improvements over PLIER, enabling the analysis of modern transcriptomic compendia such as **GTEx**, **recount2**, and **ARCHS4**.

---

## Key Features

- **Two-phase modeling**:
  - **PLIERbase**: Unsupervised factorization without priors for fast initialization.
  - **PLIERfull**: Integration of prior knowledge via regularized regression (glmnet).
- **Scalable to massive datasets** using memory-mapped matrices from `bigstatsr`.
- **Automated per-LV regularization tuning** via internal cross-validation.
- **Improved biological specificity** in LV–pathway associations.
- **7×–41× speed improvements** over PLIER.
- **Full support for very large compendia** (e.g., successful ARCHS4 modeling).

## Installation


## Quick Start

## Citation

If you use **PLIER2**, please cite:

> Subirana-Granés M, Nandi S, Zhang H, Chikina M, Pividori M.  
> *PLIER2: bigger, better and faster*. bioRxiv, 2025.  
> doi: [10.1101/2025.06.05.658122](https://doi.org/10.1101/2025.06.05.658122)

---

## License

This project is licensed under the [CC-BY 4.0 License](http://creativecommons.org/licenses/by/4.0/).

---

## Acknowledgments

Supported by the National Human Genome Research Institute,  
The Eunice Kennedy Shriver National Institute of Child Health and Human Development,  
the National Science Foundation, and the National Eye Institute.
