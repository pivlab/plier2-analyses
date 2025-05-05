library(bigstatsr)
library(here)

config <- list()

config$GENERAL <- list(
  N_CORES=max(nb_cores() %/% 4, 1),
  CHUNK_SIZE=100,
  DATA_DIR=here(file.path("data")),
  OUTPUT_DIR=here(file.path("output"))
)
  
config$ARCHS4=list(
  DATASET_NAME="archs4",
  DATASET_FOLDER=here(
    file.path(config$GENERAL$OUTPUT_DIR, "archs4")
  ),
  DATASET_ENSEMBL_VERSION=107,
  URL="https://s3.dev.maayanlab.cloud/archs4/files/human_gene_v2.5.h5",
  DATASET_FILE=here(file.path("data", "archs4", "human_gene_v2.5.h5")),
  GENES_MEAN_CUTOFF=0.5,
  GENES_VAR_CUTOFF=0.1,
  PLIER_PARAMS=list(
    RANDOM_SVD_N_CORES=1,
    RANDOM_SVD_SEED=123,
    K=942,
    MULTIPLIER=3,
    MAX_U_UPDATES=3,
    MAX_ITER=350,
    N_CORES=1
  )
)
