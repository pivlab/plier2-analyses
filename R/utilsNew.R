library(Matrix)
library(glmnet)
library(dplyr)
library(bigstatsr)
getGMT <- function(url) {
  
  tmp_file <- tempfile(fileext = ".gmt")
  
  download.file(url, tmp_file)
  
  result <- read_gmt(tmp_file)
  
  unlink(tmp_file)
  
  return(result)
}


read_gmt=function (filename) {
  gmt = list()
  lines = readLines(filename)
  for (line in lines) {
    line = gsub("\"", "", trimws(line))
    sp = unlist(strsplit(line, "\t"))
    sp[3:length(sp)] = gsub(",.*$", "", sp[3:length(sp)])
    gmt[[sp[1]]] = sort(unique(sp[3:length(sp)]))
  }
  return(gmt)
}


gmtListToSparseMat=function(gmtList){
  
  allnames=unlist(lapply(gmtList, names))
  #there are usually no duplicates
  stopifnot(all(table(allnames)==1))
  allGenes=unique(unlist(lapply(gmtList, unlist)))
  
  row_indices <- integer(0)
  col_indices <- integer(0)
  values <- integer(0)
  for (gmt in seq_along(gmtList)) {
    for (path in names(gmtList[[gmt]])) {
      thisPathGenes <- gmtList[[gmt]][[path]]
      iiGenes <- match(thisPathGenes, allGenes)
      iPath <- match(path, allnames)
      
      # Store indices and values
      row_indices <- c(row_indices, iiGenes)
      col_indices <- c(col_indices, rep(iPath, length(iiGenes)))
      values <- c(values, rep(1, length(iiGenes)))
    }
  }
  
  # Use sparseMatrix to create the matrix in one go
  
  pathMat=sparseMatrix(i = row_indices, j = col_indices, x = values,
                       dims = c(length(allGenes), length(allnames)))
  rownames(pathMat) = allGenes
  colnames(pathMat) = allnames
  pathMat
}

#' @keywords  internal
#' find the common rows of two data matrices
#' @param data1 first data matrix
#' @param data2 second data matrix
commonRows=function(data1, data2){
  intersect(rownames(data1), rownames(data2))
}


#' @keywords internal
pinv.ridge=function (m, alpha = 0) 
{
  msvd = svd(m)
  if (length(msvd$d) == 0) {
    return(array(0, dim(m)[2:1]))
  }
  else {
    if (alpha > 0) {
      ss = (msvd$d^2) + alpha^2
      msvd$d = ss/msvd$d
    }
    out = msvd$v %*% (1/msvd$d * t(msvd$u))
    rownames(out) = rownames(m)
    colnames(out) = colnames(m)
    out
  }
}

#' @keywords  internal
BH= function(pval){p.adjust(pval, method="BH")}


cleanFBM <- function(fbm, ncores = 1) {
  # 1. Block‐wise scan for max and NA
  stats <- big_apply(
    fbm,
    a.FUN = function(X, ind) {
      vals <- X[, ind, drop = FALSE]
      list(
        max  = if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE),
        na   = anyNA(vals)
      )
    },
    a.combine = function(...) {
      Reduce(function(a, b) {
        list(
          max = max(a$max, b$max, na.rm = TRUE),
          na  = a$na  || b$na
        )
      }, list(...))
    },
    ind        = cols_along(fbm),
    ncores     = ncores,
  )

  max_value <- stats$max
  has_na    <- stats$na

  # 2. Log2 transform if necessary
  if (!is.na(max_value) && max_value >= 100) {
    message("Applying log2 transformation")
    big_apply(
      fbm,
      a.FUN     = function(X, ind) { X[, ind] <- log2(X[, ind] + 1); NULL },
      ind       = cols_along(fbm),
      ncores    = ncores,
    )
  } else {
    message("Already on log scale or all NA")
  }

  # 3. Fill NAs if present
  if (has_na) {
    message("Filling NAs with 0")
    big_apply(
      fbm,
      a.FUN     = function(X, ind) { X[, ind][is.na(X[, ind])] <- 0; NULL },
      ind       = cols_along(fbm),
      ncores    = ncores,
    )
  } else {
    message("No NA values found")
  }

  return(list(max_value = max_value, had_na = has_na))
}


computeRowStatsFBM <- function(fbm, ncores = 1) {
  # Compute row sums in blocks
  row_sums <- big_apply(
    fbm,
    a.FUN     = function(X, ind) rowSums(X[, ind]),
    a.combine = "plus",
    ncores    = ncores
  )

  # Compute row sums of squares in blocks
  row_sums_sq <- big_apply(
    fbm,
    a.FUN     = function(X, ind) rowSums(X[, ind]^2),
    a.combine = "plus",
    ncores    = ncores
  )

  n_cols <- ncol(fbm)
  # Final means and variances
  row_means     <- row_sums / n_cols
  row_variances <- (row_sums_sq / n_cols) - (row_means^2)

  list(row_means = row_means, row_variances = row_variances)
}



zscoreFBM <- function(fbm, rowStats, chunk_size = 1000) {
  message("Applying Z-score transformation")
  
  row_means <- rowStats$row_means
  row_variances <- rowStats$row_variances
  
  
  # Compute standard deviations upfront
  row_sds <- sqrt(row_variances)
  
  # Iterate over columns in chunks to transform in place
  for (start_col in seq(1, ncol(fbm), by = chunk_size)) {
    end_col <- min(start_col + chunk_size - 1, ncol(fbm))
    
    # Extract the current chunk of columns
    col_chunk <- fbm[, start_col:end_col]
    
    # Compute Z-score: (X - mean) / sd
    col_chunk <- sweep(col_chunk, 1, row_means, FUN = "-") # Subtract row means
    col_chunk <- sweep(col_chunk, 1, row_sds, FUN = "/")   # Divide by row SD
    
    # Store back in FBM
    fbm[, start_col:end_col] <- col_chunk
  }
}







filterFBM<- function(fbm, rowStats, keep_samples_idx=NULL, mean_cutoff = NULL, var_cutoff = NULL, backingfile = "filtered_fbm") {
  row_means <- rowStats$row_means
  row_variances <- rowStats$row_variances
  
  # Determine rows to keep based on cutoffs
  keep_rows <- rep(TRUE, length(row_means))  # Default: keep all rows
  
  if (!is.null(mean_cutoff)) {
    keep_rows <- keep_rows & (row_means >= mean_cutoff)
  }
  
  if (!is.null(var_cutoff)) {
    keep_rows <- keep_rows & (row_variances >= var_cutoff)
  }

  if (is.null(keep_samples_idx)) {
    keep_samples_idx <- cols_along(fbm)
  }
  
  # Number of rows to keep
  n_kept <- sum(keep_rows)
  
  if (n_kept == 0) {
    stop("No rows meet the filtering criteria.")
  }
  
  # Create a new FBM with the filtered data
  fbm_filtered <- big_copy(
    X           = fbm,
    ind.row     = which(keep_rows),
    ind.col     = keep_samples_idx,
    backingfile = backingfile
  )

  
  return(list(fbm_filtered = fbm_filtered, kept_rows = which(keep_rows)))
}





randomProjection <- function(input_matrix, rd, verbose = TRUE, pos=F) {
  # Create the random projection matrix
  A <- matrix(rnorm(rd * ncol(input_matrix)), ncol = rd)
  if(pos){
  A=abs(A)
  }
  # Check if input is an FBM object
  if (inherits(input_matrix, "FBM")) {
    if (verbose) message("Processing FBM object with big_prodMat()...")
    # For FBM objects, use big_prodMat from bigstatsr package
    result <- big_prodMat(input_matrix, A)
    return(result)
  } else {
    if (verbose) message("Processing standard matrix with base R...")
    # For regular matrices, use standard matrix multiplication
    result <- input_matrix %*% A
    return(result)
  }
}

cpm_norm <- function(m) {
  # total counts per sample
  lib_sizes <- rowSums(m)
  # scaling factors in millions of reads
  sf <- lib_sizes / 1e6
  # divide each row by its scaling factor
  normalized <- m / sf
  rownames(normalized) <- rownames(m)
  colnames(normalized) <- colnames(m)
  return(normalized)
}

tpm_norm <- function(counts, gene.lengths) {
  if (!is.matrix(counts)) {
    stop("`counts` must be a matrix.")
  }
  counts <- as.matrix(counts)
  
  if (is.null(names(gene.lengths))) {
    stop("`gene.lengths` must be a named numeric vector.")
  }
  if (!all(colnames(counts) %in% names(gene.lengths))) {
    stop("All column names of `counts` must be in names(gene.lengths).")
  }
  
  # Reorder lengths to match columns
  #lengths.bp <- gene.lengths[colnames(counts)]
  
  # Convert lengths to kilobases
  lengths.kb <- gene.lengths / 1e3
  
  # 1) Divide counts by gene length in kilobases → reads per kilobase
  rpk <- sweep(counts, 2, lengths.kb, FUN = "/")
  
  # 2) Compute per-sample scaling factor: sum of RPKs
  per.sample.sum <- rowSums(rpk)
  
  # 3) Divide RPKs by the sum and multiply by 1e6 → TPM
  tpm <- sweep(rpk, 1, per.sample.sum, FUN = "/") * 1e6
  
  return(tpm)
}

