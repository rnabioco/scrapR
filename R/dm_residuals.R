library(tidyverse)
library(Seurat)
library(Matrix)
library(MatrixExtra)

dm_residuals <-
  function(seurat_object,
           pa_assay = "utrA",
           gene_pattern,
           site_prefix = NULL,
           full_sparse = FALSE) {
    # check if provided gene pattern matches >1 row
    if (length(grep(
      pattern = gene_pattern,
      rownames(seurat_object@assays[[pa_assay]]@counts),
      value = T
    )) > 1) {
      # subset matrix to only rows matching gene pattern and transpose for MGLM
      X <-
        t(seurat_object@assays[[pa_assay]]@counts[grep(pattern = gene_pattern,
                                                       rownames(seurat_object@assays[[pa_assay]]@counts),
                                                       value = T),])
      # remove cells with 0 UMIs
      X <- X[rowSums(X) > 0, ]
      # scale all cells to a total of 10 UMIs
      X <- X / rowSums(X) * 10
      # fit a Dirichlet multinomial model
      dm <-
        MGLM::MGLMfit(as.matrix(X),
                      dist = "DM")
      # calculate scaled model means
      E <- 10 * (dm@estimate / sum(dm@estimate))
      # calculate scaled model mean variances
      Var <-
        E * (1 - (dm@estimate / sum(dm@estimate))) * ((10 + sum(dm@estimate)) / (1 + sum(dm@estimate)))
      # calculate residuals for all sites in matrix
      mat1 <- t((t(X) - E) / sqrt(Var))
      
      # check if sites are to be prefixed
      if (!is.null(site_prefix)) {
        colnames(mat1) <- paste0(site_prefix, colnames(mat1))
      }
      # check if a full sparse matrix matching the original seurat cells is requested
      if (full_sparse == TRUE) {
        # replace zeros in matrix with 0.001
        mat1[mat1 == 0] <- 0.001
        # convert to sparsematrix
        mat1 <- as(mat1, "sparseMatrix")
        # create an empty dgCMatrix from "missing" cells
        mat2 <-
          MatrixExtra::emptySparse(
            nrow = length(Cells(seurat_object)[!Cells(seurat_object) %in% rownames(mat1)]),
            ncol = ncol(mat1),
            dtype = "d",
            format = "C"
          )
        # assign cell barcodes to empty matrix
        rownames(mat2) <-
          Cells(seurat_object)[!Cells(seurat_object) %in% rownames(mat1)]
        # bind matrices
        mat1 <- rbind(mat1, mat2)
        # reorder rownames to match cells in seurat object
        mat1 <- mat1[Cells(seurat_object),]
        # remove rownames to reduce matrix size
        rownames(mat1) <- NULL
        # transpose matrix to match original
        mat1 <- t(mat1)
        mat1
      } else {
        mat1
      }
    } else {
      message("Error: Only one site present")
    }
  }

rbind_csr_dm <- function(dm_matrices, seurat_obj){
  mat <- do.call("rbind_csr", dm_matrices)
  colnames(mat) <- Cells(seurat_obj)
  mat
}
