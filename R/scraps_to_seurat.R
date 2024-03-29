#' Read scraps output from umi_tools to sparseMatrix
#' 
#' @param file scraps output table containing 
#' @param n_min  minimum number of observations for a site to be kept
#' @param gene_min  minimum number of observations for gene to be kept
#' @param alt_only if TRUE, only keep genes with alternative polyA sites
#' @param filter_low if TRUE, remove low genes and cell ids, if FALSE, set to -1
#' @param cell_ids if given, use only these cell barcodes, and fill in empty ones
#' @param types filter to only these types of PA sites, set to NULL to use all
#' @param types2 another filter to only these types of PA sites, set to NULL to use all
#' @param pf data.frame with SAF first column and position factor by gene, calculated by `parse_saf_pf`
#' @import readr dplyr stringr tidyr
#' @return count matrix
#' @export
scraps_to_matrix <- function(file,
                             n_min = 5,
                             gene_min = 5,
                             alt_only = FALSE,
                             filter_low = TRUE,
                             cell_ids = NULL,
                             types = "3'UTR",
                             types2 = NULL,
                             pf = NULL) {
  # read file
  message("read file")
  dat <- read_tsv(file)
  
  # filter for certain ids
  if (!is.null(cell_ids)) {
    dat <- dat %>% 
      filter(cell %in% cell_ids)
  }
  
  # filter if total counts too low
  dat <- dat %>% group_by(gene) %>% 
    filter(sum(count) >= n_min) %>% 
    ungroup()
  
  if (!is.null(types)) {
    message("filter for ", str_c(types, collapse = ". "))
    dat <- dat %>% filter(str_detect(gene, types))
  }
  
  if (!is.null(types2)) {
    message("filter for ", str_c(types2, collapse = ". "))
    dat <- dat %>% filter(str_detect(gene, types2))
  }
  
  # if TRUE, toss observations where gene only has 1 polyA site
  if (alt_only) {
    message("filter for genes with alternative sites")
    dat2 <- dat %>% separate(gene, 
                             sep = "_", 
                             into = c("gene_name", "gene_info"), 
                             remove = FALSE,
                             extra = "merge") %>% 
      distinct(gene, gene_name, gene_info)
    
    keeps2 <- dat2$gene_name %>% table() %>%
      as.data.frame() %>%
      setNames(c("gene_name", "n")) %>%
      filter(n > 1) %>%
      pull(gene_name) %>% 
      as.vector()
    
    keeps <- dat2 %>% filter(gene_name %in% keeps2) %>% 
      pull(gene)
    
    dat <- dat %>% filter(gene %in% keeps)
  }
  
  # use psi values if pf is used
  if (!is.null(pf)) {
    message("calculate normalized PSI value")
    if (!is.null(gene_min) & gene_min > n_min) {
      if (filter_low) {
        dat <- dat %>% group_by(gene) %>% 
          mutate(total = sum(count)) %>% 
          filter(total >= gene_min) %>% 
          ungroup() %>% 
          select(-total)
        
        dat <- dat %>% inner_join(pf, by = c("gene" = "GeneID")) %>% 
          mutate(gene_name = str_remove(gene, "_.+")) %>% 
          mutate(count2 = pf * count) %>% 
          group_by(gene_name, cell) %>% 
          summarize(count = sum(count2)/sum(count)) %>% 
          ungroup() %>% 
          rename(gene = "gene_name")
      } else {
        dat <- dat %>% group_by(gene) %>% 
          mutate(total = sum(count))
        
        dat <- dat %>% inner_join(pf, by = c("gene" = "GeneID")) %>% 
          mutate(gene_name = str_remove(gene, "_.+")) %>% 
          mutate(count2 = pf * count) %>% 
          group_by(gene_name, cell) %>% 
          summarize(count = sum(count2)/sum(count), total = first(total)) %>% 
          mutate(count = ifelse(total >= gene_min, count, -1)) %>% 
          ungroup() %>% 
          rename(gene = "gene_name")
      }
    }
    

  }
  
  # transform to matrix
  message("convert to matrix")
  # mat <- dat %>% pivot_wider(names_from = cell, values_from = count, values_fill = 0) %>% 
  #   column_to_rownames("gene")
  dat <- dat %>% 
    mutate(gene = as.factor(gene),
           cell = as.factor(cell))
  mat <- Matrix::sparseMatrix(
    i = as.numeric(dat$gene),
    j = as.numeric(dat$cell),
    x = dat$count,
    dims = c(length(unique(dat$gene)),
             length(unique(dat$cell))),
    dimnames = list(
      select(dat, gene) %>%
        unique() %>%
        arrange(gene) %>%
        pull(gene),
      select(dat, cell, ) %>%
        unique() %>%
        arrange(cell) %>%
        pull(cell)
    )
  )

  # keep consistent cell ids
  if (!is.null(cell_ids)) {
    emptys <- setdiff(cell_ids, colnames(mat))
    
    if (length(emptys) > 0) {
      message("filling empty 0s")
      empty_mat <- matrix(data = 0, nrow = nrow(mat), ncol = length(emptys))
      colnames(empty_mat) <- emptys
      mat <- cbind(mat, empty_mat)
    }
    
    mat <- mat[, cell_ids]
  }
  
  return(mat)
}

#' Add scraps output to existing seurat object as new assay
#' 
#' @param file scraps output table containing 
#' @param object seurat object
#' @param n_min  minimum number of observations for a site to be kept
#' @param alt_only if TRUE, only keep genes with alternative polyA sites
#' @param types filter to only these types of PA sites, set to NULL to use all
#' @param pf data.frame with SAF first column and position factor by gene, calculated by `parse_saf_pf`
#' @import dplyr Seurat
#' @return Seurat object with inserted assay
#' @export
scraps_to_seurat <- function(file, object, 
                             assay_name = "Asite",
                             n_min = 5,
                             alt_only = TRUE,
                             types = "3'UTR",
                             pf = NULL) {
  cell_ids <- Cells(object)
  
  mat <- scraps_to_matrix(file,
                          object, 
                          n_min = n_min,
                          alt_only = alt_only,
                          cell_ids = cell_ids,
                          types = types,
                          pf = pf)

  object[[assay_name]] <- CreateAssayObject(counts = mat)
  
  return(object)
}

#' Parse SAF annotation file for position factor (pf) values
#' 
#' @param file SAF file used with Scraps
#' @param types filter to only these types of PA sites, set to NULL to use all
#' @param alt_only if TRUE, only keep genes with alternative polyA sites
#' @import readr dplyr
#' @return data.frame with SAF first column and position factor by gene
#' @export
parse_saf_pf <- function(file,
                         types = "3'UTR",
                         alt_only = TRUE) {
  bed1 <- parse_saf(file, types)
  
  # remove symbols that come from different chromosomes
  diffchrom <- bed1 %>% distinct(gene, chrom) %>% 
    group_by(gene) %>% 
    summarise(n = n()) %>% 
    filter(n > 1) %>%
    pull(gene)
  
  # position factor is 0-1
  bed2_plus <- bed1 %>% 
    filter(strand == "+") %>% 
    filter(!(gene %in% diffchrom)) %>% 
    group_by(gene) %>% 
    arrange(.by_group = TRUE, start) %>% 
    mutate(pf = row_number()) %>% 
    mutate(pf = (pf - 1) / max(pf - 1)) %>% 
    ungroup()
  
  bed2_minus <- bed1 %>% 
    filter(strand == "-") %>% 
    filter(!(gene %in% diffchrom)) %>% 
    group_by(gene) %>% 
    arrange(.by_group = TRUE, desc(start)) %>% 
    mutate(pf = row_number()) %>% 
    mutate(pf = (pf - 1) / max(pf - 1)) %>% 
    ungroup()
  
  bed3 <- bind_rows(bed2_plus, bed2_minus)
  
  if (alt_only) {
    bed3 <- bed3 %>% group_by(gene) %>% 
      filter(n() > 1) %>% 
      ungroup()
  }
  
  bed3 %>% select(GeneID, pf)
}

#' Read SAF annotation file
#' 
#' @param file SAF file used with Scraps
#' @param types filter to only these types of PA sites, set to NULL to use all
#' @param sep separator between fields for full name
#' @import readr dplyr stringr tidyr
#' @return data.frame
#' @export
parse_saf <- function(file, types = FALSE, sep = ";") {
  # read SAF and split for gene symbol
  saf <- read_tsv(file) %>% 
    mutate(GeneID2 = str_remove(GeneID, Chr)) %>% 
    separate(GeneID2, 
             sep = sep, 
             into = c("gene", 
                      "genbank", 
                      "id", 
                      "chrom", 
                      "pos",
                      "strand",
                      "class"),
             convert = TRUE,
             remove = FALSE) %>% 
    filter(!is.na(gene) | !is.na(genbank) | !is.na(id)) %>%
    mutate(pos = as.numeric(pos)) %>% 
    filter(!is.na(pos)) %>% 
    mutate(Start = Start - 1)
  
  if (types != FALSE) {
    saf <- saf %>% filter(str_detect(class, types))}
  
  # convert to bed format
  bed1 <- saf %>% select(chrom = Chr, start = Start, end = End, strand, gene, GeneID)
  bed1
}

#' Parse GTF annotation file
#' 
#' @param file SAF file used with Scraps
#' @return data.frame
#' @export
parse_gtf <- function(file) {
  gtf <- rtracklayer::import(file)
  tbl_gtf <- as.tbl_intervalcustom(gtf)
}

#' Convert GTF to tibble
#' 
#' @param x GTF from import
#' @import dplyr
#' @return data.frame
#' @export
as.tbl_intervalcustom <- function(x) {
  res <- tibble(
    chrom = as.character(x@seqnames),
    start = x@ranges@start - 1,
    end = x@ranges@start - 1 + x@ranges@width,
    name = x@elementMetadata$gene_name,
    transcript = x@elementMetadata$transcript_id,
    transcript_type =  x@elementMetadata$transcript_type,
    score = rep(".", length(x)),
    strand = as.character(x@strand),
    type = x@elementMetadata$type,
    exon_number = as.numeric(x@elementMetadata$exon_number)
  )
  res <- mutate(res, strand = ifelse(strand == "*", ".", strand))
  res
}

