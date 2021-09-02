#' @title Read the reference gene file of a given strain
#' @name read_master_gene_file
#'
#' @return a data.frame
#'
#' @export
#'
#' @examples
read_master_gene_file <- function(strain) {
  master_gene_file <- system.file("extdata", paste0("master_", strain, ".tsv"), package = "SalmoGenes")
  master_table <- read.delim(master_gene_file, comment.char = "#", header = T, stringsAsFactors = F)
  master_table
}

##
read_all_strains_genes_file <- function() {
  all_genes_list <- system.file("extdata", "all_strains_all_genes.txt", package = "SalmoGenes")
  all_genes_list <- dget(all_genes_list)
  all_genes_list
}

##
read_all_strains_proteins_file <- function() {
  all_proteins_list <- system.file("extdata", "all_strains_all_proteins.txt", package = "SalmoGenes")
  all_proteins_list <- dget(all_proteins_list)
  all_proteins_list
}
##
get_gene_id <- function(gene_names, strain) {
  # master_table <- read.delim(paste0("Strains/dict_files/master_", strain, ".tsv"), header = T, stringsAsFactors = F, sep = "\t", na.strings = c("NA", "<NA>", ""))
  master_table <- read_master_gene_file(strain)

  gene_list_by_symbol <- split(master_table, master_table$gene_name)
  gene_list_by_synonyms <- master_table %>% tidyr::separate_rows(gene_synonyms, sep = ",")
  gene_list_by_synonyms <- split(gene_list_by_synonyms, gene_list_by_synonyms$gene_synonyms)
  gene_list_by_symbol_or_synonym <- c(gene_list_by_symbol, gene_list_by_synonyms)

  convert_symbols <- function(x) {
    ifelse(!is.null(gene_list_by_symbol_or_synonym[[x]]$gene_id[1]), "",
           ifelse(!is.na(x), warning(paste0("This gene is unknown and will be converted to `NA`: ", x), call. = NA), NA))
    ifelse(!is.null(gene_list_by_symbol_or_synonym[[x]]$gene_id[1]), gene_list_by_symbol_or_synonym[[x]]$gene_id[1], NA)
  }
  gene_ids <- sapply(gene_names, FUN=convert_symbols)
  unname(gene_ids)
}

##
get_gene_id_from_protein <- function(protein_names, strain) {
  master_table <- read_master_gene_file(strain)

  protein_list <- master_table %>% tidyr::separate_rows(protein_synonyms, sep = ",")
  protein_list <- split(protein_list, protein_list$protein_synonyms)

  convert_proteins <- function(x) {
    ifelse(!is.null(protein_list[[x]]$gene_id[1]), "",
           ifelse(!is.na(x), warning(paste0("This protein is unknown and will be converted to `NA`: ", x), call. = NA), NA))
    ifelse(!is.null(protein_list[[x]]$gene_id[1]), protein_list[[x]]$gene_id[1], NA)
  }
  gene_ids <- sapply(protein_names, FUN=convert_proteins)
  unname(gene_ids)
}
##
get_gene_name <- function(gene_ids, strain) {
  # master_table <- read.delim(paste0("Strains/dict_files/master_", strain, ".tsv"), header = T, stringsAsFactors = F, sep = "\t", na.strings = c("NA", "<NA>", ""))
  master_table <- read_master_gene_file(strain)

  gene_list_by_id <- split(master_table, master_table$gene_id)
  gene_list_by_synonyms <- master_table %>% tidyr::separate_rows(gene_synonyms, sep = ",")
  gene_list_by_synonyms <- split(gene_list_by_synonyms, gene_list_by_synonyms$gene_synonyms)
  gene_list_by_symbol_or_synonym <- c(gene_list_by_id, gene_list_by_synonyms)

  convert_symbols <- function(x) {
    ifelse(!is.null(gene_list_by_symbol_or_synonym[[x]]$gene_name[1]), "",
           ifelse(!is.na(x), warning(paste0("This gene is unknown and will be converted to `NA`: ", x), call. = NA), NA))
    ifelse(!is.null(gene_list_by_symbol_or_synonym[[x]]$gene_name[1]), gene_list_by_symbol_or_synonym[[x]]$gene_name[1], NA)
  }
  gene_names <- sapply(gene_ids, FUN=convert_symbols)
  unname(gene_names)
}

##
get_strain_from_gene <- function(gene_ids) {
  all_strains <- read_all_strains_genes_file()

  check_strain <- function(x) {
    if(!is.na(x)) {
      # search <- lapply(all_strains, function(y) grep(x, y))
      search <- lapply(all_strains, function(y) x %in% y)
      # search <- search[lapply(search, length) > 0]
      res <- paste0(names(unlist(search))[which(search==TRUE)], collapse = ",")
    } else {
      res <- NA
    }
  }
  strains <- sapply(gene_ids, FUN=check_strain)
  gsub("str_", "", unname(strains))
}

##
get_strain_from_protein <- function(protein_names) {
  all_strains <- read_all_strains_proteins_file()

  check_strain <- function(x) {
    if(!is.na(x)) {
      # search <- lapply(all_strains, function(y) grep(x, y))
      search <- lapply(all_strains, function(y) x %in% y)
      # search <- search[lapply(search, length) > 0]
      res <- paste0(names(unlist(search))[which(search==TRUE)], collapse = ",")
    } else {
      res <- NA
    }
  }
  strains <- sapply(protein_names, FUN=check_strain)
  gsub("str_", "", unname(strains))
}
##
get_all_gene_ids <- function(strain) {
  master_table <- read_master_gene_file(strain)
  master_table$gene_id
}

##
get_all_gene_names <- function(strain) {
  master_table <- read_master_gene_file(strain)
  master_table$gene_name
}
