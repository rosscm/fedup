#' Returns a list of pathways from a file.
#' Currently supports the following file format: GMT, TXT, XLSX
#' @param pathway_file (char) path to file with pathway annotations
#' @param header (logical) whether pathway_file has a header (default = TRUE)
#' @param pathway_col (char) column name with pathway identifiers, one pathway per row
#'  use for non-GMT input (eg "Pathway.ID", default = NULL)
#' @param gene_col (char) column name with gene identifiers
#'  use for non-GMT input (eg "Gene.ID", default = NULL)
#' @param MIN_GENE (integer) minimum number of genes to be considered in a pathway
#'  (default = 1)
#' @param MAX_GENE (integer) maximum number of genes to be considered in a pathway
#'  (default = Inf)
#' @param GO_class (char) if GO terms are present in 'pathway_file', which class(es)
#'  should be used (i.e., BP, MF, CC) (default = NULL)
#' @return a list of vectors with pathway annotations
#' @examples
#' pathways <- writePathways(
#'  system.file("extdata", "YeastDatabase_GO_gmt.gmt", package = "FEDUP"),
#'  MIN_GENE = 10, MAX_GENE = 500, GO_class = "BP"
#')
#' pathways <- writePathways(
#'  system.file("extdata", "SAFE_terms.xlsx", package = "FEDUP"),
#'  header = TRUE, pathway_col = "Enriched.GO.names", gene_col = "Gene.ID",
#'  MIN_GENE = 10, MAX_GENE = 500
#')
#' @import openxlsx
#' @import tibble
#' @importFrom stats aggregate
#' @importFrom utils head read.delim tail
#' @export
writePathways <- function(pathway_file,
                          header = TRUE,
                          pathway_col,
                          gene_col,
                          MIN_GENE = 1L,
                          MAX_GENE = Inf,
                          GO_class = NULL) {
  # Vector of supported file types
  supported_types <- c("gmt", "txt", "xlsx")
  file_type <- sub(".*\\.", "", pathway_file) # grab pathway_file extension

  # Stop if pathway_file input is not supported
  if (!file_type %in% supported_types) {
    stop(paste0("Sorry, input pathway file type (", file_type, ") is not supported. ",
                "Please use one of the following formats: ",
                paste(supported_types, collapse = ", "), "."))
  } else {
    if (file_type == "gmt") {
      # Read in pathway_file and convert to list format
      pathway_in <- strsplit(readLines(pathway_file), "\t")
      if (header) { pathway_in <- pathway_in[-1] } # remove header if TRUE
      pathways <- lapply(pathway_in, tail, -2)
      names(pathways) <- sapply(pathway_in, head, 1)
    } else {
      # Read in pathway_file
      if (file_type == "xlsx") { pathway_in <- read.xlsx(pathway_file) }
      if (file_type == "txt")  { pathway_in <- read.delim(pathway_file, header = header) }

      # Stop if pathway_col or gene_col is missing or not found in pathway_file
      if (missing(pathway_col) || !pathway_col %in% colnames(pathway_in)) {
        stop(paste0("Pathway ID column (", pathway_col, ") not found in ", basename(pathway_file)))
      }
      if (missing(gene_col) || !gene_col %in% colnames(pathway_in)) {
        stop(paste0("Gene ID column (", gene_col, ") not found in ", basename(pathway_file)))
      }

      # Generate pathway list from dataframe
      pathway_df <- data.frame(
        pathway = pathway_in[,pathway_col],
        gene = pathway_in[,gene_col]
      )
      pathway_df <- aggregate(gene ~ pathway, data = pathway_df, FUN = paste)
      pathway_df$gene <- strsplit(pathway_df$gene, "[[:punct:]] ", perl = TRUE)
      pathways <- deframe(pathway_df)
    }

    # Subset for pathways within MIN_GENE, MAX_GENE range
    pathway_size <- lapply(pathways, length)
    pathways_sub <- pathways[which(pathway_size >= MIN_GENE & pathway_size <= MAX_GENE)]

    # Print messages to console
    cat(sprintf(
     "* Input pathway file: %s
      ** total pathways: %s
      ** pathways after size filtering (min %s, max %s): %s\n",
     basename(pathway_file), length(pathways), MIN_GENE, MAX_GENE, length(pathways_sub)
    ))

    # Subset for chosen 'GO_class' if specified
    if (!is.null(GO_class)) {
      go_keep <- paste(GO_class, collapse = "|")
      int_keep <- grep(go_keep, names(pathways_sub))
      pathways_sub <- pathways_sub[int_keep]
      cat(sprintf("      ** pathways after GO class filtering (%s): %s\n",
        paste(GO_class, sep = ", "), length(pathways_sub)))
    }

    # Finally, check for any duplicated pathway annotations
    pathways_dup <- which(duplicated(names(pathways_sub)))
    if (length(pathways_dup)) {
      pathways_sub <- pathways_sub[-pathways_dup]
      cat(sprintf("      ** pathways after filtering out duplicated annotations: %s\n",
        length(pathways_sub)))
    }

    # Stop if no pathways left after filtering
    if (!length(pathways_sub)) {
      stop(paste("No pathways left after filtering. Please alter filtering",
                 "arguments specified."))
    }
  }
  # Return subsetted list of pathways
  return(pathways_sub)
}
