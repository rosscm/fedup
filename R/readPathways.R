#' Returns a list of pathways from a file
#' Currently supports the following file format: GMT, TXT, XLSX
#' @param pathway_file (char) path to file with pathway annotations
#' @param header (logical) whether pathway_file has a header (default = TRUE)
#' @param pathway_col (char) column name with pathway identifiers, one pathway per row
#'  use for non-GMT input (eg "Pathway.ID", default = NULL)
#' @param gene_col (char) column name with gene identifiers, one gene per row
#'  use for non-GMT input (eg "Gene.ID", default = NULL)
#' @param MIN_GENE (integer) minimum number of genes to be considered in a pathway
#'  (default = 10)
#' @param MAX_GENE (integer) maximum number of genes to be considered in a pathway
#'  (default = 500)
#' @return a list of vectors with pathways
#' @examples
#' pathways <- readPathways(
#'  system.file("extdata", "YeastDatabase_GO_gmt.gmt", package = "FEDUP")
#')
#' pathways <- readPathways(
#'  system.file("extdata", "SAFE_terms.xlsx", package = "FEDUP"),
#'  header = TRUE, pathway_col = "Enriched.GO.names", gene_col = "ORF.ID"
#')
#' @export
readPathways <- function(pathway_file, header = TRUE, pathway_col, gene_col,
                         MIN_GENE = 10L, MAX_GENE = 500L) {
  supported_types <- c("gmt", "txt", "xlsx")
  file_type <- sub(".*\\.", "", pathway_file) # grab pathway_file extension
  # Stop if pathway_file type is not supported
  if (!file_type %in% supported_types) {
    stop(paste0("Sorry, input pathway file type (", file_type, ") is not supported. ",
                "Please use one of the following formats: ",
                paste(supported_types, collapse = ", "), "."))
  } else {
    if (file_type == "gmt") {
      pathway_in <- strsplit(readLines(pathway_file), "\t")
      if (header) { pathway_in <- pathway_in[-1] } # remove header if TRUE
      pathways <- lapply(pathway_in, tail, -2)
      names(pathways) <- sapply(pathway_in, head, 1)
    } else {
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
      pathways <- deframe(pathway_df)
    }
    # Subset for pathways within MIN_GENE, MAX_GENE
    pathway_size <- lapply(pathways, length)
    pathways_sub <- pathways[which(pathway_size >= MIN_GENE & pathway_size <= MAX_GENE)]
    # Print messages to console
    cat(sprintf(
     "* Input pathway file: %s
      * number of pathways: %s
      * final number after pathway size filtering (%s, %s): %s\n",
     basename(pathway_file), length(pathways), MIN_GENE, MAX_GENE, length(pathways_sub)
    ))
  }
  return(pathways_sub)
}

pathways <- readPathways("inst/extdata/YeastDatabase_GO_gmt.gmt")
pathways <- readPathways("inst/extdata/SAFE_terms.xlsx", pathway_col = "Enriched.GO.names", gene_col = "ORF.ID")
