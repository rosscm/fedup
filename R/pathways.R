#' Returns a list of pathways from various file formats.
#' Currently supports the following file format: GMT, TXT, XLSX.
#'
#' @param pathway_file (char) path to file with pathway annotations.
#' @param header (logical) whether pathway_file has a header (default TRUE).
#' @param pathway_col (char) column name with pathway identifiers.
#'  For use with for non-GMT input files (eg "Pathway.ID"; default NULL).
#' @param gene_col (char) column name with gene identifiers.
#'  For use with for non-GMT input files (eg "Gene.ID"; default NULL).
#' @param MIN_GENE (integer) minimum number of genes to be considered
#'  in a pathway (default = 1).
#' @param MAX_GENE (integer) maximum number of genes to be considered
#'  in a pathway (default = Inf).
#' @param GO_class (char) if GO terms are present in `pathway_file`, which
#'  class(es) should be subsetted (BP, MF and/or CC, default NULL).
#' @return a list of vectors with pathway annotations.
#' @examples
#' pathways <- readPathways(
#'  system.file("extdata", "Human_Reactome_November_17_2020_symbol.gmt",
#'  package = "FEDUP"), MIN_GENE = 10, MAX_GENE = 500)
#' pathways <- readPathways(
#'  system.file("extdata", "SAFE_terms.xlsx", package = "FEDUP"),
#'  header = TRUE, pathway_col = "Enriched.GO.names", gene_col = "Gene.ID",
#'  MIN_GENE = 10, MAX_GENE = 500)
#' @importFrom openxlsx read.xlsx
#' @importFrom tibble deframe
#' @importFrom stats aggregate na.omit
#' @importFrom utils head read.delim tail
#' @export
readPathways <- function(pathway_file, header = TRUE,
                        pathway_col = NULL, gene_col = NULL,
                        MIN_GENE = 1L, MAX_GENE = Inf,
                        GO_class = NULL) {
    # Vector of supported file types
    supported_types <- c("gmt", "txt", "xlsx")
    file_type <- sub(".*\\.", "", pathway_file) # grab pathway_file extension

    # Stop if pathway_file input is not supported
    if (!file_type %in% supported_types) {
    stop(paste0("Sorry, pathway file type (", file_type, ") is not supported. ",
                "Please use one of the following formats: ",
                paste(supported_types, collapse = ", "), "."))
    } else {
        message("Pathway file: ", basename(pathway_file))
        if (file_type == "gmt") {
            pathway_in <- strsplit(readLines(pathway_file), "\t")
            if (header) { pathway_in <- pathway_in[-1] } # strip header if TRUE
            pathways <- lapply(pathway_in, tail, -2)
            names(pathways) <- sapply(pathway_in, head, 1)
        } else {
            if (file_type == "xlsx") {
                pathway_in <- read.xlsx(pathway_file)
            }
            if (file_type == "txt")  {
                pathway_in <- read.delim(pathway_file, header = header)
            }
            if (missing(pathway_col) || !pathway_col %in% colnames(pathway_in)) {
                stop("Pathway ID column (", pathway_col, ") not in pathway file")
            }
            if (missing(gene_col) || !gene_col %in% colnames(pathway_in)) {
                stop("Gene ID column (", gene_col, ") not in pathway file")
            }

            # Generate pathway list from dataframe
            pathway_df <- data.frame(
                pathway = pathway_in[,pathway_col],
                gene = pathway_in[,gene_col]
            )

            # Ensure no empty entries
            pathway_df[which(pathway_df$gene == ""), "gene"] <- NA
            pathway_df <- na.omit(pathway_df)

            # Transform to list
            pathway_df <- aggregate(gene ~ pathway, data = pathway_df, FUN = paste)
            pathways <- deframe(pathway_df)
        }

        # Subset for pathways within MIN_GENE, MAX_GENE range
        pathway_size <- lapply(pathways, length)
        pathways_sub <- pathways[which(pathway_size >= MIN_GENE & pathway_size <= MAX_GENE)]
        message(" => n total pathway entries: ", length(pathways))
        message(" => n pathway entries (min ", MIN_GENE,  ", max ",
                MAX_GENE, "): ", length(pathways_sub))

        # Subset for chosen 'GO_class' if specified
        if (!is.null(GO_class)) {
            go_keep <- paste(GO_class, collapse = "|")
            int_keep <- grep(go_keep, names(pathways_sub))
            pathways_sub <- pathways_sub[int_keep]
            message(" => n pathways filtering for ", GO_class, "GO class: ",
                length(pathways_sub))
        }

        # Finally, check for any duplicated pathway annotations
        pathways_dup <- which(duplicated(names(pathways_sub)))
        if (length(pathways_dup)) {
            pathways_sub <- pathways_sub[-pathways_dup]
            message(" => n pathways after removing duplicated entries: ",
                length(pathways_sub))
        }

        # Stop if no pathways left after filtering
        if (!length(pathways_sub)) {
            stop("Oops no pathways left... try out different filtering options.")
        }
    }
    # Return named list of pathways
    return(pathways_sub)
}

#' Writes a set of pathways (list of vectors) to a GMT file.
#'
#' @param pathways (list) named list of vectors.
#' @param gmt_file (char) name of output GMT file.
#' @return GMT-formatted file. Rows represent pathways. Columns represent:
#' \itemize{
#'  \item pathway ID;
#'  \item description;
#'  \item a list of tab-delimited genes
#' }
#' @examples
#' data(pathwaysXLSX)
#' writePathways(pathwaysXLSX, tempfile("pathwaysXLSX", fileext = ".gmt"))
#' @importFrom data.table fwrite
#' @export
writePathways <- function(pathways, gmt_file) {
    tab <- data.table(
        pathway = names(pathways),
        description = gsub("\\%.*", "", names(pathways)),
        genes = unlist(lapply(pathways, paste, collapse = "\t"))
    )
    fwrite(tab, file = gmt_file, sep = "\t", col.names = FALSE, quote = FALSE)
    message("Wrote out GMT file with to ", gmt_file)
}
