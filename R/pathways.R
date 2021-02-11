#' Returns a list of pathways from various file formats.
#' Currently supports the following file format: GMT, TXT, XLSX.
#'
#' @param pathway_file (char) path to file with pathway annotations.
#' @param header (logical) whether \code{pathway_file} has a header
#'     (default FALSE).
#' @param pathway_col (char) column name with pathway identifiers.
#'     For use with non-GMT input files (eg "Pathway.ID"; default NULL).
#' @param gene_col (char) column name with gene identifiers.
#'     For use with non-GMT input files (eg "Gene.ID"; default NULL).
#' @param MIN_GENE (integer) minimum number of genes to be considered
#'     in a pathway (default = 1).
#' @param MAX_GENE (integer) maximum number of genes to be considered
#'     in a pathway (default = Inf).
#' @return a list of vectors with pathway annotations.
#' @examples
#' pathways <- readPathways(
#'     system.file("extdata", "Human_Reactome_November_17_2020_symbol.gmt",
#'     package = "FEDUP"), MIN_GENE = 10, MAX_GENE = 500)
#' pathways <- readPathways(
#'     system.file("extdata", "SAFE_terms.xlsx", package = "FEDUP"),
#'     header = TRUE, pathway_col = "Enriched.GO.names", gene_col = "Gene.ID")
#' @importFrom openxlsx read.xlsx
#' @importFrom tibble deframe
#' @importFrom stats aggregate na.omit
#' @importFrom utils head read.delim tail
#' @export
readPathways <- function(pathway_file, header = FALSE,
                        pathway_col = NULL, gene_col = NULL,
                        MIN_GENE = 1L, MAX_GENE = Inf, GO_class = NULL) {

    message("Pathway file: ", basename(pathway_file))
    s <- c("gmt", "txt", "xlsx") # supported file extensions
    f <- sub(".*\\.", "", pathway_file) # pathway_file extension
    if (!f %in% s) {
        stop(paste0("Sorry, pathway file type (", f, ") is not supported. ",
            "Supported extensions: ", paste(s, collapse = ", "), "."))
    }
    if (f == "gmt") {
        pathway_in <- strsplit(readLines(pathway_file), "\t")
        if (header) { pathway_in <- pathway_in[-1] }
        pathways <- lapply(pathway_in, tail, -2)
        names(pathways) <- vapply(pathway_in, head, n = 1, character(1))
    } else {
        if (f == "xlsx") {
            pathway_in <- read.xlsx(pathway_file)
        } else if (f == "txt") {
            pathway_in <- read.delim(pathway_file, header = header)
        }
        if (missing(pathway_col)||!pathway_col %in% colnames(pathway_in)) {
            stop("Pathway ID column (", pathway_col, ") not in file")
        } else if (missing(gene_col)||!gene_col %in% colnames(pathway_in)) {
            stop("Gene ID column (", gene_col, ") not in file")
        } else {
            pathway_df <- data.frame(
                pathway = pathway_in[,pathway_col],
                gene = pathway_in[,gene_col])
            pathway_df[which(pathway_df$gene == ""), "gene"] <- NA
            pathway_df <- na.omit(pathway_df) # ensure no NaNs
            pathway_df <- aggregate(gene ~ pathway, pathway_df, paste)
            pathways <- deframe(pathway_df) # transform df to list
        }
    }

    size <- lapply(pathways, length) # subset for pathways in [MIN:MAX] range
    pathways_s <- pathways[which(size >= MIN_GENE & size <= MAX_GENE)]
    pathways_s <- pathways_s[!duplicated(names(pathways_s))]
    message(" => n total pathways: ", length(pathways))
    message(" => n pathways (",MIN_GENE,"-",MAX_GENE, "): ", length(pathways_s))

    if (!length(pathways_s)) {
        stop("Oops, no pathways left... try different filtering options.")
    }
    return(pathways_s)
}

#' Writes a set of pathways (list of vectors) to a GMT file.
#'
#' @param pathways (list) named list of vectors.
#' @param gmt_file (char) name of output GMT file.
#' @return GMT-formatted file. Rows represent pathways. Columns represent:
#' \itemize{
#'     \item pathway ID;
#'     \item description;
#'     \item a list of tab-delimited genes
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
