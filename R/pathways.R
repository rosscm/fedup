#' Returns a list of pathways from various file formats.
#' Currently supports the following file format: gmt, txt, xlsx.
#'
#' @param pathwayFile (char) path to file with pathway annotations.
#' @param header (logical) whether \code{pathwayFile} has a header
#'     (default FALSE).
#' @param pathCol (char or int) column name or number with pathway identifiers.
#'     For use with non-GMT input files (eg "Pathway.ID" or 2; default NULL).
#' @param geneCol (char or int) column name or number with gene identifiers.
#'     For use with non-GMT input files (eg "Gene.ID" or 5; default NULL).
#' @param minGene (integer) minimum number of genes to be considered
#'     in a pathway (default 1).
#' @param maxGene (integer) maximum number of genes to be considered
#'     in a pathway (default Inf).
#' @return a list of vectors with pathway annotations.
#' @examples
#' pathways <- readPathways(
#'     system.file("extdata", "Human_Reactome_November_17_2020_symbol.gmt",
#'         package = "fedup"
#'     ),
#'     minGene = 10, maxGene = 500
#' )
#' pathways <- readPathways(
#'     system.file("extdata", "SAFE_terms.xlsx", package = "fedup"),
#'     header = TRUE, pathCol = "Enriched.GO.names", geneCol = "Gene.ID"
#' )
#' @importFrom openxlsx read.xlsx
#' @importFrom tibble deframe
#' @importFrom stats aggregate na.omit
#' @importFrom utils head read.delim tail
#' @export
readPathways <- function(pathwayFile, header = FALSE,
    pathCol = NULL, geneCol = NULL, minGene = 1L, maxGene = Inf) {
    s <- c("gmt", "txt", "xlsx")
    f <- sub(".*\\.", "", pathwayFile)
    if (!f %in% s) {
        stop(paste0(
            "Sorry, pathway file type (", f, ") is not supported. ",
            "Supported extensions: ", paste(s, collapse = ", "), "."
        ))
    }
    if (f == "gmt") {
        path_in <- strsplit(readLines(pathwayFile), "\t")
        if (header) {
            path_in <- path_in[-1]
        }
        pathways <- lapply(path_in, tail, -2)
        names(pathways) <- vapply(path_in, head, n = 1, character(1))
    } else {
        if (f == "xlsx") {
            path_in <- read.xlsx(pathwayFile)
        }
        if (f == "txt") {
            path_in <- read.delim(pathwayFile, header = header)
        }
        if (!pathCol %in% names(path_in) && !pathCol %in% seq_along(path_in)) {
            stop("Pathway ID column (", pathCol, ") not in file")
        }
        if (!geneCol %in% names(path_in) && !geneCol %in% seq_along(path_in)) {
            stop("Gene ID column (", geneCol, ") not in file")
        }
        pathway_df <- data.frame(
            pathway = path_in[, pathCol], gene = path_in[, geneCol]
        )
        pathway_df[which(pathway_df$gene == ""), "gene"] <- NA
        pathway_df <- na.omit(pathway_df)
        pathway_df <- aggregate(gene ~ pathway, pathway_df, paste)
        pathways <- deframe(pathway_df)
    }
    size <- lapply(pathways, length)
    pathways_s <- pathways[which(size >= minGene & size <= maxGene)]
    if (!length(pathways_s)) {
        stop("Oops, no pathways left... try different filtering options.")
    }
    message(
        "Pathway file: ", basename(pathwayFile),
        "\n => n total pathways: ", length(pathways),
        "\n => n pathways (", minGene, "-", maxGene, "): ", length(pathways_s)
    )
    return(pathways_s)
}

#' Writes a set of pathways (list of vectors) to a GMT file.
#'
#' @param pathways (list) named list of vectors.
#' @param gmtFile (char) name of output GMT file.
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
writePathways <- function(pathways, gmtFile) {
    tab <- data.table(
        pathway = names(pathways),
        description = gsub("\\%.*", "", names(pathways)),
        genes = unlist(lapply(pathways, paste, collapse = "\t"))
    )
    fwrite(tab, file = gmtFile, sep = "\t", col.names = FALSE, quote = FALSE)
    message("Wrote out pathway gmt file to ", gmtFile)
}
