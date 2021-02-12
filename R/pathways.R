#' Returns a list of pathways from various file formats.
#' Currently supports the following file format: gmt, txt, xlsx.
#'
#' @param pathwayFile (char) path to file with pathway annotations.
#' @param header (logical) whether \code{pathwayFile} has a header
#'     (default FALSE).
#' @param pathwayCol (char) column name with pathway identifiers.
#'     For use with non-GMT input files (eg "Pathway.ID"; default NULL).
#' @param geneCol (char) column name with gene identifiers.
#'     For use with non-GMT input files (eg "Gene.ID"; default NULL).
#' @param minGene (integer) minimum number of genes to be considered
#'     in a pathway (default 1).
#' @param maxGene (integer) maximum number of genes to be considered
#'     in a pathway (default Inf).
#' @return a list of vectors with pathway annotations.
#' @examples
#' pathways <- readPathways(
#'     system.file("extdata", "Human_Reactome_November_17_2020_symbol.gmt",
#'     package="FEDUP"), minGene=10, maxGene=500)
#' pathways <- readPathways(
#'     system.file("extdata", "SAFE_terms.xlsx", package="FEDUP"),
#'     header=TRUE, pathwayCol="Enriched.GO.names", geneCol="Gene.ID")
#' @importFrom openxlsx read.xlsx
#' @importFrom tibble deframe
#' @importFrom stats aggregate na.omit
#' @importFrom utils head read.delim tail
#' @export
readPathways <- function(pathwayFile, header=FALSE,
                        pathwayCol=NULL, geneCol=NULL,
                        minGene=1L, maxGene=Inf) {

    s <- c("gmt", "txt", "xlsx")
    f <- sub(".*\\.", "", pathwayFile)
    if (!f %in% s) {
        stop(paste0("Sorry, pathway file type (", f, ") is not supported. ",
            "Supported extensions: ", paste(s, collapse = ", "), "."))
    }
    if (f == "gmt") {
        pathway_in <- strsplit(readLines(pathwayFile), "\t")
        if (header) { pathway_in <- pathway_in[-1] }
        pathways <- lapply(pathway_in, tail, -2)
        names(pathways) <- vapply(pathway_in, head, n = 1, character(1))
    } else {
        if (f == "xlsx") {
            pathway_in <- read.xlsx(pathwayFile)
        } else if (f == "txt") {
            pathway_in <- read.delim(pathwayFile, header = header)
        }
        if (missing(pathwayCol)||!pathwayCol %in% colnames(pathway_in)) {
            stop("Pathway ID column (", pathwayCol, ") not in file")
        } else if (missing(geneCol)||!geneCol %in% colnames(pathway_in)) {
            stop("Gene ID column (", geneCol, ") not in file")
        } else {
            pathway_df <- data.frame(
                pathway = pathway_in[,pathwayCol],
                gene = pathway_in[,geneCol])
            pathway_df[which(pathway_df$gene == ""), "gene"] <- NA
            pathway_df <- na.omit(pathway_df)
            pathway_df <- aggregate(gene ~ pathway, pathway_df, paste)
            pathways <- deframe(pathway_df)
        }
    }

    size <- lapply(pathways, length)
    pathways_s <- pathways[which(size >= minGene & size <= maxGene)]
    pathways_s <- pathways_s[!duplicated(names(pathways_s))]
    if (!length(pathways_s)) {
        stop("Oops, no pathways left... try different filtering options.")
    }
    message("Pathway file: ", basename(pathwayFile),
        "\n => n total pathways: ", length(pathways),
        "\n => n pathways (",minGene,"-",maxGene, "): ", length(pathways_s))

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
#' writePathways(pathwaysXLSX, tempfile("pathwaysXLSX", fileext=".gmt"))
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
