inputObject <- function(testGene, backgroundGene, pathways) {

    pathway_genes <- unique(as.character(unlist(pathways)))
    testGene_in_pathways <- which(testGene %in% pathway_genes)
    back_gene_in_pathways <- which(backgroundGene %in% pathway_genes)

    if (is.null(testGene)) {
        stop("Oops, argument 'testGene' is empty. Supply a vector of
        genes... I promise this will work.")
    }
    if (is.null(backgroundGene)) {
        stop("Oops, argument 'backgroundGene' is empty. Supply a vector of
        genes... I promise this will work.")
    }
    if (!is.list(pathways)) {
        stop("Oops, argument 'pathways' is not in a list format...
        have you tried using readPathways() on your input pathway file?")
    }
    if (!length(testGene_in_pathways)) {
        stop("Oops, none of the genes in 'testGene' was found in 'pathways'.
        Make sure that you have some gene IDs in both inputs, otherwise how do
        you expect this works?")
    }
    if (!length(back_gene_in_pathways)) {
        stop("Oops, none of the genes in 'backgroundGenes' was found in
        'pathways'. Make sure that you have some gene IDs in both inputs,
        otherwise how do you expect this works?")
    }
    if (length(testGene) >= length(backgroundGene)) {
        stop("Oops, your test set can't have more genes than your background
        set. Have you mixed up the 'testGene' and 'backgroundGene' arguments?
        You're so close... I can feel it.")
    }

    testGene <- unique(as.character(testGene))
    backgroundGene <- unique(as.character(backgroundGene))

    list(testGene = testGene,
        backgroundGene = backgroundGene,
        pathways = pathways,
        pathways_name = names(pathways),
        pathways_size = unlist(lapply(pathways, length))
    )
}

#' Runs gene enrichment and depletion analysis for a list of pathways.
#'
#' @param testGene (char) vector of genes to use as test set.
#' @param backgroundGene (char) vector of genes to use as background set.
#' @param pathways (list) list of vectors with pathway annotations.
#' @return table of pathway enrichment and depletion results. Rows represent
#' tested pathways. Columns represent:
#' \itemize{
#'     \item pathway -- name of the pathway, corresponds to
#'         names(\code{pathways});
#'     \item size -- size of the pathway;
#'     \item real_frac -- fraction of \code{testGene} members in pathway;
#'     \item expected_frac -- fraction of \code{backgroundGene} members in
#'         pathway;
#'     \item fold_enrichment -- fold enrichment measure,
#'         evaluates as \code{real_frac} / \code{expected_frac};
#'     \item status -- indicator that pathway is enriched or depleted for
#'         \code{testGene} members;
#'     \item real_gene -- vector of \code{testGene} gene members annotated
#'         to \code{pathways};
#'     \item pvalue -- enrichment p-value calculated via Fisher's exact test;
#'     \item qvalue -- BH-adjusted p-value
#' }
#' @examples
#' data(testGene)
#' data(backgroundGene)
#' data(pathwaysGMT)
#' fedup_res <- runFedup(testGene, backgroundGene, pathwaysGMT)
#' @importFrom data.table data.table :=
#' @importFrom utils head read.delim tail
#' @importFrom stats fisher.test p.adjust
#' @export
runFedup <- function(testGene, backgroundGene, pathways) {
    inputs <- inputObject(testGene, backgroundGene, pathways)
    test <- inputs$testGene
    background <- inputs$backgroundGene
    pathways <- inputs$pathways
    pathways_name <- inputs$pathways_name
    pathways_size <- inputs$pathways_size
    message("Data input:\n => ",
        length(test), " test genes\n => ",
        length(background), " background genes\n => ",
        length(pathways), " pathawys")

    res <- data.table(pathway = pathways_name, size = pathways_size)
    res_stats <- vapply(pathways, function(x) {
        a_n <- length(test)
        b_n <- length(background)
        a <- intersect(test, x)
        b <- intersect(background, x)
        a_len <- length(a)
        b_len <- length(b)
        a_x <- (a_len / a_n) * 100
        b_x <- (b_len / b_n) * 100
        f <- a_x / b_x
        e <- ifelse(f > 1, "Enriched", "Depleted")
        m <- rbind(c(a_len, b_len), c(a_n, b_n))
        p <- fisher.test(m, alternative = "two.sided")$p.value
        return(c(
            real_frac = a_x, expected_frac = b_x, fold_enrich = f,
            status = e, pvalue = p, real_gene = paste(a, collapse = "|")))
    }, character(6))

    res[, "real_frac" := as.numeric(unlist(res_stats["real_frac",]))]
    res[, "expected_frac" := as.numeric(unlist(res_stats["expected_frac",]))]
    res[, "fold_enrichment" := as.numeric(unlist(res_stats["fold_enrich",]))]
    res[, "status" := unlist(res_stats["status",])]
    res[, "real_gene" := mapply("[", strsplit(res_stats["real_gene",], "\\|"))]
    res[, "pvalue" := as.numeric(unlist(res_stats["pvalue",]))]
    res <- res[order(res$pvalue),]
    res$qvalue <- p.adjust(res$pvalue, method = "BH")

    message("You did it! FEDUP ran successfully, feeling pretty good huh?")
    return(res)
}
