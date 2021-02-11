inputObject <- function(test_gene, background_gene, pathways) {

    pathway_genes <- unique(as.character(unlist(pathways)))
    test_gene_in_pathways <- which(test_gene %in% pathway_genes)
    back_gene_in_pathways <- which(background_gene %in% pathway_genes)

    if (is.null(test_gene)) {
        stop("Oops, argument 'test_gene' is empty. Supply a vector of
         genes ... I promise this will work.")
    } else if (is.null(background_gene)) {
        stop("Oops, argument 'background_gene' is empty. Supply a vector of
        genes ... I promise this will work")
    } else if (!is.list(pathways)) {
        stop("Oops, argument 'pathways' is not in a list format...
        have you tried using readPathways() on your input pathway file?")
    } else if (!length(test_gene_in_pathways)) {
        stop("Oops, none of the genes in 'test_gene' was found in 'pathways'.
        Make sure that you have some gene IDs in both inputs, otherwise how do
        you expect this works?")
    } else if (!length(back_gene_in_pathways)) {
        stop("Oops, none of the genes in 'background_genes' was found in
        'pathways'. Make sure that you have some gene IDs in both inputs,
        otherwise how do you expect this works?")
    } else if (length(test_gene) >= length(background_gene)) {
        stop("Oops, your test set can't have more genes than your background
        set. Have you mixed up the 'test_gene' and 'background_gene' arguments?
        You're so close... I can feel it.")
    }

    test_gene <- unique(as.character(test_gene))
    background_gene <- unique(as.character(background_gene))

    list(test_gene = test_gene,
        background_gene = background_gene,
        pathways = pathways,
        pathways_name = names(pathways),
        pathways_size = unlist(lapply(pathways, length))
    )
}

#' Runs gene enrichment and depletion analysis for a list of pathways.
#'
#' @param test_gene (char) vector of genes to use as test set.
#' @param background_gene (char) vector of genes to use as background set.
#' @param pathways (list) list of vectors with pathway annotations.
#' @return table of pathway enrichment and depletion results. Rows represent
#' tested pathways. Columns represent:
#' \itemize{
#'     \item pathway -- name of the pathway, corresponds to
#'         names(\code{pathways});
#'     \item size -- size of the pathway;
#'     \item real_frac -- fraction of \code{test_gene} members in pathway;
#'     \item expected_frac -- fraction of \code{background_gene} members in
#'         pathway;
#'     \item fold_enrichment -- fold enrichment measure,
#'         evaluates as \code{real_frac} / \code{expected_frac};
#'     \item status -- indicator that pathway is enriched or depleted for
#'         \code{test_gene} members;
#'     \item real_gene -- vector of \code{test_gene} gene members annotated
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
runFedup <- function(test_gene, background_gene, pathways) {
    inputs <- inputObject(test_gene, background_gene, pathways)
    test <- inputs$test_gene
    background <- inputs$background_gene
    pathways <- inputs$pathways
    pathways_name <- inputs$pathways_name
    pathways_size <- inputs$pathways_size
    message("Data input:\n => ",
        length(test), " test genes\n => ",
        length(background), " background genes\n => ",
        length(pathways), " pathawys")

    res <- data.table(pathway = pathways_name, size = pathways_size)
    res_stats <- vapply(pathways, function(x) {
        a_n <- length(test) # n test genes
        b_n <- length(background) # n background genes
        a <- intersect(test, x) # test genes in pathway
        b <- intersect(background, x) # background genes in pathway
        a_len <- length(a) # n test genes in pathway
        b_len <- length(b) # n background genes in pathway
        a_x <- (a_len / a_n) * 100 # fraction of test genes in pathway
        b_x <- (b_len / b_n) * 100 # fraction of background genes in pathway
        f <- a_x / b_x # fold enrichment measure
        e <- ifelse(f > 1, "Enriched", "Depleted")
        m <- rbind(c(a_len, b_len), c(a_n, b_n)) # pval contingency table
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
    res <- res[order(res$pvalue),] # BH-correct pvalues
    res$qvalue <- p.adjust(res$pvalue, method = "BH")
    message("You did it! FEDUP ran successfully, feeling pretty good huh?")
    return(res)
}
