inputObject <- function(genes, pathways) {
    if (!is.list(genes)) {
        stop("Oops, 'genes' is not in a list format... have you tried using
        using prepInput() on your input genes that you want to test?")
    }
    if (!is.list(pathways)) {
        stop("Oops, 'pathways' is not in a list format... have you tried using
        using readPathways() on your pathway annotation file?")
    }
    if (length(which(names(genes) == "background")) > 1) {
        stop("Input gene list must include only one common background set.")
    }

    pathway_gene <- unique(as.character(unlist(pathways)))
    back_genes <- genes[["background"]]
    genes[["background"]] <- NULL
    test_genes <- genes
    test_genes_all <- as.character(unlist(test_genes))

    if (!length(back_genes)) {
        stop("Oops, no background genes provided in 'genes' argument! Check
        data(geneSingle) for an example of input structure.")
    }
    if (!length(test_genes)) {
        stop("Oops, no test genes provided in 'genes' argument! Check
        data(geneSingle) for an example of input structure.")
    }
    if (length(which(lapply(test_genes, length) > length(back_genes)))) {
        n <- length(which(lapply(test_genes, length) > length(back_genes)))
        stop("Oops, ", n, " test gene sets are larger than the background set.
        Please double check your input gene list.")
    }
    if (length(which(duplicated(test_genes_all)))) {
        n <- signif((length(which(duplicated(test_genes_all))) /
            length(unique(test_genes_all))) * 100)
        warning(n, "% of genes overlap across your test gene sets (enrichment
        results may be similar across tests).")
    }

    list(
        testGenes = test_genes,
        testNames = names(test_genes),
        testLength = unlist(lapply(test_genes, length)),
        backgroundGenes = back_genes,
        pathways = pathways,
        pathways_name = names(pathways),
        pathways_size = unlist(lapply(pathways, length))
    )
}

#' Runs pathway enrichment and depletion analysis using a Fisher's
#' exact test.
#'
#' @description This function takes a list of test genes and a common background
#' set to calculate enrichment and depletion for a list of pathways. The method
#' allows for fast and efficient testing of multiple gene sets of interest.
#' @param genes (list) named list of vectors with background genes and \code{n}
#' test genes.
#' @param pathways (list) named list of vectors with pathway annotations.
#' @return List of length \code{n} with table(s) of pathway enrichment and
#' depletion results. Rows represent tested pathways. Columns represent:
#' \itemize{
#'     \item pathway -- name of the pathway, corresponds to
#'         names(\code{pathways});
#'     \item size -- size of the pathway;
#'     \item real_frac -- fraction of test gene members in pathway;
#'     \item expected_frac -- fraction of background gene members in
#'         pathway;
#'     \item fold_enrichment -- fold enrichment measure,
#'         evaluates as \code{real_frac} / \code{expected_frac};
#'     \item status -- indicator that pathway is enriched or depleted for
#'         test gene members;
#'     \item real_gene -- vector of test gene members annotated
#'         to \code{pathways};
#'     \item pvalue -- enrichment p-value calculated via Fisher's exact test;
#'     \item qvalue -- BH-adjusted p-value
#' }
#' @examples
#' data(pathwaysGMT)
#' # Run analysis with a single test set
#' data(geneSingle)
#' fedupRes <- runFedup(geneSingle, pathwaysGMT)
#' # Run analysis with two test sets
#' data(geneDouble)
#' fedupRes <- runFedup(geneDouble, pathwaysGMT)
#' # Run analysis with multiple test sets
#' data(geneMulti)
#' fedupRes <- runFedup(geneMulti, pathwaysGMT)
#' @importFrom data.table data.table :=
#' @importFrom utils head read.delim tail
#' @importFrom stats fisher.test p.adjust
#' @export
runFedup <- function(genes, pathways) {
    inputs <- inputObject(genes, pathways)
    test <- inputs$testGenes
    test_name <- inputs$testName
    test_length <- inputs$testLength
    background <- inputs$backgroundGenes
    pathways <- inputs$pathways
    pathways_name <- inputs$pathways_name
    pathways_size <- inputs$pathways_size
    message(
        "Running fedup with:",
        "\n => ", length(test), " test set(s)",
        paste0("\n  + ", test_name, ": ", test_length, " genes"),
        "\n => ", length(background), " background genes",
        "\n => ", length(pathways), " pathway annotations"
    )
    res_list <- list()
    for (i in test_name) {
        res <- data.table(pathway = pathways_name, size = pathways_size)
        res_stats <- vapply(pathways, function(x) {
            a_n <- length(test[[i]])
            b_n <- length(background)
            a <- intersect(test[[i]], x)
            b <- intersect(background, x)
            a_len <- length(a)
            b_len <- length(b)
            a_x <- (a_len / a_n) * 100
            b_x <- (b_len / b_n) * 100
            f <- ifelse(a_x == 0 && b_x == 0, 0, a_x / b_x)
            e <- ifelse(f >= 1, "enriched", "depleted")
            m <- rbind(c(a_len, b_len), c(a_n, b_n))
            p <- fisher.test(m, alternative = "two.sided")$p.value
            return(c(
                a_x = a_x, b_x = b_x, f = f, e = e, p = p,
                r = paste(a, collapse = "|")
            ))
        }, character(6))
        res[, "real_frac" := as.numeric(unlist(res_stats["a_x", ]))]
        res[, "expected_frac" := as.numeric(unlist(res_stats["b_x", ]))]
        res[, "fold_enrichment" := as.numeric(unlist(res_stats["f", ]))]
        res[, "status" := unlist(res_stats["e", ])]
        res[, "real_gene" := mapply("[", strsplit(res_stats["r", ], "\\|"))]
        res[, "pvalue" := as.numeric(unlist(res_stats["p", ]))]
        res <- res[order(res$pvalue), ]
        res$qvalue <- p.adjust(res$pvalue, method = "BH")
        res_list[[i]] <- res
    }
    message("All done!")
    return(res_list)
}
