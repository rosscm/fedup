inputObject <- function(test_gene, background_gene, pathways) {
    # Stop if test_gene is empty
    if (is.null(test_gene)) {
        stop("Oops, argument 'test_gene' is empty. Supply a vector of genes
        ...I promise this will work.")
    }

    # Stop if background_gene is empty
    if (is.null(background_gene)) {
        stop("Oops, argument 'background_gene' is empty. Supply a vector of
        genes ... I promise this will work")
    }

    # Stop if pathways is empty or not in list format
    if (!is.list(pathways)) {
        stop("Oops, argument 'pathways' is not in a list format...
        have you tried using readPathways()?")
    }

    # Stop if length of test_gene is greater than background_gene
    if (length(test_gene) >= length(background_gene)) {
        stop("Oops, your test set can't have more genes than your background
        set. Have you mixed up the 'test_gene' and 'background_gene' arguments?
        You're so close... I can feel it.")
    }

    # Vector of all unique genes in pathways
    pathway_genes <- unique(as.character(unlist(pathways)))

    # Stop if no genes from test_gene found in pathways
    test_gene_in_pathways <- which(test_gene %in% pathway_genes)
    if (!length(test_gene_in_pathways)) {
        stop("Oops, none of the genes in 'test_gene' was found in 'pathways'.
        Make sure that you have some gene IDs in both inputs, otherwise how do
        you expect this works?")
    }

    # Stop if no genes from background_gene found in pathways
    background_gene_in_pathways <- length(which(background_gene %in% pathway_genes))
    if (!length(background_gene_in_pathways)) {
        stop("Oops, none of the genes in 'background_genes' was found in
        'pathways'. Make sure that you have some gene IDs in both inputs,
        otherwise how do you expect this works?")
    }

    # Prepare input arguments
    test_gene <- unique(as.character(test_gene))
    background_gene <- unique(as.character(background_gene))

    # Store in list
    list(
        test_gene = test_gene,
        background_gene = background_gene,
        pathways = pathways,
        pathways_name = names(pathways),
        pathways_size = unlist(lapply(pathways, length))
    )
}

#' Runs gene enrichment and depletion analysis for a list of pathways.
#'
#' @param test_gene (char) vector of genes to test set.
#' @param background_gene (char) vector of genes to use as background set.
#' @param pathways (list) list of vectors with pathway annotations.
#' @return table of pathway enrichment and depletion results. Rows represent
#' tested pathways. Columns represent:
#' \itemize{
#'  \item pathway -- name of the pathway, corresponds to names(pathways);
#'  \item size -- size of the pathway;
#'  \item real_total -- number of genes in `test_gene`;
#'  \item real_pathway -- number of `test_gene` members in pathway;
#'  \item real_pathway_frac -- fraction of `test_gene` members in pathway,
#'      evaluates as (`real_pathawy` / `real_total`) * 100;
#'  \item expected_total -- number of genes in 'background_gene';
#'  \item expected_pathway -- number of `background_gene` members in pathway;
#'  \item expected_pathway_frac -- fraction of `background_gene` members in
#'      pathway, evaluates as (`expected_pathway` / `expected_total`) * 100;
#'  \item enrichment -- indicator that pathway is enriched or depleted for
#'      `test_gene` members;
#'  \item real_pathway_gene -- vector of `real_pathway` genes;
#'  \item pvalue -- enrichment p-value calculated via Fisher's exact test;
#'  \item qvalue -- BH-adjusted p-value
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
    # Get input arguments
    inputs <- inputObject(test_gene, background_gene, pathways)
    test <- inputs$test_gene
    background <- inputs$background_gene
    pathways <- inputs$pathways
    pathways_name <- inputs$pathways_name
    pathways_size <- inputs$pathways_size
    message("Data input:\n => ",
        length(test), " test genes\n => ",
        length(background), " background genes\n => ",
        length(pathways), " pathawys"
    )

    if (!length(pathways)) {
    return(data.table(
            pathway = character(),
            size = integer(),
            real_total = numeric(),
            real_pathway = numeric(),
            real_pathway_frac = numeric(),
            expected_total = numeric(),
            expected_pathway = numeric(),
            expected_pathway_frac = numeric(),
            enrichment = character(),
            real_pathway_gene = list(),
            pvalue = numeric(),
            qvalue = numeric()))
    }

    # Results table
    res <- data.table(pathway = pathways_name, size = pathways_size)
    res_stats <- sapply(pathways, function(x) {
    # Number of test/background genes
    a_n <- length(test)
    b_n <- length(background)

    # Number of test/background genes in pathway
    a <- intersect(test, x)
    b <- intersect(background, x)
    a_len <- length(a)
    b_len <- length(b)

    # Fraction of test/background genes in pathway
    a_x <- (a_len / a_n) * 100
    b_x <- (b_len / b_n) * 100

    # Pathway enriched or depleted for test genes
    e <- ifelse(a_x > b_x, "Enriched", "Depleted")

    # Set up contingency table for p value calculation
    m <- matrix(ncol = 2, nrow = 2)
    m[1,] <- c(a_len, b_len)
    m[2,] <- c(a_n, b_n)

    # Calculate p value using fishers exact test
    p <- fisher.test(m, alternative = "two.sided")$p.value

    # Return
    return(list(real_num = a_len,
                real_frac = a_x,
                expected_num = b_len,
                expected_frac = b_x,
                enrichment = e,
                pvalue = p,
                real_pathway_gene = a))
    })

    # Populate res table
    res[, "real_total" := length(test)]
    res[, "real_pathway" := unlist(res_stats["real_num",])]
    res[, "real_pathway_frac" := unlist(res_stats["real_frac",])]
    res[, "expected_total" := length(background)]
    res[, "expected_pathway" := unlist(res_stats["expected_num",])]
    res[, "expected_pathway_frac" := unlist(res_stats["expected_frac",])]
    res[, "enrichment" := unlist(res_stats["enrichment",])]
    res[, "real_pathway_gene" := mapply("[", res_stats["real_pathway_gene",])]
    res[, "pvalue" := unlist(res_stats["pvalue",])]

    # BH-correct p values
    res <- res[order(res$pvalue),]
    res$qvalue <- p.adjust(res$pvalue, method = "BH")
    message("You did it! FEDUP ran successfully, feeling pretty good huh?")

    # Return table
    return(res)
}
