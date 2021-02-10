#' Writes an enrichment dataset file for use in Cytoscape EnrichmentMap.
#'
#' @param df (data.frame) table with FEDUP enrichment results.
#'  (see runFedup() for column descriptions)
#' @param results_file (char) name of output results file.
#' @return table of gene enrichment and depletion results formatted as a
#' 'Generic results file'. Rows represent tested pathways. Columns represent:
#' \itemize{
#'     \item pathway -- pathway ID (must match pathway IDs in the GMT file
#'         provided to plotFemap());
#'     \item description -- pathway name or description;
#'     \item pvalue -- enrichment pvalue;
#'     \item qvalue -- BH-corrected pvalue;
#'     \item status -- +1 or -1, to identify enrichment in either of the two
#'         phenotypes being compared in the two-class analysis
#'         (+1 maps to red, -1 maps to blue)
#' }
#' @examples
#' data(testGene)
#' data(backgroundGene)
#' data(pathwaysGMT)
#' fedup_res <- runFedup(testGene, backgroundGene, pathwaysGMT)
#' results_file <- tempfile("fedup_res", fileext = ".txt")
#' writeFemap(fedup_res, results_file)
#' @importFrom data.table fwrite
#' @importFrom dplyr select mutate %>%
#' @export
writeFemap <- function(df, results_file) {
    df_em <- df %>%
        select("pathway", "pvalue", "qvalue", "status") %>%
        mutate("description" = gsub("\\%.*", "", df$pathway)) %>%
        mutate("status" = ifelse(df$status == "Enriched", "1", "-1")) %>%
        select("pathway", "description", "pvalue", "qvalue", "status")

    fwrite(df_em, results_file, sep = "\t", col.names = TRUE, quote = FALSE)
    message("Wrote Cytoscape-formatted FEDUP results file to ", results_file)
}

#' Draws a network representation of overlaps among enriched and depleted
#' pathways using EnrichmentMap (EM) in Cytoscape.
#'
#' @param gmt_file (char) path to GMT file (must be an absolute path).
#' @param results_file (char) path to file with FEDUP results
#'  (must be an absolute path).
#' @param pvalue (numeric) pvalue cutoff. Pathways with a higher pvalue
#'  will not be included in the EM (value between 0 and 1; default 1).
#' @param qvalue (numeric) qvalue cutoff. Pathways with a higher qvalue
#'  will not be included in the EM (value between 0 and 1; default 1).
#' @param net_name (char) name for EM in Cytoscape (default generic).
#' @param net_file (char) name of output image. Supports png, pdf, svg,
#'  jpeg image formats.
#' @return file name of image to which the network is exported. Also side
#'  effect of plotting the EM in an open session of Cytoscape.
#' @examples
#' \dontrun{
#'     # not run because Cytoscape needs to be installed and open
#'     data(testGene)
#'     data(backgroundGene)
#'     data(pathwaysGMT)
#'     gmt_file <- tempfile("pathwaysGMT", fileext = ".gmt")
#'     fedup_res <- runFedup(testGene, backgroundGene, pathwaysGMT)
#'     results_file <- tempfile("fedup_res", fileext = ".txt")
#'     net_file <- tempfile("FEDUP_EM", fileext = ".png")
#'     writePathways(pathwaysGMT, gmt_file)
#'     writeFemap(fedup_res, results_file)
#'     plotFemap(
#'         gmt_file = gmt_file,
#'         results_file = results_file,
#'         qvalue = 0.05,
#'         net_name = "FEDUP_EM",
#'         net_file = net_file)}
#' @import RCy3
#' @export
plotFemap <- function(gmt_file, results_file,
                        pvalue = 1, qvalue = 1,
                        net_name = "generic", net_file = "png") {
    # Confirm that Cytoscape is installed and opened
    cytoscapePing()
    if (net_name %in% getNetworkList()) {
        deleteNetwork(net_name)
    }

    message("Building the network")
    em_command <- paste(
        'enrichmentmap build analysisType="generic"',
        "gmtFile=", gmt_file,
        "enrichmentsDataset1=", results_file,
        "pvalue=", pvalue,
        "qvalue=", qvalue,
        "similaritycutoff=", 0.375,
        "coefficients=", "COMBINED",
        "combinedConstant=", 0.5)
    response <- commandsGET(em_command)
    renameNetwork(net_name, getNetworkSuid())

    # Node visualization (enriched = red nodes, depleted = blue nodes)
    message("Setting network chart data")
    ch_command <- paste(
        'enrichmentmap chart data="NES_VALUE"',
        "colors=", "RD_BU_9")
    response <- commandsGET(ch_command)

    # Annotate similar pathways using AutoAnnotate
    message("Annotating the network using AutoAnnotate")
    aa_command <- paste(
        "autoannotate annotate-clusterBoosted",
        "clusterAlgorithm=MCL",
        "maxWords=3",
        "network=", net_name)
    response <- commandsGET(aa_command)

    # Network layout
    message("Applying a force-directed network layout")
    ln_command <- paste(
        "layout force-directed",
        "network=", net_name)
    response <- commandsGET(ln_command)
    fitContent()

    # Draw out network to file
    message("Drawing out network to ", net_file)
    exportImage(net_file)
}
