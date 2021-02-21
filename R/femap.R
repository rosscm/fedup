#' Writes an enrichment dataset file for use in Cytoscape EnrichmentMap.
#'
#' @param df (data.frame) table with fedup enrichment results.
#'  (see runFedup() for column descriptions)
#' @param resultsFile (char) name of output results file.
#' @return table of gene enrichment and depletion results formatted as a
#' 'Generic results file'. Rows represent tested pathways. Columns represent:
#' \itemize{
#'     \item pathway -- pathway ID (must match pathway IDs in the GMT file
#'         provided to plotFemap());
#'     \item description -- pathway name or description;
#'     \item pvalue -- enrichment pvalue;
#'     \item qvalue -- BH-corrected pvalue;
#'     \item status -- +1 or -1, to identify enriched or depleted pathways
#'         (+1 maps to red, -1 maps to blue)
#' }
#' @examples
#' data(testGene)
#' data(backgroundGene)
#' data(pathwaysGMT)
#' fedupRes <- runFedup(testGene, backgroundGene, pathwaysGMT)
#' resultsFile <- tempfile("fedupRes", fileext = ".txt")
#' writeFemap(fedupRes, resultsFile)
#' @importFrom data.table fwrite
#' @importFrom dplyr select mutate %>%
#' @export
writeFemap <- function(df, resultsFile) {
    df_em <- df %>%
        select("pathway", "pvalue", "qvalue", "status") %>%
        mutate("description" = gsub("\\%.*", "", df$pathway)) %>%
        mutate("status" = ifelse(df$status == "Enriched", "1", "-1")) %>%
        select("pathway", "description", "pvalue", "qvalue", "status")

    fwrite(df_em, resultsFile, sep = "\t", col.names = TRUE, quote = FALSE)
    message("Wrote out Cytoscape-formatted fedup results file to ", resultsFile)
}

#' Draws a network representation of overlaps among enriched and depleted
#' pathways using EnrichmentMap (EM) in Cytoscape.
#'
#' @param gmtFile (char) path to GMT file (must be an absolute path).
#' @param resultsFile (char) path to file with fedup results
#'  (must be an absolute path).
#' @param pvalue (numeric) pvalue cutoff. Pathways with a higher pvalue
#'  will not be included in the EM (value between 0 and 1; default 1).
#' @param qvalue (numeric) qvalue cutoff. Pathways with a higher qvalue
#'  will not be included in the EM (value between 0 and 1; default 1).
#' @param netName (char) name for EM in Cytoscape (default generic).
#' @param netFile (char) name of output image. Supports png, pdf, svg,
#'  jpeg image formats.
#' @return file name of image to which the network is exported. Also side
#'  effect of plotting the EM in an open session of Cytoscape.
#' @examples
#' \dontrun{
#' # not run because Cytoscape needs to be installed and open
#' data(testGene)
#' data(backgroundGene)
#' data(pathwaysGMT)
#' gmtFile <- tempfile("pathwaysGMT", fileext = ".gmt")
#' fedupRes <- runFedup(testGene, backgroundGene, pathwaysGMT)
#' resultsFile <- tempfile("fedupRes", fileext = ".txt")
#' netFile <- tempfile("fedup_EM", fileext = ".png")
#' writePathways(pathwaysGMT, gmtFile)
#' writeFemap(fedupRes, resultsFile)
#' plotFemap(
#'     gmtFile = gmtFile,
#'     resultsFile = resultsFile,
#'     qvalue = 0.05,
#'     netName = "fedup_EM",
#'     netFile = netFile
#' )
#' }
#' @import RCy3
#' @export
plotFemap <- function(gmtFile, resultsFile, pvalue = 1, qvalue = 1,
    netName = "generic", netFile = "png") {

    # Confirm that Cytoscape is installed and opened
    cytoscapePing()
    if (netName %in% getNetworkList()) {
        deleteNetwork(netName)
    }

    message("Building the network")
    em_command <- paste(
        'enrichmentmap build analysisType="generic"',
        "gmtFile=", gmtFile,
        "enrichmentsDataset1=", resultsFile,
        "pvalue=", pvalue,
        "qvalue=", qvalue,
        "similaritycutoff=", 0.375,
        "coefficients=", "COMBINED",
        "combinedConstant=", 0.5
    )
    response <- commandsGET(em_command)
    renameNetwork(netName, getNetworkSuid())

    message("Setting network chart data")
    ch_command <- paste(
        'enrichmentmap chart data="NES_VALUE"',
        "colors=", "RD_BU_9"
    )
    response <- commandsGET(ch_command)

    message("Annotating the network using AutoAnnotate")
    aa_command <- paste(
        "autoannotate annotate-clusterBoosted",
        "clusterAlgorithm=MCL",
        "maxWords=3",
        "network=", netName
    )
    response <- commandsGET(aa_command)

    message("Applying a force-directed network layout")
    ln_command <- paste(
        "layout force-directed",
        "network=", netName
    )
    response <- commandsGET(ln_command)
    fitContent()

    message("Drawing out network to ", netFile)
    exportImage(netFile)
}
