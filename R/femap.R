#' Writes an enrichment dataset file for use in Cytoscape EnrichmentMap
#' @param df (data.frame) table with FEDUP enrichment results
#'  (see ?runFedup() for column descriptions)
#' @param results_file (char) name of output results file
#' @return table of gene enrichment and depletion results formatted as a
#' 'Generic results file'. Rows represent tested pathways. Columns represent:
#' \itemize{
#'    \item geneset -- gene-set ID (must match the gene-set ID in the GMT file);
#'    \item description -- gene-set name or description;
#'    \item pvalue -- enrichment pvalue;
#'    \item qvalue -- FDR correction value;
#'    \item phenotype -- +1 or -1, to identify enrichment in up- and down-regulation,
#'             or, more in general in either of the two phenotypes being compared
#'             in the two-class analysis (+1 maps to red, -1 maps to blue)
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

	# Format FEDUP results to generate generic EnrichmentMap in Cytoscape
  df_em <- df %>%
		select("pathway", "pvalue", "fdr", "enrichment") %>%
		mutate("description" = gsub("\\%.*", "", df$pathway)) %>%
    mutate("enrichment" = ifelse(df$enrichment == "enriched", "1", "-1")) %>%
		select("pathway", "description", "pvalue", "fdr", "enrichment") %>%
		`colnames<-` (c("geneset", "description", "pvalue", "qvalue", "phenotype"))

  # Write out file
	fwrite(df_em, results_file, sep = "\t", col.names = TRUE, quote = FALSE)
}

#' Plots an EnrichmentMap (EM) in Cytoscape to visualize enriched and depleted gene sets.
#' @param gmt_file (char) path to GMT file (must be absolute path)
#' @param results_file (char) path to file with enrichment results (must be absolute path)
#' @param pvalue (numeric) pvalue cutoff (value between 0 and 1). Gene set
#'  nodes with a p-value lower than the given value will not be included
#' 	in the network. (default = 1)
#' @param qvalue (numeric) FDR qvalue cutoff (value between 0 and 1). Gene set
#' 	nodes with a p-value lower than the given value will not be included
#' 	in the network. (default = 1)
#' @param net_name (char) name for network in Cytoscape
#' @return filename of image to which EM is exported. Also side effect of
#'  plotting the EM in an open session of Cytoscape.
#' @examples
#' \dontrun{
#' # not run because this function requires Cytoscape to be installed and open
#' data(testGene)
#' data(backgroundGene)
#' data(pathwaysGMT)
#' gmt_file <- tempfile("pathwaysGMT", fileext = ".gmt")
#' writePathways(pathwaysGMT, gmt_file)
#' fedup_res <- runFedup(testGene, backgroundGene, pathwaysGMT)
#' results_file <- tempfile("fedup_res", fileext = ".txt")
#' writeFemap(fedup_res, results_file)
#' plotFemap(
#'		gmt_file = gmt_file,
#'		results_file = results_file,
#'		qvalue = 0.05,
#'		net_name = "FEDUP_EM"
#' )}
#' @import RCy3
#' @export
plotFemap <- function(gmt_file,
											results_file,
											pvalue = 1,
											qvalue = 1,
							  			net_name = "generic") {

	# Confirm that Cytoscape is installed and opened
	cytoscapePing()

  # Generate EM using given parameters
	if (net_name %in% getNetworkList()) {
		deleteNetwork(net_name)
	}

	cat("* Building the network\n")
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

  cat("* Setting chart data\n")
  # Node visualization (enriched = red nodes, depleted = blue nodes)
  ch_command <- paste('enrichmentmap chart data="NES_VALUE"',
      "colors=", "RD_BU_9")
  response <- commandsGET(ch_command)

 	# Annotate similar pathways using AutoAnnotate
	cat("* Annotating the network using AutoAnnotate\n")
  aa_command <- paste(
		"autoannotate annotate-clusterBoosted",
		"clusterAlgorithm=MCL",
		"maxWords=3",
		"network=", net_name
	)
	response <- commandsGET(aa_command)

  # Network layout
  cat("* Applying network layout\n")
  ln_command <- paste(
    "layout force-directed",
    "network=", net_name
  )
  response <- commandsGET(ln_command)
	fitContent()
}
