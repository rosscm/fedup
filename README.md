# FEDUP

FEDUP is an R package for enrichment and depletion analysis on
user-defined pathways using a Fisher’s exact test. This package also
gives the option to draw a network representation of pathway overlaps
using EnrichmentMap.

## Getting started

### System prerequisites

R version &gt;= 4.0  
R packages:

-   **CRAN**: openxlsx, tibble, dplyr, data.table, ggplot2, ggthemes,
    forcats, RColorBrewer  
-   **Bioconductor**: RCy3

### Installation

Install FEDUP via devtools:

    #devtools::install_github("rosscm/FEDUP")
    devtools::load_all()

## Quick run

### Data input

Load example test genes, background genes, and pathways:

To note, the test genes comprise solely of a **muscle contraction**
pathway (Reactome ID 397014). So we would expect to see strong
*enrichment* for pathways related to muscle contraction and *depletion*
for pathways not associated with muscle contraction. Let’s see!

    data(testGene)
    data(backgroundGene)
    data(pathwaysGMT)

Take a look at the data structure:

    str(testGene)
    #>  chr [1:190] "NKX2-5" "SCN4A" "ITGB5" "SCN4B" "PAK2" "GATA4" "AKAP9" ...
    str(backgroundGene)
    #>  chr [1:10208] "PCYT1B" "PCYT1A" "PLA2G4D" "PLA2G4B" "PLA2G4C" "PLA2G4A" ...
    str(head(pathwaysGMT))
    #> List of 6
    #>  $ REGULATION OF PLK1 ACTIVITY AT G2 M TRANSITION%REACTOME%R-HSA-2565942.1          : chr [1:84] "CSNK1E" "DYNLL1" "TUBG1" "CKAP5" ...
    #>  $ GLYCEROPHOSPHOLIPID BIOSYNTHESIS%REACTOME%R-HSA-1483206.4                        : chr [1:126] "PCYT1B" "PCYT1A" "PLA2G4D" "PLA2G4B" ...
    #>  $ MITOTIC PROPHASE%REACTOME DATABASE ID RELEASE 74%68875                           : chr [1:134] "SETD8" "NUMA1" "NCAPG2" "LMNB1" ...
    #>  $ ACTIVATION OF NF-KAPPAB IN B CELLS%REACTOME%R-HSA-1169091.1                      : chr [1:67] "PSMA6" "PSMA3" "PSMA4" "PSMA1" ...
    #>  $ CD28 DEPENDENT PI3K AKT SIGNALING%REACTOME DATABASE ID RELEASE 74%389357         : chr [1:22] "CD28" "THEM4" "AKT1" "TRIB3" ...
    #>  $ UBIQUITIN-DEPENDENT DEGRADATION OF CYCLIN D%REACTOME DATABASE ID RELEASE 74%75815: chr [1:52] "PSMA6" "PSMA3" "PSMA4" "PSMA1" ...

### Pathway analysis

Now run FEDUP on sample data:

    fedup_res <- runFedup(testGene, backgroundGene, pathwaysGMT)
    #> Data input:
    #>  => 190 test genes
    #>  => 10208 background genes
    #>  => 1437 pathawys
    #> You did it! FEDUP ran successfully, feeling pretty good huh?

View output results table sorted by pvalue:

    print(fedup_res)
    #>                                                                                                                      pathway
    #>    1:                                                              MUSCLE CONTRACTION%REACTOME DATABASE ID RELEASE 74%397014
    #>    2:                                                             CARDIAC CONDUCTION%REACTOME DATABASE ID RELEASE 74%5576891
    #>    3:                                                                               ION HOMEOSTASIS%REACTOME%R-HSA-5578775.2
    #>    4:                                                       SMOOTH MUSCLE CONTRACTION%REACTOME DATABASE ID RELEASE 74%445355
    #>    5:                                                                    STRIATED MUSCLE CONTRACTION%REACTOME%R-HSA-390522.1
    #>   ---                                                                                                                       
    #> 1433:                                                                     IRAK4 DEFICIENCY (TLR2 4)%REACTOME%R-HSA-5603041.1
    #> 1434:                          ACTIVATION OF KAINATE RECEPTORS UPON GLUTAMATE BINDING%REACTOME DATABASE ID RELEASE 74%451326
    #> 1435: TNF RECEPTOR SUPERFAMILY (TNFSF) MEMBERS MEDIATING NON-CANONICAL NF-KB PATHWAY%REACTOME DATABASE ID RELEASE 74%5676594
    #> 1436:                                                                     RHO GTPASES ACTIVATE KTN1%REACTOME%R-HSA-5625970.1
    #> 1437:                                 TRANSPORT OF MATURE MRNA DERIVED FROM AN INTRONLESS TRANSCRIPT%REACTOME%R-HSA-159231.2
    #>       size real_frac expected_frac fold_enrichment   status
    #>    1:  190 100.00000     1.8612853        53.72632 Enriched
    #>    2:  124  65.26316     1.2147335        53.72632 Enriched
    #>    3:   51  26.84211     0.4996082        53.72632 Enriched
    #>    4:   37  19.47368     0.3624608        53.72632 Enriched
    #>    5:   34  17.89474     0.3330721        53.72632 Enriched
    #>   ---                                                      
    #> 1433:   11   0.00000     0.1077586         0.00000 Depleted
    #> 1434:   29   0.00000     0.2840909         0.00000 Depleted
    #> 1435:   16   0.00000     0.1567398         0.00000 Depleted
    #> 1436:   11   0.00000     0.1077586         0.00000 Depleted
    #> 1437:   37   0.00000     0.3624608         0.00000 Depleted
    #>                                       real_gene        pvalue        qvalue
    #>    1:   NKX2-5,SCN4A,ITGB5,SCN4B,PAK2,GATA4,... 1.091522e-189 1.568518e-186
    #>    2: NKX2-5,SCN4A,SCN4B,GATA4,AKAP9,KCNJ14,... 4.477692e-130 3.217222e-127
    #>    3:    SLN,STIM1,ORAI2,ORAI1,ABCC9,KCNJ11,...  1.513045e-57  7.247487e-55
    #>    4:      ITGB5,PAK2,ACTA2,VCL,MYL12B,MYL6,...  1.161897e-42  4.174116e-40
    #>    5:          VIM,TNNI3,DMD,TPM4,TPM3,TPM2,...  2.009234e-39  5.774540e-37
    #>   ---                                                                      
    #> 1433:                                            1.000000e+00  1.000000e+00
    #> 1434:                                            1.000000e+00  1.000000e+00
    #> 1435:                                            1.000000e+00  1.000000e+00
    #> 1436:                                            1.000000e+00  1.000000e+00
    #> 1437:                                            1.000000e+00  1.000000e+00

### Visualization

Plot enriched and depleted pathways (qvalue &lt; 5%) in the form of a
dot plot:

    fedup_plot <- fedup_res[which(fedup_res$qvalue < 0.05),]
    fedup_plot$log10qvalue <- -log10(fedup_plot$qvalue + 1e-10) # log10-transform qvalue for plotting
    fedup_plot$pathway <- gsub("\\%.*", "", fedup_plot$pathway) # clean pathway names
    p <- plotDotPlot(
      df = fedup_plot,
      x_var = "log10qvalue",
      y_var = "pathway",
      x_lab = "-log10(Qvalue)",
      fill_var = "status",
      fill_lab = "Enrichment\nstatus",
      size_var = "fold_enrichment",
      size_lab = "Fold enrichment"
    )
    print(p)

![](man/figures/README_FEDUP_dotplot-1.png)

As expected, we see strong enrichment for muscle-related pathways at the
top of the plot, and depletion for olfactory and amino acid metabolism
pathways at the bottom of the plot. Nice!

We can also facet the plot by enrichment status to clearly separate the
enriched and depleted pathways:

    p <- p +
      facet_grid("status", scales = "free", space = "free") +
      theme(strip.text.y = element_blank())
    print(p)

![](man/figures/README_FEDUP_dotplot_facet-1.png)

Look at all those chick… enrichments! This is a bit overwhelming, no?
What if we could see all these pathways in a summarised way that doesn’t
hurt our tired brains even more? Oh I know… let’s use an Enrichment Map!

First, make sure to have
[Cytoscape](https://cytoscape.org/download.html) downloaded and and open
on your computer. You’ll also need to install the
[EnrichmentMap](http://apps.cytoscape.org/apps/enrichmentmap) and
[AutoAnnotate](http://apps.cytoscape.org/apps/autoannotate) apps.

Then format FEDUP results for compatibility with EnrichmentMap:

    results_file <- tempfile("fedup_res", fileext = ".txt")
    writeFemap(fedup_res, results_file)
    #> Wrote Cytoscape-formatted FEDUP results file to /var/folders/mh/_0z2r5zj3k75yhtgm6l7xy3m0000gn/T//RtmpSn19kM/fedup_rese720218febc2.txt

Prepare a pathway annotation file (GMT format) from the pathway list you
passed to FEDUP (you don’t need to run this function if your pathway
annotations are already in GMT format, but it doesn’t hurt to make
sure):

    gmt_file <- tempfile("pathwaysGMT", fileext = ".gmt")
    writePathways(pathwaysGMT, gmt_file)
    #> Wrote out GMT file with to /var/folders/mh/_0z2r5zj3k75yhtgm6l7xy3m0000gn/T//RtmpSn19kM/pathwaysGMTe7206c7d0f2c.gmt

Cytoscape is open right? If so, uncomment these lines and let the magic
happen:

    #net_file <- tempfile("FEDUP_EM", fileext = ".png")
    #plotFemap(
    #  gmt_file = gmt_file,
    #  results_file = results_file,
    #  qvalue = 0.05,
    #  net_name = "FEDUP_EM",
    #  net_file = net_file
    #)

## Versioning

For the versions available, see the [tags on this
repo](https://github.com/rosscm/FEDUP/tags).

## Shoutouts

:sparkles:[**2020**](https://media.giphy.com/media/z9AUvhAEiXOqA/giphy.gif):sparkles:
