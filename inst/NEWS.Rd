\name{NEWS}
\title{News for Package 'fedup'}
\encoding{UTF-8}
\section{Changes in fedup version 0.99.7 (2021-04-24)}{
  \itemize{
    \item added info to explain example data in README and vignettes
    \item made minor changes to exported function example comments
    \item \code{plotFemap}
      + parameterized several hard-coded variables
  }
}
\section{Changes in fedup version 0.99.6 (2021-02-24)}{
  \itemize{
    \item Updated package data script paths
  }
}
\section{Changes in fedup version 0.99.5 (2021-02-24)}{
  \itemize{
    \item Trimmed external data size
  }
}
\section{Changes in fedup version 0.99.4 (2021-03-24)}{
  \itemize{
    \item updated package datasets (geneSingle, geneDouble, geneMulti)
    \item created 3 vignettes to describe package implementation using each dataset
    \item \code{runFedup}
      + runs analysis on an input list rather than single test vector
      + fold enrichment calculation evaluates 0 for 0/0 instances
      + enriched pathways defined as fold enrichment ≥ 1 (instead of > 1)
    \item \code{writeFemap}
      + writes EM-formatted tables for list of fedup results
    \item \code{plotFemap}
      + implements \code{tryCatch()} to return NULL if Cytoscape is not running locally
    \item \code{prepInput}
      + new function to prepare input gene list
  }
}
\section{Changes in fedup version 0.99.3 (2021-02-24)}{
  \itemize{
    \item Updating R version dependency to 4.1
  }
}
\section{Changes in fedup version 0.99.2 (2021-02-24)}{
  \itemize{
    \item Untracking system files
  }
}
\section{Changes in fedup version 0.99.1 (2021-02-23)}{
  \itemize{
    \item Version bump
  }
}
\section{Changes in fedup version 0.99.0 (2021-02-17)}{
  \itemize{
    \item Submitted to Bioconductor
  }
}
