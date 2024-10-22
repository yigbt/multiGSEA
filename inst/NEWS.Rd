\name{NEWS}
\title{News for Package \pkg{multiGSEA}}

\section{ Changes in version 0.99.0}{
\itemize{

  \item Package added to Bioconductor

}}


\section{ Changes in version 1.1.1 (2020-11-05)}{

\itemize{

  \item Bug fix. Prevent duplicated pathway title to cause an error.
  \item Include new function to enumerate duplicted pathway titles with trailing numbers.

}}


\section{ Changes in version 1.1.2 (2020-11-07)}{

\itemize{

  \item Speed-up the metabolite ID mapping between different ID formats

}}


\section{ Changes in version 1.1.3 (2020-11-19)}{

\itemize{

  \item Prevent warning due to import of two different select() functions.

}}


\section{ Changes in version 1.1.4 (2020-12-09)}{

\itemize{

  \item Add the correct citation of the BMC Bioinformatics article.

}}

\section{ Changes in version 1.1.5 (2021-04-19)}{

\itemize{

  \item Typos fixed.
  \item Add information and reference about metrics to calculate the feature ranks.

}}

\section{ Changes in version 1.1.7 (2021-05-10)}{

\itemize{

  \item Remove Biocarta database from vignette - it's no longer supported by graphite R package.
  \item List all available databases in error message.

}}

\section{ Changes in version 1.13.1 (2024-04-26)}{

\itemize{

  \item Utilize by default non-adjusted single-omics p-values for aggregation into a composite p-value. Adjustment for multiple testing will be done on the combined multi-omics p-values.
  \item Inlcude legacy option for p-value combination to ensure reproducibility.

}}

\section{ Changes in version 1.15.1 (2024-10-22)}{

\itemize{

  \item Remove proteomics feature with empty identifier to prevent errors in vignette building.

}}
