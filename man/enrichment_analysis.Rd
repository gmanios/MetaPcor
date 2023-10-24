\name{enrichment_analysis}
\alias{enrichment_analysis}
\encoding{UTF-8}
\title{
  A function to conduct enrichment analysis with the statistically
  significant genes that occurred from the meta-analysis
}
\description{
  A function to conduct enrichment analysis with the statistically
  significant genes that occurred from the meta-analysis.
}
\usage{
enrichment_analysis(x)
}
\arguments{
  \item{x}{
    A dataframe object containing the results of the meta_pcor function.
  }
}
\details{

}
\value{
  A list with the results of the enrichment analysis with the gProfiler tool:
  \enumerate{
    \item A dataframe with the results of the enrichment.
    \item An interactive Manhattan plot with the statistically significant terms from the sources.
  }
}
\author{
  Ioanna V. Sasilioglou, Georgios A. Manios, Polymnia V. Gkoumplia, Panagiota I. Kontou, Pantelis G. Bagos
  \email{pbagos@compgen.org}  (\emph{Correspondence})
}
\note{

}
\seealso{

}
\examples{

\dontrun{ea_results <- enrichment_analysis(meta_pcor_results)}


}
\references{
  Uku Raudvere, Liis Kolberg, Ivan Kuzmin, Tambet Arak, Priit Adler, Hedi Peterson, Jaak Vilo,
  \emph{g:Profiler: a web server for functional enrichment analysis and conversions of gene lists (2019 update)},
  Nucleic Acids Research, Volume 47, Issue W1, 02 July 2019, Pages W191–W198,
  \url{https://doi.org/10.1093/nar/gkz369}
}
\keyword{enrichment}
\keyword{meta-analysis}

