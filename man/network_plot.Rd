\name{network_plot}
\alias{network_plot}
\encoding{UTF-8}

\title{
  Plot an interactive causal association gene network with the GraphViz package
}
\description{
  A function that plots an interactive causal association gene network using the GraphViz package.
}
\usage{
  network_plot(x)
}
\arguments{
  \item{x}{
    A \code{data.table} object that contains the results of the \code{meta_pcor} function.
  }
}
\details{
  This function generates an interactive network plot of the causal association gene network using the GraphViz package.
}
\value{
  An interactive network plot of the causal association gene network.
}
\references{
  % Add references if any.
}
\author{
  Ioanna V. Sasilioglou, Georgios A. Manios, Polymnia V. Gkoumplia, Panagiota I. Kontou, Pantelis G. Bagos* (\emph{Correspondence to:} \email{pbagos@compgen.org})
}
\note{
  % Add notes if any.
}
\seealso{
  % Add related functions or packages here.
}
\examples{
 \dontrun{  # Network plot example
  network_plot(meta_pcor_results)}

}
\keyword{network-analysis}
\keyword{gene-network}
