\name{network_plot}
\alias{network_plot}
\title{
  \code{network_plot}: Plot an interactive causal association gene network with the GraphViz package
}
\description{
  A function that plots an interactive causal association gene network using the GraphViz package.
}
\usage{
  \code{network_plot(x)}
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
  % Add references to the literature or website here.
}
\author{
  Ioanna V. Sasilioglou, Georgios A. Manios, Polymnia V. Gkoumplia, Panagiota I. Kontou, Pantelis G. Bagos* (*correspondence to: \email{pbagos@compgen.org})
}
\note{
  % Add any additional notes here.
}
\seealso{
  % Add objects to see also as \code{\link{help}}.
}
\examples{
  # Network plot
  network_plot(meta_pcor_results)
}
\keyword{functions}
