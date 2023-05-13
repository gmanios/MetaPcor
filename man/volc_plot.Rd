\name{volc_plot}
\alias{volc_plot}
\title{
  \code{volc_plot}: Produce an interactive volcano plot with the Plotly package
}
\description{
  A function that produces an interactive volcano plot using the Plotly package.
}
\usage{
  \code{volc_plot_plotly(x, pval_thres, coeff_thres)}
}
\arguments{
  \item{x}{
    A \code{data.table} object that contains the results of the \code{meta_pcor} function.
  }

  \item{pval_thres}{
    A p-value threshold for the volcano plot (float).
  }

  \item{coeff_thres}{
    A coefficient threshold for the volcano plot (float).
  }
}
\details{
  This function generates an interactive volcano plot using the Plotly package. The volcano plot displays the relationship between statistical significance (p-values) and effect size (coefficients).
}
\value{
  An interactive volcano plot.
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
  # Volcano Plot
  volc_plot <- volc_plot_plotly(meta_pcor_results, pval_thres = 0.05, coeff_thres = 0.2)
}
\keyword{functions}
