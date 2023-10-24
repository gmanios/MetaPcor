\name{volc_plot}
\alias{volc_plot}
\encoding{UTF-8}

\title{
  Produce an interactive volcano plot with the Plotly package
}
\description{
  A function that produces an interactive volcano plot using the Plotly package.
}
\usage{
  volc_plot_plotly(x, pval_thres, coeff_thres)
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
  % Add references if you have any.
}
\author{
  Ioanna V. Sasilioglou, Georgios A. Manios, Polymnia V. Gkoumplia, Panagiota I. Kontou, Pantelis G. Bagos* (*correspondence to: \email{pbagos@compgen.org})
}
\note{
  % Add notes if you have any.
}
\seealso{
  % Add related functions or packages.
}
\examples{
\dontrun{  # Volcano Plot Example
  volc_plot <- volc_plot_plotly(meta_pcor_results, pval_thres = 0.05, coeff_thres = 0.2)}

}
\keyword{plotting}
\keyword{volcano-plot}
