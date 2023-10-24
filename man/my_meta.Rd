\name{my_meta}
\alias{my_meta}
\encoding{UTF-8}

\title{
  Perform meta-analysis of partial correlations
}
\description{
  A function that performs meta-analysis of partial correlations.
}
\usage{
  my_meta(correlations, method)
}
\arguments{
  \item{correlations}{
    A \code{partial correlation} \code{data.table} object, obtained from the results of the \code{partial_neighborhood} or \code{partial_shrinkage} functions.
  }

  \item{method}{
    Meta-analysis method to calculate the effect sizes. Should be either \code{"random"} or \code{"fixed"}.
  }
}
\details{
  This function performs a meta-analysis of partial correlations using the specified method.
}
\value{
  A \code{data.table} object containing the results of the meta-analysis.
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
 \dontrun{
 my_meta_results <- my_meta(correlations = sparse_partial_corr, method = 'random')
}
}
\keyword{meta-analysis}
\keyword{correlation}
