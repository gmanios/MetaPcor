\name{my_meta}
\alias{my_meta}
\title{
  \code{my_meta}: Perform meta-analysis of partial correlations
}
\description{
  A function that performs meta-analysis of partial correlations.
}
\usage{
  \code{my_meta(correlations, method)}
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
  % Add references to the literature or website here.
}
\author{
  Ioanna V. Sasilioglou, Georgios A. Manios, Polymnia V. Gkoumplia, Panagiota I. Kontou, Pantelis G. Bagos* (*correspondence to: \email{pbagos@compgen.org})
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
 my_meta_results <- my_meta(correlations = sparse_partial_corr, method = 'random')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
