\name{my_meta}
\alias{my_meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to  perform meta-analysis of partial correlations ~~
}
\description{
A function that performs meta-analysis of partial correlations
}
\usage{
  my_meta(correlations, method)
  }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{correlations}{
A partial correlation data.table object, which occurs from the results of the 
partial_neighborhood or the partial_shrinkage functions 
}

  \item{method}{
 Meta-analysis method to calculate the effect sizes ('random' or 'fixed')
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A data.table containing the results of the meta-analysis
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ioanna V. Sasilioglou, Georgios A. Manios, Polymnia V. Gkoumplia, Panagiota I. Kontou, Pantelis G. Bagos* (*correspondence to : pbagos@compgen.org)}
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
