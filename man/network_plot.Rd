\name{network_plot}
\alias{network_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to plot an interactive causal association gene netowrk with the GraphViz package... ~~
}
\description{
 A function which plots an interactive causal association gene netowrk with the GraphViz package
}
\usage{
network_plot(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A data.table object which contains the results of the meta_pcor function
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An interactive network plot of the causal association gene network
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ioanna V. Sasilioglou, Georgios A. Manios, Polymnia V. Gkoumplia, Panagiota I. Kontou, Pantelis G. Bagos* (*correspondence to : pbagos@compgen.org)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  # Network plot 
  network_plot(meta_pcor_results)
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
