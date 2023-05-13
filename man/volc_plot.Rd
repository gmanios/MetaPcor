\name{volc_plot}
\alias{volc_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~ function which produces an interactive volcano plot with the Plotly package 
} ~~
}
\description{
  A function which produces an interactive volcano plot with the Plotly package 
}
\usage{
volc_plot_plotly(x, pval_thres, coeff_thres)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A data.table object which contains the results of the meta_pcor function  
  }

  \item{pval_thres}{
   A p-value threshold for the volcano (float) 
}
  \item{coeff_thres}{
  A coefficient threshold for the volcano plot (float)  
}



}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An interactive volcano plot
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
# Volcano Plot
volc_plot<- volc_plot_plotly(meta_pcor_results, pval_thres = 0.05, 
coeff_thres = 0.2)

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
