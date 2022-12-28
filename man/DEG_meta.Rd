\name{DEG_meta}
\alias{DEG_meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
Differential expression analysis, which keeps only the statistically 
significant genes from a given set of gene expression studies 
}
\usage{
DEG_meta(folder_path, cases, controls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{folder_path}{
A folder path where the gene expression studies are stored
}
  \item{cases}{
How the cases are named in the gene expression studies
}


  \item{controls}{
How the controls are named in the gene expression studies
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

A list of data.tables containing gene expression columns of the statistically significant genes that occured from the differential epxression analysis

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
 DEG_ES <- DEG_meta(folder_path = 'studies/', cases = 'CASE', 
 controls = 'CONTROL')


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
