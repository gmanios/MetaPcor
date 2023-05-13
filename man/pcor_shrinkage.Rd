\name{pcor_shrinkage}
\alias{pcor_shrinkage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%   A function to calculate sparse partial correlation estimation with the
shrinkage method for each study.
}
\description{
 A function to calculate sparse partial correlation estimation with the
shrinkage method for each study.
}
\usage{
pcor_shrinkage(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list_of_files}{
A list of files that contain the transposed gene expression studies as data.table objects

}
  \item{significant}{
     A boolean option ("TRUE" or "FALSE") to choose only the statistically significant pairs of partial correlations, either with a p-value or a FDR or a coefficient threshold
}

  \item{pvalue_thres}{
     p-value threshold for the meta-analysis (float)

}

  \item{fdr_thres}{
     FDR value threshold for the meta-analysis (float)

}

  \item{coef_thres}{
     Coefficient threshold for the meta-analysis (float)

}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data.table object is returned with the partial correlations pairs
}

\author{
Ioanna V. Sasilioglou, Georgios A. Manios, Polymnia V. Gkoumplia, Panagiota I. Kontou, Pantelis G. Bagos* (*correspondence to : pbagos@compgen.org)}
\note{
%%  ~~further notes~~
}


\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
shrinkage_partial_corr <- pcor_shrinkage(list_of_files, significant = TRUE, pvalue_thres = 0.05)
}
