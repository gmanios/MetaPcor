\name{pcor_neighborhood}
\alias{pcor_neighborhood}
\encoding{UTF-8}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  A function to calculate sparse partial correlation estimation with
neighborhood selection approach method for each study
}
\description{
 A function to calculate sparse partial correlation estimation with
neighborhood selection approach method for each study
}
\usage{
pcor_neighborhood(list_of_files, l1, l2=0, significant, pvalue_thres=NULL,
fdr_thres= NULL, coef_thres = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list_of_files}{
A list of files that contain the transposed gene expression studies as
data.table objects

}
  \item{l1}{
     l1 Norm penalty. l1 = 0.8 is suggested (float)

}

  \item{l2}{
     l2 value for sparse partial correaltion calculation. l2 = 0.0 is suggested
     (float)

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
   \link{https://github.com/cran/space}}

\examples{
\dontrun{

# Calculate sparse partial correlation estimation with
# neighborhood selection approach method for each study

sparse_partial_corr <- pcor_neighborhood(
  list_of_files = list_of_studies,
  l1 = 0.8,
  l2 = 0.0,
  significant = TRUE,
  pvalue_thres = 0.05
)
}

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
