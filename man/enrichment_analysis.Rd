\name{enirchment_analysis}
\alias{enirchment_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
}
\description{
A function to conduct enrichment analysis with the statistically
significant genes that occured from the meta-analysis
}
\usage{
enirchment_analysis(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A dataframe object containing the results of the meta_pcor function

}
}
\details{

}
\value{
 A list with the results of the enrichment analysis with the gProfiler tool
 ([1]: A dataframe with the results of the enricment and [2] an interactive
 Manhattan plot with the statistically significant terms from the sources)
}

\author{
Ioanna V. Sasilioglou, Georgios A. Manios, Polymnia V. Gkoumplia, Panagiota I. Kontou, Pantelis G. Bagos* (*correspondence to : pbagos@compgen.org)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \link{https://biit.cs.ut.ee/gprofiler/gost}}

\examples{

ea_results <- enrichment_analysis(meta_pcor_results)

}

\references{

Uku Raudvere, Liis Kolberg,
Ivan Kuzmin, Tambet Arak,
Priit Adler,Hedi Peterson, Jaak Vilo,
g:Profiler: a web server for functional
enrichment analysis and conversions of gene lists (2019 update),
Nucleic Acids Research, Volume 47, Issue W1, 02 July 2019,
Pages W191–W198,
https://doi.org/10.1093/nar/gkz369

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
