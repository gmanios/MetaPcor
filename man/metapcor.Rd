\name{meta_pcor}
\alias{meta_pcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}


\description{
The main function of the package. All options can be run with this funciton.
}
\usage{
meta_pcor(folder_path = NULL, GEO_names = NULL, target_namespace = NULL,
option, method, meta_method = "random",
pvalue_thres = NULL, fdr_thres = NULL,
coef_thres = NULL,l1  = NULL,l2 = NULL)}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{folder_path}{
     The folder path where local gene expression studies are stored (character)
}

  \item{GEO_names}{
     Enter the GSE ID of the study (character)
}
  \item{target_namespace}{
     Enter the target namespace of gConvert for probe to gene name conversion (character).

}

  \item{option}{
     The meta-analysis option to run (integer).
      1: Pearson correlation and partial corrleation meta-analysis
      2: Partial correlation meta-analysis with thresholds
      3: Partial correlation meta-analysis without thresholds
      4: Differential expression and partial correlation meta-analysis

}

  \item{method}{
     The partial correlation option to run (character). Available options : "sparse" or "shrinkage".

}

  \item{meta_method}{
     Choose the model of the meta-analysis (character). Available options : "random" or "fixed".

}
  \item{pvalue_thres}{
     p-value threshold for the meta-analysis to use when option = 2  (float).

}

  \item{fdr_thres}{
     FDR value threshold for the meta-analysis to use when option = 2  (float).

}

  \item{coef_thres}{
     Coefficient threshold for the meta-analysis to use when option = 2 (float).

}

  \item{l1}{
     l1 value for sparse partial correaltion calculation. l1 = 0.8 is suggested. (float).

}

  \item{l2}{
     l2 value for sparse partial correaltion calculation. l2 = 0 is suggested. (float).

}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 A data.table object is returned, containing the results of the meta-analysis.
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
metapcor_result <- meta_pcor(GEO_names=c("GSE76427"),
target_namespace = c('ILLUMINA_HUMANHT_12_V4'), option = 2,
method = "sparse", meta_method = "random",
pvalue_thres = 0.01, l1 = 0.8, l2 = 0)

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
