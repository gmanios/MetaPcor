\name{meta_pcor}

\alias{meta_pcor}

\title{
Partial Correlation Meta-Analysis Function
}
\description{
This function performs meta-analysis of partial correlation coefficients using neighborhood selection approach or shrinkage methods.
The function accepts data from multiple studies as input and outputs a meta-analysis result that represents a summary of the partial correlation coefficients obtained from the individual studies.
}
\usage{
meta_pcor(file_names, option, method, meta_method= "random",
pvalue_thres = NULL, fdr_thres = NULL, 
coef_thres = NULL,l1 = NULL ,l2 = NULL, norm_data = NULL, norm_method = NULL)
}
\arguments{
\item{file_names}{a character vector containing the file names of the data to be analyzed.}

  \item{option}{
     The meta-analysis option to run (integer).
      1: Pearson correlation and partial corrleation meta-analysis
      2: Partial correlation meta-analysis with thresholds
      3: Partial correlation meta-analysis without thresholds
      4: Differential expression and partial correlation meta-analysis

}
\item{method}{a character string indicating the method to be used to estimate partial correlation coefficients. Available options are "sparse" and "shrinkage".}

\item{meta_method}{a character string indicating the method to be used for the meta-analysis of the partial correlation coefficients. Available options are "random" and "fixed". Default is "random".}

\item{pvalue_thres}{a numeric value specifying the threshold for p-values.}

\item{fdr_thres}{a numeric value specifying the threshold for false discovery rate (FDR).}

\item{coef_thres}{a numeric value specifying the threshold for correlation coefficient magnitude.}

\item{l1}{a numeric value specifying the parameter for neighborhood selection approach.}

\item{l2}{a numeric value specifying the parameter for neighborhood selection approach.}

\item{norm_data}{a character string indicating whether data normalization should be performed. Available options are "YES" and "NO". Default is "NO".}

\item{norm_method}{a character string indicating the method to be used for normalization. Available options are "RPKM" and "Z_SCOR_STAND","LOG2" and "QUANT_NORM".}
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


metapcor_result <- meta_pcor(folder_path = "/path/to/data/folder",
option = 1, method = "sparse",meta_method = "random",l1 = 0.8, l2 = 0)

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
