\name{DE_analysis}
\alias{DE_analysis}
\title{
  Differential Expression Analysis and Meta-analysis of Partial Correlation Coefficients
}
\description{
  This function performs differential expression analysis on gene expression data and performs a meta-analysis of correlation coefficients using sparse correlation and meta-correlation methods.
}
\usage{
  DE_analysis(folder_path, case, control, fold_threshold, p_value_threshold, l1)
}
\arguments{
  \item{folder_path}{
    Path to the folder containing gene expression data files.
  }
  \item{case}{
    String specifying the column name for the case samples in the gene expression data.
  }
  \item{control}{
    String specifying the column name for the control samples in the gene expression data.
  }
  \item{fold_threshold}{
    Numeric value specifying the fold-change threshold for considering genes as differentially expressed.
  }
  \item{p_value_threshold}{
    Numeric value specifying the p-value threshold for considering genes as statistically significant.
  }
  \item{l1}{
    Numeric value specifying the neighborhood size for sparse correlation calculation.
  }
}
\details{
  The function reads gene expression data files from the specified folder path and performs the following steps:
  \enumerate{
    \item Reads each data file.
    \item Splits the data into case and control samples.
    \item Applies log2 scale to the expression values.
    \item Calculates the mean expression of each group.
    \item Computes fold-change values between the case and control groups.
    \item Performs a t-test to obtain p-values for each gene.
    \item Applies fold-change and p-value thresholds to filter differentially expressed genes.
    \item Performs sparse correlation analysis using the filtered gene expression data.
    \item Performs meta-analysis of correlation coefficients using the sparse correlations.
  }
}
\value{
  Returns the meta-analysis results as a data.table.
}
\examples{
  # Example usage of DE_analysis function
  folder_path <- "/path/to/gene/expression/files"
  case <- "CASE"
  control <- "CONTROL"
  fold_threshold <- 0
  p_value_threshold <- 0.05
  l1 <- 0.8
  
  result <- DE_analysis(folder_path, case, control, fold_threshold, p_value_threshold, l1)
}
