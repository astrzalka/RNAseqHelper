% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_helpers.R
\name{make_comparison}
\alias{make_comparison}
\title{Compare expression in two strains specified by contrats using edgeR package glmTreat}
\usage{
make_comparison(
  fit,
  contrast_comp,
  design,
  genes,
  plot_md = TRUE,
  plot_volcano = TRUE,
  complete_list = FALSE,
  fold_change = 1.5,
  toptags_print = 5
)
}
\arguments{
\item{fit}{fit object made using edgeR package}

\item{contrast_comp}{strains to compare, have to present in the design}

\item{design}{design specifying which libraries belong to which strain}

\item{genes}{data_frame containing genes data}

\item{plot_md}{should md plot be shown? defaut TRUE}

\item{plot_volcano}{should volcano plot be shown? default TRUE}

\item{complete_list}{should it return complete list of genes or only those differentially expressed? default FALSE}

\item{fold_change}{fold change used in glmTreat function}

\item{toptags_print}{how many genes should be printed?}
}
\value{
result of topTags function for all or differentially expressed genes
}
\description{
Compare expression in two strains specified by contrats using edgeR package glmTreat
}
