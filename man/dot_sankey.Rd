% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dot_sankey.R
\name{dot_sankey}
\alias{dot_sankey}
\title{sankey and dot graphs showed the results of GO and KEGG analysis}
\usage{
dot_sankey(
  data,
  top = 8,
  sankey.text.size = 3,
  sankey.x.axis.text.position = 1,
  sankey.x.axis.text.size = 12,
  dot.color = "RdBu",
  dot.position = 8,
  dot.lable = "Description",
  dot.text.size = 12,
  dot.scale = 0.9,
  dot.x = 0.45,
  dot.y = 0.17,
  dot.width = 0.5,
  dot.height = 0.63,
  ...
)
}
\arguments{
\item{data}{the results after KEGG/GO analysis}

\item{top}{according to the order of p adjust value from small to large
the number of categories to show}

\item{sankey.text.size}{Sankey diagram text size}

\item{sankey.x.axis.text.position}{Location of x-axis labels in sankey plot}

\item{sankey.x.axis.text.size}{Sankey diagram x-axis text size}

\item{dot.color}{dotplot color: see "RColorBrewer::display.brewer.all()"}

\item{dot.position}{dotplot position}

\item{dot.lable}{"ID" or "Description"}

\item{dot.text.size}{dotplot text size}

\item{dot.scale}{dotplot scale}

\item{dot.x}{dotplot horizontal position}

\item{dot.y}{dotplot up and down position}

\item{dot.width}{dotplot width}

\item{dot.height}{dotplot height}

\item{...}{additional parameters}
}
\value{
figure
}
\description{
sankey and dot graphs showed the results of GO and KEGG analysis
}
\examples{
\dontrun{
data(xfbdf, package = "TCMNP")
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(tidyverse)
eg <- bitr(unique(xfbdf$target),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = "org.Hs.eg.db"
)
KK <- enrichKEGG(
  eg$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
KK <- setReadable(KK, "org.Hs.eg.db", keyType = "ENTREZID")
dot_sankey(KK)
dot_sankey(kegg.filter)
}
}
