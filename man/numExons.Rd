% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_Accessors.R
\name{numExons}
\alias{numExons}
\alias{numExons,Gene-method}
\title{Get the number of exons in a Gene}
\usage{
numExons(x)

\S4method{numExons}{Gene}(x)
}
\arguments{
\item{x}{A Gene object.}
}
\value{
The number of exons.

An integer representing the number of exons.
}
\description{
Get the number of exons in a Gene
}
\examples{
library(GenomicRanges)
gen_struct <- GeneStructure(seqnames = "chr1", start = 100, end = 500,
                            strand = "+", exon_id = 1)
gene <- RibosomalRNA(id = "gene1", symbol = "GENE1", name = "Gene 1",
                     description = "A test gene", struct = gen_struct,
                     rRNAID = "RIB1", rRNASequence = "AGGCUAG")
numExons(gene)
}
