# Load necessary packages
#' @importFrom methods new
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors Rle
#' @importFrom Biostrings AAString RNAString
NULL

# Define the GeneStructure class
setClass(
  "GeneStructure",
  contains = "GRanges"
)

#' Create a GeneStructure Object
#'
#' @param seqnames A character vector of sequence names.
#' @param start A numeric vector of start positions.
#' @param end A numeric vector of end positions.
#' @param strand A character vector of strand information.
#' @param exon_id A vector of exon IDs.
#' @examples
#' library(GenomicRanges)
#' gene_structure <- GeneStructure(seqnames = "chr1", start = c(100, 230),
#'                                 end = c(200, 240), strand = "+",
#'                                 exon_id = 1:2)
#' gene_structure
#' @return A GeneStructure object.
#' @export
GeneStructure <- function(seqnames, start, end, strand, exon_id) {
  # Perform consistency checks
  if (length(unique(seqnames)) > 1) {
    stop("All exons must be on the same chromosome.")
  }

  if (length(unique(strand)) > 1) {
    stop("All exons must be on the same strand.")
  }

  gr <- GRanges(
    seqnames = Rle(seqnames),
    ranges = IRanges(start = start, end = end),
    strand = Rle(strand),
    exon_id = exon_id
  )

  new("GeneStructure", gr)
}

# Define the virtual Gene class
#' @name Gene-class
#' @aliases Gene
#' @title Virtual Class Gene
#' @description A virtual class to represent a general gene structure.
#' @slot id A character string representing the gene ID.
#' @slot symbol A character string representing the gene symbol.
#' @slot name A character string representing the gene name.
#' @slot description A character string describing the gene.
#' @slot struct A GeneStructure object representing the gene structure.
#' @return It is not possible to instantiate a Gene Object since Gene is a virtual class
#' @seealso \code{\linkS4class{CodingGene}}, \code{\linkS4class{LncRNAGene}}, \code{\linkS4class{MicroRNAGene}}
#' @export
setClass("Gene",
         representation(
           id = "character",
           symbol = "character",
           name = "character",
           description = "character",
           struct = "GeneStructure"
         ),
         contains = "VIRTUAL"
)


# Define specific gene classes
#' @name CodingGene-class
#' @aliases CodingGene
#' @title Class CodingGene
#' @description A class to represent protein-coding genes.
#' @slot proteinID A character string representing the protein ID.
#' @slot proteinSequence A AAString representing the protein sequence.
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' gen_struct <- GeneStructure(seqnames = "chr1", start = c(100, 250),
#'                             end = c(200, 300), strand = "+", exon_id = 1:2)
#' coding_gene <- CodingGene(id = "gene1", symbol = "GENE1",
#'                    name = "Gene 1",
#'                    description = "A protein-coding gene",
#'                    struct = gen_struct, proteinID = "P12345",
#'                    proteinSequence = "MTEYKLVVVG")
#' coding_gene
#' @return A CodingGene instance
#' @seealso \code{\linkS4class{Gene}}
#' @export
setClass("CodingGene",
         contains = "Gene",
         slots = list(
           proteinID = "character",
           proteinSequence = "AAString"
         ))

#' @export
CodingGene <- function(id, symbol, name, description, struct, proteinID, proteinSequence) {
  if (!inherits(proteinSequence, "AAString")) {
    proteinSequence <- AAString(proteinSequence)
  }
  new("CodingGene", id = id, symbol = symbol, name = name, description = description, struct = struct, proteinID = proteinID, proteinSequence = proteinSequence)
}

#' @name LncRNAGene-class
#' @aliases LncRNAGene
#' @title Class LncRNAGene
#' @description A class to represent long non-coding RNA genes.
#' @slot lncRNAID A character string representing the lncRNA ID.
#' @slot RNASequence A RNAString representing the RNA sequence.
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' gen_struct <- GeneStructure(seqnames = "chr2", start = c(300, 472),
#'                             end = c(400, 900), strand = "-", exon_id = 3:4)
#' lncrna_gene <- LncRNAGene(id = "gene2", symbol = "LNC1",
#'                    name = "LncRNA 1",
#'                    description = "A long non-coding RNA gene",
#'                    struct = gen_struct, lncRNAID = "LNC123",
#'                    RNASequence = "AUGCUACG")
#' lncrna_gene
#' @return A LncRNAGene instance
#' @seealso \code{\linkS4class{Gene}}
#' @export
setClass("LncRNAGene",
         contains = "Gene",
         slots = list(
           lncRNAID = "character",
           RNASequence = "RNAString"
         ))

#' @export
LncRNAGene <- function(id, symbol, name, description, struct, lncRNAID, RNASequence) {
  if (!inherits(RNASequence, "RNAString")) {
    RNASequence <- RNAString(RNASequence)
  }
  new("LncRNAGene", id = id, symbol = symbol, name = name, description = description, struct = struct, lncRNAID = lncRNAID, RNASequence = RNASequence)
}

#' @name MicroRNAGene-class
#' @aliases MicroRNAGene
#' @title Class MicroRNAGene
#' @description A class to represent microRNA genes.
#' @slot microRNAID A character string representing the microRNA ID.
#' @slot microRNASeedSequence A RNAString representing the microRNA seed sequence.
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' gen_struct <- GeneStructure(seqnames = "chr1", start = c(100, 250),
#'                             end = c(200, 300), strand = "+", exon_id = 1:2)
#' microrna_gene <- MicroRNAGene(id = "gene3", symbol = "MIR1",
#'                      name = "MicroRNA 1", description = "A microRNA gene",
#'                      struct = gen_struct, microRNAID = "MIR123",
#'                      microRNASeedSequence = "ACGUGA")
#' microrna_gene
#' @return A MicroRNAGene instance
#' @seealso \code{\linkS4class{Gene}}
#' @export
setClass("MicroRNAGene",
         contains = "Gene",
         slots = list(
           microRNAID = "character",
           microRNASeedSequence = "RNAString"
         ))

#' @export
MicroRNAGene <- function(id, symbol, name, description, struct, microRNAID, microRNASeedSequence) {
  if (!inherits(microRNASeedSequence, "RNAString")) {
    microRNASeedSequence <- RNAString(microRNASeedSequence)
  }
  new("MicroRNAGene", id = id, symbol = symbol, name = name, description = description, struct = struct, microRNAID = microRNAID, microRNASeedSequence = microRNASeedSequence)
}

# Define additional gene classes
#' @name RibosomalRNA-class
#' @aliases RibosomalRNA
#' @title Class RibosomalRNA
#' @description A class to represent ribosomal RNA genes.
#' @slot rRNAID A character string representing the rRNA ID.
#' @slot rRNASequence A RNAString representing the rRNA sequence.
#' @return A RibosomalRNA instance
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' gen_struct <- GeneStructure(seqnames = "chr1", start = c(100, 250),
#'                             end = c(200, 300), strand = "+", exon_id = 1:2)
#' rrna_gene <- RibosomalRNA(id = "gene4", symbol = "RIB1",
#'                  name = "Ribosomal RNA 1",
#'                  description = "A ribosomal RNA gene",
#'                  struct = gen_struct, rRNAID = "RIB123",
#'                  rRNASequence = "AGGCUAG")
#' rrna_gene
#' @seealso \code{\linkS4class{Gene}}
#' @export
setClass("RibosomalRNA",
         contains = "Gene",
         slots = list(
           rRNAID = "character",
           rRNASequence = "RNAString"
         ))

#' @export
RibosomalRNA <- function(id, symbol, name, description, struct, rRNAID, rRNASequence) {
  if (!inherits(rRNASequence, "RNAString")) {
    rRNASequence <- RNAString(rRNASequence)
  }
  new("RibosomalRNA", id = id, symbol = symbol, name = name, description = description, struct = struct, rRNAID = rRNAID, rRNASequence = rRNASequence)
}

#' @name SmallNuclearRNA-class
#' @aliases SmallNuclearRNA
#' @title Class SmallNuclearRNA
#' @description A class to represent small nuclear RNA genes.
#' @slot snRNAID A character string representing the snRNA ID.
#' @slot snRNASequence A RNAString representing the snRNA sequence.
#' @return A SmallNuclearRNA instance
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' gen_struct <- GeneStructure(seqnames = "chr1", start = c(100, 250),
#'                             end = c(200, 300), strand = "+", exon_id = 1:2)
#' snrna_gene <- SmallNuclearRNA(id = "gene5", symbol = "SNR1",
#'                   name = "Small Nuclear RNA 1",
#'                   description = "A small nuclear RNA gene",
#'                   struct = gen_struct, snRNAID = "SNR123",
#'                   snRNASequence = "GCUAGCU")
#' snrna_gene
#' @seealso \code{\linkS4class{Gene}}
#' @export
setClass("SmallNuclearRNA",
         contains = "Gene",
         slots = list(
           snRNAID = "character",
           snRNASequence = "RNAString"
         ))

#' @export
SmallNuclearRNA <- function(id, symbol, name, description, struct, snRNAID, snRNASequence) {
  if (!inherits(snRNASequence, "RNAString")) {
    snRNASequence <- RNAString(snRNASequence)
  }
  new("SmallNuclearRNA", id = id, symbol = symbol, name = name, description = description, struct = struct, snRNAID = snRNAID, snRNASequence = snRNASequence)
}

#' @name PiwiInteractingRNA-class
#' @aliases PiwiInteractingRNA
#' @title Class PiwiInteractingRNA
#' @description A class to represent piwi-interacting RNA genes.
#' @slot piRNAID A character string representing the piRNA ID.
#' @slot piRNASequence A RNAString representing the piRNA sequence.
#' @examples
#' library(GenomicRanges)
#' library(Biostrings)
#' gen_struct <- GeneStructure(seqnames = "chr1", start = c(100, 250),
#'                             end = c(200, 300), strand = "+", exon_id = 1:2)
#' pirna_gene <- PiwiInteractingRNA(id = "gene6", symbol = "PIWI1",
#'                   name = "Piwi-Interacting RNA 1",
#'                   description = "A piwi-interacting RNA gene",
#'                   struct = gen_struct, piRNAID = "PIWI123",
#'                   piRNASequence = "UCGUAU")
#' pirna_gene
#' @return A PiwiInteractingRNA instance
#' @seealso \code{\linkS4class{Gene}}
#' @export
setClass("PiwiInteractingRNA",
         contains = "Gene",
         slots = list(
           piRNAID = "character",
           piRNASequence = "RNAString"
         ))

#' @export
PiwiInteractingRNA <- function(id, symbol, name, description, struct, piRNAID, piRNASequence) {
  if (!inherits(piRNASequence, "RNAString")) {
    piRNASequence <- RNAString(piRNASequence)
  }
  new("PiwiInteractingRNA", id = id, symbol = symbol, name = name, description = description, struct = struct, piRNAID = piRNAID, piRNASequence = piRNASequence)
}
