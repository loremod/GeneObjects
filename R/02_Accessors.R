# Load necessary packages
#' @importFrom methods setGeneric setMethod
NULL

# Define accessor functions for Gene class

#' Get the ID of a Gene
#'
#' @param x A Gene object.
#' @return The ID of the gene.
#' @name id
#' @aliases id
#' @export
#' @return A character string representing the gene ID.
#' @examples
#' #I create a GeneStructure Object
#' gen_struct <- GeneStructure(seqnames = "chr1", start = 100, end = 500,
#'                             strand = "+", exon_id = 1)
#'
#' #I consider the class LncRNAGene, which inherited this method from its superclass Gene
#' gene <- new("LncRNAGene", id = "gene1", symbol = "GENE1", name = "Gene 1",
#'             description = "A test gene", struct = gen_struct)
#' id(gene)
setGeneric("id", function(x) standardGeneric("id"))
#' @rdname id
setMethod("id", "Gene", function(x) x@id)

#' Get the symbol of a Gene
#'
#' @param x A Gene object.
#' @return The symbol of the gene.
#' @name symbol
#' @aliases symbol
#' @export
#' @return A character string representing the gene symbol.
#' @examples
#' #I create a GeneStructure Object
#' gen_struct <- GeneStructure(seqnames = "chr1", start = 100, end = 500,
#'                             strand = "+", exon_id = 1)
#'
#' #I consider the class CodingGene, which inherited this method from its superclass Gene
#' gene <- new("CodingGene", id = "gene1", symbol = "GENE1", name = "Gene 1",
#'             description = "A test gene", struct = gen_struct)
#' symbol(gene)
setGeneric("symbol", function(x) standardGeneric("symbol"))
#' @rdname symbol
setMethod("symbol", "Gene", function(x) x@symbol)

#' Get the name of a Gene
#'
#' @param x A Gene object.
#' @return The name of the gene.
#' @name name
#' @aliases name
#' @export
#' @return A character string representing the gene name.
#' @examples
#' #I create a GeneStructure Object
#' gen_struct <- GeneStructure(seqnames = "chr1", start = 100, end = 500, strand = "+", exon_id = 1)
#'
#' #I consider the class MicroRNAGene, which inherited this method from its superclass Gene
#' gene <- new("MicroRNAGene", id = "gene1", symbol = "GENE1", name = "Gene 1",
#'            description = "A test gene", struct = gen_struct)
#' name(gene)
setGeneric("name", function(x) standardGeneric("name"))
#' @rdname name
setMethod("name", "Gene", function(x) x@name)

#' Get the description of a Gene
#'
#' @param x A Gene object.
#' @return The description of the gene.
#' @name description
#' @aliases description
#' @export
#' @return A character string describing the gene.
#' @examples
#' #I create a GeneStructure Object
#' gen_struct <- GeneStructure(seqnames = "chr1", start = 100, end = 500,
#'                              strand = "+", exon_id = 1)
#'
#' #I consider the class CodingGene, which inherited this method from its superclass Gene
#' gene <- new("CodingGene", id = "gene1", symbol = "GENE1", name = "Gene 1",
#'             description = "A test gene", struct = gen_struct)
#' description(gene)
setGeneric("description", function(x) standardGeneric("description"))
#' @rdname description
setMethod("description", "Gene", function(x) x@description)

#' Get the structure of a Gene
#'
#' @param x A Gene object.
#' @return The structure of the gene.
#' @name struct
#' @aliases struct
#' @export
#' @return A GeneStructure object representing the gene structure.
#' @examples
#' #I create a GeneStructure Object
#' gen_struct <- GeneStructure(seqnames = "chr1", start = 100, end = 500,
#'                             strand = "+", exon_id = 1)
#'
#' #I consider the class RibosomalRNA, which inherited the method struct from its superclass Gene
#' gene <- new("RibosomalRNA", id = "gene1", symbol = "GENE1",
#'            name = "Gene 1", description = "A test gene", struct = gen_struct)
#' struct(gene)
setGeneric("struct", function(x) standardGeneric("struct"))
#' @rdname struct
setMethod("struct", "Gene", function(x) x@struct)


#' Get the number of exons in a Gene
#'
#' @param x A Gene object.
#' @return The number of exons.
#' @name numExons
#' @aliases numExons
#' @export
#' @return An integer representing the number of exons.
#' @examples
#' #I create a GeneStructure Object
#' gen_struct <- GeneStructure(seqnames = "chr1", start = 100, end = 500,
#'                            strand = "+", exon_id = 1)
#' gene <- new("RibosomalRNA", id = "gene1", symbol = "GENE1", name = "Gene 1",
#'            description = "A test gene", struct = gen_struct)
#' numExons(gene)
setGeneric("numExons", function(x) standardGeneric("numExons"))
#' @rdname numExons
setMethod("numExons", "Gene", function(x) length(unique(x@struct@ranges)))

#' Get the product sequence of a Gene
#'
#' @param x A Gene object.
#' @return The product sequence.
#' @name productSequence
#' @aliases productSequence
#' @export
#' @return A character string representing the product sequence.
#' @examples
#' #I create a GeneStructure Object
#' gen_struct <- GeneStructure(seqnames = "chr1", start = 100, end = 500,
#'                             strand = "+", exon_id = 1)
#' coding_gene <- new("CodingGene", id = "gene1", symbol = "GENE1",
#'                    name = "Gene 1", description = "A protein-coding gene",
#'                    struct = gen_struct, proteinID = "P12345",
#'                    proteinSequence = "MTEYKLVVVG")
#' productSequence(coding_gene)
setGeneric("productSequence", function(x) standardGeneric("productSequence"))

#' @rdname productSequence
setMethod("productSequence", "CodingGene", function(x) x@proteinSequence)

#' @rdname productSequence
setMethod("productSequence", "LncRNAGene", function(x) x@RNASequence)

#' @rdname productSequence
setMethod("productSequence", "MicroRNAGene", function(x) x@microRNASeedSequence)

#' @rdname productSequence
setMethod("productSequence", "RibosomalRNA", function(x) x@rRNASequence)

#' @rdname productSequence
setMethod("productSequence", "SmallNuclearRNA", function(x) x@snRNASequence)

#' @rdname productSequence
setMethod("productSequence", "PiwiInteractingRNA", function(x) x@piRNASequence)

#' Get the length of the product sequence of a Gene
#'
#' @param x A Gene object.
#' @return The length of the product sequence.
#' @name lengthProduct
#' @aliases lengthProduct
#' @export
#' @return An integer representing the length of the product sequence.
#' @examples
#' #I create a GeneStructure Object
#' gen_struct <- GeneStructure(seqnames = "chr1", start = 100, end = 500,
#'                             strand = "+", exon_id = 1)
#' coding_gene <- new("CodingGene", id = "gene1", symbol = "GENE1",
#'                      name = "Gene 1", description = "A protein-coding gene",
#'                      struct = gen_struct, proteinID = "P12345",
#'                      proteinSequence = "MTEYKLVVVG")
#' lengthProduct(coding_gene)
setGeneric("lengthProduct", function(x) standardGeneric("lengthProduct"))

#' @rdname lengthProduct
setMethod("lengthProduct", "Gene", function(x) nchar(productSequence(x)))
