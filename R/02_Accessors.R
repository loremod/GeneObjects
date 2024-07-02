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
setGeneric("struct", function(x) standardGeneric("struct"))
#' @rdname struct
setMethod("struct", "Gene", function(x) x@struct)

# Define accessor functions for GeneStructure class

#' Get the exon ID of a GeneStructure
#'
#' @param x A GeneStructure object.
#' @return The exon ID of the gene structure.
#' @name exon_id
#' @aliases exon_id
#' @export
setGeneric("exon_id", function(x) standardGeneric("exon_id"))
#' @rdname exon_id
setMethod("exon_id", "GeneStructure", function(x) x@exon_id)

#' Get the number of exons in a Gene
#'
#' @param x A Gene object.
#' @return The number of exons.
#' @name numExons
#' @aliases numExons
#' @export
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
setGeneric("lengthProduct", function(x) standardGeneric("lengthProduct"))

#' @rdname lengthProduct
setMethod("lengthProduct", "Gene", function(x) nchar(productSequence(x)))
