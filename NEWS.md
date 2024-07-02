# GeneObjects Package News

## Version 0.0.1.1000
- Initial release of the `GeneObjects` package.
- Provides a set of S4 classes to represent various types of genes, including:
  - Protein-coding genes (`CodingGene`)
  - Long non-coding RNA genes (`LncRNAGene`)
  - MicroRNA genes (`MicroRNAGene`)
  - Ribosomal RNA genes (`RibosomalRNA`)
  - Small nuclear RNA genes (`SmallNuclearRNA`)
  - Piwi-interacting RNA genes (`PiwiInteractingRNA`)
- Includes functionality to:
  - Create gene structures using the `GeneStructure` class.
  - Access gene attributes like `id`, `symbol`, `name`, `description`, `struct`, `numExons`, `productSequence`, and `lengthProduct`.
