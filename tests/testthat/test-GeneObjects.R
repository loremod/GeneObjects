library(testthat)
library(GenomicRanges)
library(GeneObjects)

test_that("CodingGene class functionality", {
  structure <- GeneStructure(seqnames = "chr13", start = c(32889611, 32973801), end = c(32973800, 32973805), strand = "-", exon_id = c("exon1", "exon2"))
  gene <- CodingGene(id = "ENSG00000139618", symbol = "BRCA2", name = "BRCA2",
                     description = "Breast cancer type 2 susceptibility protein",
                     struct = structure, proteinID = "P51587", proteinSequence = "MENQWEK")
  expect_equal(id(gene), "ENSG00000139618")
  expect_equal(symbol(gene), "BRCA2")
  expect_equal(name(gene), "BRCA2")
  expect_equal(description(gene), "Breast cancer type 2 susceptibility protein")
  expect_true(is(struct(gene), "GeneStructure"))
  expect_equal(productSequence(gene), "MENQWEK")
  expect_equal(lengthProduct(gene), nchar("MENQWEK"))
  expect_equal(numExons(gene), 2)
})

test_that("LncRNAGene class functionality", {
  structure <- GeneStructure(seqnames = "chr11", start = c(65270, 66700, 66838), end = c(66716, 66716, 66879), strand = "+", exon_id = c("exon1", "exon2", "exon3"))
  gene <- LncRNAGene(id = "ENSG00000273105", symbol = "MALAT1", name = "MALAT1",
                     description = "Metastasis Associated Lung Adenocarcinoma Transcript 1",
                     struct = structure, lncRNAID = "NONHSAG000347", RNASequence = "GATTACA")
  expect_equal(id(gene), "ENSG00000273105")
  expect_equal(symbol(gene), "MALAT1")
  expect_equal(name(gene), "MALAT1")
  expect_equal(description(gene), "Metastasis Associated Lung Adenocarcinoma Transcript 1")
  expect_true(is(struct(gene), "GeneStructure"))
  expect_equal(productSequence(gene), "GATTACA")
  expect_equal(lengthProduct(gene), nchar("GATTACA"))
  expect_equal(numExons(gene), 3)
})

test_that("MicroRNAGene class functionality", {
  structure <- GeneStructure(seqnames = "chr17", start = c(59109, 59350), end = c(59370, 59370), strand = "+", exon_id = c("exon1", "exon2"))
  gene <- MicroRNAGene(id = "ENSG00000284458", symbol = "MIR21", name = "MicroRNA 21",
                       description = "microRNA 21", struct = structure,
                       microRNAID = "MI0000077", microRNASeedSequence = "UAGCUU")
  expect_equal(id(gene), "ENSG00000284458")
  expect_equal(symbol(gene), "MIR21")
  expect_equal(name(gene), "MicroRNA 21")
  expect_equal(description(gene), "microRNA 21")
  expect_true(is(struct(gene), "GeneStructure"))
  expect_equal(productSequence(gene), "UAGCUU")
  expect_equal(lengthProduct(gene), nchar("UAGCUU"))
  expect_equal(numExons(gene), 2)
})

test_that("RibosomalRNA class functionality", {
  structure <- GeneStructure(seqnames = "chr1", start = c(69091, 69900), end = c(70008, 70008), strand = "+", exon_id = c("exon1", "exon2"))
  gene <- RibosomalRNA(id = "ENSG00000175899", symbol = "RNA5-8S5", name = "RNA, 5.8S ribosomal",
                       description = "RNA, 5.8S ribosomal RNA", struct = structure,
                       rRNAID = "NR_023363", rRNASequence = "ACCTC")
  expect_equal(id(gene), "ENSG00000175899")
  expect_equal(symbol(gene), "RNA5-8S5")
  expect_equal(name(gene), "RNA, 5.8S ribosomal")
  expect_equal(description(gene), "RNA, 5.8S ribosomal RNA")
  expect_true(is(struct(gene), "GeneStructure"))
  expect_equal(productSequence(gene), "ACCTC")
  expect_equal(lengthProduct(gene), nchar("ACCTC"))
  expect_equal(numExons(gene), 2)
})

test_that("SmallNuclearRNA class functionality", {
  structure <- GeneStructure(seqnames = "chr6", start = c(29945928, 29946000), end = c(29946086, 29946086), strand = "-", exon_id = c("exon1", "exon2"))
  gene <- SmallNuclearRNA(id = "ENSG00000283297", symbol = "RNU6-1", name = "RNA, U6 small nuclear 1",
                          description = "RNA, U6 small nuclear 1", struct = structure,
                          snRNAID = "NR_004394", snRNASequence = "GGGCC")
  expect_equal(id(gene), "ENSG00000283297")
  expect_equal(symbol(gene), "RNU6-1")
  expect_equal(name(gene), "RNA, U6 small nuclear 1")
  expect_equal(description(gene), "RNA, U6 small nuclear 1")
  expect_true(is(struct(gene), "GeneStructure"))
  expect_equal(productSequence(gene), "GGGCC")
  expect_equal(lengthProduct(gene), nchar("GGGCC"))
  expect_equal(numExons(gene), 2)
})

test_that("PiwiInteractingRNA class functionality", {
  structure <- GeneStructure(seqnames = "chr1", start = c(17368, 17400), end = c(17436, 17436), strand = "-", exon_id = c("exon1", "exon2"))
  gene <- PiwiInteractingRNA(id = "ENSG00000284332", symbol = "PIWIL4", name = "Piwi like RNA",
                             description = "Piwi like RNA", struct = structure,
                             piRNAID = "DQ571952", piRNASequence = "CTGCTG")
  expect_equal(id(gene), "ENSG00000284332")
  expect_equal(symbol(gene), "PIWIL4")
  expect_equal(name(gene), "Piwi like RNA")
  expect_equal(description(gene), "Piwi like RNA")
  expect_true(is(struct(gene), "GeneStructure"))
  expect_equal(productSequence(gene), "CTGCTG")
  expect_equal(lengthProduct(gene), nchar("CTGCTG"))
  expect_equal(numExons(gene), 2)
})
