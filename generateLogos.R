# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures", "BSgenome.Hsapiens.UCSC.hg38", "Biostrings", "ggseqlogo"))

library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(ggseqlogo)

# Download the UCSC known gene track
txdb <- makeTxDbFromUCSC(genome = "hg38", tablename = "knownGene")

# Extract splice sites
splice_sites <- bindings(txdb)

# Define functions to get splice site sequences
get_5ss <- function(splice_sites, genome, n_upstream = 3, n_downstream = 6) {
    five_ss <- GRanges(seqnames(splice_sites),
                       IRanges(start(splice_sites) - n_upstream,
                               end = start(splice_sites) + n_downstream),
                       strand = strand(splice_sites))
    five_ss_seq <- getSeq(genome, five_ss)
    return(five_ss_seq)
}

get_3ss <- function(splice_sites, genome, n_upstream = 20, n_downstream = 3) {
    three_ss <- GRanges(seqnames(splice_sites),
                        IRanges(end(splice_sites) - n_upstream,
                                end = end(splice_sites) + n_downstream),
                        strand = strand(splice_sites))
    three_ss_seq <- getSeq(genome, three_ss)
    return(three_ss_seq)
}

# Get splice site sequences
five_ss_seq <- get_5ss(splice_sites, Hsapiens)
three_ss_seq <- get_3ss(splice_sites, Hsapiens)

# Convert to character vectors and remove any sequences with 'N'
five_ss <- as.character(five_ss_seq)
five_ss <- five_ss[!grepl("N", five_ss)]

three_ss <- as.character(three_ss_seq)
three_ss <- three_ss[!grepl("N", three_ss)]

# Limit to 1000 sequences for each (for performance reasons)
five_ss <- sample(five_ss, min(1000, length(five_ss)))
three_ss <- sample(three_ss, min(1000, length(three_ss)))

# Generate sequence logos
pdf("splice_site_logos.pdf", width = 10, height = 5)

par(mfrow = c(1, 2))

# 5'ss logo
ggseqlogo(five_ss, method = "bits", seq_type = "dna") +
  ggtitle("Human 5' Splice Site (5'ss)")

# 3'ss logo
ggseqlogo(three_ss, method = "bits", seq_type = "dna") +
  ggtitle("Human 3' Splice Site (3'ss)")

dev.off()

print("Sequence logos have been generated and saved as 'splice_site_logos.pdf'")
