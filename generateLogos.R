# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures", "BSgenome.Hsapiens.UCSC.hg38", "Biostrings", "ggseqlogo","RMariaDB","txdbmaker"))

library(RMariaDB)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(ggseqlogo)
library(txdbmaker)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Download the UCSC RefSeq gene track
txdb <- makeTxDbFromUCSC(genome = "hg38", tablename = "ncbiRefSeq")

# Define canonical chromosomes
canonical_chr <- paste0("chr", c(1:22, "X", "Y", "M"))

# Get all transcripts with their lengths
all_transcripts <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
all_transcripts$tx_length <- width(all_transcripts)

# Function to identify protein-coding transcripts based on RefSeq accession
is_protein_coding <- function(tx_name) {
  startsWith(tx_name, "NM_")
}

# Filter for protein-coding genes
protein_coding_transcripts <- all_transcripts[is_protein_coding(all_transcripts$tx_name)]

# Select the longest transcript per protein-coding gene
longest_per_gene <- protein_coding_transcripts %>%
  as.data.frame() %>%
  group_by(gene_id) %>%
  slice_max(order_by = tx_length, n = 1) %>%
  pull(tx_name)

# Extract exons for the longest transcripts
exons_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
longest_tx_exons <- exons_by_tx[names(exons_by_tx) %in% longest_per_gene]

# Process internal exons more efficiently
internal_exons <- unlist(longest_tx_exons)
exon_ranks <- rle(names(internal_exons))$lengths
is_internal <- rep(FALSE, length(internal_exons))
is_internal[cumsum(exon_ranks[-length(exon_ranks)]) + 1] <- TRUE
is_internal[cumsum(exon_ranks)] <- TRUE
internal_exons <- internal_exons[!is_internal]

# Filter for canonical chromosomes
internal_exons <- internal_exons[seqnames(internal_exons) %in% canonical_chr]


# Extract splice sites
five_ss <- GRanges()
three_ss <- GRanges()

# For positive strand
pos_exons <- internal_exons[strand(internal_exons) == "+"]
five_ss <- c(five_ss, resize(pos_exons, width = 1, fix = "end"))
three_ss <- c(three_ss, resize(pos_exons, width = 1, fix = "start"))

# For negative strand
neg_exons <- internal_exons[strand(internal_exons) == "-"]
five_ss <- c(five_ss, resize(neg_exons, width = 1, fix = "start"))
three_ss <- c(three_ss, resize(neg_exons, width = 1, fix = "end"))

# Define function to get splice site sequences
get_ss_seq <- function(ss, genome, upstream, downstream) {
  pos_strand <- strand(ss) == "+"
  
  start_pos <- ifelse(pos_strand, 
                      start(ss) - upstream + 1, 
                      start(ss) - downstream)
  end_pos <- ifelse(pos_strand, 
                    end(ss) + downstream, 
                    end(ss) + upstream - 1)
  
  ss_ranges <- GRanges(seqnames(ss),
                       IRanges(start_pos, end_pos),
                       strand = strand(ss))
  
  ss_seq <- getSeq(genome, ss_ranges)
  
  # Reverse complement sequences on negative strand
  ss_seq[!pos_strand] <- reverseComplement(ss_seq[!pos_strand])
  
  return(ss_seq)
}

# Get splice site sequences
five_ss_seq <- get_ss_seq(five_ss, Hsapiens, 3, 6)
three_ss_seq <- get_ss_seq(three_ss, Hsapiens, 20, 3)

# Convert to character vectors
five_ss <- as.character(five_ss_seq)
three_ss <- as.character(three_ss_seq)

# Diagnostic step: Print the first few sequences and their lengths
cat("First few 5' splice site sequences:\n")
print(head(five_ss))
cat("\nLength of 5' splice site sequences:", unique(nchar(five_ss)), "\n")

cat("\nFirst few 3' splice site sequences:\n")
print(head(three_ss))
cat("\nLength of 3' splice site sequences:", unique(nchar(three_ss)), "\n")

# Count occurrences of different dinucleotides at splice sites
five_ss_dinucleotides <- substr(five_ss, 4, 5)
three_ss_dinucleotides <- substr(three_ss, 19, 20)  # Changed from 20, 21 to 19, 20

cat("\n5' splice site dinucleotide counts:\n")
print(table(five_ss_dinucleotides))

cat("\n3' splice site dinucleotide counts:\n")
print(table(three_ss_dinucleotides))

# Filter for canonical splice sites
five_ss_canonical <- five_ss[substr(five_ss, 4, 5) == "GT"]
three_ss_canonical <- three_ss[substr(three_ss, 18, 19) == "AG"]  # Changed from 20, 21 to 19, 20

# Generate sequence logos
if (length(five_ss_canonical) > 0 && length(three_ss_canonical) > 0) {
  pdf("canonical_splice_site_logos_protein_coding.pdf", width = 12, height = 6)
  
  # 5'ss logo
  p1 <- ggseqlogo(five_ss_canonical, method = "prob", seq_type = "dna") +
    ggtitle("Human Canonical 5' Splice Site (GT) in Protein-Coding Genes") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Position relative to splice site") +
    ylab("Probability")
  
  # 3'ss logo
  p2 <- ggseqlogo(three_ss_canonical, method = "prob", seq_type = "dna") +
    ggtitle("Human Canonical 3' Splice Site (AG) in Protein-Coding Genes") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Position relative to splice site") +
    ylab("Probability")
  
  # Arrange plots side by side
  gridExtra::grid.arrange(p1, p2, ncol = 2)
  
  dev.off()
  print("Canonical sequence logos for protein-coding genes have been generated and saved as 'canonical_splice_site_logos_protein_coding.pdf'")
} else {
  cat("\nError: Unable to generate sequence logos due to lack of canonical splice sites.\n")
}

# Print some statistics
cat("\nTotal genes:", length(unique(all_transcripts$gene_id)), "\n")
cat("Total transcripts:", length(all_transcripts), "\n")
cat("Longest transcripts selected:", length(longest_per_gene), "\n")
cat("Total internal exons:", length(internal_exons), "\n")
cat("Total 5' splice sites:", length(five_ss), "\n")
cat("Canonical 5' splice sites (GT):", length(five_ss_canonical), "\n")
cat("Percentage of canonical 5' splice sites:", round(length(five_ss_canonical) / length(five_ss) * 100, 2), "%\n\n")

cat("Total 3' splice sites:", length(three_ss), "\n")
cat("Canonical 3' splice sites (AG):", length(three_ss_canonical), "\n")
cat("Percentage of canonical 3' splice sites:", round(length(three_ss_canonical) / length(three_ss) * 100, 2), "%\n")

print("Canonical sequence logos have been generated and saved as 'canonical_splice_site_logos.pdf'")
