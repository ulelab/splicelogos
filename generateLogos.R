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

get_ss_seq <- function(ss, genome, exon_nt, intron_nt, is_five_prime = TRUE) {
  pos_strand <- strand(ss) == "+"
  
  if (is_five_prime) {
    # Correct start and end positions for the 5' splice site
    start_pos <- ifelse(pos_strand, 
                        start(ss) - exon_nt + 1, 
                        end(ss) - intron_nt + 1) 
    end_pos <- ifelse(pos_strand, 
                      start(ss) + intron_nt, 
                      end(ss) + exon_nt)          
  } else {
    # Correct start and end positions for the 3' splice site
    start_pos <- ifelse(pos_strand, 
                        start(ss) - intron_nt,    
                        end(ss) - exon_nt)    
    end_pos <- ifelse(pos_strand, 
                      start(ss) + exon_nt - 1,   
                      end(ss) + intron_nt - 1)   
  }
  
  # Ensure proper genomic ranges
  ss_ranges <- GRanges(seqnames(ss),
                       IRanges(start_pos, end_pos),
                       strand = strand(ss))
  
  ss_seq <- getSeq(genome, ss_ranges)
  
  return(ss_seq)
}

# Get splice site sequences
five_ss_seq <- get_ss_seq(five_ss, Hsapiens, 3, 6, is_five_prime = TRUE)
three_ss_seq <- get_ss_seq(three_ss, Hsapiens, 3, 8, is_five_prime = FALSE)


# Convert to character vectors
five_ss <- as.character(five_ss_seq)
three_ss <- as.character(three_ss_seq)


# Filter for canonical splice sites
five_ss_canonical <- five_ss[substr(five_ss, 4, 5) == "GT"]
three_ss_canonical <- three_ss[substr(three_ss, 7, 8) == "AG"]  


# Generate sequence logos with axis lines for both x and y axes, and y-axis limited to 1
if (length(five_ss_canonical) > 0 && length(three_ss_canonical) > 0) {
  pdf("canonical_splice_site_logos_protein_coding.pdf", width = 15, height = 5)
  
  # 5'ss logo
  p1 <- ggseqlogo(five_ss_canonical, method = "prob", seq_type = "dna") +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          panel.grid = element_blank(),                      # Remove grid lines
          panel.border = element_blank(),                    # Remove plot border
          axis.line = element_line(color = "black", size = 0.8)) +  # Add black axis lines
    xlab("Position relative to splice site") +
    ylab("Probability") +
    scale_x_continuous(breaks = 1:9, labels = c(-3:-1, 0, 1:5)) +
    geom_vline(xintercept = 3.5, linetype = "dotted", color = "black", size = 1) +  # Thicker black dotted line
    coord_cartesian(ylim = c(0, 1))  # Limit y-axis to 1
  
  # 3'ss logo
  p2 <- ggseqlogo(three_ss_canonical, method = "prob", seq_type = "dna") +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 18),
          panel.grid = element_blank(),                      # Remove grid lines
          panel.border = element_blank(),                    # Remove plot border
          axis.line = element_line(color = "black", size = 0.8)) +  # Add black axis lines
    xlab("Position relative to splice site") +
    ylab("Probability") +
    scale_x_continuous(breaks = 1:11, labels = c(-7:-1, 0, 1:3)) +  # Corrected labels for 3'ss
    geom_vline(xintercept = 8.5, linetype = "dotted", color = "black", size = 1) +  # Thicker black dotted line
    coord_cartesian(ylim = c(0, 1))  # Limit y-axis to 1
  
  # Arrange plots side by side with different widths
  grid.arrange(p1, p2, ncol = 2, widths = c(2, 3))
  
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
