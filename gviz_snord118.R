# ================================ /
# LCC Genetics /
# Author: Nils Briel, Dr., Dept. of Neurology, USZ, CH /
# Date: 2025-08-13 /
# ================================ /

# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Gviz", "biomaRt", "GenomicRanges", "VariantAnnotation"))

library(Gviz)
library(biomaRt)
library(GenomicRanges)
library(VariantAnnotation)
library(tidyverse)

# visualize the region with SNORD118 annotation and ClinVar variants

  # Define the region of interest
  chr <- "chr17"
  start <- 8173366
  end <- 8173674
  
  # Fetch SNORD118 annotation using biomaRt
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  snord118 <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "strand", 
                                   "external_gene_name",
                                   "transcript_start","transcript_end","transcription_start_site"),
                    filters = "external_gene_name",
                    values = "SNORD118",
                    mart = ensembl)

  # Create a GRanges object for SNORD118
  snord118_gr <- GRanges(seqnames = Rle(snord118$chromosome_name),
                         ranges = IRanges(start = snord118$start_position, end = snord118$end_position),
                         strand = Rle(snord118$strand),
                         gene = snord118$external_gene_name)
  
  # Fetch ClinVar variants associated with Labrune Syndrome
  clinvar_url <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
  clinvar_vcf <- readVcf(clinvar_url, "hg38", param = GRanges(seqnames = Rle(snord118$chromosome_name),
                                                              ranges = IRanges(start = start, end = end),
                                                              strand = Rle(snord118$strand)))
  clinvar_gr <- rowRanges(clinvar_vcf)
  clinvar_info <- info(clinvar_vcf)
  
  # Add combined metadata to the GRanges object
  mcols(clinvar_gr) <- cbind(mcols(clinvar_gr), clinvar_info[match(names(clinvar_gr), rownames(clinvar_info)), ])
  
  # Create annotation tracks
  itrack <- IdeogramTrack(genome = "hg38", chromosome = chr)
  genome_axis <- GenomeAxisTrack()
  
  snord118_track <- BiomartGeneRegionTrack(genome = "hg38", name = "SNORD118|TMEM107 locus", 
                                           symbol = "SNORD118", biomart = ensembl,
                                           transcriptAnnotation = "symbol")
  
  snord118_track2 <-  GeneRegionTrack(snord118_gr, chromosome = chr, 
                                      transcriptAnnotation = "symbol",
                                      collapseTranscripts = TRUE, 
                                      shape = "arrow",
                                      transcriptAnnotation = "symbol",
                                      name = "SNORD118", genome = "hg38")
  
  clinvar_filtered <- clinvar_gr[clinvar_gr$CLNVC %in% c("single_nucleotide_variant") &
             unlist(grepl("Leukoencephalopathy_with_calcifications_and_cysts|Labrune", clinvar_gr$CLNDN)),]

  # generate Lollipop plot from VCF
  library(trackViewer)
  snord118_gr$fill <- "blue3"
  clinvar_filtered$color <- as.numeric(as.factor(unlist(clinvar_filtered$CLNSIG)))
  legends <- list(labels=levels((as.factor(unlist(clinvar_filtered$CLNSIG)))), 
                  fill=c("blue",'#87CEFA',"yellow","lightblue","red3","red","grey"))
  clinvar_filtered$color <- legends$fill[clinvar_filtered$color]
  clinvar_filtered$height <- 20
  names(clinvar_filtered) <- paste0(clinvar_filtered$CLNHGVS)
  clinvar_filtered$label <- as.character(1:length(clinvar_filtered))
  clinvar_filtered$SNPsideID <- "bottom"
  pdf("../results/trackviewer_lollipop.pdf", width = 10, height = 10)
  lolliplot(clinvar_filtered, snord118_gr, ranges = snord118_gr, label_on_feature=T, yaxis=T, legend = legends)
  dev.off()

  # plot Allel-Frequency
  dTrack <- DataTrack(clinvar_filtered, name = "AF_EXAC", type = c("a", "p"))
  
  # Plot the tracks
  pdf("../results/trackviewer_gviz.pdf", width = 7, height = 3.5)
  plotTracks(list(itrack, genome_axis, snord118_track),add53 = TRUE, add35 = TRUE, from = snord118_track@start, to = snord118_track@end, chromosome = chr)  
  plotTracks(list(itrack, genome_axis, snord118_track2, dTrack), add53 = TRUE, add35 = TRUE, littleTicks = T, from = start, to = end, chromosome = chr)
  dev.off()
  
  
  # For patient-specific mutations reported as aligned to hg19 -> UCSC liftover: https://genome.ucsc.edu/cgi-bin/hgLiftOver on 2025-01-13
  # Original: chr17:8076762G>A; chr17:8076770G>C
  # Converted: chr17:8173444G>A; chr17:8173452G>C
  
  