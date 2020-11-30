library(GenomicFeatures)

#' Prepare gtf annotation
#'
#' Split GTF annotation into intron, exon, 3'UTR and 5'UTR. The exonic parts 
#' that overlap with UTRs are counted as UTRs! Intronic parts that overlap 
#' with exons or UTRs are not considered introns.
#'
#' @param gtf GRanges object
#'
#' @return GRangesList with exon, intron, 3'UTR and 5'UTR annotation
#' @export
#' 
#' @importFrom GenomicFeatures makeTxDbFromGRanges
#'
#' @examples
prepare_anno <- function(gtf){
  exon <- gtf[gtf$type == "exon"] %>% unique
  gene <- gtf[gtf$type == "gene"] %>% unique
  five_utr <- gtf[gtf$type == "five_prime_utr"] %>% unique
  three_utr <- gtf[gtf$type == "three_prime_utr"] %>% unique
  
  ## We remove all 3' and 5' UTR regions that overlap with any exons
  exon_utr <- GenomicRanges::setdiff(exon, three_utr)
  exon_unique <- GenomicRanges::setdiff(exon_utr, five_utr) %>% unique
  anno <- GRangesList(gene = gene, exon = exon_unique, three_prime_utr = three_utr, 
             five_prime_utr = five_utr)
  
  ## intron annotation
  txdb <- makeTxDbFromGRanges(gtf)
  introns <- unlist(intronsByTranscript(txdb))
  ## remove the intronic parts that overlap with exons from other transcripts
  anno[["intron"]] <- GenomicRanges::setdiff(introns, c(anno[["exon"]], 
                                         anno[["three_prime_utr"]], 
                                         anno[["five_prime_utr"]]))
  anno
}


#' Peak location barplot
#' 
#' Barplot of the location of peaks in a gene
#'
#' @param clipper GRangesList with peaks from CLIPper
#' @param omni GRangesList with peaks from omniCLIP
#' @param filepath Path to output file
#'
#' @return
#' @export
#'
#' @examples
peak_location_barplot <- function(clipper=NA, omni, filepath) {
  if (!missing(clipper)) {
    ## overlap the peaks with the annotation
    clipper_olap_len <- list()
    for (sample in names(clipper)){
      clipper_olap_len[[sample]] <- lapply(anno, function(x) 
        length(subsetByOverlaps(clipper[[sample]], x, type = "any")))
      names(clipper_olap_len[[sample]]) <- names(anno)
    } 
    names(clipper_olap_len) <- paste0("CLIPper-", names(clipper_olap_len))
  }
  
  omni_olap_len <- list()
  for (sample in names(omni)){
    omni_olap_len[[sample]] <- lapply(anno, function(x) 
      length(subsetByOverlaps(omni[[sample]], x, type = "any")))
    names(omni_olap_len[[sample]]) <- names(anno)
  } 
  names(omni_olap_len) <- paste0("omniCLIP-", names(omni_olap_len))
  
  df <- as.data.frame( t( 
    cbind( 
      if (!missing("clipper")) {sapply(clipper_olap_len, as.data.frame)},
      sapply(omni_olap_len, as.data.frame)) 
  ) )
  
  df$peaks <- rownames(df)
  df[,names(df) != "peaks"] <- apply(df[,names(df) != "peaks"], 2, as.integer)
  
  df <- df %>%gather(key = "annotation", value = "peak_number", -peaks)
  df <- df[df$annotation != "gene",]
  df$percentage <- df$peak_number / sapply(df$peaks, function(x) 
    sum(df[df$peaks == x, "peak_number"]) ) * 100
  
  ## reorder the annotation factor levels
  df$annotation <- factor(df$annotation, 
                          levels = c("exon", "five_prime_utr", 
                                     "three_prime_utr", "intron"))
  
  print(df)
  
  ## stacked barplot with the percentage of reads in the different gene regions
  print(ggplot(df, aes(x = peaks, y = percentage, fill = annotation)) +
    geom_col(position = "stack") +
    theme_bw() +
    theme(text=element_text(size=25), axis.text.y = element_text(angle = 45, 
                                                                 hjust = 1) ) +
    coord_flip() +
    theme(legend.position="bottom", legend.direction="vertical") 
  )
  # ggsave(filepath, p, width=7, height =7) 
}



#' Annotate peaks
#'
#' @param peaks GRanges object with peaks
#' @param genes GRanges object with GTF gene annotations
#'
#' @return GRanges object with new metadata columns 'gene_id', 'gene_name' and
#'   'gene_biotype'
#' @export
#'
#' @examples
add_gene_annotation <- function(peaks, genes){
  m <- match(peaks$gene_id, genes$gene_id)
  res <- peaks[, "score"][!is.na(m)]
  mcols(res) <- cbind(mcols(res), mcols(genes[m[!is.na(m)], 
                                              c("gene_id", "gene_name", 
                                                "gene_biotype")]))
  res
}



#' Gene region barplot
#' 
#' Barplot of the number of peaks at the different gene regions.
#'
#' @param peaks GRangesList with peaks 
#' @param anno GRangesList with gene annotations
#'
#' @return
#' @export
#'
#' @examples
peak_gene_region_barplot <- function(peaks, anno) {
  olap_len <- list()
  for (sample in names(peaks)){
    olap_len[[sample]] <- lapply(anno, function(x) 
      length(suppressWarnings(subsetByOverlaps(peaks[[sample]], x, type = "any"))))
    names(olap_len[[sample]]) <- names(anno)
  } 

  df <- as.data.frame( t( 
    cbind( 
      sapply(olap_len, as.data.frame)) 
  ) )
  
  df$peaks <- rownames(df)
  df[,names(df) != "peaks"] <- apply(df[,names(df) != "peaks"], 2, as.integer)
  
  df <- df %>% gather(key = "annotation", value = "peak_number", -peaks)
  df <- df[df$annotation != "gene",]
  df$percentage <- df$peak_number / sapply(df$peaks, function(x) 
    sum(df[df$peaks == x, "peak_number"]) ) * 100
  
  ## reorder the annotation factor levels
  df$annotation <- factor(df$annotation, 
                          levels = c("exon", "five_prime_utr", 
                                     "three_prime_utr", "intron"))
  print(df)
  ## stacked barplot with the percentage of reads in the different gene regions
  print(ggplot(df, aes(x = peaks, y = percentage, fill = annotation)) +
          geom_col(position = "stack") +
          theme_bw() +
          theme(text = element_text(size = 20), legend.position = "bottom", 
                legend.direction="horizontal", 
                legend.box = "horizontal")  +
          scale_x_discrete(breaks=c("unchanged", "down", "up"), 
                         labels=c("unchanged", "down", "up"),
                         limits = c("unchanged", "down", "up"),
                         name = "peaks in genes with RNA stability change") +
          guides(fill = guide_legend(nrow=2, byrow = TRUE))
  )
}
