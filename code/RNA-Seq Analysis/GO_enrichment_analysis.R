### Perform GO enrichment analysis using 'goseq' package ###

run_goseq_enrichment <- function(
  gene_list_file,
  gene_length_file,
  go_bp_file,
  go_mf_file,
  go_cc_file,
  go_term_file,
  output_dir = "GOseq_output"
)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Read input data
  gene_list <- read.table(gene_list_file, stringsAsFactors = FALSE)[,1]
  gene_length_df <- read.table(gene_length_file, sep="\t", stringsAsFactors=FALSE)
  go_bp <- read.table(go_bp_file, sep="\t", stringsAsFactors=FALSE)
  go_mf <- read.table(go_mf_file, sep="\t", stringsAsFactors=FALSE)
  go_cc <- read.table(go_cc_file, sep="\t", stringsAsFactors=FALSE)
  go_terms <- read.table(go_term_file, sep="\t", stringsAsFactors=FALSE, fill=TRUE)

  # Prepare gene vector
  assayed_genes <- gene_length_df[,1]
  gene_length <- gene_length_df[,2]
  names(gene_length) <- assayed_genes

  gene_vector <- as.integer(assayed_genes %in% gene_list)
  names(gene_vector) <- assayed_genes

  # Probability weighting function
  pwf <- nullp(gene_vector, bias.data=gene_length, plot.fit=FALSE)

  # Run enrichment
  run_single_ontology <- function(go_data, ontology_name) {

    res <- goseq(pwf,
                 gene2cat = go_data,
                 method = "Hypergeometric",
                 use_genes_without_cat = TRUE)

    res$qval <- p.adjust(res$over_represented_pvalue, method="BH")

    enriched <- subset(res, qval < 0.05)

    if (nrow(enriched) == 0) {
      message(paste("No enriched GO terms found for", ontology_name))
      return(NULL)
    }

    enriched$Ontology <- ontology_name
    return(enriched)
  }

  res_bp <- run_single_ontology(go_bp, "BP")
  res_mf <- run_single_ontology(go_mf, "MF")
  res_cc <- run_single_ontology(go_cc, "CC")

  results <- bind_rows(res_bp, res_mf, res_cc)

  if (is.null(results) || nrow(results) == 0) {
    message("No GO terms were enriched.")
    return(NULL)
  }

  # Save enrichment table
  output_table <- file.path(output_dir, "GOseq_enrichment_results.txt")
  write.table(results,
              output_table,
              sep="\t",
              quote=FALSE,
              row.names=FALSE)

  message(paste("Enrichment results saved to:", output_table))

