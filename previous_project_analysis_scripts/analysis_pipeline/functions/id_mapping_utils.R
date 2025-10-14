# =============================================================================
# ID Mapping Utilities - GTF-Based Hash Maps for Gene/Transcript Conversion
# =============================================================================
# Author: Manveer Chauhan
# Description: Fast GTF-based ID mapping for ENSG→Symbol and ENST→GeneName_TxID
#              Uses hash maps (named vectors) for O(1) lookup speed
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

#' Create ID mappings from FLAMES GTF using regex-based parsing (no rtracklayer dependency)
#'
#' @param flames_dir Path to FLAMES output directory
#' @param reference_gtf Optional path to reference GTF (not used in current implementation)
#' @param cache_file Optional path to save/load cached mappings (RDS file)
#' @return List containing gene_map, transcript_map, gene_hash, tx_hash, and metadata
#' @export
create_id_mappings_from_flames <- function(flames_dir, reference_gtf = NULL, cache_file = NULL) {

  # Check if cached version exists
  if (!is.null(cache_file) && file.exists(cache_file)) {
    log_message(sprintf("Loading cached ID mappings from %s", basename(cache_file)), "INFO")
    return(readRDS(cache_file))
  }

  # === Step 1: Find FLAMES GTF ===
  flames_gtf_path <- file.path(flames_dir, "isoform_annotated.gtf")
  if (!file.exists(flames_gtf_path)) {
    # Try unfiltered if filtered doesn't exist
    flames_gtf_path <- file.path(flames_dir, "isoform_annotated_unfiltered.gtf")
    if (!file.exists(flames_gtf_path)) {
      stop("FLAMES GTF file not found in directory: ", flames_dir)
    }
  }

  log_message(sprintf("Reading FLAMES GTF: %s", basename(flames_gtf_path)), "INFO")

  # === Step 2: Read GTF line by line (regex-based parsing) ===
  gtf_lines <- readLines(flames_gtf_path)

  log_message(sprintf("Parsing %d lines from GTF", length(gtf_lines)), "DEBUG")

  # Initialize lists for hash maps
  gene_id_to_name <- list()
  transcript_id_to_gene <- list()
  transcript_id_to_name <- list()

  # Regex patterns (adapted from feature_ID_converter.py)
  gene_id_pattern <- 'gene_id\\s+"([^"]+)"'
  gene_name_pattern <- 'gene_name\\s+"([^"]+)"'
  transcript_id_pattern <- 'transcript_id\\s+"([^"]+)"'

  # === Step 3: Parse line by line ===
  for (line in gtf_lines) {
    # Skip comment lines
    if (startsWith(line, "#")) next

    # Extract gene_id, gene_name, transcript_id using regex
    gene_id_match <- str_match(line, gene_id_pattern)
    gene_name_match <- str_match(line, gene_name_pattern)
    transcript_id_match <- str_match(line, transcript_id_pattern)

    # Build gene_id → gene_name hash map
    if (!is.na(gene_id_match[2]) && !is.na(gene_name_match[2])) {
      gene_id_clean <- gsub("\\..*", "", gene_id_match[2])  # Remove version
      gene_name <- gene_name_match[2]

      # Only store if gene_name is not empty
      if (gene_name != "") {
        gene_id_to_name[[gene_id_clean]] <- gene_name
      }
    }

    # Build transcript_id → gene_name hash map
    if (!is.na(transcript_id_match[2])) {
      transcript_id_clean <- gsub("\\..*", "", transcript_id_match[2])

      # Determine gene symbol for this transcript
      gene_symbol <- NA

      # Handle novel isoforms (contains "-")
      # Example: ENSG00000104435-79611117-79665011-1 → parent gene is ENSG00000104435
      if (grepl("-", transcript_id_clean)) {
        # Extract parent gene ID
        parent_gene <- strsplit(transcript_id_clean, "-", fixed = TRUE)[[1]][1]

        # Look up parent gene symbol from our hash map
        if (!is.null(gene_id_to_name[[parent_gene]])) {
          gene_symbol <- gene_id_to_name[[parent_gene]]
        } else {
          gene_symbol <- parent_gene  # Fallback to ENSG
        }
      } else {
        # Not a novel isoform - use gene_name from this line
        if (!is.na(gene_name_match[2]) && gene_name_match[2] != "") {
          gene_symbol <- gene_name_match[2]
        } else if (!is.na(gene_id_match[2])) {
          # Fallback: use gene_id if gene_name not available
          gene_id_clean <- gsub("\\..*", "", gene_id_match[2])
          gene_symbol <- gene_id_clean
        }
      }

      # Create composite ID: GeneName_TranscriptID
      if (!is.na(gene_symbol)) {
        composite_id <- paste0(gene_symbol, "_", transcript_id_clean)
        transcript_id_to_name[[transcript_id_clean]] <- composite_id
        transcript_id_to_gene[[transcript_id_clean]] <- gene_symbol
      }
    }
  }

  log_message(sprintf("Parsed %d genes, %d transcripts from GTF",
                      length(gene_id_to_name), length(transcript_id_to_name)), "DEBUG")

  # === Step 4: Convert lists to data frames ===
  gene_map <- data.frame(
    ensembl_gene_id = names(gene_id_to_name),
    gene_symbol = unlist(gene_id_to_name, use.names = FALSE),
    stringsAsFactors = FALSE
  )

  transcript_map <- data.frame(
    ensembl_transcript_id = names(transcript_id_to_name),
    gene_symbol = unlist(transcript_id_to_gene, use.names = FALSE),
    gene_symbol_tx_id = unlist(transcript_id_to_name, use.names = FALSE),
    stringsAsFactors = FALSE
  )

  # === Step 5: Create hash maps (named vectors for O(1) lookup) ===
  log_message("Creating hash maps for fast lookup", "DEBUG")

  gene_hash <- setNames(gene_map$gene_symbol, gene_map$ensembl_gene_id)
  tx_hash <- setNames(transcript_map$gene_symbol_tx_id, transcript_map$ensembl_transcript_id)

  # === Step 6: Package results ===
  result <- list(
    gene_map = gene_map,              # Data frame: ENSG → Symbol
    transcript_map = transcript_map,   # Data frame: ENST → GeneName_TxID
    gene_hash = gene_hash,             # Named vector for fast lookup
    tx_hash = tx_hash,                 # Named vector for fast lookup
    metadata = list(
      n_genes = nrow(gene_map),
      n_transcripts = nrow(transcript_map),
      flames_gtf = flames_gtf_path,
      reference_gtf = reference_gtf,
      created = Sys.time()
    )
  )

  # === Step 7: Cache results if requested ===
  if (!is.null(cache_file)) {
    # Create directory if it doesn't exist
    cache_dir <- dirname(cache_file)
    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE)
    }

    log_message(sprintf("Caching ID mappings to %s", basename(cache_file)), "DEBUG")
    saveRDS(result, cache_file)
  }

  log_message(sprintf("Created ID mappings: %d genes, %d transcripts",
                      result$metadata$n_genes, result$metadata$n_transcripts), "INFO")

  return(result)
}

#' Apply gene symbol conversion to count matrix
#'
#' @param count_matrix Count matrix with ENSG IDs as rownames
#' @param gene_hash Named vector (ENSG → Gene Symbol)
#' @param handle_duplicates How to handle duplicate gene symbols ("sum", "keep_first", "keep_highest")
#' @return Count matrix with gene symbols as rownames
#' @export
apply_gene_symbol_conversion <- function(count_matrix, gene_hash, handle_duplicates = "sum") {

  ensg_ids <- rownames(count_matrix)

  # Map ENSG to gene symbols
  gene_symbols <- gene_hash[ensg_ids]

  # Keep original ENSG for unmapped features
  gene_symbols[is.na(gene_symbols)] <- ensg_ids[is.na(gene_symbols)]

  # Add symbols as rownames
  rownames(count_matrix) <- gene_symbols

  # Handle duplicate gene symbols
  if (any(duplicated(gene_symbols))) {
    n_dup <- sum(duplicated(gene_symbols))
    log_message(sprintf("Found %d duplicate gene symbols, handling with: %s", n_dup, handle_duplicates), "DEBUG")

    if (handle_duplicates == "sum") {
      # Sum counts for duplicate symbols
      count_matrix <- rowsum(count_matrix, gene_symbols)
    } else if (handle_duplicates == "keep_first") {
      # Keep first occurrence
      count_matrix <- count_matrix[!duplicated(gene_symbols), ]
    } else if (handle_duplicates == "keep_highest") {
      # Keep the feature with highest total counts
      total_counts <- Matrix::rowSums(count_matrix)
      count_df <- data.frame(
        gene_symbol = gene_symbols,
        total_counts = total_counts,
        row_idx = 1:length(gene_symbols)
      )
      keep_idx <- count_df %>%
        group_by(gene_symbol) %>%
        slice_max(total_counts, n = 1, with_ties = FALSE) %>%
        pull(row_idx)
      count_matrix <- count_matrix[keep_idx, ]
    }
  }

  return(count_matrix)
}

#' Apply transcript label conversion to count matrix
#'
#' @param count_matrix Count matrix with ENST IDs as rownames
#' @param tx_hash Named vector (ENST → GeneName_TxID)
#' @return Count matrix with GeneName_TxID as rownames
#' @export
apply_transcript_label_conversion <- function(count_matrix, tx_hash) {

  enst_ids <- rownames(count_matrix)

  # Map ENST to GeneName_TxID
  tx_labels <- tx_hash[enst_ids]

  # Keep original ENST for unmapped features
  tx_labels[is.na(tx_labels)] <- enst_ids[is.na(tx_labels)]

  # Add labels as rownames
  rownames(count_matrix) <- tx_labels

  # No deduplication needed for transcripts (each ENST is unique)

  return(count_matrix)
}

log_message("ID mapping utilities loaded successfully", "INFO")
