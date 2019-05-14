biomartGeneDownload <- function(mart) {
  df <- getBM(attributes = c(
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "gene_biotype",
    "ensembl_gene_id",
    "external_gene_name",
    "description"
  ), mart=mart)
  df <- dplyr::rename(df,
    chr = chromosome_name,
    start = start_position,
    end = end_position,
    gene_id = ensembl_gene_id,
    gene_name = external_gene_name
  )
}

biomartTranscriptDownload <- function(mart) {
  df <- getBM(attributes = c(
    "chromosome_name",
    "transcript_start",
    "transcript_end",
    "strand",
    "transcription_start_site",
    "ensembl_transcript_id",
    "external_transcript_name",
    "ensembl_gene_id",
    "external_gene_name",
    "description"
  ), mart=mart)
  df <- dplyr::rename(df,
    chr = chromosome_name,
    start = transcript_start,
    end = transcript_end,
    transcript_id = ensembl_transcript_id,
    gene_id = ensembl_gene_id,
    gene_name = external_gene_name,
    transcript_name = external_transcript_name
  )
}

biomartGODownload <- function(mart, gene_ids) {
  df <- getBM(
    attributes = c("ensembl_gene_id", "go_id"),
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = mart
  )
  df <- dplyr::rename(df,
    gene_id = ensembl_gene_id
  )
  df <- df[df$go_id != "", ]
}

biomartGODescriptions <- function(mart) {
  df <- getBM(attributes = c(
    "go_id",
    "name_1006",
    "namespace_1003"
  ), mart=mart)
  df <- dplyr::rename(df,
    go_name = name_1006,
    go_domain = namespace_1003
  )
}

