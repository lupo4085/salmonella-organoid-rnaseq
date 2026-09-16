#!/usr/bin/env Rscript

# Frozen mouse KEGG GSEA used for Figure 1C-D of the final MSci report.
# Usage: Rscript analysis/GSEA_Mouse_fixed.R /path/to/coldata_ms.csv

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(patchwork)
  library(BiocParallel)
})

analysis_seed <- 20260915L
gsea_cutoff <- 0.05
min_gene_set_size <- 10L
max_gene_set_size <- 500L
top_per_direction <- 8L
expected_counts <- c(73L, 16L, 29L, 2L)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(args) >= 1L) args[[1L]] else file.path(repo_root, "data", "coldata_ms.csv")
input_csv <- normalizePath(input_csv, mustWork = TRUE)
kegg_snapshot_file <- file.path(repo_root, "reference", "kegg_mmu_snapshot_2026-09-15.rds")
if (!file.exists(kegg_snapshot_file)) stop("Missing frozen KEGG mapping: ", kegg_snapshot_file)

results_dir <- file.path(repo_root, "results")
figures_dir <- file.path(repo_root, "figures")
environment_dir <- file.path(repo_root, "environment")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(environment_dir, recursive = TRUE, showWarnings = FALSE)

# Serial execution avoids platform-dependent Windows worker behaviour.
register(SerialParam())

kegg_snapshot <- readRDS(kegg_snapshot_file)
term2gene <- kegg_snapshot$KEGGPATHID2EXTID[, c("from", "to")]
term2name <- kegg_snapshot$KEGGPATHID2NAME[, c("from", "to")]
colnames(term2gene) <- c("term", "gene")
colnames(term2name) <- c("term", "name")

raw <- read.csv(input_csv, check.names = FALSE)
count_matrix <- raw |>
  tibble::column_to_rownames("gene_id") |>
  as.matrix()

expected_columns <- c(
  "m411_NT1", "m411_NT2", "m411_NT3",
  "m411_MBCD1", "m411_MBCD2", "m411_MBCD3",
  "m411_STm1", "m411_STm2", "m411_STm3",
  "m411_STm+MBCD1", "m411_STm+MBCD2", "m411_STm+MBCD3"
)
stopifnot(identical(colnames(count_matrix), expected_columns))

condition <- factor(
  c(rep("NT", 3), rep("MBCD", 3), rep("STm", 3), rep("STm_MBCD", 3)),
  levels = c("NT", "MBCD", "STm", "STm_MBCD")
)
coldata <- data.frame(row.names = colnames(count_matrix), condition_ms = condition)

dds <- DESeqDataSetFromMatrix(round(count_matrix), coldata, ~ condition_ms)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds, quiet = TRUE)

make_ranked_gene_list <- function(result_table) {
  ranks <- result_table$log2FoldChange
  names(ranks) <- sub("\\..*$", "", rownames(result_table))
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  id_map <- bitr(names(ranks), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  id_map <- id_map[order(match(id_map$ENSEMBL, names(ranks)), id_map$ENTREZID), ]
  mapped_ranks <- ranks[match(id_map$ENSEMBL, names(ranks))]
  names(mapped_ranks) <- id_map$ENTREZID
  mapped_ranks <- mapped_ranks[!is.na(names(mapped_ranks))]
  mapped_ranks <- mapped_ranks[!duplicated(names(mapped_ranks))]
  sort(mapped_ranks, decreasing = TRUE)
}

run_fixed_gsea <- function(numerator, denominator) {
  result_table <- as.data.frame(results(
    dds,
    contrast = c("condition_ms", numerator, denominator),
    alpha = gsea_cutoff
  ))
  gene_list <- make_ranked_gene_list(result_table)
  set.seed(analysis_seed)
  gsea <- GSEA(
    geneList = gene_list,
    exponent = 1,
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    eps = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea",
    TERM2GENE = term2gene,
    TERM2NAME = term2name
  )
  list(gsea = gsea, gene_list = gene_list)
}

prepare_results <- function(run, comparison) {
  out <- as.data.frame(run$gsea@result)
  out <- out[is.finite(out$NES) & !is.na(out$p.adjust), ]
  out$direction <- ifelse(out$NES > 0, "Positive NES", "Negative NES")
  out$core_gene_count <- vapply(strsplit(out$core_enrichment, "/", fixed = TRUE), length, integer(1))
  out$comparison <- comparison
  out
}

intra <- run_fixed_gsea("STm", "NT")
extra <- run_fixed_gsea("STm_MBCD", "MBCD")
intra_all <- prepare_results(intra, "STm vs NT")
extra_all <- prepare_results(extra, "STm_MBCD vs MBCD")
intra_sig <- intra_all[intra_all$p.adjust < gsea_cutoff, ]
extra_sig <- extra_all[extra_all$p.adjust < gsea_cutoff, ]

write.csv(intra_all, file.path(results_dir, "GSEA_STm_vs_NT_all.csv"), row.names = FALSE)
write.csv(extra_all, file.path(results_dir, "GSEA_STm_MBCD_vs_MBCD_all.csv"), row.names = FALSE)
write.csv(intra_sig, file.path(results_dir, "GSEA_STm_vs_NT_padj005.csv"), row.names = FALSE)
write.csv(extra_sig, file.path(results_dir, "GSEA_STm_MBCD_vs_MBCD_padj005.csv"), row.names = FALSE)

summary_table <- bind_rows(intra_sig, extra_sig) |>
  count(comparison, direction, name = "pathways")
write.csv(summary_table, file.path(results_dir, "GSEA_summary.csv"), row.names = FALSE)

observed_counts <- c(
  sum(intra_sig$NES > 0), sum(intra_sig$NES < 0),
  sum(extra_sig$NES > 0), sum(extra_sig$NES < 0)
)
if (!identical(as.integer(observed_counts), expected_counts)) {
  stop("Frozen-result validation failed. Expected ", paste(expected_counts, collapse = "/"),
       "; observed ", paste(observed_counts, collapse = "/"), ".")
}

select_for_plot <- function(result_table) {
  result_table |>
    group_by(direction) |>
    arrange(p.adjust, desc(abs(NES)), .by_group = TRUE) |>
    slice_head(n = top_per_direction) |>
    ungroup()
}

wrap_label <- function(x, width = 31) {
  vapply(x, function(label) paste(strwrap(label, width = width), collapse = "\n"), character(1))
}

intra_plot_data <- select_for_plot(intra_sig)
extra_plot_data <- select_for_plot(extra_sig)
common_nes_limit <- ceiling(max(abs(c(intra_plot_data$NES, extra_plot_data$NES))) * 10) / 10 + 0.1
common_p_adjust_limits <- range(c(intra_plot_data$p.adjust, extra_plot_data$p.adjust), finite = TRUE)
common_core_gene_limits <- range(c(intra_plot_data$core_gene_count, extra_plot_data$core_gene_count), finite = TRUE)

make_nes_plot <- function(plot_data, title) {
  plot_data <- plot_data |>
    arrange(NES) |>
    mutate(
      display_label = wrap_label(Description),
      display_label = factor(display_label, levels = unique(display_label))
    )
  ggplot(plot_data, aes(x = NES, y = display_label)) +
    geom_vline(xintercept = 0, linewidth = 0.55, colour = "grey35") +
    geom_segment(aes(x = 0, xend = NES, yend = display_label), linewidth = 0.75, colour = "grey72") +
    geom_point(aes(size = core_gene_count, colour = p.adjust), alpha = 0.95) +
    scale_x_continuous(
      limits = c(-common_nes_limit, common_nes_limit),
      breaks = seq(-floor(common_nes_limit), floor(common_nes_limit), by = 1),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_colour_gradient(low = "#D73027", high = "#4575B4", limits = common_p_adjust_limits) +
    scale_size_continuous(range = c(3.2, 8.5), limits = common_core_gene_limits) +
    labs(title = title, x = "Normalised enrichment score (NES)", y = NULL,
         colour = "BH-adjusted p", size = "Core genes") +
    theme_classic(base_size = 12, base_family = "Times New Roman") +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
      axis.text.y = element_text(size = 9.5, lineheight = 0.92),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11, margin = margin(t = 8)),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(8, 8, 8, 8)
    )
}

mbcd_label <- paste0("M", intToUtf8(0x03B2), "CD")
plot_c <- make_nes_plot(intra_plot_data, "STm vs NT\n(invasion-permitting)")
plot_d <- make_nes_plot(extra_plot_data,
                        paste0("STm + ", mbcd_label, " vs ", mbcd_label, "\n(invasion-restricted)"))
combined <- plot_c + plot_d + plot_layout(guides = "collect", widths = c(1, 1)) &
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "Figure1_C_bidirectional_GSEA.png"), plot_c,
       width = 7.1, height = 8.2, units = "in", dpi = 300, bg = "white")
ggsave(file.path(figures_dir, "Figure1_D_bidirectional_GSEA.png"), plot_d,
       width = 7.1, height = 8.2, units = "in", dpi = 300, bg = "white")
ggsave(file.path(figures_dir, "Figure1_CD_bidirectional_GSEA_preview.png"), combined,
       width = 14.5, height = 8.2, units = "in", dpi = 300, bg = "white")

capture.output(sessionInfo(), file = file.path(environment_dir, "sessionInfo.txt"))
writeLines(c(
  paste0("analysis_seed=", analysis_seed),
  paste0("gsea_cutoff_adjusted_p=", gsea_cutoff),
  "pAdjustMethod=BH", "eps=0",
  paste0("minGSSize=", min_gene_set_size),
  paste0("maxGSSize=", max_gene_set_size),
  "execution=serial",
  paste0("input_md5=", unname(tools::md5sum(input_csv))),
  paste0("kegg_snapshot_md5=", unname(tools::md5sum(kegg_snapshot_file)))
), file.path(environment_dir, "analysis_manifest.txt"))

print(summary_table)
message("Frozen GSEA validation passed.")
