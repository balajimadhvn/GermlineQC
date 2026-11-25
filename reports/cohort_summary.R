#!/usr/bin/env Rscript
# Run-level summary straight from the SQLite variant store.
# Usage: Rscript reports/cohort_summary.R results_validation/database/variants.db reports/figures

suppressPackageStartupMessages({
  library(DBI); library(RSQLite); library(dplyr); library(ggplot2); library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
db_path <- if (length(args) >= 1) args[1] else "results_validation/database/variants.db"
out_dir <- if (length(args) >= 2) args[2] else "reports/figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

con <- dbConnect(SQLite(), db_path)
on.exit(dbDisconnect(con))

# QC metrics in wide form, one row per sample
qc <- dbGetQuery(con, "
  SELECT sample_id, metric, value
  FROM qc_metrics
  WHERE metric IN ('mean_target_depth', 'pct_target_20x', 'duplication_rate', 'ts_tv', 'q30_rate')
") %>% pivot_wider(names_from = metric, values_from = value)

# PASS calls per sample by type, het/hom ratio
calls <- dbGetQuery(con, "
  SELECT sample_id,
         SUM(var_type = 'SNV')                          AS snvs,
         SUM(var_type = 'INDEL')                        AS indels,
         ROUND(1.0 * SUM(genotype = 'het') /
               NULLIF(SUM(genotype = 'hom_alt'), 0), 2) AS het_hom_ratio
  FROM variants
  WHERE filter = 'PASS'
  GROUP BY sample_id
")

summary <- qc %>% inner_join(calls, by = "sample_id") %>% arrange(sample_id)
write.table(summary, file.path(dirname(out_dir), "cohort_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(as.data.frame(summary), row.names = FALSE)

# Consequence spectrum of PASS calls
effects <- dbGetQuery(con, "
  SELECT impact, effect, COUNT(DISTINCT chrom || pos || ref || alt) AS n
  FROM variants
  WHERE filter = 'PASS' AND effect IS NOT NULL AND effect != ''
  GROUP BY impact, effect
  ORDER BY n DESC
  LIMIT 12
")

p1 <- ggplot(effects, aes(x = n, y = reorder(effect, n), fill = impact)) +
  geom_col() +
  scale_fill_manual(values = c(HIGH = "#b2182b", MODERATE = "#ef8a62", LOW = "#67a9cf", MODIFIER = "#bdbdbd")) +
  labs(x = "Distinct PASS variants (all samples)", y = NULL, fill = "SnpEff impact",
       title = "Most frequent predicted consequences") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "consequences.png"), p1, width = 8, height = 4.5, dpi = 120, bg = "white")

# Depth vs number of SNVs called (titration samples show the saturation curve)
p2 <- ggplot(summary, aes(x = mean_target_depth, y = snvs)) +
  geom_point(size = 2.5, colour = "#2c7fb8") +
  scale_x_continuous(trans = "log2") +
  labs(x = "Mean target depth (log2)", y = "PASS SNVs", title = "PASS SNVs vs depth") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "depth_vs_calls.png"), p2, width = 6, height = 4, dpi = 120, bg = "white")
cat("figures written to", out_dir, "\n")
