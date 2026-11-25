-- 1. Sample QC overview
SELECT s.sample_id,
       s.qc_status,
       MAX(CASE WHEN q.metric = 'mean_target_depth'  THEN q.value END) AS mean_depth,
       MAX(CASE WHEN q.metric = 'pct_target_20x'     THEN q.value END) AS pct_20x,
       MAX(CASE WHEN q.metric = 'duplication_rate'   THEN q.value END) AS dup_rate,
       MAX(CASE WHEN q.metric = 'ts_tv'              THEN q.value END) AS ts_tv
FROM samples s
JOIN qc_metrics q USING (sample_id)
GROUP BY s.sample_id
ORDER BY s.sample_id;

-- 2. Any QC metric outside its acceptance range
SELECT sample_id, metric, value, min_allowed, max_allowed, status
FROM qc_metrics
WHERE status IN ('FAIL', 'WARN');

-- 3. PASS variant counts by type and genotype per sample
SELECT sample_id, var_type, genotype, COUNT(*) AS n
FROM variants
WHERE filter = 'PASS'
GROUP BY sample_id, var_type, genotype
ORDER BY sample_id, var_type, genotype;

-- 4. Coding / splice variants per gene (PASS only)
SELECT gene, impact, COUNT(DISTINCT chrom || ':' || pos || ref || '>' || alt) AS distinct_variants,
       COUNT(DISTINCT sample_id) AS samples
FROM variants
WHERE filter = 'PASS' AND impact IN ('HIGH', 'MODERATE', 'LOW')
GROUP BY gene, impact
ORDER BY distinct_variants DESC
LIMIT 15;

-- 5. Variants with a ClinVar assertion
SELECT DISTINCT chrom, pos, ref, alt, gene, hgvs_p, clinvar_sig, clinvar_review
FROM variants
WHERE clinvar_sig IS NOT NULL AND filter = 'PASS'
ORDER BY pos
LIMIT 20;

-- 6. Low-confidence calls that still passed hard filters (candidates for manual review / orthogonal confirmation)
SELECT sample_id, chrom, pos, ref, alt, genotype, depth, gq, vaf
FROM variants
WHERE filter = 'PASS' AND (gq < 20 OR depth < 15 OR (genotype = 'het' AND (vaf < 0.25 OR vaf > 0.75)))
ORDER BY sample_id, pos
LIMIT 20;

-- 7. Variants seen in every sample vs. private to one sample
WITH per_variant AS (
    SELECT chrom, pos, ref, alt, COUNT(DISTINCT sample_id) AS n_samples
    FROM variants WHERE filter = 'PASS' GROUP BY chrom, pos, ref, alt
)
SELECT CASE WHEN n_samples = (SELECT COUNT(*) FROM samples) THEN 'all samples'
            WHEN n_samples = 1 THEN 'private' ELSE 'shared (some)' END AS category,
       COUNT(*) AS variants
FROM per_variant
GROUP BY category;

-- 8. ACMG classes per sample (only when the ACMG module ran)
SELECT sample_id, acmg_class, COUNT(*) AS variants
FROM acmg
GROUP BY sample_id, acmg_class
ORDER BY sample_id, CASE acmg_class WHEN 'Pathogenic' THEN 1 WHEN 'Likely pathogenic' THEN 2
         WHEN 'Uncertain significance' THEN 3 WHEN 'Likely benign' THEN 4 ELSE 5 END;

-- 9. Variants needing review: P/LP or VUS with coding impact, joined to call quality
SELECT a.sample_id, a.gene, a.hgvs_c, a.hgvs_p, a.criteria, a.acmg_class, v.depth, v.gq, v.vaf
FROM acmg a
JOIN variants v USING (sample_id, chrom, pos, ref, alt)
WHERE a.acmg_class IN ('Pathogenic', 'Likely pathogenic')
   OR (a.acmg_class = 'Uncertain significance' AND a.consequence NOT LIKE '%intron%'
       AND a.consequence NOT LIKE '%UTR%' AND a.consequence != '')
ORDER BY a.acmg_class, a.gene
LIMIT 25;
