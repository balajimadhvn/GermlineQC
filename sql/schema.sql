-- GermlineQC variant store
PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS variants;
DROP TABLE IF EXISTS qc_metrics;
DROP TABLE IF EXISTS samples;

CREATE TABLE samples (
    sample_id   TEXT PRIMARY KEY,
    qc_status   TEXT CHECK (qc_status IN ('PASS', 'FAIL')),
    loaded_at   TEXT DEFAULT (datetime('now'))
);

CREATE TABLE qc_metrics (
    sample_id   TEXT NOT NULL REFERENCES samples(sample_id),
    metric      TEXT NOT NULL,
    value       REAL,
    min_allowed REAL,
    max_allowed REAL,
    status      TEXT CHECK (status IN ('PASS', 'WARN', 'FAIL', 'INFO')),
    PRIMARY KEY (sample_id, metric)
);

CREATE TABLE variants (
    variant_id      INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id       TEXT NOT NULL REFERENCES samples(sample_id),
    chrom           TEXT NOT NULL,
    pos             INTEGER NOT NULL,
    ref             TEXT NOT NULL,
    alt             TEXT NOT NULL,
    var_type        TEXT CHECK (var_type IN ('SNV', 'INDEL')),
    filter          TEXT,
    qual            REAL,
    genotype        TEXT,
    depth           INTEGER,
    gq              INTEGER,
    vaf             REAL,
    gene            TEXT,
    effect          TEXT,
    impact          TEXT,
    transcript      TEXT,
    hgvs_c          TEXT,
    hgvs_p          TEXT,
    clinvar_id      TEXT,
    clinvar_sig     TEXT,
    clinvar_review  TEXT,
    clinvar_disease TEXT,
    UNIQUE (sample_id, chrom, pos, ref, alt)
);

CREATE INDEX idx_variants_locus  ON variants(chrom, pos);
CREATE INDEX idx_variants_gene   ON variants(gene);
CREATE INDEX idx_variants_impact ON variants(impact);

-- convenience view: PASS variants with a reportable consequence
CREATE VIEW reportable AS
SELECT *
FROM variants
WHERE filter = 'PASS'
  AND genotype IN ('het', 'hom_alt')
  AND (impact IN ('HIGH', 'MODERATE') OR clinvar_sig LIKE '%athogenic%');
