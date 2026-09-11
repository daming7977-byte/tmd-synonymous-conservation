from pathlib import Path
import os
import csv
import gzip
import json
import hashlib
from collections import defaultdict, Counter

ROOT = Path(os.environ["P5_DATA_PARENT"]) / "paper5_evolution_spatial_bridge"
P3 = Path(os.environ["P5_DATA_PARENT"]) / "paper3_GSE300977_raw"

REQ = Path("P5_STAGE1A_FULL_ALIGNMENT_REQUEST_AUDIT.tsv")
TARGET = Path("P5_STAGE1A_FULL_TARGET_PROTEIN_UNIVERSE.tsv")

TMD = (
    ROOT / "03_TMD_MAPPING_FEASIBILITY" /
    "P3_TO_P5_PRIMARY_PLUS75_HUMAN_TMD_COORDINATES_FROZEN.tsv"
)

QC = (
    P3 / "tmd_anchor_preregister/bridge" /
    "TMD_BRIDGE_TRANSCRIPT_STRUCTURAL_QC_v1p1.tsv"
)

CACHE = Path("ENSEMBL_FULL_ALIGNMENT_CACHE")

for p in (REQ, TARGET, TMD, QC, CACHE):
    if not p.exists():
        raise FileNotFoundError(p)

# ============================================================
# BLINDING
# NO P3 spatial outcome
# NO P4 evolutionary result
#
# Frozen Stage0 projection rules:
#   target core coverage >= 0.80
#   N-boundary support = >=1 target nongap in first 3 TMD residues
#   C-boundary support = >=1 target nongap in last 3 TMD residues
#   projected length ratio = 0.65–1.50
# ============================================================

def norm(x):
    return (
        str(x)
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
        .upper()
    )

def ungap(x):
    return x.replace("-", "").replace(".", "")

def seq_sha(x):
    return hashlib.sha256(
        ungap(x).encode("ascii")
    ).hexdigest()

def load_gz_json(path):
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)

def source_position_to_column(aln):
    out = {}
    pos = 0
    for col, aa in enumerate(aln):
        if aa not in "-.":
            out[pos] = col
            pos += 1
    return out

def target_column_to_position(aln):
    out = {}
    pos = 0
    for col, aa in enumerate(aln):
        if aa not in "-.":
            out[col] = pos
            pos += 1
    return out

def write_tsv(path, rows):
    if not rows:
        raise RuntimeError(f"no rows for {path}")

    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=list(rows[0].keys()),
            delimiter="\t",
            lineterminator="\n"
        )
        w.writeheader()
        w.writerows(rows)

# ============================================================
# 1. Exact-source formal human genes
# ============================================================

with open(REQ, encoding="utf-8") as fh:
    req_rows = list(
        csv.DictReader(fh, delimiter="\t")
    )

genes = sorted({
    r["human_ensembl_gene"]
    for r in req_rows
    if r["source_exact_frozen_P3"] == "1"
})

if len(genes) != 1921:
    raise RuntimeError(
        f"expected 1921 exact-source genes, got {len(genes)}"
    )

gene_set = set(genes)

# ============================================================
# 2. Exact frozen target gene/protein relationships
# ============================================================

expected = defaultdict(dict)

with open(TARGET, encoding="utf-8") as fh:
    target_rows = list(
        csv.DictReader(fh, delimiter="\t")
    )

if len(target_rows) != 21397:
    raise RuntimeError(
        f"expected 21397 target relationships, got {len(target_rows)}"
    )

for r in target_rows:
    gene = r["human_ensembl_gene"]
    sp = r["target_species"]

    if gene not in gene_set:
        raise RuntimeError(
            f"non-formal gene in target universe: {gene}"
        )

    if sp in expected[gene]:
        raise RuntimeError(
            f"duplicate target relationship: {gene} {sp}"
        )

    expected[gene][sp] = {
        "gene": r["target_ensembl_gene"],
        "protein": r["target_protein_id"],
    }

if sum(len(x) for x in expected.values()) != 21397:
    raise RuntimeError("formal relationship count mismatch")

# ============================================================
# 3. Frozen human TMDs
# ============================================================

tmds = defaultdict(list)
gene_tx = defaultdict(set)

with open(TMD, encoding="utf-8") as fh:
    r = csv.DictReader(fh, delimiter="\t")

    for row in r:
        gene = row["human_ensembl_gene"]

        if gene not in gene_set:
            continue

        tx = row["human_ensembl_transcript"]
        gene_tx[gene].add(tx)

        tmds[gene].append({
            "tmd_index": int(row["tmd_index"]),
            "is_internal": int(row["is_internal"]),
            "start0": int(row["human_tmd_aa_start0"]),
            "end0": int(row["human_tmd_aa_end0_exclusive"]),
            "length": int(row["human_tmd_length_aa"]),
        })

for gene in genes:
    if len(gene_tx[gene]) != 1:
        raise RuntimeError(
            f"non-unique transcript for {gene}: {gene_tx[gene]}"
        )

    tmds[gene].sort(
        key=lambda x: x["tmd_index"]
    )

n_tmd = sum(len(tmds[g]) for g in genes)
n_internal = sum(
    t["is_internal"]
    for g in genes
    for t in tmds[g]
)

if n_tmd != 9827:
    raise RuntimeError(
        f"expected 9827 formal TMDs, got {n_tmd}"
    )

if n_internal != 8420:
    raise RuntimeError(
        f"expected 8420 internal TMDs, got {n_internal}"
    )

# ============================================================
# 4. Frozen P3 protein SHA
# ============================================================

qc = {}

with open(QC, encoding="utf-8") as fh:
    r = csv.DictReader(fh, delimiter="\t")

    for row in r:
        qc[row["target_transcript"]] = (
            row["target_protein_sha256"]
        )

# ============================================================
# 5. Full TMD projection
# ============================================================

detail = []
gene_audit = []

for gi, gene in enumerate(genes, start=1):

    cache = CACHE / f"{gene}.protein_alignment.json.gz"

    if not cache.exists():
        raise RuntimeError(
            f"missing alignment cache: {gene}"
        )

    data = load_gz_json(cache)

    entries = data.get("data", [])

    if len(entries) != 1:
        raise RuntimeError(
            f"unexpected alignment schema for {gene}: {len(entries)}"
        )

    homologies = entries[0].get("homologies", [])

    by_species = defaultdict(list)

    for h in homologies:

        if h.get("type") != "ortholog_one2one":
            continue

        sp = h.get("target", {}).get("species", "")

        if sp in expected[gene]:
            by_species[sp].append(h)

    tx = next(iter(gene_tx[gene]))

    if tx not in qc:
        raise RuntimeError(
            f"P3 frozen protein SHA missing: {tx}"
        )

    frozen_sha = qc[tx]

    n_rel = 0
    n_gene_match = 0
    n_protein_match = 0
    n_source_exact = 0

    for sp in sorted(expected[gene]):

        hits = by_species.get(sp, [])

        if len(hits) != 1:
            raise RuntimeError(
                f"expected unique one2one missing/duplicated: "
                f"{gene} {sp} n={len(hits)}"
            )

        n_rel += 1

        h = hits[0]
        src = h.get("source", {})
        tgt = h.get("target", {})

        target_gene = str(tgt.get("id", ""))
        target_protein = str(tgt.get("protein_id", ""))

        gene_match = (
            target_gene == expected[gene][sp]["gene"]
        )

        protein_match = (
            target_protein == expected[gene][sp]["protein"]
        )

        if gene_match:
            n_gene_match += 1

        if protein_match:
            n_protein_match += 1

        src_aln = norm(src.get("align_seq", ""))
        tgt_aln = norm(tgt.get("align_seq", ""))

        src_sha = (
            seq_sha(src_aln)
            if src_aln else ""
        )

        source_exact = (
            bool(src_aln)
            and src_sha == frozen_sha
        )

        if source_exact:
            n_source_exact += 1

        same_len = (
            bool(src_aln)
            and bool(tgt_aln)
            and len(src_aln) == len(tgt_aln)
        )

        if same_len:
            source_map = source_position_to_column(src_aln)
            target_map = target_column_to_position(tgt_aln)
        else:
            source_map = {}
            target_map = {}

        for t in tmds[gene]:

            start0 = t["start0"]
            end0 = t["end0"]
            L = t["length"]

            failure = []

            if not gene_match:
                failure.append("TARGET_GENE_MISMATCH")

            if not protein_match:
                failure.append("TARGET_PROTEIN_MISMATCH")

            if not source_exact:
                failure.append("COMPARA_SOURCE_NOT_EXACT_P3")

            if not same_len:
                failure.append(
                    "ALIGNMENT_LENGTH_MISMATCH_OR_MISSING"
                )

            positions = list(range(start0, end0))

            complete = (
                same_len
                and all(p in source_map for p in positions)
            )

            if not complete:

                failure.append(
                    "HUMAN_TMD_PROJECTION_INCOMPLETE"
                )

                core_cov = None
                n_support = False
                c_support = False
                target_len = None
                ratio = None
                target_start = None
                target_end = None

            else:

                cols = [
                    source_map[p]
                    for p in positions
                ]

                core_nongap = sum(
                    tgt_aln[c] not in "-."
                    for c in cols
                )

                core_cov = core_nongap / L

                first_cols = cols[:min(3, len(cols))]
                last_cols = cols[-min(3, len(cols)):]

                n_support = any(
                    tgt_aln[c] not in "-."
                    for c in first_cols
                )

                c_support = any(
                    tgt_aln[c] not in "-."
                    for c in last_cols
                )

                span_cols = range(
                    min(cols),
                    max(cols) + 1
                )

                target_positions = [
                    target_map[c]
                    for c in span_cols
                    if c in target_map
                ]

                target_len = len(target_positions)

                ratio = (
                    target_len / L
                    if L else None
                )

                if target_positions:
                    target_start = min(target_positions) + 1
                    target_end = max(target_positions) + 1
                else:
                    target_start = None
                    target_end = None

                if core_cov < 0.80:
                    failure.append(
                        "TARGET_CORE_COVERAGE_LT_0.80"
                    )

                if not n_support:
                    failure.append(
                        "N_BOUNDARY_UNSUPPORTED"
                    )

                if not c_support:
                    failure.append(
                        "C_BOUNDARY_UNSUPPORTED"
                    )

                if (
                    ratio is None
                    or ratio < 0.65
                    or ratio > 1.50
                ):
                    failure.append(
                        "PROJECTED_LENGTH_RATIO_OUTSIDE_0.65_1.50"
                    )

            mapping_pass = (
                len(failure) == 0
            )

            detail.append({
                "human_ensembl_gene": gene,
                "human_ensembl_transcript": tx,
                "target_species": sp,
                "expected_target_gene":
                    expected[gene][sp]["gene"],
                "returned_target_gene":
                    target_gene,
                "target_protein_id":
                    target_protein,
                "tmd_index":
                    t["tmd_index"],
                "is_internal":
                    t["is_internal"],
                "human_tmd_length_aa":
                    L,
                "source_exact_p3_hash":
                    int(source_exact),
                "target_gene_matches_frozen":
                    int(gene_match),
                "target_protein_matches_formal":
                    int(protein_match),
                "alignment_length_match":
                    int(same_len),
                "human_tmd_projection_complete":
                    int(complete),
                "target_core_coverage":
                    "" if core_cov is None else core_cov,
                "n_boundary_support":
                    int(n_support),
                "c_boundary_support":
                    int(c_support),
                "target_projected_length_aa":
                    "" if target_len is None else target_len,
                "target_projected_length_ratio":
                    "" if ratio is None else ratio,
                "target_projected_start_1based":
                    "" if target_start is None else target_start,
                "target_projected_end_1based":
                    "" if target_end is None else target_end,
                "mapping_pass":
                    int(mapping_pass),
                "failure_reason":
                    "PASS"
                    if mapping_pass
                    else ";".join(failure),
            })

    gene_audit.append({
        "human_ensembl_gene":
            gene,
        "human_ensembl_transcript":
            tx,
        "expected_one2one_species":
            len(expected[gene]),
        "returned_one2one_species":
            n_rel,
        "target_gene_match_species":
            n_gene_match,
        "target_protein_match_species":
            n_protein_match,
        "source_exact_frozen_P3_species":
            n_source_exact,
        "gene_integrity_pass":
            int(
                n_rel == len(expected[gene])
                and n_gene_match == len(expected[gene])
                and n_protein_match == len(expected[gene])
                and n_source_exact == len(expected[gene])
            ),
    })

    if gi % 100 == 0 or gi == len(genes):
        print(
            f"PROJECTED_GENES {gi}/{len(genes)}",
            flush=True
        )

# ============================================================
# 6. Relationship integrity
# ============================================================

if sum(x["gene_integrity_pass"] for x in gene_audit) != 1921:
    raise RuntimeError(
        "full relationship/source integrity did not remain exact"
    )

# ============================================================
# 7. TMD summaries
# ============================================================

by_tmd = defaultdict(list)

for r in detail:
    key = (
        r["human_ensembl_gene"],
        r["human_ensembl_transcript"],
        int(r["tmd_index"]),
        int(r["is_internal"]),
    )
    by_tmd[key].append(r)

summary = []

for key, rows in sorted(by_tmd.items()):

    gene, tx, idx, internal = key

    mapped = sum(
        int(r["mapping_pass"]) == 1
        for r in rows
    )

    summary.append({
        "human_ensembl_gene":
            gene,
        "human_ensembl_transcript":
            tx,
        "tmd_index":
            idx,
        "is_internal":
            internal,
        "expected_one2one_species":
            len(rows),
        "mapping_pass_species":
            mapped,
        "mapping_ge6_species":
            int(mapped >= 6),
    })

if len(summary) != 9827:
    raise RuntimeError(
        f"expected 9827 TMD summaries, got {len(summary)}"
    )

internal = [
    x for x in summary
    if x["is_internal"] == 1
]

if len(internal) != 8420:
    raise RuntimeError(
        f"expected 8420 internal TMDs, got {len(internal)}"
    )

internal_ge6 = [
    x for x in internal
    if x["mapping_ge6_species"] == 1
]

# ============================================================
# 8. Failure census
# ============================================================

fail = Counter()

for r in detail:
    if r["failure_reason"] == "PASS":
        continue

    for x in r["failure_reason"].split(";"):
        fail[x] += 1

pass_pairs = sum(
    r["mapping_pass"] == 1
    for r in detail
)

# ============================================================
# 9. Write
# ============================================================

detail.sort(
    key=lambda x: (
        x["human_ensembl_gene"],
        int(x["tmd_index"]),
        x["target_species"],
    )
)

gene_audit.sort(
    key=lambda x:
        x["human_ensembl_gene"]
)

summary.sort(
    key=lambda x: (
        x["human_ensembl_gene"],
        int(x["tmd_index"]),
    )
)

write_tsv(
    "P5_STAGE1A_FULL_TMD_PROJECTION_DETAIL.tsv",
    detail
)

write_tsv(
    "P5_STAGE1A_FULL_TMD_PROJECTION_REQUEST_AUDIT.tsv",
    gene_audit
)

write_tsv(
    "P5_STAGE1A_FULL_TMD_PROJECTION_SUMMARY.tsv",
    summary
)

fraction = (
    len(internal_ge6)
    / len(internal)
)

lines = [
    "P5 STAGE1A FULL TMD PROJECTION",
    "=" * 72,
    "",
    f"EXACT_SOURCE_GENES = {len(genes)}",
    "FORMAL_ONE2ONE_RELATIONSHIPS = 21397",
    "",
    f"FORMAL_PRIMARY_PLUS75_TMD = {len(summary)}",
    f"FORMAL_INTERNAL_TMD = {len(internal)}",
    "",
    f"TMD_SPECIES_PAIRS = {len(detail)}",
    f"PROJECTION_PASS_PAIRS = {pass_pairs}",
    f"PROJECTION_FAIL_PAIRS = {len(detail)-pass_pairs}",
    "",
    f"INTERNAL_TMD_WITH_GE6_MAPPED_SPECIES = {len(internal_ge6)}",
    f"FRACTION_INTERNAL_TMD_WITH_GE6_MAPPED_SPECIES = {fraction:.6f}",
    "",
    "HUMAN_SOURCE_EXACT_FROZEN_P3_GENES = 1921",
    "TARGET_GENE_RELATIONSHIP_MATCHES = 21397",
    "TARGET_PROTEIN_RELATIONSHIP_MATCHES = 21397",
    "",
    "PROJECTION_FAILURE_REASON_COUNTS:"
]

if fail:
    for k, v in sorted(fail.items()):
        lines.append(
            f"  {k}\t{v}"
        )
else:
    lines.append("  NONE")

lines += [
    "",
    "FULL_PROJECTION_INTEGRITY = PASS",
    "",
    "PROJECTION_RULES_CHANGED_FROM_STAGE0_PILOT = FALSE",
    "P3_SPATIAL_OUTCOME_LOADED = FALSE",
    "P4_EVOLUTIONARY_RESULT_LOADED = FALSE",
    "PRIMARY_EVOLUTIONARY_SCORE_DEFINITION_CHANGED = FALSE",
    "PRIMARY_WINDOW_CHANGED = FALSE",
]

Path(
    "P5_STAGE1A_FULL_TMD_PROJECTION_CENSUS.txt"
).write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8"
)

print()
print("\n".join(lines))
