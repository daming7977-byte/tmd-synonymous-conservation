from pathlib import Path
import os
import csv
import gzip
import json
import math
from collections import defaultdict, Counter
from statistics import median, mean

ROOT = Path(os.environ["P5_DATA_PARENT"]) / "paper5_evolution_spatial_bridge"

TMD = (
    ROOT / "03_TMD_MAPPING_FEASIBILITY" /
    "P3_TO_P5_PRIMARY_PLUS75_HUMAN_TMD_COORDINATES_FROZEN.tsv"
)

PROJECTION = Path(
    "P5_STAGE1A_FULL_TMD_PROJECTION_DETAIL.tsv"
)

PROJECTION_SUMMARY = Path(
    "P5_STAGE1A_FULL_TMD_PROJECTION_SUMMARY.tsv"
)

ROUTE = Path(
    "P5_STAGE1A_FULL_TARGET_PROTEIN_CDS_ROUTE.tsv"
)

HUMAN_CDS = Path(
    "P5_STAGE1A_FULL_HUMAN_RELEASE104_CDS_FROZEN.fa"
)

ALIGN_CACHE = Path(
    "ENSEMBL_FULL_ALIGNMENT_CACHE"
)

CDS_CACHE = Path(
    "ENSEMBL_FULL_CDS_CACHE"
)

for p in (
    TMD,
    PROJECTION,
    PROJECTION_SUMMARY,
    ROUTE,
    HUMAN_CDS,
    ALIGN_CACHE,
    CDS_CACHE,
):
    if not p.exists():
        raise FileNotFoundError(p)

# ============================================================
# BLINDING
# ============================================================
# NO P3 spatial outcome.
# NO P4 evolutionary result.
# NO bridge analysis.
#
# Primary predictor and window are frozen from pilot.
# ============================================================

CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',

    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',

    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',

    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

FOURFOLD_AA = {"A", "G", "P", "T", "V"}
STOP = {"TAA", "TAG", "TGA"}

WINDOW_PLUS_START = 30
WINDOW_PLUS_END = 90
WINDOW_LENGTH = 61


def clean_seq(x):
    return (
        str(x)
        .replace("-", "")
        .replace(".", "")
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
        .upper()
    )


def clean_alignment(x):
    return (
        str(x)
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
        .upper()
    )


def load_gz(path):
    with gzip.open(
        path,
        "rt",
        encoding="utf-8"
    ) as fh:
        return json.load(fh)


def raw_codons(cds):
    cds = clean_seq(cds)

    if len(cds) % 3 != 0:
        return None

    return [
        cds[i:i+3]
        for i in range(0, len(cds), 3)
    ]


def codons_for_protein_length(cds, protein_length):
    xs = raw_codons(cds)

    if xs is None:
        return None

    if len(xs) == protein_length:
        return xs

    if (
        len(xs) == protein_length + 1
        and xs[-1] in STOP
    ):
        return xs[:-1]

    return None


def parse_human_fasta(path):
    by_tx = {}

    header = None
    parts = []

    def flush():
        nonlocal header, parts

        if header is None:
            return

        tx = header.split("|")[0]

        if tx in by_tx:
            raise RuntimeError(
                f"duplicate human CDS transcript: {tx}"
            )

        by_tx[tx] = "".join(parts).upper()

    with open(path, encoding="utf-8") as fh:

        for line in fh:

            line = line.rstrip()

            if line.startswith(">"):

                flush()

                header = line[1:]
                parts = []

            else:

                parts.append(
                    line.strip()
                )

        flush()

    return by_tx


# ============================================================
# 1. Exact-source formal TMD universe
# ============================================================

with open(
    PROJECTION_SUMMARY,
    encoding="utf-8"
) as fh:

    ps = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )

if len(ps) != 9827:
    raise RuntimeError(
        f"expected 9827 formal TMDs, got {len(ps)}"
    )

formal_keys = {
    (
        r["human_ensembl_gene"],
        r["human_ensembl_transcript"],
        int(r["tmd_index"])
    )
    for r in ps
}

formal_genes = {
    r["human_ensembl_gene"]
    for r in ps
}

if len(formal_genes) != 1921:
    raise RuntimeError(
        f"expected 1921 formal genes, got {len(formal_genes)}"
    )


# ============================================================
# 2. Human TMD coordinates
# ============================================================

tmd_by_gene = defaultdict(list)

with open(
    TMD,
    encoding="utf-8"
) as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:

        key = (
            row["human_ensembl_gene"],
            row["human_ensembl_transcript"],
            int(row["tmd_index"])
        )

        if key not in formal_keys:
            continue

        tmd_by_gene[
            row["human_ensembl_gene"]
        ].append({
            "gene":
                row["human_ensembl_gene"],

            "transcript":
                row["human_ensembl_transcript"],

            "tmd_index":
                int(row["tmd_index"]),

            "is_first":
                int(row["is_first"]),

            "is_internal":
                int(row["is_internal"]),

            "tmd_start0":
                int(row["human_tmd_aa_start0"]),

            "tmd_end0_exclusive":
                int(
                    row[
                        "human_tmd_aa_end0_exclusive"
                    ]
                ),
        })

for gene in tmd_by_gene:
    tmd_by_gene[gene].sort(
        key=lambda x:x["tmd_index"]
    )

if sum(
    len(x)
    for x in tmd_by_gene.values()
) != 9827:

    raise RuntimeError(
        "formal TMD coordinate count mismatch"
    )


# ============================================================
# 3. Human frozen CDS
# ============================================================

human_cds_by_tx = parse_human_fasta(
    HUMAN_CDS
)


# ============================================================
# 4. Frozen projection rows
# ============================================================

projection_by_gene_tmd = defaultdict(list)

with open(
    PROJECTION,
    encoding="utf-8"
) as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:

        gene = row[
            "human_ensembl_gene"
        ]

        if gene not in formal_genes:
            continue

        key = (
            gene,
            int(row["tmd_index"])
        )

        projection_by_gene_tmd[
            key
        ].append(row)


# ============================================================
# 5. Frozen target protein/CDS routes
# ============================================================

route = {}

with open(
    ROUTE,
    encoding="utf-8"
) as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:

        pid = row[
            "target_protein_id"
        ]

        if pid in route:
            raise RuntimeError(
                f"duplicate route protein: {pid}"
            )

        route[pid] = row

if len(route) != 21397:
    raise RuntimeError(
        f"expected 21397 route proteins, got {len(route)}"
    )


# ============================================================
# 6. Main scoring, processed one human gene at a time
# ============================================================

species_detail = []
tmd_summary = []

status_counts = Counter()
pair_failure_counts = Counter()

genes = sorted(
    formal_genes
)


for gi, gene in enumerate(
    genes,
    start=1
):

    gene_tmds = tmd_by_gene[
        gene
    ]

    if not gene_tmds:
        raise RuntimeError(
            f"no TMDs for {gene}"
        )

    txs = {
        x["transcript"]
        for x in gene_tmds
    }

    if len(txs) != 1:
        raise RuntimeError(
            f"non-unique transcript for {gene}"
        )

    tx = next(
        iter(txs)
    )

    if tx not in human_cds_by_tx:
        raise RuntimeError(
            f"human frozen CDS missing: {tx}"
        )

    # --------------------------------------------------------
    # Load Compara protein alignment once for this gene.
    # --------------------------------------------------------

    aln_path = (
        ALIGN_CACHE /
        f"{gene}.protein_alignment.json.gz"
    )

    obj = load_gz(
        aln_path
    )

    entries = obj.get(
        "data",
        []
    )

    if len(entries) != 1:
        raise RuntimeError(
            f"unexpected alignment schema: {gene}"
        )

    homologies = entries[
        0
    ].get(
        "homologies",
        []
    )

    homology_lookup = {}

    human_protein_candidates = set()

    for h in homologies:

        if h.get(
            "type"
        ) != "ortholog_one2one":
            continue

        src = h.get(
            "source",
            {}
        )

        tgt = h.get(
            "target",
            {}
        )

        sp = tgt.get(
            "species",
            ""
        )

        pid = str(
            tgt.get(
                "protein_id",
                ""
            )
        )

        if not sp or not pid:
            continue

        src_aln = clean_alignment(
            src.get(
                "align_seq",
                ""
            )
        )

        tgt_aln = clean_alignment(
            tgt.get(
                "align_seq",
                ""
            )
        )

        if (
            not src_aln
            or not tgt_aln
            or len(src_aln) != len(tgt_aln)
        ):
            continue

        human_protein = clean_seq(
            src_aln
        )

        human_protein_candidates.add(
            human_protein
        )

        homology_lookup[
            (sp, pid)
        ] = (
            src_aln,
            tgt_aln
        )

    if len(
        human_protein_candidates
    ) != 1:

        raise RuntimeError(
            f"non-unique human source protein for {gene}"
        )

    human_protein = next(
        iter(
            human_protein_candidates
        )
    )

    human_codons = codons_for_protein_length(
        human_cds_by_tx[
            tx
        ],
        len(
            human_protein
        )
    )

    if human_codons is None:
        raise RuntimeError(
            f"human CDS/protein length mismatch: {gene} {tx}"
        )

    # --------------------------------------------------------
    # Build each gene×species evolutionary background once.
    # --------------------------------------------------------

    species_to_pid = {}

    for t in gene_tmds:

        rows = projection_by_gene_tmd[
            (
                gene,
                t["tmd_index"]
            )
        ]

        for prow in rows:

            if prow[
                "mapping_pass"
            ] != "1":
                continue

            sp = prow[
                "target_species"
            ]

            pid = prow[
                "target_protein_id"
            ]

            if not pid:
                continue

            if sp in species_to_pid:

                if species_to_pid[
                    sp
                ] != pid:

                    raise RuntimeError(
                        f"inconsistent target protein "
                        f"{gene} {sp}"
                    )

            else:

                species_to_pid[
                    sp
                ] = pid


    pair_data = {}


    for sp,pid in sorted(
        species_to_pid.items()
    ):

        rt = route.get(
            pid
        )

        if rt is None:

            pair_failure_counts[
                "TARGET_ROUTE_MISSING"
            ] += 1

            continue

        if rt[
            "route_status"
        ] != "PASS":

            pair_failure_counts[
                "TARGET_ROUTE_NOT_PASS:"
                + rt[
                    "route_status"
                ]
            ] += 1

            continue

        parent = rt[
            "parent_transcript"
        ]

        if not parent:

            pair_failure_counts[
                "TARGET_PARENT_MISSING"
            ] += 1

            continue

        cds_path = (
            CDS_CACHE /
            f"{parent}.cds.json.gz"
        )

        if not cds_path.exists():

            pair_failure_counts[
                "TARGET_CDS_CACHE_MISSING"
            ] += 1

            continue

        aln_pair = homology_lookup.get(
            (
                sp,
                pid
            )
        )

        if aln_pair is None:

            pair_failure_counts[
                "HOMOLOGY_ALIGNMENT_MISSING"
            ] += 1

            continue

        src_aln,tgt_aln = aln_pair

        if clean_seq(
            src_aln
        ) != human_protein:

            pair_failure_counts[
                "SOURCE_NOT_FROZEN_HUMAN_PROTEIN"
            ] += 1

            continue

        target_protein = clean_seq(
            tgt_aln
        )

        cds_obj = load_gz(
            cds_path
        )

        target_codons = codons_for_protein_length(
            cds_obj.get(
                "seq",
                ""
            ),
            len(
                target_protein
            )
        )

        if target_codons is None:

            pair_failure_counts[
                "TARGET_CDS_PROTEIN_LENGTH_MISMATCH"
            ] += 1

            continue

        source_pos = -1
        target_pos = -1

        by_human_position = {}

        background = defaultdict(
            lambda: {
                "N":0,
                "M":0
            }
        )

        for col in range(
            len(
                src_aln
            )
        ):

            s = src_aln[
                col
            ]

            ta = tgt_aln[
                col
            ]

            if s not in "-.":
                source_pos += 1

            if ta not in "-.":
                target_pos += 1

            if s in "-.":
                continue

            target_aa = (
                None
                if ta in "-."
                else ta
            )

            rec = {
                "target_aa":
                    target_aa,

                "target_pos":
                    (
                        None
                        if ta in "-."
                        else target_pos
                    ),
            }

            by_human_position[
                source_pos
            ] = rec

            if ta in "-.":
                continue

            if s != ta:
                continue

            if s not in FOURFOLD_AA:
                continue

            if not (
                0 <= source_pos
                < len(
                    human_codons
                )
            ):
                continue

            if not (
                0 <= target_pos
                < len(
                    target_codons
                )
            ):
                continue

            hc = human_codons[
                source_pos
            ]

            tc = target_codons[
                target_pos
            ]

            if CODE.get(
                hc
            ) != s:

                continue

            if CODE.get(
                tc
            ) != ta:

                continue

            third_match = int(
                hc[2]
                ==
                tc[2]
            )

            rec.update({
                "fourfold_informative":
                    True,

                "aa":
                    s,

                "human_codon":
                    hc,

                "target_codon":
                    tc,

                "third_match":
                    third_match,
            })

            background[
                s
            ][
                "N"
            ] += 1

            background[
                s
            ][
                "M"
            ] += third_match


        pair_data[
            sp
        ] = {
            "pid":
                pid,

            "by_human_position":
                by_human_position,

            "background":
                dict(
                    background
                ),
        }


    # --------------------------------------------------------
    # Score TMDs for this gene
    # --------------------------------------------------------

    for tmd in gene_tmds:

        idx = tmd[
            "tmd_index"
        ]

        protein_len = len(
            human_protein
        )

        window_start0 = (
            tmd[
                "tmd_end0_exclusive"
            ]
            +
            WINDOW_PLUS_START
            - 1
        )

        window_end0_exclusive = (
            tmd[
                "tmd_end0_exclusive"
            ]
            +
            WINDOW_PLUS_END
        )

        geometry_ok = (
            window_start0 >= 0
            and
            window_end0_exclusive
            <= protein_len
            and
            (
                window_end0_exclusive
                -
                window_start0
            )
            == WINDOW_LENGTH
        )

        if geometry_ok:

            human_window_codons = (
                human_codons[
                    window_start0:
                    window_end0_exclusive
                ]
            )

            human_nt = "".join(
                human_window_codons
            )

            gc3 = (
                sum(
                    c[2] in {
                        "G",
                        "C"
                    }
                    for c in human_window_codons
                )
                /
                len(
                    human_window_codons
                )
            )

            cpg = (
                sum(
                    human_nt[
                        i:i+2
                    ] == "CG"
                    for i in range(
                        len(
                            human_nt
                        ) - 1
                    )
                )
                /
                max(
                    len(
                        human_nt
                    ) - 1,
                    1
                )
            )

        else:

            gc3 = None
            cpg = None


        valid_z = []
        informative_valid = []
        aa_identity_valid = []


        if geometry_ok:

            for prow in projection_by_gene_tmd[
                (
                    gene,
                    idx
                )
            ]:

                if prow[
                    "mapping_pass"
                ] != "1":
                    continue

                sp = prow[
                    "target_species"
                ]

                pid = prow[
                    "target_protein_id"
                ]

                pair = pair_data.get(
                    sp
                )

                if (
                    pair is None
                    or
                    pair[
                        "pid"
                    ]
                    != pid
                ):
                    continue

                obs = 0.0
                exp = 0.0
                variance = 0.0
                informative = 0

                aa_k = Counter()
                aa_x = Counter()

                aa_same = 0
                aa_non_gap = 0


                for hp in range(
                    window_start0,
                    window_end0_exclusive
                ):

                    rec = pair[
                        "by_human_position"
                    ].get(
                        hp
                    )

                    if rec is None:
                        continue

                    taa = rec.get(
                        "target_aa"
                    )

                    if taa is not None:

                        aa_non_gap += 1

                        if (
                            taa
                            ==
                            human_protein[
                                hp
                            ]
                        ):

                            aa_same += 1


                    if not rec.get(
                        "fourfold_informative",
                        False
                    ):
                        continue

                    aa = rec[
                        "aa"
                    ]

                    informative += 1
                    aa_k[
                        aa
                    ] += 1

                    aa_x[
                        aa
                    ] += int(
                        rec[
                            "third_match"
                        ]
                    )


                if informative < 5:

                    status = (
                        "LT5_INFORMATIVE_4D_SITES"
                    )

                    z = None

                else:

                    for aa,k in aa_k.items():

                        bg = pair[
                            "background"
                        ].get(
                            aa,
                            {
                                "N":0,
                                "M":0
                            }
                        )

                        N = bg[
                            "N"
                        ]

                        M = bg[
                            "M"
                        ]

                        x = aa_x[
                            aa
                        ]

                        if (
                            N < k
                            or
                            N <= 0
                        ):
                            continue

                        obs += x

                        p = (
                            M / N
                        )

                        exp += (
                            k * p
                        )

                        if N > 1:

                            variance += (
                                k
                                *
                                p
                                *
                                (
                                    1-p
                                )
                                *
                                (
                                    (
                                        N-k
                                    )
                                    /
                                    (
                                        N-1
                                    )
                                )
                            )


                    if variance <= 0:

                        status = (
                            "NULL_VARIANCE_ZERO"
                        )

                        z = None

                    else:

                        z = (
                            obs-exp
                        ) / math.sqrt(
                            variance
                        )

                        status = "PASS"


                status_counts[
                    status
                ] += 1


                row = {
                    "human_ensembl_gene":
                        gene,

                    "human_ensembl_transcript":
                        tx,

                    "tmd_index":
                        idx,

                    "is_first":
                        tmd[
                            "is_first"
                        ],

                    "is_internal":
                        tmd[
                            "is_internal"
                        ],

                    "target_species":
                        sp,

                    "target_protein_id":
                        pid,

                    "window_start0":
                        window_start0,

                    "window_end0_exclusive":
                        window_end0_exclusive,

                    "window_length_codons":
                        WINDOW_LENGTH,

                    "informative_4d_sites":
                        informative,

                    "observed_exact_third_base_matches":
                        (
                            obs
                            if z is not None
                            else ""
                        ),

                    "expected_matches":
                        (
                            exp
                            if z is not None
                            else ""
                        ),

                    "null_variance":
                        (
                            variance
                            if z is not None
                            else ""
                        ),

                    "excess_conservation_z":
                        (
                            z
                            if z is not None
                            else ""
                        ),

                    "window_aa_non_gap_fraction":
                        (
                            aa_non_gap
                            /
                            WINDOW_LENGTH
                        ),

                    "window_aa_identity_fraction":
                        (
                            aa_same
                            /
                            WINDOW_LENGTH
                        ),

                    "human_window_GC3":
                        gc3,

                    "human_window_CpG_density":
                        cpg,

                    "status":
                        status,
                }

                species_detail.append(
                    row
                )

                if status == "PASS":

                    valid_z.append(
                        z
                    )

                    informative_valid.append(
                        informative
                    )

                    aa_identity_valid.append(
                        aa_same
                        /
                        WINDOW_LENGTH
                    )


        score_available = (
            geometry_ok
            and
            len(
                valid_z
            ) >= 6
        )


        tmd_summary.append({
            "human_ensembl_gene":
                gene,

            "human_ensembl_transcript":
                tx,

            "tmd_index":
                idx,

            "is_first":
                tmd[
                    "is_first"
                ],

            "is_internal":
                tmd[
                    "is_internal"
                ],

            "primary_window_geometry_ok":
                int(
                    geometry_ok
                ),

            "valid_species_z":
                len(
                    valid_z
                ),

            "primary_score_available":
                int(
                    score_available
                ),

            "POST30_90_4D_EXCESS_CONSERVATION_Z_MEDIAN":
                (
                    median(
                        valid_z
                    )
                    if score_available
                    else ""
                ),

            "POST30_90_4D_EXCESS_CONSERVATION_Z_MEAN":
                (
                    mean(
                        valid_z
                    )
                    if score_available
                    else ""
                ),

            "median_informative_4d_sites":
                (
                    median(
                        informative_valid
                    )
                    if informative_valid
                    else ""
                ),

            "median_window_aa_identity":
                (
                    median(
                        aa_identity_valid
                    )
                    if aa_identity_valid
                    else ""
                ),

            "human_window_GC3":
                (
                    gc3
                    if gc3 is not None
                    else ""
                ),

            "human_window_CpG_density":
                (
                    cpg
                    if cpg is not None
                    else ""
                ),
        })


    if (
        gi % 100 == 0
        or
        gi == len(
            genes
        )
    ):

        print(
            f"SCORED_GENES "
            f"{gi}/{len(genes)}",
            flush=True
        )


# ============================================================
# 7. Full census
# ============================================================

if len(
    tmd_summary
) != 9827:

    raise RuntimeError(
        f"expected 9827 TMD score rows, "
        f"got {len(tmd_summary)}"
    )


internal = [
    x
    for x in tmd_summary
    if x[
        "is_internal"
    ] == 1
]


if len(
    internal
) != 8420:

    raise RuntimeError(
        f"expected 8420 internal TMDs, "
        f"got {len(internal)}"
    )


geometry_internal = [
    x
    for x in internal
    if x[
        "primary_window_geometry_ok"
    ] == 1
]


if len(
    geometry_internal
) != 6523:

    raise RuntimeError(
        f"expected 6523 full-window internal TMDs, "
        f"got {len(geometry_internal)}"
    )


score_internal = [
    x
    for x in internal
    if x[
        "primary_score_available"
    ] == 1
]


fraction_all = (
    len(
        score_internal
    )
    /
    len(
        internal
    )
)


fraction_geometry = (
    len(
        score_internal
    )
    /
    len(
        geometry_internal
    )
)


# ============================================================
# 8. Output
# ============================================================

species_detail.sort(
    key=lambda x:(
        x[
            "human_ensembl_gene"
        ],
        int(
            x[
                "tmd_index"
            ]
        ),
        x[
            "target_species"
        ]
    )
)


tmd_summary.sort(
    key=lambda x:(
        x[
            "human_ensembl_gene"
        ],
        int(
            x[
                "tmd_index"
            ]
        )
    )
)


with open(
    "P5_STAGE1A_FULL_GENE_SPECIES_TMD_EVOLUTION.tsv",
    "w",
    newline="",
    encoding="utf-8"
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            species_detail[
                0
            ].keys()
        ),
        delimiter="\t",
        lineterminator="\n"
    )

    w.writeheader()
    w.writerows(
        species_detail
    )


with open(
    "P5_STAGE1A_FULL_TMD_EVOLUTION_SCORE.tsv",
    "w",
    newline="",
    encoding="utf-8"
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            tmd_summary[
                0
            ].keys()
        ),
        delimiter="\t",
        lineterminator="\n"
    )

    w.writeheader()
    w.writerows(
        tmd_summary
    )


lines = [
    "P5 STAGE1A FULL EVOLUTIONARY SCORE",
    "=" * 72,
    "",
    f"FORMAL_GENES = {len(genes)}",
    f"FORMAL_TMD = {len(tmd_summary)}",
    f"FORMAL_INTERNAL_TMD = {len(internal)}",
    "",
    "PRIMARY_WINDOW = TMD_end+30..+90",
    "PRIMARY_WINDOW_LENGTH_CODONS = 61",
    "",
    f"INTERNAL_TMD_WITH_FULL_PRIMARY_WINDOW = {len(geometry_internal)}",
    f"INTERNAL_TMD_WITH_PRIMARY_SCORE = {len(score_internal)}",
    f"FRACTION_ALL_INTERNAL_TMD_WITH_PRIMARY_SCORE = {fraction_all:.6f}",
    f"FRACTION_GEOMETRY_ELIGIBLE_INTERNAL_TMD_WITH_PRIMARY_SCORE = {fraction_geometry:.6f}",
    "",
    "GENE_SPECIES_TMD_STATUS_COUNTS:"
]


for k,v in sorted(
    status_counts.items()
):

    lines.append(
        f"  {k}\t{v}"
    )


lines += [
    "",
    "PAIR_BACKGROUND_FAILURE_COUNTS:"
]


if pair_failure_counts:

    for k,v in sorted(
        pair_failure_counts.items()
    ):

        lines.append(
            f"  {k}\t{v}"
        )

else:

    lines.append(
        "  NONE"
    )


lines += [
    "",
    "PRIMARY_SCORE = POST30_90_4D_EXCESS_CONSERVATION_Z_MEDIAN",
    "MIN_VALID_SPECIES_PER_TMD = 6",
    "MIN_INFORMATIVE_4D_SITES_PER_GENE_SPECIES_TMD = 5",
    "",
    "PRIMARY_SCORE_DEFINITION_CHANGED_FROM_PILOT = FALSE",
    "PRIMARY_WINDOW_CHANGED_FROM_PILOT = FALSE",
    "P3_SPATIAL_OUTCOME_LOADED = FALSE",
    "P4_EVOLUTIONARY_RESULT_LOADED = FALSE",
]


Path(
    "P5_STAGE1A_FULL_EVOLUTION_SCORE_CENSUS.txt"
).write_text(
    "\n".join(
        lines
    ) + "\n",
    encoding="utf-8"
)


print()
print(
    "\n".join(
        lines
    ),
    flush=True
)
